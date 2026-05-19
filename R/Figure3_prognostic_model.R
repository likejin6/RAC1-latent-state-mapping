############################################################
# RAC1-latent-state-mapping
# Figure3_prognostic_model.R
#   Figure 3: TCGA-LUAD prognostic modeling and validation
#             Fig.3A (LASSO), Fig.3B (KM-curve), Fig.3C (ROC)
#
# INPUT：
#   - RAC1_TopDRgenes.csv
#   - 下载 TCGA-LUAD 数据
#
# OUTPUT：
#   - TCGA_LUAD_raw.rds
#   - Figure3A_LASSO.pdf
#   - Figure3B_KM.pdf
#   - Time_dependent_ROC.pdf
#   - 相关CSV表格
############################################################



library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(survival)
library(glmnet)
library(survminer)
library(SummarizedExperiment)
library(survcomp)
library(timeROC)

# 检查必要文件
if(!file.exists("RAC1_TopDRgenes.csv")) stop("请先运行 0_preprocessing_and_vKO.R 生成 RAC1_TopDRgenes.csv")

######
## Step 26：下载TCGA-LUAD表达数据
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp)
luad_exp <- GDCprepare(query_exp)

saveRDS(
  luad_exp,
  "TCGA_LUAD_raw.rds"
)

######
## Step 27：处理矩阵数据
# expression matrix
luad_exp <- readRDS("TCGA_LUAD_raw.rds")
exp_mat <- assay(luad_exp,"fpkm_uq_unstrand")
# 去除NA
gene_names <- rowData(luad_exp)$gene_name
keep <- !is.na(gene_names)
exp_mat <- exp_mat[keep, ]
gene_names <- gene_names[keep]
# 设置行名为基因名
rownames(exp_mat) <- gene_names
# 去除重复基因
exp_mat <- exp_mat[!duplicated(rownames(exp_mat)), ]
# 标准化
exp_mat <- log2(exp_mat + 1)
# 矩阵转置
exp_df <- as.data.frame(t(exp_mat))

######
## Step 28：临床数据整合+生存变量构建
# 提取患者ID
exp_df$patient <- substr(rownames(exp_df), 1, 12)
# 提取临床信息
clinical <- as.data.frame(colData(luad_exp))
clinical2 <- clinical[, c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")]
# 构建总生存时间(OS.time)和生存状态(OS)
clinical2$OS.time <- ifelse(is.na(clinical2$days_to_death), clinical2$days_to_last_follow_up, clinical2$days_to_death)
clinical2$OS <- ifelse(clinical2$vital_status == "Dead", 1, 0)
# 合并表达矩阵和临床数据
merged_df <- merge(exp_df, clinical2, by.x = "patient", by.y = "submitter_id")
# 去除NA
merged_df <- merged_df[!is.na(merged_df$OS.time), ]

######
## Step 29：筛选 RAC1 DRgenes
# extract DRgenes
rac1_dr <- read.csv("RAC1_TopDRgenes.csv")
rac1_sig <- rac1_dr %>%
  dplyr::filter(
    p.adj < 0.05
  ) %>%
  dplyr::arrange(desc(distance))
candidate_genes <- rac1_sig$gene
# 去除RAC1自身
candidate_genes <- setdiff(candidate_genes,"RAC1")
candidate_genes <- candidate_genes[candidate_genes %in% colnames(merged_df)]
# 过滤低表达基因
expr_keep <- sapply(candidate_genes, function(g){
  mean(merged_df[[g]] > 1, na.rm = TRUE) > 0.3
})
filtered_genes <- candidate_genes[expr_keep]
# 取Top20作为核心标签
rac1_sig2 <- rac1_sig[rac1_sig$gene %in% filtered_genes, ]
sig_genes <- head(rac1_sig2$gene, 20)

######
## Step 30：构建建模数据集
sig_data <- merged_df[, c("patient", "OS.time", "OS", sig_genes)]
# 去除NA
sig_data <- na.omit(sig_data) 
table(sig_data$OS)

######
## Step 31：单因素Cox回归筛选预后基因
uni_results <- data.frame()
for(g in sig_genes){
  formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", g))
  fit <- coxph(formula, data = sig_data)
  s <- summary(fit)
  uni_results <- rbind(uni_results, data.frame(gene = g, HR = s$coef[1,2], p.value = s$coef[1,5]))
}

# 筛选P<0.1的基因用于LASSO
uni_results <- uni_results[order(uni_results$p.value), ]
write.csv(uni_results, "RAC1_signature_univariate_cox.csv", row.names = FALSE)

cox_genes <- uni_results %>% filter(p.value < 0.1) %>% pull(gene)

######
## Step 32：LASSO-Cox回归（risk model 构建）
# 去除OS.time≤0的样本
sig_data2 <- sig_data %>% filter(OS.time > 0)
# 构建LASSO输入矩阵
x <- as.matrix(sig_data2[, cox_genes])
y <- Surv(sig_data2$OS.time, sig_data2$OS)

# LASSO建模
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 5, cox.ties = "breslow")

pdf("Figure3A_LASSO.pdf", width = 7, height = 6)
plot(cvfit)
dev.off()

# 提取lambda.min对应的非零系数基因
coef_min <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(coef_min)[coef_min[,1] != 0]
write.csv(data.frame(gene = lasso_genes), "RAC1_LASSO_genes.csv", row.names = FALSE)

######
## Step 33：多因素Cox回归构建最终模型（risk model 构建）
formula_str <- paste(lasso_genes, collapse = " + ")
formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", formula_str))
multi_cox <- coxph(formula, data = sig_data2)
multi_summary <- summary(multi_cox)

multi_results <- data.frame(gene = rownames(multi_summary$coef), HR = multi_summary$coef[,2], p.value = multi_summary$coef[,5])
write.csv(multi_results, "RAC1_multivariate_cox.csv", row.names = FALSE)

######
## Step 34：risk score + KM curve
# 计算risk score
risk_score <- predict(multi_cox, type = "risk")
sig_data2$risk_score <- risk_score
# 按中位数分为高低风险组
sig_data2$risk_group <- ifelse(sig_data2$risk_score > median(sig_data2$risk_score), "High-risk", "Low-risk")
table(sig_data2$risk_group)
# KM curve
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = sig_data2)
km_plot <- ggsurvplot(
  fit, data = sig_data2, pval = TRUE, risk.table = TRUE,
  conf.int = FALSE, surv.median.line = "hv",
  legend.title = "Risk group", legend.labs = c("High-risk", "Low-risk"),
  palette = c("#E41A1C", "#377EB8"),
  xlab = "Time (days)", ylab = "Overall survival probability"
)

pdf("Figure3B_KM.pdf", width = 7, height = 6)
print(km_plot)
dev.off()

write.csv(sig_data2, "RAC1_tumor_state_risk_model.csv", row.names = FALSE)
# 风险评分表
risk_table <- sig_data2[, c("patient", "OS.time", "OS", "risk_score", "risk_group")]
write.csv(risk_table, "RAC1_risk_scores.csv", row.names = FALSE)

####################################
# 验证模型
####################################


######
## Step 35：提取临床变量
clinical_info <- clinical[, c(
  "submitter_id",
  "age_at_diagnosis",
  "gender",
  "ajcc_pathologic_stage"
)]
#去除重复患者ID
clinical_info <- clinical_info[
  !duplicated(clinical_info$submitter_id),
]

######
## Step 36：合并临床信息与风险评分表
clinical_risk <- merge(
  risk_table,
  clinical_info,
  by.x = "patient",
  by.y = "submitter_id"
)

dim(clinical_risk)
head(clinical_risk)

######
## Step 37：将age换算为year
clinical_risk$age <-
  clinical_risk$age_at_diagnosis / 365

summary(clinical_risk$age)

######
## Step 38：简化病理分期
clinical_risk$stage_simple <-
  gsub("A|B", "", clinical_risk$ajcc_pathologic_stage)
#去除stage
clinical_risk$stage_simple <-
  gsub("Stage ", "", clinical_risk$stage_simple)

######
## Step 39：变量转换为临床影响预后的因素
clinical_risk$stage_simple <-
  factor(clinical_risk$stage_simple)

clinical_risk$gender <-
  factor(clinical_risk$gender)

table(clinical_risk$stage_simple)
#去除NA
clinical_risk <- clinical_risk %>%
  dplyr::filter(
    !is.na(age),
    !is.na(gender),
    !is.na(stage_simple)
  )

clinical_risk <- clinical_risk %>%
  dplyr::filter(OS.time > 0)

######
## Step 40：单因素Cox检验risk score的预后价值
uni_cox <- coxph(
  Surv(OS.time, OS) ~ risk_score,
  data = clinical_risk
)

summary(uni_cox)

######
## Step 41：多因素Cox进行 risk score 的独立性检验
multi_cox <- coxph(
  Surv(OS.time, OS) ~
    age +
    gender +
    stage_simple +
    risk_score,
  data = clinical_risk
)

summary(multi_cox)
multi_summary <- summary(multi_cox)

multi_results <- data.frame(
  variable = rownames(multi_summary$coef),
  HR = multi_summary$coef[,2],
  lower95 = multi_summary$conf.int[,3],
  upper95 = multi_summary$conf.int[,4],
  p.value = multi_summary$coef[,5]
)

multi_results
# 保存结果
write.csv(
  multi_results,
  "Independent_Cox_results.csv",
  row.names = FALSE
)

write.csv(
  clinical_risk,
  "Clinical_Risk_Integrated_Table.csv",
  row.names = FALSE
)

######
## Step 42：Time-dependent ROC analysis
##构图：Figure3C_ROC
roc_res <- timeROC(
  T = clinical_risk$OS.time,
  delta = clinical_risk$OS,
  marker = clinical_risk$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(365, 1095, 1825),
  iid = TRUE
)

pdf("Time_dependent_ROC.pdf",width = 6,height = 6)
plot(
  roc_res,
  time = 365,
  col = "red",
  lwd = 2,
  title = FALSE
)

plot(
  roc_res,
  time = 1095,
  col = "blue",
  lwd = 2,
  add = TRUE
)

plot(
  roc_res,
  time = 1825,
  col = "darkgreen",
  lwd = 2,
  add = TRUE
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "1-year AUC = ",
      round(roc_res$AUC[1], 3)
    ),
    
    paste0(
      "3-year AUC = ",
      round(roc_res$AUC[2], 3)
    ),
    
    paste0(
      "5-year AUC = ",
      round(roc_res$AUC[3], 3)
    )
  ),
  
  col = c(
    "red",
    "blue",
    "darkgreen"
  ),
  
  lwd = 2
)
dev.off()

# save AUC table
auc_table <- data.frame(
  Time = c(
    "1-year",
    "3-year",
    "5-year"
  ),
  AUC = roc_res$AUC
)

write.csv(
  auc_table,
  "TimeROC_AUC_results.csv",
  row.names = FALSE
)

roc_res$AUC

######
## Step 43：C-index analysis
cindex_res <- concordance.index(
  x = clinical_risk$risk_score,
  surv.time = clinical_risk$OS.time,
  surv.event = clinical_risk$OS,
  method = "noether"
)

cindex_res$c.index
# save C-index 
cindex_table <- data.frame(
  C_index = cindex_res$c.index,
  lower95 = cindex_res$lower,
  upper95 = cindex_res$upper
)

write.csv(cindex_table,"C_index_results.csv",row.names = FALSE)

cat("预后模型分析完成。\n")