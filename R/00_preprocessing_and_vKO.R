############################################################
# RAC1-latent-state-mapping
# Step 00: preprocessing, candidate gene screening, and virtual KO runs
#
# INPUT：
#   - GSE131907_UMI_matrix.rds (原始 UMI 矩阵)
#   - GSE131907_annotation.txt.gz (细胞注释)
# OUTPUT：
#   - expr_sub.rds (下采样并过滤后的表达矩阵_5Kcells)
#   - {gene}_vko.rds (每个候选基因的 scTenifoldKnk 结果)
#   - {gene}_TopDRgenes.csv (每个基因的 top 100 DR genes)
#   - WT_GRN_matrix.rds, RAC1_vKO_GRN_matrix.rds (由后续Fig1脚本生成，但这里也保存)
############################################################



library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(igraph)
library(scTenifoldKnk)

# 检查必要文件
if(!file.exists("GSE131907_UMI_matrix.rds")){
  stop("请将 GSE131907_UMI_matrix.rds 放在工作目录下。")
}
if(!file.exists("GSE131907_annotation.txt.gz")){
  stop("请将 GSE131907_annotation.txt.gz 放在工作目录下。")
}

######
## Step 1
## 读入矩阵和注释,用作者annotation筛选 tLUNG & Epithelial cells
anno <- fread("GSE131907_annotation.txt.gz")

tumor_cells <- anno %>%
  filter(
    Sample_Origin == "tLung",
    Cell_type.refined == "Epithelial cells"
  )

cells_use <- tumor_cells$Index

counts <- readRDS("GSE131907_UMI_matrix.rds")
cells_use <- intersect(cells_use, colnames(counts))
expr <- counts[, cells_use]

dim(expr)

######
## Step 2
## 下采样抽取5000cells（减少运算负担）,去低表达基因（至少在25个细胞中表达）
set.seed(123)
cells_sub <- sample(colnames(expr), 5000)
expr_sub <- expr[, cells_sub]

# 保存小数据（成为以后读入数据的文件）
saveRDS(expr_sub, file = "GSE131907_tumor_epithelial_5k.rds")

# 去低表达基因
gene_filter <- rowSums(expr_sub > 0) >= 25
expr_sub <- expr_sub[gene_filter, ]

######
## Step 3：Seurat标准化，找高变基因（HVG）2000个
seu <- CreateSeuratObject(counts = expr_sub)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)

hvg <- VariableFeatures(seu)

######
## Step 4：计算所有基因的平均表达，取高表达基因top2000
avg_exp <- rowMeans(expr_sub)
high_exp <- names(sort(avg_exp, decreasing = TRUE))[1:2000]

######
## Step 5： 取hvg和high_exp的交集筛选candidate genes
candidate_genes <- intersect(hvg, high_exp)

length(candidate_genes)

######
## Step 6：准备vKO矩阵（即scTenifoldKnk中用于构建GRN的矩阵）
expr_vko <- expr_sub[candidate_genes, ]
dim(expr_vko)

######
## Step 7：构建基因相关性网络
# 控制数量
candidate_genes2 <- intersect(hvg, high_exp)
if(length(candidate_genes2) > 500){
  candidate_genes2 <- candidate_genes2[1:500]
}
expr_mat <- as.matrix(expr_sub[candidate_genes2, ])
# 计算相关性
cor_mat <- abs(cor(t(expr_mat)))
diag(cor_mat) <- 0

g <- graph_from_adjacency_matrix(cor_mat, weighted = TRUE, mode = "undirected")

###### 
## Step 8 （本步骤用于辅助确定后续 multi-gene virtual knockout screening的候选扰动基因。）
## 这里的 hub ranking 仅作为参考，不作为最终 vKO genes 的唯一自动筛选依据
# 计算degree&betweenness，筛选top30 hub genes（作为候选vko基因）
deg <- degree(g)
btw <- betweenness(g, weights = 1/E(g)$weight)
# 初步查看候选基因的网络中心性排序
score <- deg + btw
# 提取 top-ranked genes 作为后续人工筛选和生物学解释的参考集合
hub_genes <- names(sort(score, decreasing = TRUE))[1:30]
hub_genes

###### 
## Step 8.5：进一步进行人工生物学功能筛选
## 以下 genes 代表一组经过 network-informed prioritization和 biological relevance 综合选择后的 candidate perturbation targets，
## 而不是由单一 hub score 完全自动筛选得到的 top genes。
genes <- c("RAC1","LDHA","AREG","ANXA2","S100A6",
          "LPCAT1","TMSB10","TAGLN2","CTSD","PON2")

###### 
## Step 9：运行vKO,获得DRgenes
genes <- intersect(genes, rownames(expr_vko))
for (g in genes) {
  
  cat("Running vKO for:", g, "\n")
  
  res <- tryCatch({
    scTenifoldKnk(
      countMatrix = expr_vko,
      gKO = g,
      qc = TRUE,
      nc_nNet = 8,
      nc_nCells = 500,
      nCores = 2
    )
  }, error = function(e){
    cat("❌ Error in", g, ":", e$message, "\n")
    return(NULL)
  })

  if(!is.null(res)){
    saveRDS(res, paste0(g, "_vko.rds"))
    cat("✔ Saved:", g, "\n")
  }
  
  rm(res)
  gc()
}

######
## Step 10：提取DRgenes
dr_list <- list()
for(g in genes){
  cat("Processing:", g, "\n")
  res <- readRDS(paste0(g, "_vko.rds"))
  # 正确对象
  dr <- res$diffRegulation
  # 按distance排序
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  # Top 100
  top_dr <- head(dr, 100)
  dr_list[[g]] <- top_dr
  write.csv(
    top_dr,
    paste0(g, "_TopDRgenes.csv"),
    row.names = FALSE
  )
}
###### 
##Step 11：比较vKO基因的 DRgene 强度
summary_df <- data.frame()

for(g in genes){
  dr <- read.csv(paste0(g, "_TopDRgenes.csv"))
  summary_df <- rbind(
    summary_df,
    data.frame(
      gene = g,
      mean_distance = mean(dr$distance),
      max_distance = max(dr$distance),
      sig_genes = sum(dr$p.adj < 0.05)
    )
  )
}

summary_df <- summary_df[
  order(summary_df$mean_distance, decreasing = TRUE),
]

summary_df
######
## Step 12：根据summary_df结果,选取强度较高的基因进行后续分析
top_vko <- c(
  "RAC1",
  "AREG",
  "LPCAT1",
  "TMSB10"
)
backup <- c(
  "CTSD",
  "LDHA"
)
##########################
##经过生物学功能分析,选取strongest vKO 基因 RAC1 进行分析
##########################

######
## Step 13：提取 RAC1 的 DRgenes
rac1_ko <- readRDS("RAC1_vko.rds")

# 提取并保存 WT 和 KO 网络矩阵（供 Figure1 和 Python 验证使用）
wt_net <- as.matrix(rac1_ko$tensorNetworks$WT)
ko_net <- as.matrix(rac1_ko$tensorNetworks$KO)
saveRDS(wt_net, "WT_GRN_matrix.rds")
saveRDS(ko_net, "RAC1_vKO_GRN_matrix.rds")

# 保存处理后的表达矩阵供后续使用
saveRDS(expr_sub, "expr_sub.rds")

cat("预处理及 vKO 全部完成。RAC1 被选为 strongest vKO gene 用于后续分析。\n")
