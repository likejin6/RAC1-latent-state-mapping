############################################################
# RAC1-latent-state-mapping
# Complete pipeline
############################################################

############################################################
###使用GEO公开数据集：GSE131907
############################################################



############################################################
#PHASE I：筛选vKO基因&确定后续分析基因——RAC1
############################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(igraph)
library(scTenifoldKnk)

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
## Step 8
## 计算degree&betweenness，筛选top30 hub genes（作为候选vko基因）
deg <- degree(g)
btw <- betweenness(g, weights = 1/E(g)$weight)
score <- deg + btw

hub_genes <- names(sort(score, decreasing = TRUE))[1:30]

hub_genes
###### 
## Step 8.5：对得到的top30 hub genes 进行人工生物学功能筛选
# 筛选扰动可能性较高的基因
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



############################################################
#PHASE II：latent tumor-state analysis
############################################################

library(igraph)
library(Matrix)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(dplyr)
rac1_ko <- readRDS("RAC1_vko.rds")

######
## Step 14：extract sciTenfoldKnk创建的 WT & vKO GRN
wt_net <- rac1_ko$tensorNetworks$WT
ko_net <- rac1_ko$tensorNetworks$KO
wt_net <- as.matrix(wt_net)
ko_net <- as.matrix(ko_net)

stopifnot(all(rownames(wt_net) == rownames(ko_net)))
genes <- rownames(wt_net)

saveRDS(wt_net,"WT_GRN_matrix.rds")
saveRDS(ko_net,"RAC1_vKO_GRN_matrix.rds")

######
## Step 15：GRN rewiring quantification
delta_net <- ko_net - wt_net
# symmetrize
delta_sym <- (delta_net + t(delta_net)) / 2
diag(delta_sym) <- 0
summary(as.vector(delta_sym))

quantile(
  abs(delta_sym[abs(delta_sym) > 0]),
  probs = c(0.5, 0.75, 0.9, 0.95, 0.99)
)

rac1_edges <- abs(delta_sym["RAC1", ])
rac1_cut <- quantile(rac1_edges[rac1_edges > 0],0.9)
delta_sym2 <- delta_sym
delta_sym2["RAC1",abs(delta_sym2["RAC1", ]) > rac1_cut] <- 0
delta_sym2[abs(delta_sym2[, "RAC1"]) > rac1_cut,"RAC1"] <- 0

cutoff <- 0.03
delta_thr <- delta_sym2
delta_thr[abs(delta_thr) < cutoff] <- 0

diag(delta_thr) <- 0

#######
## STEP 16：构图：Rewiring graph 
#Figure1A_GRN_rewiring
g_rewire <- graph_from_adjacency_matrix(
  delta_thr,
  weighted = TRUE,
  mode = "undirected",
  diag = FALSE
)
g_rewire <- delete_vertices(g_rewire,degree(g_rewire) == 0)
wt_sym <- (wt_net + t(wt_net)) / 2
ko_sym <- (ko_net + t(ko_net)) / 2

diag(wt_sym) <- 0
diag(ko_sym) <- 0

wt_cut <- quantile(abs(wt_sym[abs(wt_sym) > 0]),0.95)
ko_cut <- quantile(abs(ko_sym[abs(ko_sym) > 0]),0.95)

wt_thr <- wt_sym
ko_thr <- ko_sym

wt_thr[abs(wt_thr) < wt_cut] <- 0
ko_thr[abs(ko_thr) < ko_cut] <- 0

g_wt <- graph_from_adjacency_matrix(
  abs(wt_thr),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g_ko <- graph_from_adjacency_matrix(
  abs(ko_thr),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

deg_wt <- degree(g_wt)
deg_ko <- degree(g_ko)
deg_delta <- deg_ko - deg_wt

V(g_rewire)$delta_degree <- deg_delta[V(g_rewire)$name]

hub_thresh <- quantile(deg_ko,0.75)
V(g_rewire)$hub_status <- "unchanged"
V(g_rewire)$hub_status[
  V(g_rewire)$delta_degree > 0 &
    deg_ko[V(g_rewire)$name] >= hub_thresh
] <- "new_hub"
V(g_rewire)$hub_status[
  V(g_rewire)$delta_degree < 0 &
    deg_wt[V(g_rewire)$name] >= hub_thresh
] <- "lost_hub"

tg <- as_tbl_graph(g_rewire)

label_genes_hub <- V(g_rewire)$name[
  V(g_rewire)$hub_status %in% c("new_hub", "lost_hub")
]

label_genes_delta <- names(
  sort(
    abs(deg_delta[V(g_rewire)$name]),
    decreasing = TRUE
  )
)[1:min(12, vcount(g_rewire))]

hub_genes <- unique(
  c(label_genes_hub,label_genes_delta))

set.seed(123)
p_rewire <- ggraph(tg, layout = "fr", weights = sqrt(abs(E(g_rewire)$weight))) +
  geom_edge_link(aes(color = weight, width = abs(weight)), alpha = 0.75, lineend = "round") +
  scale_edge_color_gradient2(low = "steelblue3", mid = "grey85", high = "firebrick3",
                             midpoint = 0, name = "delta regulation") +
  scale_edge_width(range = c(0.5, 2.5), guide = "none") +
  geom_node_point(aes(size = abs(delta_degree),color = hub_status,alpha = hub_status)) +
  scale_size_continuous(range = c(4, 14),name = "delta degree") +
  scale_color_manual(values = c("new_hub" = "firebrick3","lost_hub" = "steelblue3","unchanged" = "grey80")) +
  scale_alpha_manual(values = c("new_hub" = 0.95,"lost_hub" = 0.95,"unchanged" = 0.45),guide = "none") +
  geom_node_text(aes(label = ifelse(name %in% hub_genes,name,"")),
                 repel = TRUE,
                 size = 3.5,
                 family = "sans"
  ) +
  theme_void(base_size = 16) +
  ggtitle("RAC1-centered GRN rewiring after virtual perturbation") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
        legend.position = "right"
  )

ggsave("Figure1A_GRN_rewiring.pdf", p_rewire, width = 12, height = 9)

edge_table <- as_data_frame(
  g_rewire,
  what = "edges"
)
write.csv(
  edge_table,
  "Figure1A_Rewiring_edges_Optimized.csv",
  row.names = FALSE
)

node_table <- data.frame(
  gene = names(V(g_rewire)),
  delta_degree = V(g_rewire)$delta_degree,
  hub_status = V(g_rewire)$hub_status
)
write.csv(
  node_table,
  "Figure1A_Rewiring_nodes_Optimized.csv",
  row.names = FALSE
)

#######
## STEP 17：Figure1B_coordination_network
expr_sub <- readRDS("expr_sub.rds")
rac1_dr <- read.csv("RAC1_TopDRgenes.csv")
# 去除RAC1自身
rac1_dr <- rac1_dr %>%
  filter(gene != "RAC1")
#定义 programs
#（program的特征基因来源于TopDRgenes和program的经典markergene）
alveolar_program <- c("SFTPA1", "SFTPA2", "SFTPB", "SFTPD", "NAPSA", "CTSH")
plasticity_program <- c("S100A4", "MARCKS", "KRT17", "SNCG")
epithelial_program <- c("PSCA", "S100P", "MUC4", "CD99")
stress_program <- c("GADD45B", "DDIT3", "HEXIM1", "ASS1", "SLC39A8", "LRRK2", "CDKN2A")
oncogenic_program <- c("HMGA1", "MACC1", "PRSS23", "CD24")
remodeling_program <- c("SPINT2", "C19orf33")

program_genes <- unique(
  c(
    alveolar_program,
    plasticity_program,
    epithelial_program,
    stress_program,
    oncogenic_program,
    remodeling_program
  )
)

plot_genes <- unique(
  c(
    head(rac1_dr$gene, 25),
    program_genes
  )
)

plot_genes <- intersect(
  plot_genes,
  rownames(expr_sub)
)


sub_expr <- expr_sub[rownames(expr_sub) %in% plot_genes,]
# 标准化
sub_expr_scaled <- t(scale(t(sub_expr)))
# 计算相关性
cor_mat2 <- cor(t(sub_expr_scaled),method = "pearson")
# 只保留强关系
cor_cutoff <- 0.5
cor_mat2[
  abs(cor_mat2) < cor_cutoff
] <- 0

diag(cor_mat2) <- 0

g2 <- graph_from_adjacency_matrix(
  cor_mat2,
  weighted = TRUE,
  mode = "undirected",
  diag = FALSE
)
# 删除孤立节点
g2 <- delete_vertices(g2, degree(g2) == 0)
# 定义 programs
alveolar_program <- c("SFTPA1","SFTPA2","SFTPB","SFTPD","NAPSA", "CTSH")
plasticity_program <- c("S100A4","MARCKS","KRT17","SNCG")
epithelial_program <- c("PSCA", "S100P", "MUC4", "CD99")
stress_program <- c("GADD45B","DDIT3","HEXIM1","ASS1","SLC39A8", "LRRK2", "CDKN2A")
oncogenic_program <- c("HMGA1","MACC1","PRSS23","CD24")
remodeling_program <- c("SPINT2", "C19orf33")

V(g2)$program <- ifelse(V(g2)$name %in% alveolar_program, "alveolar",
                        ifelse(V(g2)$name %in% plasticity_program, "plasticity",
                               ifelse(V(g2)$name %in% epithelial_program,"epithelial",
                                      ifelse(V(g2)$name %in% stress_program, "stress",
                                             ifelse(V(g2)$name %in% oncogenic_program, "oncogenic", 
                                                    ifelse(V(g2)$name %in% remodeling_program,"remodeling","other"))))))


# 画图颜色
program_colors <- c(
  "alveolar" = "gold",
  "plasticity" = "skyblue",
  "epithelial" = "orange",
  "stress" = "forestgreen",
  "oncogenic" = "firebrick3",
  "remodeling" = "mediumpurple",
  "other" = "grey70"
)
V(g2)$color <- program_colors[V(g2)$program]
# edge styling
E(g2)$color <- ifelse(E(g2)$weight > 0, "firebrick3", "steelblue3")
E(g2)$width <- abs(E(g2)$weight) * 4
# node size by degree
deg2 <- degree(g2)
V(g2)$size <- deg2 * 2 + 8
V(g2)$label.cex <- 0.8

set.seed(123)
layout1 <- layout_with_fr(g2, niter = 10000)
program_colors_used <- program_colors[
  names(program_colors) %in% unique(V(g2)$program)
]

pdf(
  "Figure1B_coordination_network.pdf",
  width = 10,
  height = 9
)
plot(
  g2,
  layout = layout1,
  vertex.label.color = "black",
  vertex.label.family = "sans",
  vertex.frame.color = NA,
  edge.curved = 0.15,
  margin = 0.2,
  main = "RAC1-responsive coordination architecture in WT tumor cells"
)

legend(
  "topright",
  legend = names(program_colors),
  col = program_colors,
  pch = 19,
  pt.cex = 1.5,
  bty = "n",
  title = "program"
)
dev.off()

# export edge & node tables 
edge_table <- as_data_frame(g2, what = "edges")
write.csv(edge_table, "Figure1B_coordination_network_edges.csv", row.names = FALSE)

node_table <- data.frame(
  gene = names(deg2),
  degree = deg2,
  program = V(g2)$program
)
write.csv(node_table, "Figure1B_coordination_network_nodes.csv", row.names = FALSE)



############################################################
#PHASE III: Perturbation-informed Latent Tumor-State Manifold
############################################################
library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
######
## STEP 18：(MANIFOLD)WT Cells Projected into RAC1 Perturbation-Responsive Space
expr_sub <- readRDS("expr_sub.rds")
rac1_dr <- read.csv("RAC1_TopDRgenes.csv")
# 选取 perturbation-responsive genes
manifold_genes <- rac1_dr$gene[rac1_dr$gene != "RAC1"][1:300]
manifold_genes <- intersect( manifold_genes,rownames(expr_sub))

length(manifold_genes)

sub_expr_manifold <- expr_sub[rownames(expr_sub) %in% manifold_genes,]
# create Seurat object
plot_seu <- CreateSeuratObject(counts = sub_expr_manifold)
plot_seu <- NormalizeData(plot_seu)
plot_seu <- ScaleData(plot_seu,features = rownames(plot_seu))
plot_seu <- RunPCA(plot_seu,features = rownames(plot_seu),npcs = 10)
plot_seu <- RunUMAP(plot_seu,dims = 1:10)# UMAP manifold
plot_seu <- FindNeighbors(plot_seu,dims = 1:10)
plot_seu <- FindClusters(plot_seu,resolution = 0.4)

saveRDS(
  plot_seu,
  "RAC1_manifold_seurat.rds"
)

#######
## STEP 19：(UCELL) UCell-based program scoring
alveolar_program <- c("SFTPA1","SFTPA2","SFTPB","SFTPD","NAPSA","CTSH")
plasticity_program <- c("S100A4","MARCKS","KRT17","SNCG")
epithelial_program <- c("PSCA","S100P","MUC4","CD99")
stress_program <- c("GADD45B","DDIT3","HEXIM1","ASS1","SLC39A8","LRRK2","CDKN2A")
oncogenic_program <- c("HMGA1","MACC1","PRSS23","CD24")
remodeling_program <- c("SPINT2","C19orf33")

alveolar_program <- intersect(alveolar_program,rownames(plot_seu))
plasticity_program <- intersect(plasticity_program,rownames(plot_seu))
epithelial_program <- intersect(epithelial_program, rownames(plot_seu))
stress_program <- intersect(stress_program,rownames(plot_seu))
oncogenic_program <- intersect(oncogenic_program,rownames(plot_seu))
remodeling_program <- intersect(remodeling_program,rownames(plot_seu))

program_list <- list(
  alveolar_program = alveolar_program,
  plasticity_program = plasticity_program,
  epithelial_program = epithelial_program,
  stress_program = stress_program,
  oncogenic_program = oncogenic_program,
  remodeling_program = remodeling_program
)

plot_seu <- AddModuleScore_UCell(plot_seu,features = program_list)

saveRDS(
  plot_seu,
  "RAC1_UCell_landscape_seurat.rds"
)

######
## STEP 20：Cluster-program annotation
program_cols <- c(
  "alveolar_program_UCell",
  "plasticity_program_UCell",
  "epithelial_program_UCell",
  "remodeling_program_UCell",
  "stress_program_UCell",
  "oncogenic_program_UCell"
)

meta_df <- plot_seu@meta.data
meta_df$cluster <- Idents(plot_seu)

cluster_program_summary <- meta_df %>%
  group_by(cluster) %>%
  summarise(
    alveolar = mean(alveolar_program_UCell),
    plasticity = mean(plasticity_program_UCell),
    epithelial = mean(epithelial_program_UCell),
    remodeling = mean(remodeling_program_UCell),
    stress = mean(stress_program_UCell),
    oncogenic = mean(oncogenic_program_UCell),
    .groups = "drop"
  )

cluster_program_summary$dominant_program <- apply(
  cluster_program_summary[, c(
    "alveolar",
    "plasticity",
    "epithelial",
    "remodeling",
    "stress",
    "oncogenic"
  )],
  1,
  function(x){
    names(x)[which.max(x)]
  }
)

plot_seu$dominant_program <- cluster_program_summary$dominant_program[
  match(
    as.character(Idents(plot_seu)),
    as.character(cluster_program_summary$cluster)
  )
]

#######
## STEP 21：构图：Figure2A_latent_manifold
program_mat <- as.data.frame(
  cluster_program_summary[, c(
    "alveolar",
    "plasticity",
    "epithelial",
    "remodeling",
    "stress",
    "oncogenic"
  )]
)
rownames(program_mat) <- as.character(cluster_program_summary$cluster)

program_z <- scale(program_mat)
program_z <- as.data.frame(program_z)
program_z$cluster <- rownames(program_z)

auto_state <- apply(
  program_z[, c(
    "alveolar",
    "plasticity",
    "epithelial",
    "remodeling",
    "stress",
    "oncogenic"
  )],
  1,
  function(x){
    top2 <- names(
      sort(
        x,
        decreasing = TRUE
      )
    )[1:2]
    paste0(
      top2[1],
      "-",
      top2[2],
      " state"
    )
  }
)

auto_state_table <- data.frame(
  cluster = rownames(program_z),
  state_annotation = auto_state,
  program_z[, c(
    "alveolar",
    "plasticity",
    "epithelial",
    "remodeling",
    "stress",
    "oncogenic"
  )],
  row.names = NULL
)
# semantic refinement
refined_labels <- c(
  "stress-alveolar state" = "alveolar-stress",
  "alveolar-stress state" = "alveolar-stress",
  
  "oncogenic-alveolar state" = "alveolar-oncogenic",
  "alveolar-oncogenic state" = "alveolar-oncogenic",
  
  "oncogenic-plasticity state" = "oncogenic-plasticity",
  "plasticity-oncogenic state" = "oncogenic-plasticity",
  
  "remodeling-oncogenic state" = "oncogenic-remodeling",
  "oncogenic-remodeling state" = "oncogenic-remodeling",
  
  "epithelial-plasticity state" = "epithelial-plasticity",
  "plasticity-epithelial state" = "epithelial-plasticity",
  
  "plasticity-remodeling state" = "plasticity-remodeling",
  "remodeling-plasticity state" = "plasticity-remodeling",
  
  "remodeling-epithelial state" = "epithelial-remodeling",
  "epithelial-remodeling state" = "epithelial-remodeling"
)

auto_state_table$refined_state_annotation <- refined_labels[
  auto_state_table$state_annotation
]

auto_state_table$refined_state_annotation[
  is.na(auto_state_table$refined_state_annotation)
] <- auto_state_table$state_annotation[
  is.na(auto_state_table$refined_state_annotation)
]


plot_seu$state_annotation <- auto_state_table$refined_state_annotation[
  match(
    as.character(Idents(plot_seu)),
    auto_state_table$cluster
  )
]

write.csv(
  auto_state_table,
  "Figure2A_Automatic_State_Annotation_Refined_Table.csv",
  row.names = FALSE
)

table(plot_seu$state_annotation)

p_1 <- DimPlot(
  plot_seu,
  reduction = "umap",
  group.by = "state_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.35
) +
  ggtitle(
    "RAC1-responsive latent tumor-state architecture"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18,
      margin = margin(b = 12)
    ),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 20, 15, 20)
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(size = 3)
    )
  )

ggsave(
  "Figure2A_latent_manifold.pdf",
  p_1,
  width = 11,
  height = 8
)

######
## STEP 22：构图：Figure2B_program_landscape
program_plot_info <- list(
  list(
    feature = "alveolar_program_UCell",
    title = "Alveolar identity",
    color = "goldenrod3"
  ),
  list(
    feature = "plasticity_program_UCell",
    title = "Plasticity",
    color = "steelblue4"
  ),
  list(
    feature = "epithelial_program_UCell",
    title = "Epithelial remodeling",
    color = "orange2"
  ),
  list(
    feature = "remodeling_program_UCell",
    title = "Remodeling",
    color = "purple3"
  ),
  list(
    feature = "stress_program_UCell",
    title = "Stress adaptation",
    color = "forestgreen"
  ),
  list(
    feature = "oncogenic_program_UCell",
    title = "Oncogenic activation",
    color = "firebrick4"
  )
)

p_list_program <- lapply(
  program_plot_info,
  function(x){
    FeaturePlot(
      plot_seu,
      features = x$feature,
      reduction = "umap",
      pt.size = 0.22,
      order = TRUE,
      cols = c("grey95", x$color)
    ) +
      ggtitle(x$title) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        ),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(5, 12, 5, 5)
      )
  }
)

p_2 <- wrap_plots(
  p_list_program,
  ncol = 2
) +
  plot_annotation(
    title = "Program-resolved functional landscape across the RAC1-responsive latent manifold",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18
      )
    )
  )

ggsave(
  "Figure2B_program_landscape.pdf",
  p_2,
  width = 11,
  height = 9
)

#####
## STEP 23：构图：Figure2C_cluster_program_heatmap
cluster_program_summary$state_annotation <-
  auto_state_table$refined_state_annotation[
    match(
      as.character(cluster_program_summary$cluster),
      auto_state_table$cluster
    )
  ]

cluster_program_summary$cluster_state <- paste0(
  "C",
  cluster_program_summary$cluster,
  " | ",
  cluster_program_summary$state_annotation
)

cluster_program_summary$cluster_state <- factor(
  cluster_program_summary$cluster_state,
  levels = rev(cluster_program_summary$cluster_state)
)

write.csv(
  cluster_program_summary,
  "Figure2C_Cluster_Program_Summary.csv",
  row.names = FALSE
)

cluster_program_long <- cluster_program_summary %>%
  pivot_longer(
    cols = c(
      alveolar,
      epithelial,
      plasticity,
      remodeling,
      stress,
      oncogenic
    ),
    names_to = "program",
    values_to = "score"
  )

cluster_program_long$program <- factor(
  cluster_program_long$program,
  levels = c(
    "alveolar",
    "epithelial",
    "plasticity",
    "remodeling",
    "stress",
    "oncogenic"
  )
)

p_3_heat <- ggplot(
  cluster_program_long,
  aes(
    x = program,
    y = cluster_state,
    fill = score
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "firebrick3",
    name = "Mean UCell"
  ) +
  labs(
    title = "Cluster-program architecture",
    subtitle = "RAC1-responsive tumor-state space",
    x = "Program",
    y = "Cluster | state"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 13,
      margin = margin(b = 12)
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 20, 15, 25)
  )

ggsave(
  "Figure2C_cluster_program_heatmap.pdf",
  p_3_heat,
  width = 11,
  height = 7
)

######
## STEP 24：构图：Figure2D_functional_gradients
plot_seu$gradient_plasticity_alveolar <-
  plot_seu$plasticity_program_UCell -
  plot_seu$alveolar_program_UCell

plot_seu$gradient_stress_oncogenic <-
  plot_seu$stress_program_UCell -
  plot_seu$oncogenic_program_UCell

plot_seu$gradient_remodeling_epithelial <-
  plot_seu$remodeling_program_UCell -
  plot_seu$epithelial_program_UCell

gradient_features <- c(
  "gradient_plasticity_alveolar",
  "gradient_stress_oncogenic",
  "gradient_remodeling_epithelial"
)

gradient_titles <- c(
  "Plasticity-alveolar gradient",
  "Stress-oncogenic gradient",
  "Remodeling-epithelial gradient"
)

p_list_gradient <- lapply(
  seq_along(gradient_features),
  function(i){
    
    df_i <- FetchData(
      plot_seu,
      vars = c(
        "umap_1",
        "umap_2",
        gradient_features[i]
      )
    )
  
    colnames(df_i) <- c(
      "UMAP_1",
      "UMAP_2",
      "gradient"
    )
    
    ggplot(
      df_i,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = gradient
      )
    ) +
      geom_point(
        size = 0.28,
        alpha = 0.95
      ) +
      scale_color_gradient2(
        low = "navy",
        mid = "grey92",
        high = "firebrick3",
        midpoint = 0,
        name = "gradient"
      ) +
      ggtitle(
        gradient_titles[i]
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(
          size = 15,
          face = "bold",
          hjust = 0.5
        ),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.margin = margin(10, 35, 10, 10)
      )
  }
)

p_4_gradient <- wrap_plots(
  p_list_gradient,
  ncol = 1
) +
  plot_annotation(
    title ="Functional program gradients\nin RAC1-responsive latent state space",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18,
        lineheight = 1.1
      )
    )
  )

ggsave(
  "Figure2D_functional_gradients.pdf",
  p_4_gradient,
  width = 7.8,
  height = 13.3
)

######
## STEP 25：export UMAP coordinates with programs

umap_df <- as.data.frame(
  Embeddings(
    plot_seu,
    reduction = "umap"
  )
)

umap_df$cell <- rownames(umap_df)
umap_df$cluster <- Idents(plot_seu)
umap_df$dominant_program <- plot_seu$dominant_program
umap_df$state_annotation <- plot_seu$state_annotation

for(col in program_cols){
  umap_df[[col]] <- plot_seu@meta.data[[col]]
}

umap_df$gradient_plasticity_alveolar <-
  plot_seu$gradient_plasticity_alveolar

umap_df$gradient_stress_oncogenic <-
  plot_seu$gradient_stress_oncogenic

umap_df$gradient_remodeling_epithelial <-
  plot_seu$gradient_remodeling_epithelial

write.csv(
  umap_df,
  "RAC1_UMAP_coordinates_with_programs.csv",
  row.names = FALSE
)

saveRDS(
  plot_seu,
  "RAC1_UCell_landscape_program_annotated.rds"
)



############################################################
#PHASE Ⅳ: clinical translation：Prognostic modeling
############################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

library(dplyr)
library(survival)
library(glmnet)
library(survminer)
library(SummarizedExperiment)
library(survcomp)

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
candidate_genes <- setdiff(candidate_genes,"RAC1")# remove RAC1 itself
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
# summary
cat("========== 最终多因素Cox模型总结 ==========\n")
print(multi_summary)
cat("================ 分析完成 ================\n")


####################################
# 验证模型
####################################

library(dplyr)
library(survival)
library(survminer)
library(timeROC)
library(survcomp)

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

# 加legend
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



############################################################
# THE END
############################################################
