############################################################
# RAC1-latent-state-mapping
# Figure1_GRN_rewiring.R
#   Figure 1: RAC1 GRN rewiring(Fig.1A) and coordination network(Fig.1B)
#
# INPUT：
#   - RAC1_vko.rds (包含由scTenifoldKnk生成的tensorNetworks)
#   - expr_sub.rds (WT 表达矩阵)
#   - RAC1_TopDRgenes.csv
#
# OUTPUT：
#   - Figure1A_GRN_rewiring.pdf
#   - Figure1A_Rewiring_edges_Optimized.csv
#   - Figure1A_Rewiring_nodes_Optimized.csv
#   - Figure1B_coordination_network.pdf
#   - Figure1B_coordination_network_edges.csv
#   - Figure1B_coordination_network_nodes.csv
############################################################



library(igraph)
library(Matrix)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(dplyr)

# 检查必要文件
if(!file.exists("RAC1_vko.rds")) stop("请先运行 00_preprocessing_and_vKO.R 生成 RAC1_vko.rds")
if(!file.exists("expr_sub.rds")) stop("请先运行 00_preprocessing_and_vKO.R 生成 expr_sub.rds")
if(!file.exists("RAC1_TopDRgenes.csv")) stop("请先运行 00_preprocessing_and_vKO.R 生成 RAC1_TopDRgenes.csv")

# 读取数据
rac1_ko <- readRDS("RAC1_vko.rds")

wt_net <- rac1_ko$tensorNetworks$WT
ko_net <- rac1_ko$tensorNetworks$KO
wt_net <- as.matrix(wt_net)
ko_net <- as.matrix(ko_net)

stopifnot(all(rownames(wt_net) == rownames(ko_net)))
genes <- rownames(wt_net)

saveRDS(wt_net,"WT_GRN_matrix.rds")
saveRDS(ko_net,"RAC1_vKO_GRN_matrix.rds")


# ---------- Figure 1A: GRN rewiring ----------
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

# ---------- Figure 1B: Coordination network in WT ----------
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
cor_mat2[abs(cor_mat2) < cor_cutoff] <- 0
diag(cor_mat2) <- 0

g2 <- graph_from_adjacency_matrix(
  cor_mat2,
  weighted = TRUE,
  mode = "undirected",
  diag = FALSE
)
# 删除孤立节点
g2 <- delete_vertices(g2, degree(g2) == 0)

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

cat("Figure1 完成。\n")