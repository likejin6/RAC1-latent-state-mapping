############################################################
# RAC1-latent-state-mapping
# Figure2_latent_manifold.R
#   Figure 2: perturbation-informed latent manifold(Fig.2A), program landscape( Fig.2B), 
#             cluster heatmap(Fig.2C),and gradients(Fig.2D)
#
# INPUT：
#   - expr_sub.rds
#   - RAC1_TopDRgenes.csv
#
# OUTPUT：
#   - RAC1_manifold_seurat.rds
#   - RAC1_UCell_landscape_seurat.rds
#   - RAC1_UCell_landscape_program_annotated.rds
#   - Figure2A_latent_manifold.pdf
#   - Figure2B_program_landscape.pdf
#   - Figure2C_cluster_program_heatmap.pdf
#   - Figure2D_functional_gradients.pdf
#   - 相关CSV表格
############################################################



library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 检查必要文件
if(!file.exists("expr_sub.rds")) stop("请先运行 00_preprocessing_and_vKO.R 生成 expr_sub.rds")
if(!file.exists("RAC1_TopDRgenes.csv")) stop("请先运行 00_preprocessing_and_vKO.R 生成 RAC1_TopDRgenes.csv")

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
## STEP 19：(UCell 评分) UCell-based program scoring
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

cat("Figure2 全部完成。\n")