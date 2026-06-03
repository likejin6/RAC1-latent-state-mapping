# Python Graph Representation Validation of RAC1-induced GRN Rewiring

This folder contains a supplementary Python-based graph representation analysis for the RAC1-LUAD perturbation-informed tumor-state mapping project.

The goal of this analysis is to evaluate whether RAC1-responsive differentially rewired genes (DR genes), identified from the main scTenifoldKnk-based workflow, also show larger network-role or graph-representation shifts between WT and RAC1-vKO gene regulatory networks (GRNs).

## Important note

The WT and RAC1-vKO GRN adjacency matrices used in this Python analysis were exported from the main R-based scTenifoldKnk workflow. Therefore, this analysis is **not** an independent data-source validation and is **not** an independent GRN inference workflow.

Instead, it serves as a secondary graph-based computational validation and representation-learning extension based on scTenifoldKnk-derived WT/vKO networks. The independent aspect lies in the graph analysis perspective and quantitative representation methods, including topology metrics, adjacency-vector PCA, and aligned node2vec embeddings.

## Input files

The input files are located in the `python_inputs/` folder:

```text
python_inputs/
├── WT_GRN_matrix.csv
├── RAC1_vKO_GRN_matrix.csv
├── RAC1_TopDRgenes.csv
└── GRN_gene_list.csv
```

These files were exported from the main R workflow:

* `WT_GRN_matrix.csv`: WT GRN adjacency matrix
* `RAC1_vKO_GRN_matrix.csv`: RAC1-vKO GRN adjacency matrix
* `RAC1_TopDRgenes.csv`: RAC1-responsive DR gene ranking from scTenifoldKnk
* `GRN_gene_list.csv`: gene list used in the GRN matrices

## Analysis overview

The notebook performs three complementary graph-based analyses:

### 1. Network topology shift validation

WT and RAC1-vKO adjacency matrices were converted into weighted NetworkX graphs. The top 5% strongest absolute-weight edges were retained to construct comparable sparse graph representations.

For each gene, the following topology metrics were calculated in WT and RAC1-vKO graphs:

* degree centrality
* weighted strength
* betweenness centrality
* weighted clustering coefficient

A composite topology shift score was then calculated for each gene.

### 2. Adjacency-vector PCA graph representation shift

Each gene was represented by its full GRN adjacency vector, meaning that each gene was described by its weighted connection profile to all other genes.

WT and RAC1-vKO gene vectors were jointly embedded using PCA. For each gene, the Euclidean distance between its WT and RAC1-vKO positions in PCA space was calculated as the graph representation shift distance.

### 3. Aligned node2vec embedding shift validation

Node2Vec was used to learn graph embeddings from WT and RAC1-vKO GRNs.

Because node2vec embeddings trained separately on two graphs can differ by arbitrary rotation or reflection, RAC1-vKO embeddings were aligned to WT embeddings using orthogonal Procrustes alignment before calculating node-level embedding shift distances.

This provides a random-walk-based graph representation sensitivity check.

## Main results

Across all three validation layers, RAC1-responsive DR genes showed larger WT-to-RAC1-vKO network-role or graph-representation shifts than non-DR genes.

| Validation layer                     | DR median shift | non-DR median shift | Mann–Whitney U p value | Top20 DR fraction |
| ------------------------------------ | --------------: | ------------------: | ---------------------: | ----------------: |
| Network topology shift               |        0.007286 |            0.000000 |               1.10e-10 |              0.70 |
| Adjacency-vector PCA embedding shift |        0.856612 |            0.148545 |               1.30e-37 |              0.95 |
| Aligned node2vec embedding shift     |        0.513932 |            0.232779 |               8.06e-10 |              0.55 |

These results support that RAC1-responsive DR genes are not only highly ranked by scTenifoldKnk perturbation distance, but also occupy more shifted positions in WT versus RAC1-vKO graph representation spaces.

## Key figures

The main figures are located in the `figures/` folder.

Recommended figures for presentation:

```text
figures/
├── Figure_topology_shift_boxplot.pdf
├── Figure_top_shifted_genes_barplot.pdf
├── Figure_PCA_embedding_shift_boxplot.pdf
├── Figure_top_embedding_shifted_genes_barplot.pdf
├── Figure_node2vec_embedding_shift_boxplot.pdf
└── Figure_top_node2vec_shifted_genes_barplot.pdf
```

The PCA gene representation scatter plot is also provided as a diagnostic visualization:

```text
Figure_PCA_gene_embedding_space.pdf
```

## Result tables

The result tables are located in the `results/` folder.

The most important summary table is:

```text
results/RAC1_python_graph_validation_final_summary.csv
```

Additional result tables include node-level shift scores, statistical tests, and top-20 DR gene enrichment analyses for topology shift, PCA embedding shift, and node2vec shift.

## Notebook

The full analysis is provided in:

```text
RAC1_graph_representation_validation.ipynb
```

The notebook includes:

1. input loading and inspection
2. matrix consistency checks
3. weighted GRN graph construction
4. topology shift validation
5. adjacency-vector PCA embedding shift validation
6. aligned node2vec embedding shift validation
7. final validation summary and output export

## Python packages

Main Python packages used:

```text
pandas
numpy
networkx
matplotlib
scikit-learn
scipy
node2vec
```

## How to reproduce

Open the notebook in Google Colab.

If running in Google Colab, upload the input files or `python_inputs.zip` when prompted. Then run all cells sequentially.

The notebook will generate:

```text
figures/
results/
RAC1_python_graph_validation_results_figures.zip
```

## Interpretation

This supplementary Python analysis should be interpreted as a graph-based re-analysis of scTenifoldKnk-derived WT and RAC1-vKO GRNs. It does not replace the main R-based perturbation workflow. Instead, it provides additional computational evidence that RAC1-responsive DR genes show stronger network-role and graph-representation shifts after RAC1 virtual knockout.
