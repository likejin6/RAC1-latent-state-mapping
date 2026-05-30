# RAC1-latent-state-mapping

**Perturbation-informed latent tumor-state mapping reveals RAC1-responsive regulatory network rewiring in lung adenocarcinoma**

This repository contains the analysis code, processed results, and figures for a computational biology project investigating how RAC1-responsive regulatory network rewiring informs latent tumor-state organization in lung adenocarcinoma (LUAD).

Rather than focusing only on RAC1 expression or conventional prognostic modeling, this project uses RAC1 as a **virtual perturbation anchor**. RAC1 virtual knockout-derived regulatory rewiring genes are converted into **perturbation-informed features** to reconstruct tumor-state architecture from single-cell transcriptomic data.

---

## Project overview

Tumor cells in LUAD exhibit heterogeneous and plastic transcriptional states. However, conventional differential expression analysis may not fully capture perturbation-induced changes in regulatory network structure.

This project integrates:

- single-cell RNA-seq preprocessing
- virtual knockout analysis using `scTenifoldKnk`
- gene regulatory network (GRN) rewiring analysis
- RAC1-responsive feature selection
- perturbation-informed latent manifold reconstruction
- UCell-based functional program scoring
- TCGA-LUAD prognostic modeling

The central question is:

> Can RAC1 virtual perturbation-derived regulatory rewiring genes be used to reconstruct and interpret latent tumor-state organization in LUAD?

---
## Quick links

- [One-page project brief](docs/1_page_Project_Brief_Likejin_RAC1_LUAD.pdf)
- [Full English project report](docs/Li kejin_RAC1_LUAD_Perturbation_Mapping_Project Report.pdf)
- [Main figures](figures/)
- [Analysis scripts](R/)
  
---

## Workflow

1. LUAD single-cell RNA-seq preprocessing
2. Tumor epithelial cell extraction
3. Candidate perturbation gene prioritization
4. Multi-gene virtual knockout screening
5. RAC1 selection as the main perturbation target
6. RAC1-vKO GRN rewiring analysis
7. RAC1-responsive gene identification
8. Perturbation-informed latent manifold construction
9. Functional program scoring and hybrid state annotation
10. TCGA-LUAD prognostic relevance analysis

---

## Key findings

- RAC1 showed the strongest perturbation response among successfully evaluated candidate genes.
- RAC1 virtual knockout induced regulatory network rewiring and hub-status reorganization.
- 99 RAC1-responsive genes reconstructed a perturbation-informed latent tumor-state manifold from WT tumor cells.
- The latent manifold contained hybrid functional states, including alveolar-stress, oncogenic-plasticity, oncogenic-remodeling, epithelial-plasticity, plasticity-remodeling, and epithelial-remodeling states.
- Functional program gradients suggested continuous organization across the RAC1-responsive latent state space.
- An eight-gene RAC1-responsive prognostic signature retained prognostic relevance in the TCGA-LUAD cohort.

---

## Main figures

### Figure 1. RAC1-responsive regulatory network rewiring

- **Figure 1A:** RAC1-centered GRN rewiring after virtual perturbation
- **Figure 1B:** RAC1-responsive coordination architecture in WT tumor cells

Files are available in:

- `figures/fig1/`
- `results/results_figure1/`

### Figure 2. Perturbation-informed latent tumor-state architecture

- **Figure 2A:** RAC1-responsive latent tumor-state manifold
- **Figure 2B:** Program-resolved functional landscape
- **Figure 2C:** Cluster-program architecture heatmap
- **Figure 2D:** Functional program gradients

Files are available in:

- `figures/fig2/`
- `results/results_figure2/`

### Figure 3. TCGA-LUAD prognostic relevance

- **Figure 3A:** LASSO-Cox feature selection
- **Figure 3B:** Kaplan-Meier survival analysis
- **Figure 3C:** Time-dependent ROC analysis

Files are available in:

- `figures/fig3/`
- `results/results_figure3/`

---

## Repository structure

- `R/`  
  R scripts for preprocessing, virtual knockout, GRN rewiring, latent manifold construction, and prognostic modeling.

- `data/`  
  Core processed input files and data description.

- `figures/`  
  Manuscript-quality figures in PDF and PNG formats.

- `results/`  
  Processed result tables supporting each figure.

- `supplementary/`  
  Supplementary tables supporting the research project report.

- `docs/`  
  Research project report and project brief.

- `environment/`  
  R session information used for the analysis.

---

## How to reproduce the analysis

Run the R scripts in the following order:

```r
source("R/00_preprocessing_and_vKO.R")
source("R/Figure1_GRN_rewiring.R")
source("R/Figure2_latent_manifold.R")
source("R/Figure3_prognostic_model.R")
```

For a complete workflow, run:

```r
source("R/complete_analysis_pipeline.R")
```

---

## Data availability

The single-cell RNA-seq dataset used in this project is:

- **GSE131907**: LUAD single-cell RNA-seq dataset

The raw UMI count matrix and annotation file are not included in this repository because of file-size limitations. Please download the raw data from GEO and place them in the appropriate local directory before running the preprocessing script.

TCGA-LUAD transcriptomic and clinical data were downloaded using the `TCGAbiolinks` R package.

Processed intermediate files required for downstream analyses are described in:

- `data/README.md`
- `results/README.md`
- `supplementary/README.md`

---

## Supplementary materials

Supplementary tables are provided in:

- `supplementary/`

These include:

- virtual knockout screening summary
- RAC1 top differentially rewired genes
- cluster-level UCell program scores and hybrid state annotations
- prognostic model summary tables
- marker-based functional program gene sets

---

## Environment

The R session information used for this analysis is provided in:

- `sessionInfo.txt`

---

## Project positioning

This project is intended as a **computational biology research project report** for demonstrating a perturbation-informed single-cell analysis framework.

The main conceptual contribution is:

> converting virtual perturbation-derived regulatory network response information into graph-informed features for latent tumor-state mapping.

---

## License

This repository is released under the MIT License.
