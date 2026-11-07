# Endothelial cells as sources of immune checkpoint ligands and candidate regulators of PD-1 signalling in non-small cell lung cancer 

## Repository layout

```
config/paths.R            input/output locations (set with environment variables)
R/utils.R                 shared helpers (Seurat workflow, QC, inferCNV, marker plots)
analysis/
  01_scrna_processing/    QC and initial processing of public and internal sc/snRNA-seq datasets
  02_endothelial_cells/   EC identification, filtering, batch correction, and subtyping
  03_t_cells/             T-cell lineages, CD4 subtypes, CD8 trajectories
  04_icl_expression/      Statistical comparisons of ICL expression
  05_cellchat/            Scripts for running CellChat and generating figures
  06_spatial/             CosMx, Xenium and CODEX processing and proximity statistics
figures.R    code for every computational figure panel
```
## Data

| Dataset | Source |
|---|---|
| Internal sc/snRNA-seq (6 fresh-frozen, 28 FFPE LUSC) | See manuscript for data accession numbers |
| NSCLC and normal lung scRNA-seq atlas | Salcher et al., 2022, CZ CELLxGENE |
| Cell-type annotation reference | Guimarães et al., 2024, CZ CELLxGENE |
| Human Primary Cell Atlas (SingleR) | `celldex::HumanPrimaryCellAtlasData()` |
| CosMx SMI NSCLC FFPE dataset | Bruker Spatial Biology |
| CODEX mIF, 35 NSCLC patients | Zheng et al., 2025, Dryad |
| Xenium Prime 5K human LUAD | 10x Genomics |

Data files are not stored in this repository.

## Software

R 4.2.2 with Seurat 5.2.1, harmony 1.2.4, SingleR 2.0.0, celldex, scDblFinder 1.12.0, scater, inferCNV 1.3.3, slingshot 2.6.0, CellChat 2.1.2, glmmTMB 1.1.11, emmeans, edgeR, limma, speckle, terra 1.7-39, sf 1.0-21, sfheaders, RANN, FNN, zellkonverter, biomaRt, qs2, tidyverse, data.table, ggpubr, ggsci, ggforce, ggsignif, paletteer.