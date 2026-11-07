# =============================================================================
# R/utils.R
#
# Shared helper functions used across the analysis and figure scripts.
# Source from the repository root:
#   source("R/utils.R")
# This also sources config/paths.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(qs2)
})

source("config/paths.R")

options(bitmapType = "cairo")

# -----------------------------------------------------------------------------
# Standard Seurat RNA workflow: log-normalise, HVGs, scale, PCA, UMAP
# -----------------------------------------------------------------------------
seurat_rna_pipeline <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  seu %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:10, reduction.name = "rna.umap")
}

# -----------------------------------------------------------------------------
# Per-sample QC on a SingleCellExperiment (scater + scDblFinder)
# Returns list(filtered sce, per-cell metrics data.frame)
# -----------------------------------------------------------------------------
qc_filter <- function(sce, mito_list) {
  metrics <- scater::perCellQCMetrics(sce, subsets = list(Mito = mito_list))
  metrics_df <- as.data.frame(metrics) %>%
    rownames_to_column(var = "barcode")

  lib   <- scater::isOutlier(metrics_df$sum, log = TRUE, type = "lower")
  genes <- scater::isOutlier(metrics_df$detected, log = TRUE, type = "lower")
  # mito <- metrics_df$subsets_Mito_percent > 10
  mito  <- metrics_df$subsets_Mito_percent >
    median(metrics_df$subsets_Mito_percent) + 3 * sd(metrics_df$subsets_Mito_percent)

  discard <- lib | genes | mito
  metrics_df$discard <- discard
  sce <- sce[, !discard]

  # doublet removal
  sce <- scDblFinder::scDblFinder(sce)
  doublet_df <- data.frame(barcode = colnames(sce),
                           doublet = sce$scDblFinder.class)
  sce <- sce[, which(sce$scDblFinder.class == "singlet")]

  metrics_df <- metrics_df %>% left_join(doublet_df, by = "barcode")
  list(sce, metrics_df)
}

# -----------------------------------------------------------------------------
# Cell-cycle regression (S - G2M difference), used for T cells
# -----------------------------------------------------------------------------
cellCycleReg <- function(seu) {
  seu <- JoinLayers(seu, assay = "RNA")
  seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  ScaleData(seu, vars.to.regress = "CC.Difference", features = rownames(seu))
}

# -----------------------------------------------------------------------------
# inferCNV on one sample (immune cells as reference)
# -----------------------------------------------------------------------------
runInferCNV <- function(seu, sampleAnnot, save_dir, ref = NULL,
                        gene_order_file = file.path(REF_DIR, "gencode_v21_gen_pos.complete.txt")) {
  counts_matrix <- seu[["RNA"]]$counts
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file  = sampleAnnot,
    delim             = "\t",
    gene_order_file   = gene_order_file,
    ref_group_names   = ref,
    chr_exclude       = c("chrY", "chrM"))

  infercnv::run(infercnv_obj,
                cutoff             = 0.1,   # 0.1 for 10x Genomics
                out_dir            = save_dir,
                cluster_by_groups  = FALSE,
                cluster_references = FALSE,
                denoise            = FALSE,
                HMM                = FALSE,
                output_format      = NA,
                num_threads        = 10,
                no_prelim_plot     = TRUE,
                no_plot            = TRUE)
}

# -----------------------------------------------------------------------------
# Violin/box plot of mean marker-set expression per cluster
# geneList: named list of character vectors
# -----------------------------------------------------------------------------
plotCluster_knownMarkers <- function(seu, geneList) {
  meanExpr_sub <- data.frame(barcode = character(), cluster = numeric(),
                             category = character(), meanExpr = numeric())

  for (i in seq_along(geneList)) {
    genes <- intersect(Features(seu), geneList[[i]])
    if (length(genes) <= 1) next
    meanExpr <- seu[["RNA"]]$data[genes, ] %>% colMeans(na.rm = TRUE)
    meanExpr_sub <- rbind(meanExpr_sub, data.frame(
      barcode  = colnames(seu),
      cluster  = seu$seurat_clusters,
      category = names(geneList)[i],
      meanExpr = meanExpr))
  }

  meanExpr_sub %>%
    ggplot(aes(x = cluster, y = meanExpr)) +
    geom_violin(width = 1.2, scale = "width") +
    geom_boxplot(width = 0.1) +
    facet_wrap(~category, ncol = 2) +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text   = element_text(size = 16),
          axis.text    = element_text(size = 16),
          axis.title   = element_text(size = 20),
          legend.text  = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Cluster", y = "Mean Expression of Marker Genes Per Cell")
}