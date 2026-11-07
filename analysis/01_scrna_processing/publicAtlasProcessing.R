# =============================================================================
# External scRNA-seq cohort (Salcher et al., 2022; CZ CELLxGENE).
# Converts the atlas h5ad into one Seurat object per donor
# broad cell-type labels from the pan-cancer reference (Guimaraes et al., 2024).
#
# Run twice: once with DISEASE_SET = "tumour" (LUAD + LUSC) and once with
# DISEASE_SET = "normal".
# =============================================================================

source("R/utils.R")
library(zellkonverter)
library(biomaRt)

DISEASE_SET <- Sys.getenv("IEI_DISEASE_SET", "tumour")   # "tumour" or "normal"
diseases <- if (DISEASE_SET == "normal") "normal" else
  c("lung adenocarcinoma", "squamous cell lung carcinoma")
suffix <- if (DISEASE_SET == "normal") "_normal_annotatedSeu.RDS" else "_annotatedSeu.RDS"
dir.create(PUBLIC_OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

salker <- readH5AD(dpath("salker_cellxgene.h5ad"))
salker_counts <- assay(salker, "count")
salker_meta <- colData(salker) %>% as.data.frame()
salker_meta_lung <- salker_meta %>%
  dplyr::filter(grepl("10x 3", assay), disease %in% diseases) %>%
  group_by(donor_id) %>%
  mutate(numCells = n()) %>%
  ungroup() %>%
  dplyr::filter(numCells >= 1000)

# Ensembl IDs -> HGNC symbols
ensembl <- readRDS(file.path(REF_DIR, "ensembl.RDS"))
annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = rownames(salker_counts),
                     mart = ensembl)
annotations_filtered <- dplyr::filter(annotations, !duplicated(hgnc_symbol))
cellIndices <- which(colData(salker)$donor_id %in% salker_meta_lung$donor_id)

salker_counts_filtered <- salker_counts[annotations_filtered$ensembl_gene_id, cellIndices]
rownames(salker_counts_filtered) <- annotations_filtered$hgnc_symbol
salker_counts_filtered <- salker_counts_filtered[rownames(salker_counts_filtered) != "", ]

salker_seu <- CreateSeuratObject(counts = salker_counts_filtered)
salker_seu <- AddMetaData(salker_seu, salker_meta)

# One annotated object per donor
BoroniReference <- readRDS(file.path(REF_DIR, "boroniCellTypeAnnotRef.RDS"))

for (di in unique(salker_seu$donor_id)) {
  seuSubset <- subset(salker_seu, donor_id == di) %>% seurat_rna_pipeline()
  seuSubset <- RenameCells(seuSubset, new.names = paste(colnames(seuSubset), gsub("-", "_", di), sep = "-"))

  anchors <- FindTransferAnchors(reference = BoroniReference, query = seuSubset,
                                 dims = 1:30, reference.reduction = "rpca")
  predictions <- TransferData(anchorset = anchors, refdata = BoroniReference$cellType, dims = 1:30)
  seuSubset <- AddMetaData(seuSubset, predictions)
  seuSubset <- SCTransform(seuSubset) %>% RunPCA()
  seuSubset$sample <- gsub("-", "_", di)

  saveRDS(seuSubset, file.path(PUBLIC_OBJ_DIR, paste0("salcher_", di, suffix)))
}