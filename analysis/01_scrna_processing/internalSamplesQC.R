# =============================================================================
# Internal scRNA/snRNA-seq cohort 
# =============================================================================

source("R/utils.R")

# QC from Cell Ranger output


# cell type annotation
sample_meta <- fread(dpath("patientMetadata_publicData.csv"))
internal_samples <- sample_meta$sample[file.exists(
  file.path(INTERNAL_OBJ_DIR, sample_meta$sample, "rna_ProcessedSeurat.RDS"))]

BoroniReference <- readRDS(file.path(REF_DIR, "boroniCellTypeAnnotRef.RDS"))

for (f in internal_samples) {
  obj_path <- file.path(INTERNAL_OBJ_DIR, f, "rna_ProcessedSeurat.RDS")
  seu <- readRDS(obj_path)

  anchors <- FindTransferAnchors(reference = BoroniReference, query = seu,
                                 dims = 1:30, reference.reduction = "rpca")
  predictions <- TransferData(anchorset = anchors, refdata = BoroniReference$cellType, dims = 1:30)
  seu <- AddMetaData(seu, predictions)

  saveRDS(seu, obj_path)
}