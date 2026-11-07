# =============================================================================
# Copy-number inference to identify malignant cells (inferCNV v1.3.3)
# =============================================================================

source("R/utils.R")

publicObj <- list.files(PUBLIC_OBJ_DIR)

for (f in publicObj) {
  seu <- readRDS(file.path(PUBLIC_OBJ_DIR, f))
  sample <- seu$sample[[1]]
  out_dir <- file.path(INFERCNV_DIR, paste0(sample, "_onlyMaligObs"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sample_annot <- data.frame(sample_id = colnames(seu)) %>%
    mutate(label = ifelse(seu$predicted.id %in% c("Malignant"), "Observation", "Reference"))
  annot_file <- file.path(out_dir, "rna_SampleAnnot_onlyMaligObs.tsv")
  write_tsv(sample_annot, annot_file, col_names = FALSE)

  tryCatch(runInferCNV(seu, sampleAnnot = annot_file, save_dir = out_dir, ref = "Reference"),
           error = function(e) message("inferCNV error: ", f))
}
