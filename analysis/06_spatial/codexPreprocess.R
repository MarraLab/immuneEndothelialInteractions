# =============================================================================
# CODEX multiplex immunofluorescence (Zheng et al., 2025; Dryad).
# =============================================================================

source("R/utils.R")

mat <- readMM(dpath("counts_nsclcCODEX.mtx"))
barcodes <- fread(dpath("barcodes_nsclcCODEX.tsv"))
features <- read_tsv(dpath("features_nsclcCODEX.tsv"), col_names = FALSE)
meta <- fread(dpath("metadata_nsclcCODEX.tsv")) %>% mutate(V1 = paste0("cell_", V1))

rownames(mat) <- paste0("cell_", barcodes$V1)
colnames(mat) <- features$X1
meta <- column_to_rownames(as.data.frame(meta), var = "V1")
mif <- CreateSeuratObject(t(mat), meta.data = meta)

mif$ctype <- case_when(
  mif$phenotype == "Vascular" ~ "Endothelial",
  mif$phenotype == "Th"       ~ "CD4",
  mif$phenotype == "Tc"       ~ "CD8",
  .default = "Other")
mif$x <- mif$`Centroid X µm`
mif$y <- mif$`Centroid Y µm`
mif$section <- "Section 1"
mif$patient <- mif$patient_id

ptMeta <- fread(dpath("deident_clinical.csv"))
luad <- ptMeta$patient_id[which(ptMeta$`Histology class (0: LUAD, 1: LUSC, 2: Other)` == 0)]
lusc <- ptMeta$patient_id[which(ptMeta$`Histology class (0: LUAD, 1: LUSC, 2: Other)` == 1)]
mif$histology <- case_when(
  mif$patient_id %in% luad ~ "LUAD",
  mif$patient_id %in% lusc ~ "LUSC",
  .default = "Other")

qs_save(mif, dpath("codex_allSamples.qs2"))