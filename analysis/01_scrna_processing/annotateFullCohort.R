# =============================================================================
# Merge every internal and public sample (LUAD, LUSC, normal lung) into one
# Seurat object, then normalise, PCA, Harmony (by sample) and UMAP.
# =============================================================================

source("R/utils.R")

publicObj  <- list.files(PUBLIC_OBJ_DIR)
ecMeta     <- fread(dpath("patientMetadata_publicData.csv"))
lungSamples <- ecMeta$sample

seuList <- list()
for (f in lungSamples) {
  internal_path <- file.path(INTERNAL_OBJ_DIR, f, "rna_ProcessedSeurat.RDS")
  if (file.exists(internal_path)) {
    seu <- readRDS(internal_path)
  } else {
    e <- publicObj[grepl(f, publicObj)]
    seu <- readRDS(file.path(PUBLIC_OBJ_DIR, e))
  }

  DefaultAssay(seu) <- "RNA"
  seu <- DietSeurat(seu, assay = "RNA")
  seuList[[f]] <- CreateSeuratObject(counts = seu[["RNA"]]$counts, meta.data = seu[[]])
}

seuAll <- Reduce(merge, seuList)
seuAll <- JoinLayers(seuAll, assay = "RNA")
DefaultAssay(seuAll) <- "RNA"
seuAll <- seuAll %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "sample", reduction.save = "harmony", verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:15)

qs_save(seuAll, dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))

# =============================================================================
# Broad cell-type labels for the merged cohort.
#   1. Sketch 75,000 cells (leverage score), Harmony, cluster, and project the
#      clustering back onto all cells.
#   2. Assign celltype_broad by cluster.
#   3. SingleR (Human Primary Cell Atlas) as the orthogonal EC annotation.
#   4. Add the sample-level metadata used by the statistics and CellChat
#      scripts (cancerSubtype, platform, histology, tissue, dataset).
# =============================================================================

library(SingleR)
library(celldex)

allCells <- qs_read(dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))

# ---- 1. sketch-based clustering ---------------------------------------------
allCells <- SketchData(object = allCells, ncells = 75000,
                       method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(allCells) <- "sketch"
allCells <- FindVariableFeatures(allCells) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "sample", reduction.save = "harmony", verbose = FALSE)
allCells <- RunUMAP(allCells, reduction = "harmony", dims = 1:15, return.model = TRUE)
allCells <- FindNeighbors(allCells, reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

allCells <- ProjectIntegration(object = allCells, sketched.assay = "sketch",
                               assay = "RNA", reduction = "harmony")
allCells <- ProjectData(object = allCells, assay = "RNA",
                        full.reduction = "harmony.full", sketched.assay = "sketch",
                        sketched.reduction = "harmony", umap.model = "umap",
                        dims = 1:50, refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(allCells) <- "RNA"

# ---- 2. broad cell types by cluster -----------------------------------------
allCells$celltype_broad <- case_when(
  allCells$cluster_full == 0  ~ "Malignant",
  allCells$cluster_full == 1  ~ "Myeloid",
  allCells$cluster_full == 2  ~ "Fibroblast",
  allCells$cluster_full == 3  ~ "Epithelial",
  allCells$cluster_full == 4  ~ "T",
  allCells$cluster_full == 5  ~ "Malignant",
  allCells$cluster_full == 6  ~ "Epithelial",
  allCells$cluster_full == 7  ~ "Endothelial",
  allCells$cluster_full == 8  ~ "T",
  allCells$cluster_full == 9  ~ "B",
  allCells$cluster_full == 10 ~ "Malignant",
  allCells$cluster_full == 11 ~ "Malignant",
  allCells$cluster_full == 12 ~ "Pericyte",
  allCells$cluster_full == 13 ~ "Epithelial",
  allCells$cluster_full == 14 ~ "Malignant",
  allCells$cluster_full == 15 ~ "T",
  allCells$cluster_full == 16 ~ "Endothelial",
  allCells$cluster_full == 17 ~ "T",
  allCells$cluster_full == 18 ~ "T",
  allCells$cluster_full == 19 ~ "T",
  allCells$cluster_full == 20 ~ "T"
)

# ---- 3. SingleR (orthogonal EC call) ----------------------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
out <- SingleR(test = allCells[["RNA"]]$data, ref = ref, labels = ref$label.main) %>%
  as.data.frame()
allCells <- AddMetaData(allCells, out)
write_csv(rownames_to_column(out, var = "barcode"), dpath("allCells_singleRAnnot.csv"))

# ---- 4. sample-level metadata -----------------------------------------------
ecMeta <- fread(dpath("patientMetadata_publicData.csv"))
luadSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "LUAD")]
luscSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "LUSC")]
normSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "normal")]

allCells$cancerSubtype <- case_when(
  allCells$sample %in% luadSamples ~ "LUAD",
  allCells$sample %in% luscSamples ~ "LUSC",
  allCells$sample %in% normSamples ~ "Normal")

# "FFPE" = internal GEM-X Flex (JN_*) samples
allCells$platform <- case_when(grepl("JN", allCells$sample) ~ "FFPE", .default = "Other")

# RECONSTRUCTED (not in the uploaded scripts): columns read by
# 04_icl_expression/01_icl_expression_models.R and 05_cellchat/01_run_cellchat.R.
# Confirm these match how the object used for the paper was annotated.
allCells$histology <- ifelse(allCells$cancerSubtype %in% c("LUAD", "LUSC"), allCells$cancerSubtype, NA)
allCells$tissue    <- ifelse(allCells$cancerSubtype == "Normal", "normal", "tumour")
if (!"dataset" %in% colnames(allCells@meta.data)) allCells$dataset <- NA_character_
allCells$dataset[is.na(allCells$dataset)] <- "Internal_GSC"

qs_save(allCells, dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))