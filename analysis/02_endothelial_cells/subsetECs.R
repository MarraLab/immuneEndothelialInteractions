# =============================================================================
# Endothelial-cell identification, contaminant removal and batch correction.
# =============================================================================

source("R/utils.R")

allCells <- qs_read(dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))
ecMeta <- fread(dpath("patientMetadata_publicData.csv"))
luadSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "LUAD")]
luscSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "LUSC")]
normSamples <- ecMeta$sample[which(ecMeta$cancerSubtype == "normal")]

# ---- 1. ECs by either method ------------------------------------------------
ecInd <- which(allCells$pruned.labels == "Endothelial_cells" | allCells$predicted.id == "Endothelial")
endoAll <- subset(allCells, cells = colnames(allCells)[ecInd]) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "sample", reduction.save = "harmony", verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:15)
qs_save(endoAll, dpath("endoAllCells_unfiltered.qs2"))

# ---- 2. cluster and find contaminants ----------------------------------------
endoAll <- FindNeighbors(endoAll, reduction = "harmony", assay = "rna", dims = 1:30) %>%
  FindClusters(resolution = 0.5)
endoAll_markers <- FindAllMarkers(endoAll, logfc.threshold = 1, min.pct = 0.25)
write_tsv(endoAll_markers, dpath("allEndo_filtering_markers.tsv"))
qs_save(endoAll, dpath("endoAllCells_unfiltered.qs2"))

# Inspected with:
# FeaturePlot(endoAll, reduction = "umap",
#             features = c("PTPRC", "CD3D", "LYZ", "EPCAM", "PECAM1", "COL1A1"))
contaminant_clusters <- c(7, 11, 12, 13, 15, 16, 17)
endoFiltered <- subset(endoAll, cells = colnames(endoAll)[
  which(!(endoAll$seurat_clusters %in% contaminant_clusters))]) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

endoFiltered$cancerSubtype <- case_when(
  endoFiltered$sample %in% luadSamples ~ "LUAD",
  endoFiltered$sample %in% luscSamples ~ "LUSC",
  endoFiltered$sample %in% normSamples ~ "Normal")
endoFiltered$platform <- case_when(grepl("JN", endoFiltered$sample) ~ "FFPE", .default = "Other")

# ---- 3. batch correction by platform and re-clustering -----------------------
endoFiltered <- endoFiltered %>%
  RunHarmony(group.by.vars = c("platform"), reduction.save = "harmony", verbose = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:15)
endoFiltered <- FindNeighbors(endoFiltered, reduction = "harmony", assay = "rna",
                              dims = 1:15, k.param = 20)
endoFiltered <- FindClusters(endoFiltered, resolution = 0.2)

qs_save(endoFiltered, dpath("endoAllCells_filtered.qs2"))