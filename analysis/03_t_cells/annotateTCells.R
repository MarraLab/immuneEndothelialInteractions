# =============================================================================
# T CELL ANNOTATIONS
# =============================================================================

source("R/utils.R")

tCellSpecific_list <- list()
for (f in lungSamples) {
  internal_path <- file.path(INTERNAL_OBJ_DIR, f, "rna_ProcessedSeurat.RDS")
  seu <- if (file.exists(internal_path)) readRDS(internal_path) else readRDS(file.path(PUBLIC_OBJ_DIR, f))

  if (sum(seu$predicted.id == "NK/T_cells") <= 1) next
  seu <- subset(seu, predicted.id == "NK/T_cells")

  dat <- seu[["RNA"]]$data
  cd8_indices    <- which(dat["CD8A", ] > 0 & dat["CD4", ] == 0)
  cd4_indices    <- which(dat["CD8A", ] == 0 & dat["CD4", ] > 0)
  dp_indices     <- which(dat["CD8A", ] > 0 & dat["CD4", ] > 0)
  dn_cd3_indices <- which(colSums(dat[grepl("CD3[A-Z]", rownames(dat)), ]) > 0)
  nk_indices     <- which(dat["NCAM1", ] > 0)

  tCell_idents <- rep("notT", ncol(seu))
  tCell_idents[dn_cd3_indices] <- "DoubleNeg"
  tCell_idents[cd4_indices]    <- "CD4"
  tCell_idents[cd8_indices]    <- "CD8"
  tCell_idents[dp_indices]     <- "DoublePos"
  tCell_idents[nk_indices]     <- "NK"
  seu$tIdent <- tCell_idents

  for (cellType in unique(seu$tIdent)) {
    filteredSeu <- subset(seu, tIdent == cellType)
    tCellSpecific_list[[cellType]] <- c(tCellSpecific_list[[cellType]], list(filteredSeu))
  }
}
saveRDS(tCellSpecific_list, dpath("LungCohort", "allSamples_tCellsSubsetList_publicData_wNormal.RDS"))

# lineage lists, each with the double negatives appended
doubleNeg <- lapply(tCellSpecific_list[["DoubleNeg"]], function(x)
  AddMetaData(x, data.frame(barcode = colnames(x), subtype = "doubleNeg", row.names = colnames(x))))

for (cellType in c("CD4", "CD8", "NK", "DoublePos")) {
  cellTypeList <- lapply(tCellSpecific_list[[cellType]], function(x)
    AddMetaData(x, data.frame(barcode = colnames(x), subtype = cellType, row.names = colnames(x))))
  saveRDS(c(cellTypeList, doubleNeg),
          dpath("LungCohort", paste0("allLung_", cellType, "_list_wNormal_wDoubleNegs.RDS")))
}

library(slingshot)
library(SingleCellExperiment)

lc <- function(...) dpath("LungCohort", ...)

# ---- A. refine lineages -------------------------------------------------------
doublePos <- readRDS(lc("allLung_DoublePos_wNormal_wDoubleNegs_batchCorrected.RDS"))
nkt       <- readRDS(lc("allLung_NK_wNormal_wDoubleNegs_batchCorrected.RDS"))
cd4       <- readRDS(lc("allLung_CD4_wNormal_wDoubleNegs_batchCorrected.RDS"))
cd8       <- readRDS(lc("allLung_CD8_wNormal_wDoubleNegs_batchCorrected.RDS"))

dnBarcodes_all <- colnames(doublePos)[which(doublePos$subtype == "doubleNeg")]
write_csv(data.frame(barcode = dnBarcodes_all), dpath("dnBarcodes_all.csv"))

# NOTE: the uploaded version of this function first subset each object to the
# normal-lung samples only (a re-run dated "Sept 3", saved to
# tCellSubtyping_inProgress_normals.qs2). That subset is removed here because
# the downstream steps read tCellSubtyping_inProgress.qs2 (all samples).
reprocT <- function(seu) {
  seu %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = c("sample"), reduction.save = "harmony", verbose = FALSE, max_iter = 100) %>%
    RunUMAP(reduction = "harmony", dims = 1:15) %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 0.2)
}

tCell_list <- lapply(list(doublePos, nkt, cd4, cd8), reprocT)   # order: DP, NK/T, CD4, CD8

# keep double negatives only in clusters where they co-cluster with the lineage
# (second iteration of manual cluster review)
clusterKeep_list <- list(
  doublePos = c(100),
  nkt       = c(100),
  cd4       = c(0, 1, 2, 3, 4, 6, 7, 8),
  cd8       = c(0, 1, 3, 4, 5, 6, 7, 8))

tCell_list2 <- lapply(seq_along(tCell_list), function(i) {
  seu <- tCell_list[[i]]
  toKeep <- which(seu$subtype != "doubleNeg" | seu$seurat_clusters %in% clusterKeep_list[[i]])
  subset(seu, cells = colnames(seu)[toKeep])
})
tCell_list2 <- lapply(tCell_list2, reprocT)

# NOTE: the original saves tCell_list (before the double-negative cleanup),
# not tCell_list2. Kept as in the original; confirm which was intended.
qs_save(tCell_list, dpath("tCellSubtyping_inProgress.qs2"))

# ---- B. CD4+ subtypes -----------------------------------------------------------
cd4 <- tCell_list[[3]]
cd4 <- FindClusters(cd4, resolution = 0.25)
# markers inspected: IFNG, GZMB, SELL, CCR7, CXCR6, CXCL13, FOXP3, IL2RA, MKI67
cd4$identity <- case_when(
  cd4$seurat_clusters == 0  ~ "SELL+CCR7+",
  cd4$seurat_clusters == 1  ~ "IFNG+",
  cd4$seurat_clusters == 2  ~ "Tregs",
  cd4$seurat_clusters == 3  ~ "Tph",
  cd4$seurat_clusters == 4  ~ "Trm",
  cd4$seurat_clusters == 5  ~ "Cycling",
  cd4$seurat_clusters == 6  ~ "SELL+CCR7+",
  cd4$seurat_clusters == 7  ~ "Cytotoxic",
  cd4$seurat_clusters == 8  ~ "SELL+CCR7+",
  cd4$seurat_clusters == 9  ~ "SELL+CCR7+",
  cd4$seurat_clusters == 10 ~ "SELL+CCR7+")
cd4Idents <- data.frame(barcode = colnames(cd4), identity = cd4$identity) %>% drop_na()

# ---- C. CD8+ subset and trajectory ---------------------------------------------
cd8_allCells <- tCell_list[[4]]
cd8_allCells <- FindClusters(cd8_allCells, resolution = 0.3)

# clusters 5, 7, 10 are NKT / cycling cells and are excluded from the trajectory
cd8_allCells$identity <- case_when(
  cd8_allCells$seurat_clusters == 5  ~ "NKT",
  cd8_allCells$seurat_clusters == 7  ~ "Cycling",
  cd8_allCells$seurat_clusters == 10 ~ "NKT",
  .default = paste("CD8", cd8_allCells$seurat_clusters, sep = "_"))

indices <- which(!(cd8_allCells$seurat_clusters %in% c(5, 7, 10)))
cd8_subset <- subset(cd8_allCells, cells = colnames(cd8_allCells)[indices])
cd8_subset <- FindNeighbors(cd8_subset, reduction = "harmony", k = 10) %>%
  FindClusters(resolution = 0.8)

# merged clusters used as Slingshot nodes; "1" = naive (CCR7/SELL/LEF1) start
cd8_subset$clustRev <- case_when(
  cd8_subset$seurat_clusters %in% c(13)        ~ "1",
  cd8_subset$seurat_clusters %in% c(7, 9, 14)  ~ "2",
  cd8_subset$seurat_clusters %in% c(4, 5)      ~ "3",
  cd8_subset$seurat_clusters %in% c(1, 12)     ~ "4",
  cd8_subset$seurat_clusters %in% c(8, 11, 15) ~ "5",
  cd8_subset$seurat_clusters %in% c(0)         ~ "6",
  cd8_subset$seurat_clusters %in% c(3, 10)     ~ "7",
  cd8_subset$seurat_clusters %in% c(6)         ~ "8",
  cd8_subset$seurat_clusters %in% c(2)         ~ "9")
cd8_subset <- RunUMAP(cd8_subset, reduction = "harmony", dims = 1:15)
qs_save(cd8_subset, dpath("cd8TCells.qs2"))

cd8_joined <- JoinLayers(cd8_subset, overwrite = TRUE)
sce <- as.SingleCellExperiment(cd8_joined)
reducedDims(sce)$UMAP <- Embeddings(cd8_subset, reduction = "umap")
reducedDims(sce)$PCA  <- Embeddings(cd8_subset, reduction = "harmony")
sce$cluster <- cd8_joined$clustRev

set.seed(123)
sce <- slingshot(sce,
                 clusterLabels = "cluster",
                 reducedDim    = "PCA",
                 start.clus    = "1",
                 end.clus      = c("4", "9"),   # activated or exhausted
                 extend        = "n")
qs_save(sce, lc("allCD8_slingshotRes.qs2"))

embedded <- embedCurves(sce, "UMAP")
qs_save(embedded, lc("allCD8_slingshotRes_umapEmbed.qs2"))

# pseudotime bins (width 2.5) per lineage
pt <- slingPseudotime(sce)
breaks <- seq(0, 37.5, by = 2.5)
lineage1_ID <- cut(as.numeric(pt[, 1]), breaks = breaks, labels = paste("Lin1", breaks, sep = "_")[-1])
breaks <- seq(0, 50, by = 2.5)
lineage2_ID <- cut(as.numeric(pt[, 3]), breaks = breaks, labels = paste("Lin2", breaks, sep = "_")[-1])

ptIdents <- data.frame(barcode = names(pt[, 1]), lin1Ident = lineage1_ID, lin2Ident = lineage2_ID) %>%
  mutate(identity = paste(lin1Ident, lin2Ident, sep = "_")) %>%
  dplyr::select(barcode, identity)
ptIdents$identity <- gsub("NA_", "", ptIdents$identity)
ptIdents$identity <- gsub("_NA", "", ptIdents$identity)

ptIdents_filtered <- dplyr::filter(ptIdents, identity != "NA") %>%
  group_by(identity) %>%
  mutate(numCells = n()) %>%
  dplyr::select(-numCells) %>%
  ungroup()

cd8Idents <- rbind(
  data.frame(barcode = colnames(cd8_allCells), identity = cd8_allCells$identity) %>%
    dplyr::filter(!(barcode %in% ptIdents_filtered$barcode)),
  ptIdents_filtered)
write_csv(cd8Idents, dpath("cd8Idents.csv"))

# ---- D. compile identities ---------------------------------------------------------
overlap <- intersect(cd8Idents$barcode, cd4Idents$barcode)   
cd4Idents <- dplyr::filter(cd4Idents, !(barcode %in% overlap))

nk <- tCell_list[[2]]
cd3Neg_indices <- which(colSums(nk[["RNA"]]$data[grepl("CD3[A-Z]", rownames(nk[["RNA"]]$data)), ]) == 0)
nkIdents <- rep("NKT", ncol(nk))
nkIdents[cd3Neg_indices] <- "NK"           # NCAM1+ CD3-
nkIdents_df <- data.frame(barcode = colnames(nk), identity = nkIdents)

dpSubset <- subset(tCell_list[[1]], subtype == "DoublePos")
dpIdents <- data.frame(barcode = colnames(dpSubset), identity = "DoublePos")

allTCells <- data.frame(
  barcode = c(unlist(lapply(tCell_list, colnames)), dnBarcodes_all) %>% unique()) %>%
  left_join(Reduce(rbind, list(cd4Idents, cd8Idents, nkIdents_df, dpIdents)), by = "barcode")
allTCells$identity[is.na(allTCells$identity)] <- "DoubleNeg"
write_csv(allTCells, lc("tCellIdentities.csv"))