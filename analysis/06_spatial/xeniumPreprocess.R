# =============================================================================
# Xenium Prime 5K (10x Genomics) single LUAD section.
# =============================================================================

source("R/utils.R")
library(SingleR)
library(celldex)
library(FNN)

# ---- 1. load, normalise, SingleR ----------------------------------------------------
cells <- fread(dpath("nsclc_xenium", "cells.csv.gz"))
mat <- Read10X_h5(dpath("nsclc_xenium", "cell_feature_matrix.h5"))
obj <- CreateSeuratObject(counts = mat[["Gene Expression"]])
cells <- as.data.frame(cells)
rownames(cells) <- cells$cell_id
obj <- AddMetaData(obj, metadata = cells)

obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
obj <- FindNeighbors(obj, dims = 1:20) %>% FindClusters()
obj <- RunUMAP(obj, dims = 1:20)

ref <- celldex::HumanPrimaryCellAtlasData()
out <- SingleR(test = obj[["RNA"]]$data, ref = ref, labels = ref$label.main) %>% as.data.frame()
obj <- AddMetaData(obj, out)
write_csv(rownames_to_column(out, var = "barcode"), dpath("xenium_singleRAnnot.csv"))
allCells <- obj

# ---- 2a. endothelial cells -------------------------------------------------------------
endoX <- subset(allCells, pruned.labels == "Endothelial_cells") %>% seurat_rna_pipeline()
endoX <- FindNeighbors(endoX) %>% FindClusters(resolution = 0.5)
# cluster 8 removed after inspecting PECAM1 / COL2A1 / CD3E / PTPRC
endoX2 <- subset(endoX, seurat_clusters != "8") %>% seurat_rna_pipeline()
endoX2 <- FindNeighbors(endoX2) %>% FindClusters(resolution = 0.5)
# subtypes from plotCluster_knownMarkers(endoX2, geneSets_goveia)
endoX2$subtype <- case_when(
  endoX2$seurat_clusters == "0"  ~ "CapillaryTypeII",
  endoX2$seurat_clusters == "1"  ~ "CapillaryTypeI",
  endoX2$seurat_clusters == "2"  ~ "Arterial",
  endoX2$seurat_clusters == "3"  ~ "Tip",
  endoX2$seurat_clusters == "4"  ~ "CapillaryTypeII",
  endoX2$seurat_clusters == "5"  ~ "Stalk",
  endoX2$seurat_clusters == "6"  ~ "CapillaryTypeII",
  endoX2$seurat_clusters == "7"  ~ "Stalk",
  endoX2$seurat_clusters == "8"  ~ "Venous/PCV",
  endoX2$seurat_clusters == "9"  ~ "Lymphatic",
  endoX2$seurat_clusters == "10" ~ "EndMT")
qs_save(endoX2, dpath("endo_xenium.qs2"))

# ---- 2b. T cells ----------------------------------------------------------------------------
tX <- subset(allCells, pruned.labels == "T_cells") %>% seurat_rna_pipeline()
tX <- FindNeighbors(tX) %>% FindClusters(resolution = 0.5)
# from CD3E / CD8A / CD4 / NCAM1 / FOXP3
tX$subtype <- case_when(
  tX$seurat_clusters == "0"  ~ "CD4",
  tX$seurat_clusters == "1"  ~ "CD4",
  tX$seurat_clusters == "2"  ~ "CD4",
  tX$seurat_clusters == "3"  ~ "CD8",
  tX$seurat_clusters == "4"  ~ "",
  tX$seurat_clusters == "5"  ~ "CD8",
  tX$seurat_clusters == "6"  ~ "CD4",
  tX$seurat_clusters == "7"  ~ "",
  tX$seurat_clusters == "8"  ~ "",
  tX$seurat_clusters == "9"  ~ "",
  tX$seurat_clusters == "10" ~ "CD4")
qs_save(tX, dpath("tCells_xenium.qs2"))

tX2 <- subset(tX, seurat_clusters %in% c("3", "5")) %>% seurat_rna_pipeline()
tX2 <- FindNeighbors(tX2) %>% FindClusters(resolution = 0.5)
# qs_save(tX2, dpath("cd8T_xenium.qs2")) 

# ---- 3. ctype for the proximity analysis --------------------------------------------------
xen_seu <- allCells
xen_seu$ctype <- xen_seu$pruned.labels
xen_seu$ctype <- ifelse(colnames(xen_seu) %in% colnames(endoX2), "Endothelial", xen_seu$ctype)
xen_seu <- AddMetaData(xen_seu, tX[[]][, "subtype", drop = FALSE])
xen_seu$ctype <- ifelse(is.na(xen_seu$subtype), xen_seu$ctype, xen_seu$subtype)
xen_seu$patient <- "Xenium"
xen_seu$section <- "Section 1"
# qs_save(xen_seu, dpath("annotated_xeniumObj.qs2"))

# ---- 4. candidate NECTIN2+ EC / TIGIT+ CD8+ pairs (for locating example regions) --------------
endoNectin <- WhichCells(endoX2, expression = NECTIN2 > 0)
cd8Tigit   <- WhichCells(tX2, expression = TIGIT > 0)

coords <- allCells[[]][, c("x_centroid", "y_centroid")]
rownames(coords) <- allCells$cell_id
nectin_coords <- coords[endoNectin, ]
tigit_coords  <- coords[cd8Tigit, ]

nn <- get.knnx(data = tigit_coords, query = nectin_coords)
pairs <- do.call(rbind, lapply(1:10, function(i) data.frame(
  nectin_cell = rownames(nectin_coords),
  tigit_cell  = rownames(tigit_coords)[nn$nn.index[, i]],
  distance    = nn$nn.dist[, i],
  x1 = nectin_coords[, 1], y1 = nectin_coords[, 2],
  x2 = tigit_coords[nn$nn.index[, i], 1], y2 = tigit_coords[nn$nn.index[, i], 2])))
pairs <- unique(pairs[pairs$distance <= 10, ])
write_csv(pairs, rpath("xenium_nectin2_tigit_candidate_pairs.csv"))