# =============================================================================
# Harmony batch correction (by sample) of one T-cell lineage list 
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(harmony)
})

args <- commandArgs(trailingOnly = TRUE)
seuList_path <- args[1]
outPath <- args[2]
seuList <- readRDS(seuList_path)
message("Loaded list object")

# cell cycle regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cellCycleReg <- function(seu){
  seu <- JoinLayers(seu, assay = "RNA")
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  seu <- ScaleData(seu, vars.to.regress = "CC.Difference", features = rownames(seu))
  return(seu)
}

# prepping for merge
preprocess_list <- function(obj_list, assay = "RNA") {
  lapply(obj_list, function(x) {
    DefaultAssay(x) <- assay
    x <- JoinLayers(x, assay = assay)
    x <- DietSeurat(x, assay = assay, counts = TRUE, data = TRUE, scale.data = FALSE)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x
  })
}

# separate out samples with fewer than 50 cells
seuList_small <- list()
small_sampleIDs <- c()
smallIndices <- c()

for(i in seq(1, length(seuList))){
    seu <- seuList[[i]]
    if(ncol(seu) == 1){
        smallIndices <- c(smallIndices, i)
    } else{
        if(ncol(seu) <= 50){
        seuList_small <- c(seuList_small, list(seu[["RNA"]]$counts))
        small_sampleIDs <- c(small_sampleIDs, seu$sample)
        smallIndices <- c(smallIndices, i)
        }
    }
}

if(length(seuList_small) == 0 & length(smallIndices) > 0){
    seuList[smallIndices] <- NULL
    seuListFull <- seuList
    } else{
    if(length(seuList_small) > 0){
        seuList[smallIndices] <- NULL
        sharedGenes <- Reduce(intersect, lapply(seuList_small, rownames))
        seuList_small <- lapply(seuList_small, function(x) x[sharedGenes,])
        
        small_mat <- Reduce(cbind, seuList_small)
        
        rownames(small_mat) <- sharedGenes
        dupInd <- which(duplicated(rownames(small_mat)))
        
        if(ncol(small_mat) <= 50){
        seuListFull <- seuList
        } else{
        small <- CreateSeuratObject(counts = small_mat) %>%
            NormalizeData() %>%
            FindVariableFeatures() %>%
            ScaleData() %>%
            # SCTransform() %>%
            RunPCA()
        small$sample <- small_sampleIDs
        
        seuListFull <- c(seuList, small)
        }
    } else{
        seuListFull <- seuList
    }
}

# regress out cell cycle
seuListFull <- preprocess_list(seuListFull)
merged <- Reduce(merge, seuListFull)
merged <- cellCycleReg(merged)
merged <- RunPCA(merged)

batchCorrected <- RunHarmony(merged, group.by.vars = "sample", reduction.save = "harmony", verbose = FALSE)
batchCorrected <- RunUMAP(batchCorrected, reduction = "harmony", dims = 1:15)

batchCorrected <- FindNeighbors(batchCorrected, reduction = "harmony") %>%
    FindClusters(resolution = 0.8)
batchCorrected <- JoinLayers(batchCorrected, assay = "RNA")

saveRDS(batchCorrected, outPath)