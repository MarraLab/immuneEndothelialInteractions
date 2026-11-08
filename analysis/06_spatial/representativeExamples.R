# =============================================================================
# 06_representative_examples.R
#
# Finds candidate ligand+ EC / receptor+ T-cell neighbours used ONLY to pick
# the representative regions drawn in Figures 3B-E and 4D-G. The statistics
# are in spatialProximityStats.R.
# =============================================================================

source("R/utils.R")
library(RANN)
library(sf)
library(sfheaders)
sf::sf_use_s2(FALSE)

spatialSamples <- list.files(COSMX_RAW_DIR)
sample_dir <- function(f) dpath("SpatialCohort", f)
flat_dir   <- function(f) file.path(COSMX_RAW_DIR, f, paste0(f, "-Flat_files_and_images"))

# =============================================================================
# A. CODEX
# =============================================================================
mat <- readMM(dpath("counts_nsclcCODEX.mtx"))
barcodes <- fread(dpath("barcodes_nsclcCODEX.tsv"))
features <- read_tsv(dpath("features_nsclcCODEX.tsv"), col_names = FALSE)
meta <- fread(dpath("metadata_nsclcCODEX.tsv")) %>% mutate(V1 = paste0("cell_", V1))
rownames(mat) <- paste0("cell_", barcodes$V1)
colnames(mat) <- features$X1


# ---- Helper: find unique cross-set pairs within a radius ----
pairs_within_radius <- function(
  A_xy, A_ids,
  B_xy, B_ids,
  radius_um
) {
  if (nrow(A_xy) == 0L || nrow(B_xy) == 0L) {
    return(tibble(A_id = character(), B_id = character(), dist = numeric()))
  }

  k_cap <- min(50L, nrow(A_xy))

  nn <- RANN::nn2(
    data  = A_xy,
    query = B_xy,
    k     = k_cap
  )

  out <- purrr::map_dfr(seq_len(nrow(B_xy)), function(i) {

    idxs  <- nn$nn.idx[i, ]
    valid <- idxs[idxs > 0L]
    if (!length(valid)) return(NULL)

    dists <- sqrt(
      (A_xy[valid, 1] - B_xy[i, 1])^2 +
      (A_xy[valid, 2] - B_xy[i, 2])^2
    )

    keep <- dists <= radius_um
    if (!any(keep)) return(NULL)

    tibble(
      A_id = A_ids[valid[keep]],
      B_id = B_ids[i],
      dist = dists[keep]
    )
  })
  
  if (nrow(out) == 0L) {
    return(tibble(A_id = character(), B_id = character(), dist = numeric()))
  }

  dplyr::distinct(out, A_id, B_id, .keep_all = TRUE)
}

# ---- Main function ----
compute_PD1_PDL1_adjacencies <- function(
  df,
  q_pos = 0.95,
  radius_um = 10,
  verbose = TRUE
) {

  req_cols <- c(
    "CellID", "Centroid X µm", "Centroid Y µm",
    "phenotype", "patient_id", "core_id",
    "PD-L1", "PD-1"
  )
  miss <- setdiff(req_cols, names(df))
  if (length(miss) > 0){
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }
    
  df <- df |>
    dplyr::mutate(
      `Centroid X µm` = as.numeric(`Centroid X µm`),
      `Centroid Y µm` = as.numeric(`Centroid Y µm`),
      CellID = as.character(CellID)
    )

  res <- df |>
    dplyr::group_by(patient_id, core_id) |>
    dplyr::group_modify(~{

      g <- .x

      # Positivity thresholds (per patient/core)
      thr_PDL1 <- stats::quantile(g[["PD-L1"]], q_pos, na.rm = TRUE)
      thr_PD1  <- stats::quantile(g[["PD-1"]],  q_pos, na.rm = TRUE)

      PDL1_pos <- g[["PD-L1"]] >= thr_PDL1 & !is.na(g[["PD-L1"]])
      PD1_pos  <- g[["PD-1"]]  >= thr_PD1  & !is.na(g[["PD-1"]])

      # T vs non‑T
      is_T     <- g[["phenotype"]] %in% c("Tc", "Th", "Treg", "NK Cell")
      is_nonT  <- !is_T

      T_mask    <- is_T    & PD1_pos
      nonT_mask <- is_nonT & PDL1_pos

      if (!any(T_mask) || !any(nonT_mask))
        return(tibble())

      # Helper
      get_xy_ids <- function(mask) {
        sub <- g[mask, ]
        list(
          xy   = as.matrix(sub[, c("Centroid X µm", "Centroid Y µm")]),
          ids  = sub[["CellID"]],
          type = sub[["phenotype"]]
        )
      }

      T    <- get_xy_ids(T_mask)
      nonT <- get_xy_ids(nonT_mask)

      # Spatial pairs: non‑T ↔ T
      pairs <- pairs_within_radius(
        A_xy = nonT$xy, A_ids = nonT$ids,
        B_xy = T$xy,    B_ids = T$ids,
        radius_um = radius_um
      )

      if (nrow(pairs) == 0L)
        return(tibble())

      # Map phenotypes back onto pairs
      tibble(
        patient_id   = unique(g$patient_id),
        core_id      = unique(g$core_id),
        T_cell_id    = pairs$B_id,
        adj_cell_id  = pairs$A_id
      ) |>
        dplyr::left_join(
          dplyr::select(g, CellID, T_cell_type = phenotype),
          by = c("T_cell_id" = "CellID")
        ) |>
        dplyr::left_join(
          dplyr::select(g, CellID, adj_cell_type = phenotype),
          by = c("adj_cell_id" = "CellID")
        )
    }) |>
    dplyr::ungroup()

  if (verbose) {
    message(
      "Done. Returned ",
      nrow(res),
      " PD‑1⁺ T ↔ PD‑L1⁺ non‑T adjacencies ",
      "(radius = ", radius_um, " µm)."
    )
  }

  res
}



markers <- c("PD-L1", "PD-1", "CD8", "CD4", "HLA-A", "HLA-DR")

res_list <- list()
for (p in unique(meta$patient_id)) {
  subMeta <- dplyr::filter(meta, patient_id == p)

  for (c in unique(subMeta$core_id)) {

    cells <- dplyr::filter(subMeta, core_id == c) |> pull(V1)

    inputDat <- dplyr::filter(subMeta, core_id == c) |>
      dplyr::select(
        CellID, `Centroid X µm`, `Centroid Y µm`,
        phenotype, patient_id, core_id
      ) |>
      cbind(
        mat[cells, c("PD-L1","PD-1")] |>
          as.data.frame()
      )

    idx <- paste(p, c, sep = "-")

    res_list[[idx]] <- compute_PD1_PDL1_adjacencies(
      inputDat,
      q_pos = 0.90,
      radius_um = 10,
      verbose = FALSE
    )
  }
}

res <- dplyr::bind_rows(res_list)
write_csv(res, dpath("lung_codexRes_10um.csv"))

# =============================================================================
# B. CosMx polygon neighbours
# =============================================================================
polygon_adjacency <- function(f, d_thresh = 5) {
  samplePath <- flat_dir(f)
  seg <- fread(file.path(samplePath, paste0(f, "_polygons2.csv")))
  seu <- readRDS(file.path(sample_dir(f), "rna_ProcessedSeurat.RDS"))

  # Precompute ligand/receptor expression for all cells
  ligands <- intersect(Features(seu), genes_df$ligand)
  receptors <- intersect(Features(seu), genes_df$receptor)
  lig_rec_expr <- as.data.frame(t(seu[["SCT"]]$data[unique(c(ligands,receptors)), ]))
  lig_rec_expr$cell <- rownames(lig_rec_expr)
  setDT(lig_rec_expr)

  # Build polygons once
  setorder(seg, fov, cellID, vertex_index)
  seg[, poly_order := .GRP, by = .(fov, cellID)]
  poly_groups <- unique(seg[, .(poly_order, fov, cellID)])

  sfc <- sfheaders::sfc_polygon(seg, x = "x_global_px", y = "y_global_px", polygon_id = "poly_order", close = TRUE)
  polys_sf <- st_sf(poly_groups, geometry = sfc)

  # Attach cell type
  meta <- seu[[]]
  polys_sf$newID <- paste0(polys_sf$cellID, "_", polys_sf$fov)
  polys_sf$class <- meta[polys_sf$newID, "predicted.id"]

  # Subset T cells
  oc_sf <- polys_sf[polys_sf$class != "NK/T_cells", ]
  tc_sf <- polys_sf[polys_sf$class == "NK/T_cells", ]
  oc_by <- split(oc_sf, oc_sf$fov, drop = TRUE)
  tc_by <- split(tc_sf, tc_sf$fov, drop = TRUE)
  common_fovs <- intersect(names(oc_by), names(tc_by))

  pair_list <- vector("list", length(common_fovs))
  names(pair_list) <- common_fovs

  for (i in seq_along(common_fovs)) {
    fv <- common_fovs[i]
    oc_fv <- oc_by[[fv]]
    tc_fv <- tc_by[[fv]]

    nbrs <- sf::st_is_within_distance(tc_fv, oc_fv, dist = d_thresh)

    if (all(lengths(nbrs) == 0L)) {
      pair_list[[fv]] <- NULL
      next
    }

    # Build (t_cell, o_cell) pairs from neighbor lists
    tc_idx_rep <- rep.int(seq_along(nbrs), lengths(nbrs))
    oc_idx_vec <- unlist(nbrs, use.names = FALSE)
    dist_px <- as.numeric(sf::st_distance(tc_fv[tc_idx_rep, ], oc_fv[oc_idx_vec, ], by_element = TRUE))

    pair_list[[fv]] <- data.table(
      fov     = fv,
      o_cell  = oc_fv$newID[oc_idx_vec],
      t_cell  = tc_fv$newID[tc_idx_rep],
      dist_px = dist_px
    )
  }

  adj_pairs <- rbindlist(pair_list, use.names = TRUE, fill = TRUE)
  if (nrow(adj_pairs) == 0L) {
    return(invisible(NULL))
  }

  adj_pairs <- merge(adj_pairs, lig_rec_expr, by.x = "o_cell", by.y = "cell", all.x = TRUE)
  adj_pairs <- merge(adj_pairs, lig_rec_expr, by.x = "t_cell", by.y = "cell", all.x = TRUE, suffixes = c(".oc", ".tc"))

  adj_pairs[, sample := f]

  fwrite(adj_pairs, file.path(samplePath, paste0(f, "_adjacency.csv")))
}

for (f in spatialSamples) polygon_adjacency(f)

seuList <- readRDS(dpath("SpatialCohort", "allSamples_list.RDS"))
tSpatialCells <- readRDS(dpath("SpatialCohort", "tCells_batchCorrected.RDS"))

adjList <- list()
for (f in spatialSamples) {
  adj_file <- file.path(flat_dir(f), paste0(f, "_adjacency.csv"))
  if (!file.exists(adj_file)) next
  adjList[[f]] <- fread(adj_file) %>% mutate(sample = f)
}
adjacentOtherTs <- Reduce(rbind, adjList)

# assigning identity
cellIdents_list <- lapply(seuList, function(x) data.frame(barcode = colnames(x), cellType = x$predicted.id))
allCellIdents <- Reduce(rbind, cellIdents_list) %>%
  separate(barcode, "\\.", into = c("barcode", "sample"))

tCellIdents <- data.frame(
  barcode = colnames(tSpatialCells),
  subtype = tSpatialCells$subtype
) %>%
  separate(barcode, "\\.", into = c("barcode", "sample"))

# which pairs
# 1) Build one logical vector per pair (receptor.tc > 0 & ligand.oc > 0)
pair_labels <- paste(genes_df$ligand, genes_df$receptor, sep = "-")

conds <- map2(genes_df$ligand, genes_df$receptor, function(ligand, receptor) {
  rec_tc <- paste0(receptor, ".tc")
  lig_oc <- paste0(ligand, ".oc")
  if (!(rec_tc %in% names(adjacentOtherTs)) || !(lig_oc %in% names(adjacentOtherTs))) {
    return(rep(FALSE, nrow(adjacentOtherTs)))
  }
  rec_vals <- suppressWarnings(as.numeric(adjacentOtherTs[[rec_tc]]))
  lig_vals <- suppressWarnings(as.numeric(adjacentOtherTs[[lig_oc]]))
  dplyr::coalesce(rec_vals, 0) > 0 & dplyr::coalesce(lig_vals, 0) > 0
})

# 2) Combine into a logical matrix (rows = observations, cols = pairs)
matches_mat <- if (length(conds)) do.call(cbind, conds) else matrix(FALSE, nrow(adjacentOtherTs), 0)
colnames(matches_mat) <- pair_labels

long_matches <- as_tibble(matches_mat) %>%
  mutate(.row_id = row_number()) %>%              # track source row
  pivot_longer(-.row_id, names_to = "pair", values_to = "matched") %>%
  dplyr::filter(matched) %>%
  dplyr::select(-matched)

# Join back to the original data to keep metadata (t_cell, o_cell, fov, sample, etc.)
adjacentOtherTs_long <- long_matches %>%
  left_join(
    adjacentOtherTs %>% mutate(.row_id = row_number()),
    by = ".row_id"
  ) %>%
  dplyr::select(-.row_id)

adjacentOtherTs_filtered_cellTyped <- inner_join(adjacentOtherTs_long, allCellIdents, by = c("o_cell" = "barcode", "sample"), suffix = c("", ".oc")) %>%
  inner_join(tCellIdents, by = c("t_cell" = "barcode", "sample"), suffix = c("", ".tc"))

print(table(adjacentOtherTs_filtered_cellTyped$cellType, adjacentOtherTs_filtered_cellTyped$pair))
write_csv(adjacentOtherTs_filtered_cellTyped, dpath("SpatialCohort", "allTCell_ligInteractions.csv"))
write_csv(allCellIdents, dpath("SpatialCohort", "allCellIdents.csv"))