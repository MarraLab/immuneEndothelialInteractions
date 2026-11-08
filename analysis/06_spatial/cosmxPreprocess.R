# =============================================================================
# CosMx SMI NSCLC FFPE dataset (Bruker; 8 sections from 5 LUAD patients).
# =============================================================================

source("R/utils.R")
library(sf)
library(terra)
sf::sf_use_s2(FALSE)

spatialSamples <- list.files(COSMX_RAW_DIR)
sample_dir <- function(f) dpath("SpatialCohort", f)
flat_dir   <- function(f) file.path(COSMX_RAW_DIR, f, paste0(f, "-Flat_files_and_images"))

# ---- cell-boundary polygons and Nanostring objects ------------------------------
get_fov_id_from_path <- function(fp, return_padded = FALSE, pad_width = 4) {
  base <- basename(fp)

  # 1) FOV[_-]?00...<id>
  m <- regexpr("(?i)FOV[_-]?0*([A-Za-z0-9]+)", base, perl = TRUE)
  if (m[1] != -1) {
    id <- sub("(?i)^.*FOV[_-]?0*", "", regmatches(base, m), perl = TRUE)

  } else {
    # 2) F00...<id>  (e.g., F001)
    m <- regexpr("(?i)F[_-]?0*([A-Za-z0-9]+)", base, perl = TRUE)
    if (m[1] != -1) {
      id <- sub("(?i)^.*F[_-]?0*", "", regmatches(base, m), perl = TRUE)

    } else {
      # 3) Last number before extension
      id <- sub("^.*?(\\d+)(?:\\.[^.]+)?$", "\\1", base, perl = TRUE)
    }
  }

  # Normalize: prefer digits; fall back to alphanum string if needed
  id_num <- suppressWarnings(as.integer(gsub("\\D", "", id)))
  if (!is.na(id_num)) {
    if (return_padded) sprintf(paste0("%0", pad_width, "d"), id_num) else as.character(id_num)
  } else {
    id
  }
}

# Read a Cell Labels TIF into an integer matrix [rows = height, cols = width]
read_label_matrix <- function(tif_path) {
  arr <- tiff::readTIFF(tif_path, as.is = TRUE)
  if (length(dim(arr)) == 3L) arr <- arr[,,1]
  matrix(as.integer(round(arr)), nrow = nrow(arr), ncol = ncol(arr))
}

labels_to_polygons <- function(label_mat, fov_id) {

  # terra expects raster with origin at top-left 
  r <- rast(label_mat)

  # Define extent so that:
  # x: 0 → width
  # y: 0 → height (TOP-LEFT origin preserved)
  ext(r) <- c(0, ncol(label_mat), 0, nrow(label_mat))

  p <- as.polygons(r, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  sf_polys <- st_as_sf(p)
  
  colnames(sf_polys)[colnames(sf_polys) == "lyr.1"] <- "cell_ID"
  sf_polys$fov <- as.character(fov_id)
  sf_polys <- sf_polys[!is.na(sf_polys$cell_ID) & sf_polys$cell_ID != 0, ]
  sf_polys
}

extract_vertices <- function(geometry, fov, cell_ID, x0, y0, fov_h, fov_w) {
  crd <- sf::st_coordinates(geometry)
  dt  <- as.data.table(crd)

  # Keep outer ring only (L3 == 1)
  if ("L3" %in% names(dt)) dt <- dt[L3 == 1]

  dt[, .(
    fov = fov,
    cell_ID = as.integer(cell_ID),
    x_global_px = as.numeric(X) + x0,
    y_global_px = as.numeric(Y) + y0,
    vertex_index = seq_len(.N) - 1L
  )]
}

for (f in spatialSamples) {
  samplePath <- paste0(flat_dir(f), "/")
  tx <- fread(paste0(samplePath, f, "_tx_file.csv"))

  if ("FOV" %in% names(tx)) setnames(tx, "FOV", "fov")
  tx[, fov := as.character(fov)]

  # sanity: required columns
  stopifnot(all(c("fov","x_global_px","y_global_px","x_local_px","y_local_px") %in% names(tx)))

  # Derive UL origins per FOV from molecules
  fov_from_tx <- tx[, .(
    x0 = median(x_global_px - x_local_px, na.rm = TRUE),
    y0 = median(y_global_px - y_local_px, na.rm = TRUE)
  ), by = fov]

  labels_dir <- paste0(samplePath, "CellLabels")
  lab_files  <- list.files(labels_dir, pattern = "\\.tif(f)?$", full.names = TRUE, recursive = TRUE)
  stopifnot(length(lab_files) > 0)

  polys_local <- list()

  for (lab_path in lab_files) {

    fov_id <- as.character(get_fov_id_from_path(lab_path))
    lab <- read_label_matrix(lab_path)
    fov_h <- nrow(lab)
    fov_w <- ncol(lab)

    sf_polys <- labels_to_polygons(lab, fov_id)
    polys_local[[fov_id]] <- sf_polys
  }

  polys_local_sf <- rbindlist(polys_local, fill = TRUE)
  polys_local_sf <- merge(polys_local_sf, fov_from_tx, by = "fov")

  coords_list <- lapply(seq_len(nrow(polys_local_sf)), function(i) {

    extract_vertices(
      geometry = polys_local_sf$geometry[i],
      fov = polys_local_sf$fov[i],
      cell_ID = polys_local_sf$cell_ID[i],
      x0 = polys_local_sf$x0[i],
      y0 = polys_local_sf$y0[i],
      fov_h = fov_h,
      fov_w = fov_w
    )
  })

  seg_out <- rbindlist(coords_list)
  setorder(seg_out, fov, cell_ID, vertex_index)

  seg_out <- seg_out[, {
    r <- .SD

    # check if first and last vertex are identical
    first_x <- r$x_global_px[1]
    first_y <- r$y_global_px[1]

    last_x  <- r$x_global_px[.N]
    last_y  <- r$y_global_px[.N]

    if (first_x != last_x || first_y != last_y) {
      r <- rbind(
        r,
        data.table(
          fov = r$fov[1],
          cell_ID = r$cell_ID[1],
          x_global_px = first_x,
          y_global_px = first_y,
          vertex_index = max(r$vertex_index) + 1L
        )
      )
    }

    r
  }, by = .(fov, cell_ID)]

  seg_out_final <- seg_out[, .(
    fov = as.integer(fov),
    cellID = as.integer(cell_ID),   # IMPORTANT: numeric only
    x_global_px = x_global_px,
    y_global_px = y_global_px,
    vertex_index = as.integer(vertex_index)
  )]

  write_csv(seg_out_final, paste0(samplePath, f, "_polygons2.csv"))

  nano.obj <- LoadNanostring(data.dir = samplePath, fov = f)
  saveRDS(nano.obj, file.path(sample_dir(f), "nanostrObj.RDS"))
}