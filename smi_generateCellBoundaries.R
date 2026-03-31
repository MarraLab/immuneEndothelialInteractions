source("../../scFunctions.R")

library(tiff)
library(sf)
sf::sf_use_s2(FALSE)

library(sfheaders)
library(terra)

# Root of the extracted sample bundle
root_dir <- "PATH"
spatialSamples <- list.files(root_dir)
savePath <- "sPATH"

# Extract FOV ID from file path
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

# Read a TIF with cell labels into an integer matrix [rows = height, cols = width]
read_label_matrix <- function(tif_path) {
  arr <- tiff::readTIFF(tif_path, as.is = TRUE)
  if (length(dim(arr)) == 3L) arr <- arr[,,1]
  matrix(as.integer(round(arr)), nrow = nrow(arr), ncol = ncol(arr))
}

labels_to_polygons <- function(label_mat, fov_id) {
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
  
  return(sf_polys)
}

# Fast point→label lookup using direct matrix indexing.
# Assumes transcripts x=column, y=row in pixel space with 0-based or 1-based—handled robustly.
lookup_labels_matrix <- function(label_mat, x, y) {
  h <- nrow(label_mat); w <- ncol(label_mat)
  # Convert to 1-based integer indices. If coords already 1-based, floor(x) works too.
  cx <- pmin(pmax(floor(x) + 1L, 1L), w)
  cy <- pmin(pmax(floor(y) + 1L, 1L), h)
  # R matrices are column-major: element at [row, col]
  lab <- label_mat[cbind(cy, cx)]
  as.integer(lab)
}

# Pixel areas per cell id from the label matrix
cell_areas_from_matrix <- function(label_mat) {
  v <- as.integer(label_mat)
  v <- v[v > 0L]
  if (length(v) == 0) return(data.table(cell_id = integer(), area_px = numeric()))
  tab <- tabulate(v, nbins = max(v))
  data.table(cell_id = which(tab > 0L), area_px = as.numeric(tab[tab > 0L]))
}

align_rows <- function(m, genes) {
    miss <- setdiff(genes, colnames(m))
    if (length(miss)) {
        add <- Matrix(0, ncol = ncol(m), nrow = length(miss),
                      dimnames = list(miss, colnames(m)))
        m <- rbind(m, add)
    }
    m[genes, , drop = FALSE]
}

# This creates a light JSON with vertices as [x,y] lists per cell
to_vertices_list <- function(geom) {
  # Works for POLYGON or MULTIPOLYGON; flattens to a single outer ring per part
  cls <- unique(sf::st_geometry_type(geom))
  coords_list <- list()
  i <- 1L
  for (g in geom) {
    m <- sf::st_coordinates(g)
    # st_coordinates returns [x,y,L1,L2,...]; split by polygon id (L1)
    parts <- split(as.data.frame(m)[,1:2], f = m[,"L1"])
    # Use the outer ring (first) of each part
    verts <- lapply(parts, function(df) unname(split(df, seq_len(nrow(df)))))
    # Flatten one level
    coords_list[[i]] <- verts
    i <- i + 1L
  }
  coords_list
}

# Extract a per-cell data.table of vertices with part & per-part vertex_index
# - geometry: an sfc_POLYGON or sfc_MULTIPOLYGON (single row from your sf)
# - fov, cell_ID: identifiers you want carried forward
# - keep_holes: set FALSE to keep only outer ring(s)
# - add_global: set TRUE if you pass x0, y0 and want x_global_px/y_global_px added
extract_vertices_with_parts <- function(geometry, fov, cell_ID,
                                        keep_holes = FALSE,
                                        add_global = FALSE,
                                        x0 = NA_real_, y0 = NA_real_) {
  # st_coordinates returns: X, Y and indexing columns (L1, L2, L3, ...)
  crd <- sf::st_coordinates(geometry)
  dt  <- as.data.table(crd)

  # If present, L3 is ring index: 1 = outer, >1 = holes. Optionally drop holes.
  if (!keep_holes && "L3" %in% names(dt)) {
    dt <- dt[L3 == 1]
  }

  # Define 'part':
  #  - If L2 exists, use it (it separates sub-polygons inside a MULTIPOLYGON).
  #  - If L2 doesn't exist (simple POLYGON), set part = 1 for all rows.
  dt[, part := if ("L2" %in% names(dt)) as.integer(L2) else 1L]

  # Build output vertex table in local pixels
  out <- dt[, .(fov         = as.character(fov),
                cell_ID     = as.integer(cell_ID),
                part        = part,
                x_local_px  = as.numeric(X),
                y_local_px  = as.numeric(Y))]

  # Per-(fov, cell_ID, part) 0-based vertex_index
  setorder(out, fov, cell_ID, part)
  out[, vertex_index := seq_len(.N) - 1L, by = .(fov, cell_ID, part)]

  # Optionally add global pixel coordinates from x0, y0 (per-FOV origin)
  if (isTRUE(add_global)) {
    stopifnot(is.finite(x0), is.finite(y0))
    out[, `:=`(
      x_global_px = x_local_px + x0,
      y_global_px = y_local_px + y0
    )]
  }

  out[]
}

# Optional: close each ring by repeating its first vertex at the end.
# This is sometimes helpful for downstream tools that expect explicitly closed paths.
close_rings <- function(vertices_dt) {
  setorder(vertices_dt, fov, cellID, vertex_index)
  vertices_dt[, {
    r <- .SD
    # Compare first and last vertices
    if (!all(r[1, .(x_global_px, y_global_px)] == r[.N, .(x_global_px, y_global_px)])) {
      # Append first row to close the ring
      r <- rbind(r, r[1])
      r[, vertex_index := seq_len(.N) - 1L]
    }
    r
  }, by = .(fov, cellID)]
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

for(f in spatialSamples){
  samplePath <- paste0(root_dir, f, "/", f, "-Flat_files_and_images/")
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
    cellID = as.integer(cell_ID),
    x_global_px = x_global_px,
    y_global_px = y_global_px,
    vertex_index = as.integer(vertex_index)
  )]

  write_csv(seg_out_final, paste0(samplePath, f, "_polygons2.csv"))

  segs <- fread(paste0(samplePath, f, "_polygons2.csv"))
  nano.obj <- LoadNanostring(data.dir = samplePath, fov = f)

  saveRDS(nano.obj, paste0(savePath, f, "/nanostrObj.RDS"))
}
