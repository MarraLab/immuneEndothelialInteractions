# =============================================================================
# Spatial proximity / enrichment of ligand+ ECs and receptor+ T cells for
# CosMx, Xenium and CODEX (Supplementary Tables 5-7).
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(RANN); library(glmmTMB); library(qs2)
})

ENDO <- "Endothelial"
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------------------------------------------------------------------------
## 1. FEATURE ALIASES AND SPEC RESOLUTION
## ---------------------------------------------------------------------------
as_alias_list <- function(x) {
  if (is.null(x)) list() else if (is.list(x)) x else setNames(as.list(x), x)
}
norm_nm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))    # "PD-L1" == "PDL1" == "pd_l1"

find_feature <- function(aliases, available) {
  hit <- match(norm_nm(aliases), norm_nm(available))
  hit <- hit[!is.na(hit)]
  if (length(hit)) available[hit[1]] else NA_character_
}
resolve_group <- function(group, available) {
  group <- as_alias_list(group)
  if (!length(group)) return(character(0))
  res <- vapply(group, find_feature, character(1), available = available)
  res[!is.na(res)]
}

resolve_spec <- function(spec, available) {
  axes <- list()
  for (nm in names(spec$axes)) {
    ax  <- spec$axes[[nm]]
    lig <- resolve_group(ax$lig, available)
    rec <- resolve_group(ax$rec, available)
    if (!length(lig) || !length(rec)) {
      message("[", spec$name, "] axis '", nm, "' SKIPPED: ",
              if (!length(lig)) "no ligand feature measured" else "receptor not measured")
      next
    }
    miss <- setdiff(names(as_alias_list(ax$lig)), names(lig))
    if (length(miss))
      message("[", spec$name, "] axis '", nm, "': ligand(s) not measured, ignored: ",
              paste(miss, collapse = ", "))
    # One axis per T-cell subset when spec$t_sets is given (otherwise the axis's own t_types)
    sets <- if (!is.null(spec$t_sets)) spec$t_sets else list(all = ax$t_types)
    for (ts in names(sets)) {
      key <- if (is.null(spec$t_sets)) nm else paste(nm, ts, sep = "_")
      axes[[key]] <- list(lig = lig, rec = rec, t_types = sets[[ts]], base_axis = nm, t_set = ts)
    }
  }
  s_all <- resolve_group(spec$ec_strict$all, available)
  s_any <- resolve_group(spec$ec_strict$any, available)
  strict_ok <- !is.null(spec$ec_strict) &&
    length(s_all) == length(as_alias_list(spec$ec_strict$all)) && length(s_any) >= 1
  if (!strict_ok) message("[", spec$name, "] strict endothelial definition unavailable -> 'strict' scenarios dropped")
  mac    <- resolve_group(spec$mac, available)
  mac_ok <- identical(spec$type, "rna") && length(mac) >= 1
  if (!mac_ok) message("[", spec$name, "] macrophage-contamination filter unavailable -> 'purity' scenarios dropped")

  spec$resolved <- list(axes = axes, strict_all = s_all, strict_any = s_any,
                        strict_ok = strict_ok, mac = mac, mac_ok = mac_ok)
  spec
}

report_spec <- function(spec) {
  rs <- spec$resolved
  ax <- imap_dfr(rs$axes, function(a, nm) {
    tibble(modality = spec$name, axis = nm,
           ligands  = paste(names(a$lig), a$lig, sep = "=", collapse = ", "),
           receptor = paste(a$rec, collapse = ", "),
           t_types  = paste(a$t_types, collapse = "/"))
  })
  print(ax)
  message("strict endothelial definition: ", rs$strict_ok, " | macrophage-contamination filter: ", rs$mac_ok)
  invisible(ax)
}

## ---------------------------------------------------------------------------
## 2. MODALITY SPECS
## ---------------------------------------------------------------------------
# type          "rna" (counts) or "protein" (intensities)
# x_col, y_col  coordinate columns in metadata; um_per_unit converts to microns
# ctype_map     optional named vector
# t_sets        T-cell subsets to analyse separately
# pos           positivity rule for ligand/receptor:
#                 method "count"    -> transcript count >= value
#                 method "quantile" -> value >= the `value` quantile of that feature
#                                      within the section (use for protein)
# strict_pos    same idea, for the strict endothelial definition
# primary       pre-specified scenario/distance 
# dists         distance grid in microns. Centroid distances start at 10 um 
# tiles         tile sizes (um) for the permutation strata; results carry `tile_um`

spec_cosmx <- list(
  name = "CosMx", type = "rna", x_col = "x_px", y_col = "y_px", um_per_unit = 0.12,
  ctype_map = NULL,
  t_sets = list(CD4 = "CD4", CD8 = "CD8", T_all = c("CD4", "CD8")),
  axes = list(                              # t_types are ignored when t_sets is given
    PD1_PDL1      = list(lig = list(PDL1 = "CD274"),
                         rec = list(PD1 = "PDCD1")),
    PD1_PDL2      = list(lig = list(PDL2 = "PDCD1LG2"),
                         rec = list(PD1 = "PDCD1")),
    TIGIT_NECTIN2 = list(lig = list(NECTIN2 = "NECTIN2"),
                         rec = list(TIGIT = "TIGIT")),
    TIGIT_PVR     = list(lig = list(PVR = "PVR"),
                         rec = list(TIGIT = "TIGIT"))),
  ec_strict = list(all = list(PECAM1 = "PECAM1"), any = list(CLDN5 = "CLDN5", VWF = "VWF")),
  mac       = list(CD68 = "CD68", CD163 = "CD163"),
  pos       = list(method = "count", values = 1:3),
  strict_pos = list(method = "count", value = 1),
  primary   = list(pos_value = 1, ec_def = "base", purity = FALSE, dist_um = 10, tile_um = 500,
                   t_sets = c("CD4", "CD8")),    
  dists     = c(10, 15, 20, 25, 30),     
  tiles     = c(250, 500, 1000))        

spec_xenium <- modifyList(spec_cosmx, list(
  name = "Xenium", x_col = "x_centroid", y_col = "y_centroid", um_per_unit = 1))   # already um

spec_codex <- list(
  name = "CODEX", type = "protein", x_col = "x", y_col = "y", um_per_unit = 1,
  ctype_map = NULL,
  t_sets = list(CD4 = "CD4", CD8 = "CD8", T_all = c("CD4", "CD8")),
  axes = list(
    PD1_PDL1 = list(lig = list(PDL1 = c("PD-L1", "PDL1", "CD274")),              
                    rec = list(PD1  = c("PD-1", "PD1", "PDCD1")))),
  ec_strict = list(all = list(CD31 = c("CD31", "PECAM1")), any = list(CD34 = "CD34", VWF = "VWF")),
  mac       = NULL,                                                            
  pos       = list(method = "quantile", values = c(0.75, 0.90, 0.95)),
  strict_pos = list(method = "quantile", value = 0.90),
  primary   = list(pos_value = 0.90, ec_def = "base", purity = FALSE, dist_um = 10, tile_um = 500,
                   t_sets = c("CD4", "CD8")),
  dists     = c(10, 15, 20, 25, 30),    
  tiles     = c(250, 500, 1000))     

## ---------------------------------------------------------------------------
## 3. BUILD THE CELL TABLE FROM ANY MODALITY
## ---------------------------------------------------------------------------
# mat : features x cells (counts for RNA, intensities for protein; sparse ok)
# meta: data.frame, rownames = cell ids, with columns
#         patient, section, ctype, and the coordinate columns named in the spec
is_pos <- function(v, method = c("count", "quantile"), value, ref = v) {
  method <- match.arg(method)
  if (method == "count") v >= value else v >= quantile(ref, value, na.rm = TRUE)
}

build_cells <- function(mat, meta, spec) {
  stopifnot(!is.null(spec$resolved))
  common <- intersect(colnames(mat), rownames(meta))
  mat  <- mat[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]
  rs   <- spec$resolved
  feats <- unique(c(unlist(lapply(rs$axes, function(a) c(a$lig, a$rec))),
                    rs$strict_all, rs$strict_any, rs$mac))
  vals <- t(as.matrix(mat[feats, , drop = FALSE]))              # cells x features

  ct <- as.character(meta$ctype)
  if (!is.null(spec$ctype_map))
    ct <- ifelse(ct %in% names(spec$ctype_map), spec$ctype_map[ct], ct)

  cells <- tibble(cell = rownames(meta), patient = meta$patient, section = meta$section,
                  ctype = ct,
                  x_um = meta[[spec$x_col]] * spec$um_per_unit,
                  y_um = meta[[spec$y_col]] * spec$um_per_unit)
  vdf <- as_tibble(vals); names(vdf) <- paste0("val__", colnames(vals))
  cells <- bind_cols(cells, vdf)
  if ("histology" %in% names(meta)) cells$histology <- as.character(meta$histology)   # LUAD / LUSC

  sec_id <- paste(cells$patient, cells$section, sep = "|")
  if (rs$strict_ok) {
    cells$ec_strict <- NA
    for (ix in split(seq_len(nrow(cells)), sec_id)) {
      pf <- function(f) is_pos(vals[ix, f], spec$strict_pos$method, spec$strict_pos$value,
                               ref = vals[ix, f])
      a_ok <- Reduce(`&`, lapply(rs$strict_all, pf))
      o_ok <- Reduce(`|`, lapply(rs$strict_any, pf))
      cells$ec_strict[ix] <- a_ok & o_ok
    }
  }
  if (rs$mac_ok) cells$mac_marker <- rowSums(vals[, rs$mac, drop = FALSE])
  cells
}

# Strata for a given tile size
tile_strata <- function(df, tile_um) {
  paste(df$patient, df$section, floor(df$x_um / tile_um), floor(df$y_um / tile_um), sep = "|")
}

# Spatial strata for permutations 
add_tiles <- function(cells, tile_um = 500) {
  cells %>% mutate(stratum = paste(patient, section, floor(x_um / tile_um),
                                   floor(y_um / tile_um), sep = "|"))
}

## ---------------------------------------------------------------------------
## 3b. ADAPTER: CosMx Seurat object (from Seurat::LoadNanostring) -> mat + meta
## ---------------------------------------------------------------------------
# LoadNanostring gives a Seurat object with one FOV image per field of view
#   tissue_col : metadata column identifying the tissue section / sample
#   ctype_col  : cell-type annotation column
#   patient_from_tissue : maps tissue -> patient. 
#   fov_as_section : TRUE if coordinates are LOCAL to each FOV (see coord_check)
cosmx_inputs_from_seurat <- function(obj, tissue_col = NULL, ctype_col,
                                     patient_from_tissue = function(x) sub("_Rep[0-9]+$", "", x),
                                     assay = SeuratObject::DefaultAssay(obj),
                                     fov_as_section = FALSE) {
  fovs <- SeuratObject::Images(obj)
  coords <- map_dfr(fovs, function(f) {
    xy <- SeuratObject::GetTissueCoordinates(obj[[f]], which = "centroids")
    if (!"cell" %in% names(xy)) xy$cell <- rownames(xy)
    xy$fov_name <- f
    xy[, c("cell", "x", "y", "fov_name")]
  })
  coords <- coords[!duplicated(coords$cell), ]
  meta <- obj@meta.data
  meta$cell <- rownames(meta)
  meta <- inner_join(meta, coords, by = "cell")
  rownames(meta) <- meta$cell
  meta$x_px <- meta$x; meta$y_px <- meta$y            # pixels; spec_cosmx converts (0.12 um/px)
  # If each image is one tissue (e.g. Lung5_Rep1, Lung12, ...), the image name IS the tissue
  tissue <- if (is.null(tissue_col)) as.character(meta$fov_name) else as.character(meta[[tissue_col]])
  meta$patient <- patient_from_tissue(tissue)
  meta$section <- if (fov_as_section) paste(tissue, meta$fov_name, sep = "|") else tissue
  meta$ctype   <- as.character(meta[[ctype_col]])
  mat <- SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
  list(mat = mat, meta = meta)
}

# Are coordinates global (FOVs tile one slide) or local (each FOV starts near 0)?
coord_check <- function(meta) {
  meta %>% group_by(section, fov_name) %>%
    summarise(x_min = min(x_px), x_max = max(x_px), y_min = min(y_px), y_max = max(y_px),
              n = n(), .groups = "drop") %>% print(n = 50)
}

## ---------------------------------------------------------------------------
## 4. GEOMETRY AND PERMUTATION HELPERS
## ---------------------------------------------------------------------------
# positive if ANY of the axis features is positive (OR across PD-L1 / PD-L2, etc.)
any_pos <- function(df, cols, method, value, ref_df = df) {
  m <- vapply(cols, function(cn) is_pos(df[[cn]], method, value, ref = ref_df[[cn]]),
              logical(nrow(df)))
  if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
  rowSums(m) > 0
}

# centroid distance from each row of q to the nearest row of d (microns)
nearest_default <- function(q, d) {
  if (nrow(d) == 0) return(rep(Inf, nrow(q)))
  drop(RANN::nn2(as.matrix(d[, c("x_um", "y_um")]),
                 as.matrix(q[, c("x_um", "y_um")]), k = 1)$nn.dists)
}
# OPTIONAL boundary-to-boundary distance (CosMx / Xenium cell polygons, in microns)
make_nearest_poly <- function(polys) {
  function(q, d) {
    if (nrow(d) == 0) return(rep(Inf, nrow(q)))
    qg <- polys[match(q$cell, names(polys))]
    dg <- polys[match(d$cell, names(polys))]
    j  <- sf::st_nearest_feature(qg, dg)
    as.numeric(sf::st_distance(qg, dg[j], by_element = TRUE))
  }
}

# Extract cell polygons from a Seurat/Nanostring object (segmentation boundaries)
cosmx_polygons_from_seurat <- function(obj, cells, px_um = 0.12) {
  polys <- map(SeuratObject::Images(obj), function(f) {
    fov <- obj[[f]]
    if (!"segmentation" %in% SeuratObject::Boundaries(fov)) {
      message("no segmentation boundary in ", f, " (reload with segmentations, or use centroids)")
      return(NULL)
    }
    v <- SeuratObject::GetTissueCoordinates(fov, which = "segmentation")   # vertices: x, y, cell
    v <- v[v$cell %in% cells, c("x", "y", "cell")]
    if (!nrow(v)) return(NULL)
    v$x <- v$x * px_um; v$y <- v$y * px_um
    sfheaders::sf_polygon(v, x = "x", y = "y", polygon_id = "cell", close = TRUE)
  })
  polys <- do.call(rbind, compact(polys))
  setNames(sf::st_geometry(polys), polys$cell)
}
radius_nn <- function(query_xy, data_xy, radius, kmax = 50) {
  if (nrow(data_xy) == 0 || nrow(query_xy) == 0)
    return(list(idx = matrix(0L, nrow(query_xy), 1), dist = matrix(Inf, nrow(query_xy), 1)))
  r <- RANN::nn2(data_xy, query_xy, k = min(kmax, nrow(data_xy)),
                 searchtype = "radius", radius = radius)
  list(idx = r$nn.idx, dist = r$nn.dists)
}

# One-sided permutation test: do `group` cells have a higher rate of `outcome`
# than expected? Group labels are permuted within `strata` (vectorised)
perm_test <- function(outcome, group, strata, n_perm = 1000, eps = 1e-3, keep_null = TRUE) {
  n <- length(group)
  s <- as.integer(factor(strata)); base <- order(s)
  obs  <- mean(outcome[group])
  null <- vapply(seq_len(n_perm), function(i) {
    shuf <- order(s, runif(n))                   # random order WITHIN each stratum block
    g <- logical(n); g[base] <- group[shuf]
    mean(outcome[g])
  }, numeric(1))
  ng  <- sum(group)
  out <- tibble(n_group = ng, n_total = n, obs_rate = obs, exp_rate = mean(null),
                obs_events = round(obs * ng), exp_events = mean(null) * ng,
                log2OE    = log2((obs + eps) / (mean(null) + eps)),
                log2OE_pc = log2((obs * ng + 0.5) / (mean(null) * ng + 0.5)),
                p_perm = (1 + sum(null >= obs)) / (n_perm + 1))
  if (keep_null) out$null_events <- list(round(null * ng))
  out
}

## ---------------------------------------------------------------------------
## 5. ONE SECTION x ONE AXIS x ONE SCENARIO x DISTANCE GRID
## ---------------------------------------------------------------------------
prox_section <- function(sec, ax, dists, pos_method, pos_value,
                         ec_def = "base", purity = FALSE,
                         nearest_fun = nearest_default, n_perm = 1000, min_pos = 5,
                         keep_null = TRUE, tile_um = NULL) {
  sec <- drop_na(sec)
  lig_cols <- paste0("val__", ax$lig); rec_cols <- paste0("val__", ax$rec)
  ec_ok <- if (ec_def == "strict") sec$ec_strict else rep(TRUE, nrow(sec))
  ec <- sec[sec$ctype == ENDO & ec_ok, ]
  tc <- sec[sec$ctype %in% ax$t_types, ]
  if (nrow(ec) < 2 * min_pos || nrow(tc) < min_pos) return(NULL)

  lig <- any_pos(ec, lig_cols, pos_method, pos_value, ref_df = sec)
  if (purity) {                                  # drop (not relabel) contaminated ligand+ ECs
    bad <- lig & ec$mac_marker >= 1
    ec <- ec[!bad, ]; lig <- lig[!bad]
  }
  pd1 <- any_pos(tc, rec_cols, pos_method, pos_value, ref_df = sec)
  if (sum(lig) < min_pos || sum(pd1) < min_pos || sum(!lig) < min_pos) return(NULL)

  s_ec <- if (is.null(tile_um)) ec$stratum else tile_strata(ec, tile_um)   # permutation strata
  s_tc <- if (is.null(tile_um)) tc$stratum else tile_strata(tc, tile_um)
  d_ec <- nearest_default(ec, tc[pd1, ])             # each endothelial cell -> nearest receptor+ T cell
  d_t  <- nearest_default(tc, ec)                    # each T  -> nearest endothelial cell

  map_dfr(dists, function(d) {
    bind_rows(
      perm_test(d_ec <= d, lig, s_ec, n_perm, keep_null = keep_null) %>%
        mutate(test = "ligand+ vs ligand- endothelial cells: receptor+ T within d"),
      perm_test(d_t  <= d, pd1, s_tc, n_perm, keep_null = keep_null) %>%
        mutate(test = "receptor+ vs receptor- T: endothelial cell within d")) %>%
      mutate(dist_um = d)
  })
}

make_scenarios <- function(spec) {
  rs <- spec$resolved
  expand_grid(pos_method = spec$pos$method, pos_value = spec$pos$values,
              ec_def = c("base", if (rs$strict_ok) "strict"),
              purity = c(FALSE, if (rs$mac_ok) TRUE),
              tile_um = spec$tiles %||% 500)
}

run_grid <- function(cells, spec, axis_name, scen, n_perm = 200, ...) {
  ax   <- spec$resolved$axes[[axis_name]]
  secs <- cells %>% group_split(patient, section)
  map_dfr(secs, function(sec) {
    map_dfr(seq_len(nrow(scen)), function(i) {
      sc <- scen[i, ]
      r  <- prox_section(sec, ax, spec$dists, sc$pos_method, sc$pos_value,
                         sc$ec_def, sc$purity, n_perm = n_perm, tile_um = sc$tile_um, ...)
      if (is.null(r)) return(NULL)
      bind_cols(tibble(modality = spec$name, axis = axis_name,
                       t_set = ax$t_set %||% "all", base_axis = ax$base_axis %||% axis_name,
                       receptor = paste(ax$rec, collapse = "+"),
                       patient = sec$patient[1], section = sec$section[1]), sc, r)
    })
  })
}

# All measurable axes of one modality. If the cell table has a `histology` column, each histology
# is analysed SEPARATELY (scenarios, permutations and summaries never mix them) and the result
# carries a `histology` column. `scen` overrides the default scenario grid (e.g. primary only).
run_modality <- function(cells, spec, n_perm = 200, axes = names(spec$resolved$axes),
                         by = "histology", scen = NULL, ...) {
  if (!is.null(by) && by %in% names(cells)) {
    grp <- as.character(cells[[by]])
    if (anyNA(grp)) message("[", spec$name, "] ", sum(is.na(grp)), " cells with missing ", by, " excluded")
    levs <- sort(unique(grp[!is.na(grp)]))
    message("[", spec$name, "] running separately by ", by, ": ", paste(levs, collapse = ", "))
    return(map_dfr(levs, function(h) {
      sub <- cells[which(grp == h), ]
      if (!ENDO %in% sub$ctype) {
        message("[", spec$name, "] ", h, ": no '", ENDO, "' cells, skipped"); return(NULL)
      }
      r <- run_modality(sub, spec, n_perm, axes, by = NULL, scen = scen, ...)
      if (nrow(r)) r[[by]] <- h
      r
    }))
  }
  stopifnot(length(axes) > 0)
  present <- unique(cells$ctype)
  if (!ENDO %in% present)
    stop("no '", ENDO, "' cells in cells$ctype (found: ", paste(present, collapse = ", "),
         "); fix ctype / ctype_map")
  keep_ax <- vapply(axes, function(a) any(spec$resolved$axes[[a]]$t_types %in% present), logical(1))
  if (any(!keep_ax))
    message("[", spec$name, "] skipped (no cells of the T-cell type in ctype): ",
            paste(axes[!keep_ax], collapse = ", "))
  axes <- axes[keep_ax]
  stopifnot(length(axes) > 0)
  scen <- scen %||% make_scenarios(spec)
  res  <- map_dfr(axes, function(a) run_grid(cells, spec, a, scen, n_perm, ...))
  if (!nrow(res)) message("[", spec$name, "] no section had enough positive cells; check thresholds")
  res
}

## ---------------------------------------------------------------------------
## 6. PATIENT-LEVEL SUMMARIES, CURVES, CROSS-MODALITY CONSISTENCY
## ---------------------------------------------------------------------------
scen_keys <- c("pos_method", "pos_value", "ec_def", "purity", "tile_um")

summarise_patients <- function(res) {
  sk <- intersect(scen_keys, names(res)); hk <- intersect("histology", names(res))
  by_pat <- res %>%
    group_by(across(all_of(c("modality", hk, "axis", "patient", "test", "dist_um", sk)))) %>%
    summarise(log2OE = mean(log2OE), .groups = "drop")          # sections averaged within patient
  overall <- by_pat %>%
    group_by(across(all_of(c("modality", hk, "axis", "test", "dist_um", sk)))) %>%
    summarise(n_patients = n(), mean_log2OE = mean(log2OE), n_patients_pos = sum(log2OE > 0),
              # with n = 5 patients the smallest possible one-sided p is 1/32 (0.031)
              p_wilcox = tryCatch(wilcox.test(log2OE, mu = 0, alternative = "greater")$p.value,
                                  error = function(e) NA_real_),
              .groups = "drop")
  list(by_patient = by_pat, overall = overall)
}

plot_curves <- function(summ, prefix = "prox") {
  d0 <- summ$by_patient
  d0$scen <- paste0("thr=", d0$pos_value, " | ", d0$ec_def, " | purity=", d0$purity,
                    if ("tile_um" %in% names(d0)) paste0(" | tile=", d0$tile_um) else "")
  by <- intersect(c("modality", "histology", "axis"), names(d0))
  d0 %>% group_by(across(all_of(by))) %>% group_split() %>% walk(function(d) {
    h <- if ("histology" %in% names(d)) d$histology[1] else "all"
    p <- ggplot(d, aes(dist_um, log2OE, colour = patient, group = interaction(patient, scen))) +
      geom_hline(yintercept = 0, linetype = 2) + geom_line(alpha = 0.6) + geom_point(size = 1) +
      facet_grid(test ~ scen) +
      labs(title = paste(d$modality[1], h, d$axis[1]), x = "Distance threshold (um)",
           y = "log2(observed / expected)") + theme_bw()
    ggsave(sprintf("%s_%s_%s_%s.pdf", prefix, d$modality[1], h, d$axis[1]), p, width = 16, height = 6)
  })
}

# Keep only the T-cell subsets named in spec$primary$t_sets 
filter_tset <- function(d, pr) {
  if (is.null(pr$t_sets) || !"t_set" %in% names(d)) d else dplyr::filter(d, t_set %in% pr$t_sets)
}

# Keep only the pre-specified primary tile size (other tile sizes are sensitivity analysis)
filter_tile <- function(d, pr) {
  if (is.null(pr$tile_um) || !"tile_um" %in% names(d)) d else dplyr::filter(d, near(tile_um, pr$tile_um))
}

# One row per modality x axis x test at each modality's primary scenario
primary_summary <- function(res, specs) {
  map_dfr(specs, function(sp) {
    pr <- sp$primary
    d <- res %>%
      filter(modality == sp$name, near(pos_value, pr$pos_value), ec_def == pr$ec_def,
             purity == pr$purity, dist_um == pr$dist_um) %>%
      filter_tset(pr) %>% filter_tile(pr)
    hk <- intersect("histology", names(d))
    d %>%
      group_by(across(all_of(c("modality", hk, "axis", "test", "patient")))) %>%
      summarise(log2OE = mean(log2OE), .groups = "drop") %>%
      group_by(across(all_of(c("modality", hk, "axis", "test")))) %>%
      summarise(n_patients = n(), mean_log2OE = mean(log2OE),
                n_patients_pos = sum(log2OE > 0), .groups = "drop")
  })
}

## ---------------------------------------------------------------------------
## 6b. POOLED (STRATIFIED) PERMUTATION TEST ACROSS SECTIONS
## ---------------------------------------------------------------------------
# Labels were permuted within tiles, independently in each section, so null draw i
# from every section can be ADDED to give draw i of a pooled statistic
#   A "events"    : total observed events summed over sections. Sections count in
#                   proportion to how many events they can contribute.
#   B "patient_z" : each patient's events (replicate sections summed first) are
#                   standardised against that patient's own null, then averaged, so
#                   every patient counts equally.
pool_keys <- c("modality", "histology", "axis", "t_set", "test", "dist_um", "pos_method",
               "pos_value", "ec_def", "purity", "tile_um")

pooled_one <- function(x) {
  L <- lengths(x$null_events)
  if (length(unique(L)) != 1) return(tibble())      # all sections must share n_perm
  n_perm <- L[1]
  obs <- x$obs_events
  N   <- do.call(rbind, x$null_events)              # sections x n_perm
  tol <- 1e-9

  # A: pooled events
  TA_obs  <- sum(obs)
  TA_null <- colSums(N)
  exp_ev  <- mean(TA_null)

  # B: patient-equal z
  pat   <- as.integer(factor(x$patient))
  Pobs  <- as.numeric(rowsum(obs, pat))
  Pnull <- rowsum(N, pat)                           # patients x n_perm
  mu    <- rowMeans(Pnull)
  sdv   <- apply(Pnull, 1, sd)
  ok    <- sdv > 0                                 
  zobs  <- (Pobs - mu) / pmax(sdv, 1e-12); zobs[!ok] <- 0
  znull <- (Pnull - mu) / pmax(sdv, 1e-12); znull[!ok, ] <- 0
  TB_obs  <- mean(zobs)
  TB_null <- colMeans(znull)

  tibble(n_sections = nrow(x), n_patients = length(unique(x$patient)),
         obs_events = TA_obs, exp_events = exp_ev,
         log2_pooled_OE = log2((TA_obs + 0.5) / (exp_ev + 0.5)),
         p_enriched_events   = (1 + sum(TA_null >= TA_obs - tol)) / (n_perm + 1),
         p_depleted_events   = (1 + sum(TA_null <= TA_obs + tol)) / (n_perm + 1),
         mean_patient_z      = TB_obs,
         p_enriched_patientz = (1 + sum(TB_null >= TB_obs - tol)) / (n_perm + 1),
         p_depleted_patientz = (1 + sum(TB_null <= TB_obs + tol)) / (n_perm + 1),
         n_perm = n_perm)
}

# min_sections = 1 so a single-section dataset (e.g. one Xenium slide) still returns rows;
# with one section the pooled test is simply that section's own permutation test
dedupe_test2 <- function(res) {
  if (!all(c("receptor", "t_set") %in% names(res))) return(res)  
  is_t2 <- grepl("^receptor", res$test)
  key <- paste(res$modality, res$receptor, res$t_set, res$patient, res$section, res$pos_method,
               res$pos_value, res$ec_def, res$purity, res$tile_um, res$dist_um, sep = "|")
  dup <- rep(FALSE, nrow(res)); dup[is_t2] <- duplicated(key[is_t2])
  res <- res[!dup, ]
  is_t2 <- grepl("^receptor", res$test)
  res$axis[is_t2] <- paste0(res$receptor[is_t2], "_", res$t_set[is_t2], " (receptor test)")
  res
}

pooled_perm_test <- function(res, min_sections = 1, dedupe = TRUE, drop_test2_axes = character(0)) {
  stopifnot("null_events" %in% names(res))
  if (dedupe) res <- dedupe_test2(res)
  if (length(drop_test2_axes))
    res <- res %>% filter(!(grepl("^receptor", test) & axis %in% drop_test2_axes))
  keys <- intersect(pool_keys, names(res))
  out <- res %>%
    group_by(across(all_of(keys))) %>%
    group_modify(~ pooled_one(.x)) %>%
    ungroup()
  if (!"n_sections" %in% names(out)) return(out)
  filter(out, n_sections >= min_sections)
}

# The pre-specified primary scenario per modality, with Holm adjustment across those rows
primary_pooled <- function(pooled, specs) {
  map_dfr(specs, function(sp) {
    pr <- sp$primary
    pooled %>%
      filter(modality == sp$name, near(pos_value, pr$pos_value), ec_def == pr$ec_def,
             purity == pr$purity, dist_um == pr$dist_um) %>%
      filter_tset(pr) %>% filter_tile(pr)
  }) %>%
    group_by(across(any_of("histology"))) %>%            # Holm within each histology
    mutate(p_holm_events   = p.adjust(p_enriched_events,   "holm"),
           p_holm_patientz = p.adjust(p_enriched_patientz, "holm")) %>%
    ungroup()
}

# The primary scenario across ALL tile sizes
tile_sensitivity <- function(pooled, sp) {
  pr <- sp$primary
  pooled %>%
    filter(modality == sp$name, near(pos_value, pr$pos_value), ec_def == pr$ec_def,
           purity == pr$purity, dist_um == pr$dist_um) %>%
    filter_tset(pr) %>%
    select(any_of(c("histology", "axis", "test", "tile_um", "obs_events", "exp_events",
                    "log2_pooled_OE", "p_enriched_events", "p_depleted_events",
                    "mean_patient_z", "p_enriched_patientz")))
}

# Pooled observed/expected across the distance grid (filled = p_enriched_events < 0.05, uncorrected)
plot_pooled <- function(pooled, file = "pooled_curves.pdf", pos_value = NULL, ec_def = "base",
                        purity = FALSE, tile = NULL) {
  d <- pooled %>% filter(ec_def == !!ec_def, purity == !!purity)
  if (!is.null(pos_value)) d <- d %>% filter(near(pos_value, !!pos_value))
  if ("tile_um" %in% names(d)) { tl <- tile %||% 500; d <- d %>% filter(near(tile_um, !!tl)) }
  d$panel <- paste0(d$modality, if ("histology" %in% names(d)) paste0(" ", d$histology) else "",
                    " thr=", d$pos_value)
  p <- ggplot(d, aes(dist_um, log2_pooled_OE, colour = axis)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line() +
    geom_point(aes(shape = p_enriched_events < 0.05), size = 2) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "p < 0.05") +
    facet_grid(test ~ panel) +
    labs(x = "Distance threshold (um)", y = "log2((obs + 0.5) / (exp + 0.5)) pooled events") +
    theme_bw()
  ggsave(file, p, width = 10, height = 6); invisible(p)
}

strip_nulls <- function(res) select(res, -any_of("null_events"))

## ---------------------------------------------------------------------------
## 7. MIXED MODEL ADJUSTING FOR LOCAL T-CELL DENSITY
## ---------------------------------------------------------------------------
ec_level_table <- function(sec, ax, d, pos_method, pos_value, local_r = 50,
                           nearest_fun = nearest_default) {
  ec <- sec[sec$ctype == ENDO, ]; tc <- sec[sec$ctype %in% ax$t_types, ]
  lig <- any_pos(ec, paste0("val__", ax$lig), pos_method, pos_value, ref_df = sec)
  pd1 <- any_pos(tc, paste0("val__", ax$rec), pos_method, pos_value, ref_df = sec)
  nnT <- radius_nn(as.matrix(ec[, c("x_um", "y_um")]),
                   as.matrix(tc[, c("x_um", "y_um")]), local_r, kmax = 500)
  tibble(patient = sec$patient[1], section = sec$section[1], stratum = ec$stratum, lig_pos = lig,
         has_nb  = nearest_fun(ec, tc[pd1, ]) <= d,
         local_T = rowSums(nnT$idx > 0))
}

fit_ec_model <- function(cells, spec, axis_name, d = spec$primary$dist_um,
                         pos_method = spec$pos$method, pos_value = spec$primary$pos_value) {
  ax  <- spec$resolved$axes[[axis_name]]
  tab <- cells %>% group_split(patient, section) %>%
    map_dfr(ec_level_table, ax = ax, d = d, pos_method = pos_method, pos_value = pos_value)
  f <- if (length(unique(tab$patient)) >= 2) {
    has_nb ~ lig_pos + scale(log1p(local_T)) + (1 | patient/section)
  } else {
    has_nb ~ lig_pos + scale(log1p(local_T)) + (1 | stratum)
  }
  m <- glmmTMB(f, family = binomial, data = tab)
  list(model = m, odds_ratio = exp(fixef(m)$cond["lig_posTRUE"]), data = tab)
}

## ---------------------------------------------------------------------------
## 8. CONTEXT: who supplies the ligand+ neighbours of receptor+ T cells?
##    Null = shuffle cell-type labels among ligand+ non-T cells within strata.
##    log2OE ~ 0 for ECs = neither over- nor under-represented given their share
##    of ligand+ cells (report share_of_lig_pos alongside).
## ---------------------------------------------------------------------------
sender_context <- function(sec, ax, d, pos_method, pos_value,
                           senders = c(ENDO, "Myeloid", "Malignant", "Fibroblast"),
                           n_perm = 1000, kmax = 50) {
  senders <- intersect(senders, unique(sec$ctype))
  tc  <- sec[sec$ctype %in% ax$t_types, ]
  pd1 <- any_pos(tc, paste0("val__", ax$rec), pos_method, pos_value, ref_df = sec)
  src <- sec[sec$ctype %in% senders, ]
  src <- src[any_pos(src, paste0("val__", ax$lig), pos_method, pos_value, ref_df = sec), ]
  if (sum(pd1) < 5 || nrow(src) < 5 || !length(senders)) return(NULL)

  nn  <- radius_nn(as.matrix(tc[pd1, c("x_um", "y_um")]),
                   as.matrix(src[, c("x_um", "y_um")]), d, kmax)
  ok  <- nn$idx > 0
  idx <- pmax(nn$idx, 1L)
  rate <- function(ct, s) mean(rowSums((matrix(ct[idx], nrow(idx)) == s) & ok) > 0)

  s_int <- as.integer(factor(src$stratum)); base <- order(s_int)
  ct0   <- src$ctype
  obs   <- vapply(senders, function(s) rate(ct0, s), numeric(1))
  null  <- replicate(n_perm, {
    shuf <- order(s_int, runif(length(s_int)))
    ct <- character(length(ct0)); ct[base] <- ct0[shuf]
    vapply(senders, function(s) rate(ct, s), numeric(1))
  })
  if (is.null(dim(null))) null <- matrix(null, nrow = length(senders))
  tibble(patient = sec$patient[1], section = sec$section[1], dist_um = d, sender = senders,
         share_of_lig_pos = vapply(senders, function(s) mean(ct0 == s), numeric(1)),
         obs_rate = obs, exp_rate = rowMeans(null),
         log2OE = log2((obs + 1e-3) / (rowMeans(null) + 1e-3)),
         p_enriched = (1 + rowSums(null >= obs)) / (n_perm + 1),
         p_depleted = (1 + rowSums(null <= obs)) / (n_perm + 1))
}

## ---------------------------------------------------------------------------
## 9. EXAMPLE DRIVERS (edit and run)
## ---------------------------------------------------------------------------

cosmx_seu <- qs_read(dpath("cosmxObj_allCells.qs2"))
xen_seu <- qs_read(dpath("annotated_xeniumObj.qs2"))
mif <- qs_read(dpath("codex_allSamples.qs2"))
out_dir <- dpath("SpatialAdjOutputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (TRUE) {

  # ---- CosMx (Seurat object; 7 sections / 5 patients) ----------------------
  names(cosmx_seu@meta.data); SeuratObject::Images(cosmx_seu)     # 8 images = 8 tissue sections
  inp  <- cosmx_inputs_from_seurat(cosmx_seu, ctype_col = "predicted.id") 
  coord_check(inp$meta)   # global coordinates: ranges are large and plausible for a slide;
                          # each FOV starting near 0 means LOCAL coordinates -> fov_as_section = TRUE
  print(table(inp$meta$section, inp$meta$patient))   # confirm 8 sections -> 5 patients
  mat <- inp$mat; meta <- inp$meta
  spec_cosmx <- resolve_spec(spec_cosmx, rownames(mat)); report_spec(spec_cosmx)
  cosmx     <- build_cells(mat, meta, spec_cosmx) %>% add_tiles(500)

  set.seed(1)                                                         
  # (1) Sensitivity grid: thresholds x endothelial definitions x macrophage filter x TILE SIZES
  #     (250 / 500 / 1000 um) x distances (10-30 um).
  res_cosmx_grid <- run_modality(cosmx, spec_cosmx, n_perm = 2500, keep_null = FALSE)

  # (2) Pre-specified primary scenario at every tile size
  pr <- spec_cosmx$primary
  scen_primary <- tibble(pos_method = spec_cosmx$pos$method, pos_value = pr$pos_value,
                         ec_def = pr$ec_def, purity = pr$purity, tile_um = spec_cosmx$tiles)
  res_cosmx <- run_modality(cosmx, spec_cosmx, n_perm = 2500, scen = scen_primary)

  saveRDS(res_cosmx, file.path(out_dir, "cosmx_res_with_nulls.rds"))                      # needed for pooled tests
  readr::write_tsv(strip_nulls(res_cosmx), file.path(out_dir, "cosmx_adjRes_primary_tiles.tsv"))
  readr::write_tsv(res_cosmx_grid, file.path(out_dir, "cosmx_adjRes_grid.tsv"))

  # Pooled stratified permutation test, by axis AND T-cell subset (CD4 / CD8 / T_all)
  pooled <- pooled_perm_test(res_cosmx)             # Test 2 duplicates removed automatically
  readr::write_tsv(pooled, file.path(out_dir, "cosmx_adjRes_pooled.tsv"))

  # primary_pooled(pooled, list(spec_cosmx))          # primary tile only, Holm across CD4/CD8 rows
  # tile_sensitivity(pooled, spec_cosmx)              # same scenario at 250 / 500 / 1000 um tiles
#   plot_pooled(pooled, "cosmx_pooled_curves.pdf", pos_value = 1)

  # Contact-scale analysis with boundary-to-boundary distances (polygons)
#   keep    <- cosmx$cell[cosmx$ctype %in% c(ENDO, "CD4", "CD8")]        # only cells used in Tests 1-2
#   polys   <- cosmx_polygons_from_seurat(cosmx_seu, keep)               # named sfc, microns
#   cosmx_b <- cosmx %>% filter(cell %in% names(polys))
#   spec_contact <- modifyList(spec_cosmx, list(dists = c(0.6, 2, 5, 10, 20), tiles = 500))
#   res_contact  <- run_modality(cosmx_b, spec_contact, n_perm = 200,
#                                nearest_fun = make_nearest_poly(polys))

  # ---- Xenium (single section, n = 1; TIGIT axes) --------------------------
  mat  <- Seurat::GetAssayData(xen_seu, layer = "counts")
  meta <- xen_seu@meta.data             # needs patient, section, ctype, x_centroid, y_centroid
  meta$patient <- "X1"; meta$section <- "X1"            # one tissue: constant labels
  spec_xenium <- resolve_spec(spec_xenium, rownames(mat)); report_spec(spec_xenium)
  xen     <- build_cells(mat, meta, spec_xenium) %>% add_tiles(500)
  
  set.seed(1)
  res_xen <- run_modality(xen, spec_xenium, n_perm = 2500)
  pooled_xen <- pooled_perm_test(res_xen)               # one section: equals the section-level test
  
  saveRDS(res_xen, file.path(out_dir, "xenium_res_with_nulls.rds"))                      # needed for pooled tests
  readr::write_tsv(strip_nulls(res_xen), file.path(out_dir, "xenium_adjRes_primary_tiles.tsv"))
  readr::write_tsv(pooled_xen, file.path(out_dir, "xenium_adjRes_pooled.tsv"))

#   primary_pooled(pooled_xen, list(spec_xenium))
  # What this can and cannot say: p-values describe THIS section only (no replicate, so no
  # between-tissue variability; CosMx replicates of one patient disagreed). Wilcoxon and
  # patient-level summaries are not meaningful with n = 1: use p_perm, obs/exp events and
  # log2OE_pc. Optional: repeat the primary scenario in large spatial blocks (e.g. 2 x 2
  # quadrants) to see whether an effect is localised; blocks are NOT independent replicates.

  # ---- CODEX (protein; h5ad -> R via zellkonverter::readH5AD or anndata) ----
  # X = cells x markers intensity matrix; obs = centroids, cell type, sample
  # mat <- t(as.matrix(assay(sce, "X")))   or   t(as.matrix(intensity_df))   # -> markers x cells
  # BEFORE running, edit spec_codex: ctype_map (your labels -> Endothelial/CD4/CD8/Myeloid/...),
  # marker aliases, x_col / y_col and um_per_unit (pixels -> microns; often ~0.38 um/px at 20x).
  # meta needs patient (one per patient), section (core/ROI/sample), ctype, coordinates AND
  # histology ("LUAD" / "LUSC"): LUAD and LUSC live in one object, so run_modality() splits on
  # `histology` and analyses each separately (results carry a `histology` column).
  mat <- Seurat::GetAssayData(mif, layer = "counts")
  meta <- mif@meta.data

  spec_codex <- resolve_spec(spec_codex, rownames(mat)); report_spec(spec_codex)
  codex <- build_cells(mat, meta, spec_codex) %>% add_tiles(500)
  table(codex$ctype, codex$histology)                    # Endothelial / CD4 / CD8 must exist in BOTH

  # Large cohort: run the FULL GRID cheaply (no null draws), then the PRIMARY scenario only
  # (primary distance, all tile sizes) with 2,500 permutations and null draws.
  set.seed(1)
  res_codex_grid <- run_modality(codex, spec_codex, n_perm = 2500, keep_null = FALSE)
  readr::write_tsv(res_codex_grid, file.path(out_dir, "codex_adjRes_grid.tsv"))

  pr <- spec_codex$primary
  spec_primary <- modifyList(spec_codex, list(dists = pr$dist_um))
  scen_primary <- tibble(pos_method = spec_codex$pos$method, pos_value = pr$pos_value,
                         ec_def = pr$ec_def, purity = pr$purity, tile_um = spec_codex$tiles)
  res_codex <- run_modality(codex, spec_primary, n_perm = 2500, scen = scen_primary)
  
  pooled_codex <- pooled_perm_test(res_codex)            # grouped by histology automatically
  primary_pooled(pooled_codex, list(spec_codex))         # Holm within each histology
  tile_sensitivity(pooled_codex, spec_codex)
  qs_save(res_codex, file.path(out_dir, "codex_res_with_nulls.qs2"))                      # needed for pooled tests

  readr::write_tsv(strip_nulls(res_codex), file.path(out_dir, "codex_adjRes_primary_tiles.tsv"))
  readr::write_tsv(pooled_codex, file.path(out_dir, "codex_adjRes_pooled.tsv"))

  # ---- Combine (never pool cells across modalities) -------------------------
  res_all <- bind_rows(res_cosmx, res_xen, res_codex)
  readr::write_tsv(res_all, file.path(out_dir, "allModes_adjRes_primary_unpooled.tsv"))

  summ    <- summarise_patients(res_all)
  # plot_curves(summ, "prox")
  primary_summary(res_all, list(spec_cosmx, spec_xenium, spec_codex))

  # ---- Primary analysis: pre-specified scenario, more permutations ----------
  # (edit the spec$primary values first; then re-run that single scenario with n_perm = 5000)

  # ---- Mixed model and sender context, per modality x axis ------------------
  # fit_cd8 <- fit_ec_model(cosmx, spec_cosmx, "PD1_PDL1_CD8"); summary(fit_cd8$model); fit_cd8$odds_ratio
  # fit_cd4 <- fit_ec_model(cosmx, spec_cosmx, "PD1_PDL1_CD4"); summary(fit_cd4$model); fit_cd4$odds_ratio
  # ctx <- cosmx %>% group_split(patient, section) %>%
  #   map_dfr(sender_context, ax = spec_cosmx$resolved$axes$PD1_PDL1_CD8, d = 10,
  #           pos_method = "count", pos_value = 1)
  # ctx %>% group_by(sender) %>%
  #   summarise(mean_log2OE = mean(log2OE), mean_share = mean(share_of_lig_pos))
}