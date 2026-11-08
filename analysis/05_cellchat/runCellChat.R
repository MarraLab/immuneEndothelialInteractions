# =============================================================================
# Per-sample CellChat (v2) ligand-receptor inference between ECs and T cells
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(Matrix)
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(glmmTMB); library(qs2)
})

## ---------------------------------------------------------------------------
## 0. INPUT AND SETTINGS
## ---------------------------------------------------------------------------
# Merged Seurat object with a log-normalised "data" layer AND a "counts" layer.
# Required metadata columns:
#   sample, dataset, histology ("LUAD"/"LUSC"), platform ("Flex"/Other) 
#   celltype_broad
#   ec_subtype (ECs only), t_state (T cells only, used when celltype_broad=="T")
seu <- qs_read(dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)                       

## ---------------------------------------------------------------------------
## 0a. RESOLVE HYBRID CD8 LINEAGE LABELS (Lin1_x_Lin2_y -> two separate cells)
## ---------------------------------------------------------------------------
hybrid_pattern <- "^(Lin1_[0-9]+(?:\\.[0-9]+)?)_(Lin2_[0-9]+(?:\\.[0-9]+)?)$"

meta0            <- seu@meta.data
t_state_raw      <- meta0$t_state
looks_hybrid     <- !is.na(t_state_raw) & grepl("Lin1.*Lin2|Lin2.*Lin1", t_state_raw)
matches_pattern  <- !is.na(t_state_raw) & grepl(hybrid_pattern, t_state_raw)

hybrid_labels_found     <- sort(unique(t_state_raw[looks_hybrid & matches_pattern]))
unmatched_hybrid_labels <- sort(unique(t_state_raw[looks_hybrid & !matches_pattern]))
message(length(hybrid_labels_found), " distinct hybrid Lin1/Lin2 t_state label(s) matched and will be split: ",
        paste(hybrid_labels_found, collapse = ", "))
if (length(unmatched_hybrid_labels) > 0)
  warning(length(unmatched_hybrid_labels), " t_state value(s) look hybrid (mention both Lin1 and Lin2) but ",
          "did NOT match hybrid_pattern, so they were left UNCHANGED rather than guessed at: ",
          paste(unmatched_hybrid_labels, collapse = ", "),
          " -- widen hybrid_pattern above if these should also be split.")

hybrid_cells <- rownames(meta0)[matches_pattern]

if (length(hybrid_cells) == 0) {
  message("No cells matched hybrid_pattern -- seu left unchanged.")
} else {
  lin1_label <- sub(hybrid_pattern, "\\1", t_state_raw[matches_pattern])
  lin2_label <- sub(hybrid_pattern, "\\2", t_state_raw[matches_pattern])

  # Copy 1: the ORIGINAL cells, kept under their original barcodes, t_state
  # overwritten to the Lin1 half of the hybrid label.
  seu_lin1 <- seu[, hybrid_cells]
  seu_lin1$t_state <- lin1_label

  # Copy 2: NEW duplicate cells (unique barcodes via RenameCells(), so Seurat
  # doesn't collide them with copy 1), t_state set to the Lin2 half.
  seu_lin2 <- seu[, hybrid_cells]
  seu_lin2 <- RenameCells(seu_lin2, new.names = paste0(colnames(seu_lin2), "_dupLin2"))
  seu_lin2$t_state <- lin2_label

  # Every other cell 
  seu_rest <- seu[, setdiff(colnames(seu), hybrid_cells)]

  n_before <- ncol(seu)
  seu <- merge(seu_rest, list(seu_lin1, seu_lin2))
  seu <- JoinLayers(seu)

  message("Resolved ", length(hybrid_cells), " hybrid cell(s) into ",
          2 * length(hybrid_cells), " singly-labelled cells (",
          length(hybrid_cells), " kept under their original barcode, now Lin1; ",
          length(hybrid_cells), " newly created duplicate(s), now Lin2). seu: ",
          n_before, " -> ", ncol(seu), " cells (expect +", length(hybrid_cells), ").")
  stopifnot("seu cell count after hybrid-label duplication doesn't match expectation" =
              ncol(seu) == n_before + length(hybrid_cells))
}
rm(meta0, t_state_raw, looks_hybrid, matches_pattern)

cfg <- list(
  cap_ec       = 150,   
  cap_other    = 500,   
  min_ec       = 20,    
  min_t        = 20,   
  min_cells_cc = 10,   
  nboot        = 100,
  seed         = 1,
  alpha        = 0.05,
  cores        = max(1L, parallel::detectCores(logical = FALSE) * 0.4),
  cache_dir     = dpath("cellchat_cache"),
  cache_version = 1,
  use_cache     = TRUE,
  ec_label      = "Endothelial"
)
dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## 0b. T-CELL/EC GROUPING COLUMNS
## ---------------------------------------------------------------------------
clean <- function(x) gsub("[^A-Za-z0-9_]", "_", as.character(x))
meta  <- seu@meta.data
b     <- meta$celltype_broad

#   cc_group     : EC aggregated;   T cells split by state; others broad
#   cc_group_sub : EC split by subtype (used for pooled, subtype-level runs)
seu$cc_group <- clean(case_when(
  b == cfg$ec_label ~ cfg$ec_label,
  b == "T"           ~ meta$t_state,
  TRUE               ~ as.character(b)))
seu$cc_group_sub <- clean(case_when(
  b == cfg$ec_label ~ meta$ec_subtype,
  b == "T"           ~ meta$t_state,
  TRUE               ~ as.character(b)))
meta <- seu@meta.data

ec_subtype_labels   <- clean(unique(na.omit(meta$ec_subtype)))
is_ec_subtype_group <- function(g) g %in% ec_subtype_labels

is_cd4_subtype_group <- function(g) grepl("^CD4", g) | grepl("^Treg", g)
is_cd8_subtype_group <- function(g) grepl("^Lin", g)

is_T_group <- function(g) {
  is_cd4_subtype_group(g) | is_cd8_subtype_group(g) |
    grepl("^CD[48]", g) | g %in% c("NKT", "DoubleNeg", "NK", "T")
}

all_groups_check <- sort(unique(seu$cc_group))
message("is_T_group() split -- VERIFY this before trusting anything downstream:")
message("  Matched as T:     ", paste(all_groups_check[is_T_group(all_groups_check)], collapse = ", "))
message("  NOT matched as T: ", paste(all_groups_check[!is_T_group(all_groups_check)], collapse = ", "))

sender_types  <- c(cfg$ec_label, "Myeloid", "Malignant", "Epithelial", "Fibroblast")
db_int <- CellChatDB.human$interaction
poi    <- db_int$interaction_name[grepl("PDCD1|TIGIT|CD226|PVRIG", db_int$interaction_name)]
print(poi)
targets_main <- c("CD274_PDCD1", "PDCD1LG2_PDCD1", "NECTIN2_TIGIT", "PVR_TIGIT")

## Positive-control interactions
.lig_col0 <- if ("ligand.symbol"   %in% names(db_int)) "ligand.symbol"   else "ligand"
.rec_col0 <- if ("receptor.symbol" %in% names(db_int)) "receptor.symbol" else "receptor"
pos_ctrl_cd4 <- db_int$interaction_name[
  grepl("^HLA-DR", db_int[[.lig_col0]]) & grepl("^CD4$", db_int[[.rec_col0]])]
pos_ctrl_cd8 <- db_int$interaction_name[
  grepl("^HLA-A$", db_int[[.lig_col0]]) & grepl("^CD8", db_int[[.rec_col0]])]
pos_ctrl_pairs <- c(pos_ctrl_cd4, pos_ctrl_cd8)
print(list(pos_ctrl_cd4 = pos_ctrl_cd4, pos_ctrl_cd8 = pos_ctrl_cd8))
if (length(pos_ctrl_cd4) == 0)
  warning("No interaction matched ligand ~ 'HLA-DR*' & receptor == 'CD4' in db_int -- ",
          "run db_int[grepl('HLA-D', db_int[[.lig_col0]]), ] to find the actual spelling ",
          "and set pos_ctrl_cd4 by hand.")
if (length(pos_ctrl_cd8) == 0)
  warning("No interaction matched ligand == 'HLA-A' & receptor ~ 'CD8*' in db_int -- ",
          "run db_int[grepl('^HLA-A', db_int[[.lig_col0]]), ] to find the actual spelling ",
          "and set pos_ctrl_cd8 by hand.")

## ---------------------------------------------------------------------------
## 1. HELPERS
## ---------------------------------------------------------------------------
# Cap cells per sample x broad cell type (ECs get cfg$cap_ec, others cfg$cap_other)
cap_cells <- function(meta, cap_ec, cap_other, seed, ec_label = cfg$ec_label) {
  set.seed(seed)
  meta %>%
    rownames_to_column("cell") %>%
    group_by(sample, celltype_broad) %>%
    filter(row_number() %in% sample.int(
      n(), min(if (first(celltype_broad) == ec_label) cap_ec else cap_other, n()))) %>%
    ungroup() %>%
    pull(cell)
}

# One CellChat run. `type`/`trim` control how group-average expression is computed
run_cc <- function(obj, group_col = "cc_group", type = "truncatedMean", trim = 0.1,
                   min_cells = cfg$min_cells_cc, nboot = cfg$nboot, seed = cfg$seed,
                   de_filter = TRUE, db = CellChatDB.human) {
  stopifnot("sample" %in% colnames(obj@meta.data))
  n_samp <- n_distinct(obj@meta.data$sample)
  if (n_samp != 1) stop("run_cc() received ", n_samp, " samples -- CellChat must be run ",
                        "per sample only (no pooling). Subset to one sample before calling.")
  m <- obj@meta.data
  keep_groups <- names(which(table(m[[group_col]]) >= min_cells))
  cells <- rownames(m)[m[[group_col]] %in% keep_groups]
  obj <- obj[, cells]
  m   <- obj@meta.data
  m[[group_col]] <- factor(m[[group_col]])
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  cc <- createCellChat(object = data.input, meta = m, group.by = group_col)
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, thresh.p = if (de_filter) 0.05 else 1)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- suppressMessages(computeCommunProb(
    cc, type = type, trim = trim, population.size = FALSE,
    nboot = nboot, seed.use = seed))
  filterCommunication(cc, min.cells = min_cells)
}

# Group sizes for ONE sample's cells against min_cells
check_group_sizes <- function(obj, group_col, min_cells = cfg$min_cells_cc) {
  m <- obj@meta.data
  tibble(group = names(table(m[[group_col]])), n_cells = as.integer(table(m[[group_col]]))) %>%
    mutate(kept = n_cells >= min_cells,
           reason = if_else(kept, NA_character_,
                            paste0("< min_cells (", min_cells, ")")))
}

# Run CellChat
run_cc_one_sample <- function(obj, group_col = "cc_group", min_cells = cfg$min_cells_cc,
                              min_groups = 2, ...) {
  stopifnot("sample" %in% colnames(obj@meta.data))
  s <- unique(obj@meta.data$sample)
  if (length(s) != 1) stop("run_cc_one_sample() received ", length(s),
                           " samples -- pass exactly one sample's cells.")
  sizes <- check_group_sizes(obj, group_col, min_cells) %>% mutate(sample = s, group_col = group_col)
  n_ok  <- sum(sizes$kept)
  if (n_ok < min_groups) {
    message("[", s, "] skipped entirely: only ", n_ok, " group(s) meet min_cells (",
            min_cells, ") for '", group_col, "' -- need >= ", min_groups)
    return(list(sample = s, cc = NULL, flagged = sizes %>% mutate(run_attempted = FALSE)))
  }
  cc <- tryCatch(run_cc(obj, group_col = group_col, min_cells = min_cells, ...),
                error = function(e) { message("[", s, "] run_cc() error: ", conditionMessage(e)); NULL })
  list(sample = s, cc = cc, flagged = sizes %>% mutate(run_attempted = TRUE))
}

cache_path <- function(sample, group_col, cap_label, type, trim, min_cells) {
  key <- paste(sample, group_col, cap_label, type, trim, min_cells, sep = "__")
  file.path(cfg$cache_dir, paste0(clean(key), ".rds"))
}

run_cc_one_sample_cached <- function(obj, group_col = "cc_group", cap_label = "na",
                                     type = "truncatedMean", trim = 0.1,
                                     min_cells = cfg$min_cells_cc, ...) {
  s <- unique(obj@meta.data$sample)
  stopifnot(length(s) == 1)
  f <- cache_path(s, group_col, cap_label, type, trim, min_cells)

  if (isTRUE(cfg$use_cache) && file.exists(f)) {
    cached <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$cache_version, cfg$cache_version)) {
      message("[", s, "] (", cap_label, "/", group_col, ") loaded from cache: ", f)
      return(cached$result)
    }
    if (!is.null(cached))
      message("[", s, "] cache at ", f, " is from cache_version ", cached$cache_version,
              " (current: ", cfg$cache_version, ") -- re-running.")
  }

  r <- run_cc_one_sample(obj, group_col = group_col, min_cells = min_cells,
                        type = type, trim = trim, ...)
  tryCatch(
    saveRDS(list(cache_version = cfg$cache_version, result = r), f),
    error = function(e) warning("[", s, "] could not write cache file ", f, ": ", conditionMessage(e))
  )
  r
}

# table of probability and p-value for chosen sources x targets x pairs.
# sources = NULL -> all non-T groups; targets = NULL -> all T groups.
extract_lr <- function(cc, sources = NULL, targets = NULL, pairs = NULL) {
  P  <- cc@net$prob;  PV <- cc@net$pval
  g  <- dimnames(P)[[1]]; lr_all <- dimnames(P)[[3]]
  s  <- if (is.null(sources)) g[!is_T_group(g)] else intersect(sources, g)
  t  <- if (is.null(targets)) g[is_T_group(g)]  else intersect(targets, g)
  lr <- if (is.null(pairs)) lr_all else intersect(pairs, lr_all)
  if (!length(s) || !length(t) || !length(lr)) return(NULL)
  out <- expand.grid(source = s, target = t, interaction = lr, stringsAsFactors = FALSE)
  ix  <- cbind(match(out$source, g), match(out$target, g), match(out$interaction, lr_all))
  out$prob <- P[ix]; out$pval <- PV[ix]
  as_tibble(out)
}

## ---------------------------------------------------------------------------
## 2. CAP CELLS (used everywhere below)
## ---------------------------------------------------------------------------
keep_cells <- cap_cells(meta, cfg$cap_ec, cfg$cap_other, cfg$seed)
seu_c  <- seu[, keep_cells]
meta_c <- seu_c@meta.data
seu_c$stratum <- paste(meta_c$histology, meta_c$platform, sep = "_")
meta_c <- seu_c@meta.data

sample_meta <- meta_c %>% distinct(sample, dataset, histology, platform)

seu_full <- seu
seu_full$stratum <- paste(seu_full@meta.data$histology, seu_full@meta.data$platform, sep = "_")

## ---------------------------------------------------------------------------
## 3. (A) ROBUSTNESS TO EXPRESSION `type`, PER SAMPLE (no pooling)
##    Question: do PD-1 / TIGIT calls survive when sparse genes are not zeroed?
## ---------------------------------------------------------------------------
type_grid <- tribble(
  ~label,     ~type,           ~trim,
  "triMean",  "triMean",       0.10,
  "trunc10",  "truncatedMean", 0.10,
  "trunc05",  "truncatedMean", 0.05)

primary_type_label <- "trunc10"
primary_row <- type_grid %>% filter(label == primary_type_label)
stopifnot("primary_type_label must match exactly one row of type_grid" = nrow(primary_row) == 1)

run_type_grid <- function(seu_input, cap_label_prefix, grid = type_grid, pairs = poi) {
  meta_input <- seu_input@meta.data
  raw <- pmap(grid, function(label, type, trim) {
    pmap_par(unique(meta_input$sample), function(s) {
      cells <- rownames(meta_input)[meta_input$sample == s]
      r <- run_cc_one_sample_cached(seu_input[, cells], group_col = "cc_group_sub",
                                    cap_label = cap_label_prefix, type = type, trim = trim)
      list(label = label, sample = s, cc = r$cc, flagged = r$flagged %>% mutate(type_label = label))
    })
  }) %>% flatten()

  flagged <- map_dfr(raw, "flagged") %>% mutate(cap_label = cap_label_prefix)

  tab <- map_dfr(raw, function(r) {
    if (is.null(r$cc)) return(NULL)
    lr <- extract_lr(r$cc, sources = NULL, targets = NULL, pairs = pairs)
    if (is.null(lr)) return(NULL)
    lr %>% mutate(type_label = r$label, sample = r$sample)
  }) %>%
    left_join(meta_input %>% distinct(sample, dataset, histology, platform), by = "sample") %>%
    mutate(cap_label = cap_label_prefix)

  list(raw = raw, flagged = flagged, table = tab)
}

## Capped arm, full sensitivity grid 
tg_capped          <- run_type_grid(seu_c, "typecheck", type_grid)
type_check_raw     <- tg_capped$raw
flagged_typecheck  <- tg_capped$flagged
type_check         <- tg_capped$table

type_summary <- type_check %>%
  filter(is_ec_subtype_group(source)) %>%
  group_by(histology, platform, type_label, interaction) %>%
  summarise(n_samples = n_distinct(sample), max_prob = max(prob),
            n_sig_pairs = sum(pval < cfg$alpha & prob > 0), .groups = "drop")

## Baseline arm
tg_baseline           <- run_type_grid(seu_full, "baseline_typecheck", primary_row)
flagged_typecheck_bl  <- tg_baseline$flagged
type_check_baseline   <- tg_baseline$table
print(type_summary %>% filter(interaction %in% targets_main), n = 100)

## Positive-control 
tg_baseline_posctrl         <- run_type_grid(seu_full, "baseline_typecheck", primary_row, pairs = pos_ctrl_pairs)
type_check_posctrl_baseline <- tg_baseline_posctrl$table

message(sum(!flagged_typecheck$kept), " group-in-sample exclusions across the type-robustness grid; ",
        "see flagged_typecheck for exactly which groups/samples and why.")

## ---------------------------------------------------------------------------
## 4. (B) PER-SAMPLE CELLCHAT WITH BENCHMARK SENDERS 
## ---------------------------------------------------------------------------
type_use <- primary_row$type; trim_use <- primary_row$trim

run_main_pipeline <- function(seu_input, cap_label, type_use = "truncatedMean",
                              trim_use = 0.10, min_cells = cfg$min_cells_cc) {
  meta_input   <- seu_input@meta.data
  sample_meta_i <- meta_input %>% distinct(sample, dataset, histology, platform)
  samples <- unique(meta_input$sample)

  per_sample_raw <- pmap_par(samples, function(s) {
    cells_s <- rownames(meta_input)[meta_input$sample == s]
    run_cc_one_sample_cached(seu_input[, cells_s], group_col = "cc_group", cap_label = cap_label,
                            min_cells = min_cells, type = type_use, trim = trim_use)
  })

  flagged_main <- map_dfr(per_sample_raw, "flagged") %>% mutate(cap_label = !!cap_label)

  per_sample <- map(per_sample_raw, function(r) {
    if (is.null(r$cc)) return(NULL)
    main <- extract_lr(r$cc, sources = sender_types, pairs = poi)            # benchmark senders
    ctrl <- extract_lr(r$cc, sources = cfg$ec_label, pairs = NULL)           # all L-R, EC -> T
    if (is.null(main) || is.null(ctrl)) return(NULL)
    list(sample = r$sample, main = main,
         ctrl = ctrl %>% filter(prob > 0, pval < cfg$alpha) %>% distinct(interaction),
         ctrl_full = ctrl %>% group_by(interaction) %>%
           summarise(prob_max = max(prob), sig = any(pval < cfg$alpha & prob > 0), .groups = "drop"),
         ctrl_raw = ctrl %>% select(target, interaction, prob, pval))
  })
  per_sample <- compact(per_sample)
  saveRDS(per_sample, dpath(paste0("cellchat_per_sample_", cap_label, ".rds")))

  eligible <- tibble(sample = map_chr(per_sample, "sample")) %>% left_join(sample_meta_i, by = "sample")
  n_elig   <- nrow(eligible)
  n_skipped <- length(samples) - n_elig
  message("[", cap_label, "] ", n_elig, " / ", length(samples), " samples produced a CellChat run; ",
          n_skipped, " skipped entirely (fewer than 2 groups meeting min_cells = ", min_cells,
          " for cc_group) -- see the returned $flagged for the full per-group breakdown.")
  print(table(eligible$histology, eligible$platform))

  flagged_main_summary <- flagged_main %>%
    group_by(sample) %>%
    summarise(run_attempted = any(run_attempted),
              n_groups_present = n(), n_groups_kept = sum(kept),
              groups_dropped = paste(group[!kept], collapse = ", "), .groups = "drop") %>%
    left_join(sample_meta_i, by = "sample") %>%
    mutate(cap_label = !!cap_label) %>%
    arrange(run_attempted, n_groups_kept)
  print(flagged_main_summary, n = 200)

  main_res <- map_dfr(per_sample, ~ mutate(.x$main, sample = .x$sample)) %>%
    left_join(sample_meta_i, by = "sample") %>% mutate(cap_label = !!cap_label)

  bench <- main_res %>%
    group_by(cap_label, sample, dataset, histology, platform, source, interaction) %>%
    summarise(prob_max = max(prob), sig = any(pval < cfg$alpha & prob > 0), .groups = "drop")

  bench_tab <- bench %>%
    filter(interaction %in% targets_main) %>%
    group_by(cap_label, histology, source, interaction) %>%
    summarise(n = n(), frac_sig = mean(sig), median_prob = median(prob_max), .groups = "drop")
  print(bench_tab, n = 100)

  list(cap_label = cap_label, per_sample = per_sample, flagged_main = flagged_main,
       flagged_main_summary = flagged_main_summary, eligible = eligible,
       main_res = main_res, bench = bench, bench_tab = bench_tab)
}

## ---- Arm 1: capped / downsampled ----
res_capped <- run_main_pipeline(seu_c, cap_label = "capped")
list2env(setNames(res_capped[c("per_sample", "flagged_main", "eligible", "main_res", "bench", "bench_tab")],
                  c("per_sample", "flagged_main", "eligible", "main_res", "bench", "bench_tab")),
         envir = environment())  

n_elig <- nrow(eligible)

p_bench <- res_capped$bench %>%
  filter(interaction %in% targets_main) %>%
  ggplot(aes(source, prob_max, fill = histology)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ interaction, scales = "free_y") +
  labs(x = "Sender (-> T cells)", y = "Max communication probability per sample") +
  theme_bw()
ggsave(rpath("diagnostics", "cellchat_benchmark_senders_capped.jpg"), p_bench, width = 8, height = 4)

## ---------------------------------------------------------------------------
## 4b. BASELINE ARM: same pipeline, NO downsampling 
## ---------------------------------------------------------------------------
res_baseline <- run_main_pipeline(seu_full, cap_label = "baseline_uncapped")

p_bench_baseline <- res_baseline$bench %>%
  filter(interaction %in% targets_main) %>%
  ggplot(aes(source, prob_max, fill = histology)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ interaction, scales = "free_y") +
  labs(x = "Sender (-> T cells)", y = "Max communication probability per sample") +
  theme_bw()
ggsave(rpath("diagnostics", "cellchat_benchmark_senders_baseline.jpg"), p_bench_baseline, width = 8, height = 4)

saveRDS(list(capped = res_capped, baseline_uncapped = res_baseline), dpath("cellchat_flagged_exclusions.rds"))

## ---------------------------------------------------------------------------
## 4c. CAPPED vs BASELINE COMPARISON
## ---------------------------------------------------------------------------
bench_compare <- bind_rows(res_capped$bench, res_baseline$bench) %>%
  filter(interaction %in% targets_main)

samples_capped   <- unique(res_capped$eligible$sample)
samples_baseline <- unique(res_baseline$eligible$sample)
message("Samples eligible in capped only: ",  paste(setdiff(samples_capped, samples_baseline), collapse = ", "))
message("Samples eligible in baseline only: ", paste(setdiff(samples_baseline, samples_capped), collapse = ", "))

call_compare <- bench_compare %>%
  select(cap_label, sample, source, interaction, sig, prob_max) %>%
  tidyr::pivot_wider(names_from = cap_label, values_from = c(sig, prob_max)) %>%
  filter(!is.na(sig_capped), !is.na(sig_baseline_uncapped)) %>%       # only samples present in BOTH arms
  mutate(call_flip = sig_capped != sig_baseline_uncapped)

message(sum(call_compare$call_flip, na.rm = TRUE), " / ", nrow(call_compare),
        " sample x sender x interaction significance calls differ between capped and baseline ",
        "(present in both arms) -- see call_compare for exactly which.")
print(call_compare %>% filter(call_flip), n = 100)

p_compare <- ggplot(call_compare, aes(prob_max_capped, prob_max_baseline_uncapped, colour = call_flip)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), name = "Call flipped") +
  facet_wrap(~ interaction) +
  labs(x = "Max probability (capped)", y = "Max probability (baseline, uncapped)") + theme_bw()
ggsave(rpath("diagnostics", "cellchat_capped_vs_baseline.jpg"), p_compare, width = 8, height = 4)

## ---------------------------------------------------------------------------
## 5. (C) EXPRESSION-MATCHED NEGATIVE CONTROLS
##    Compare PD-1 / TIGIT pairs against pairs whose ligand and
##    receptor are detected at similar frequencies (within a factor of 2^tol).
## ---------------------------------------------------------------------------
ctrl_freq <- map_dfr(per_sample, ~ mutate(.x$ctrl, sample = .x$sample)) %>%
  count(interaction, name = "n_sig") %>%
  mutate(freq = n_sig / n_elig)

lig_col <- if ("ligand.symbol"   %in% names(db_int)) "ligand.symbol"   else "ligand"
rec_col <- if ("receptor.symbol" %in% names(db_int)) "receptor.symbol" else "receptor"

split_genes <- function(x) strsplit(x, ",\\s*")
all_genes <- unique(unlist(split_genes(c(db_int[[lig_col]], db_int[[rec_col]]))))
all_genes <- intersect(all_genes[!is.na(all_genes)], rownames(seu_c))

det_frac <- function(obj, cells, genes) {
  m <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes, cells, drop = FALSE]
  setNames(Matrix::rowMeans(m > 0), genes)
}
cplx_frac <- function(sym, fr) {           
  vapply(split_genes(sym), function(g) {
    g <- g[g %in% names(fr)]
    if (!length(g)) NA_real_ else min(fr[g])
  }, numeric(1))
}

ec_cells <- rownames(meta_c)[meta_c$cc_group == cfg$ec_label]
t_cells  <- rownames(meta_c)[is_T_group(meta_c$cc_group)]
fr_ec <- det_frac(seu_c, ec_cells, all_genes)
fr_t  <- det_frac(seu_c, t_cells,  all_genes)

lr_tab <- tibble(interaction = db_int$interaction_name,
                 lig_fr = cplx_frac(db_int[[lig_col]], fr_ec),
                 rec_fr = cplx_frac(db_int[[rec_col]], fr_t)) %>%
  left_join(ctrl_freq, by = "interaction") %>%
  mutate(freq = replace_na(freq, 0))

match_controls <- function(target, lr_tab, tol_log2 = 1, eps = 1e-3, exclude = poi) {
  t <- lr_tab[lr_tab$interaction == target, ]
  lr_tab %>% filter(
    !interaction %in% exclude, !is.na(lig_fr), !is.na(rec_fr),
    abs(log2((lig_fr + eps) / (t$lig_fr + eps))) <= tol_log2,
    abs(log2((rec_fr + eps) / (t$rec_fr + eps))) <= tol_log2)
}

ctrl_out <- map_dfr(intersect(targets_main, lr_tab$interaction), function(tg) {
  ctl <- match_controls(tg, lr_tab)
  tf  <- lr_tab$freq[lr_tab$interaction == tg]
  tibble(target = tg, target_freq = tf, n_controls = nrow(ctl),
         control_median = median(ctl$freq),
         control_q95    = quantile(ctl$freq, 0.95),
         empirical_p    = (1 + sum(ctl$freq >= tf)) / (nrow(ctl) + 1))
})
print(ctrl_out)  

p_ctrl <- map_dfr(intersect(targets_main, lr_tab$interaction), function(tg) {
  match_controls(tg, lr_tab) %>% mutate(target = tg)
}) %>%
  ggplot(aes(freq)) + geom_histogram(bins = 20) +
  geom_vline(data = ctrl_out, aes(xintercept = target_freq), colour = "red") +
  facet_wrap(~ target, scales = "free_y") +
  labs(x = "Fraction of samples with significant EC -> T call",
       y = "Expression-matched control L-R pairs") + theme_bw()
ggsave(rpath("diagnostics", "cellchat_negative_controls.jpg"), p_ctrl, width = 9, height = 4)

## ---------------------------------------------------------------------------
## 6. (A) STRENGTH OF EC -> T vs. (i) OTHER CELL TYPES, (ii) MATCHED NEGATIVE
## ---------------------------------------------------------------------------

## (i) vs. other broad cell types: EC aggregated, not split by subtype
p_vs_celltypes_baseline <- res_baseline$bench %>%
  filter(interaction %in% targets_main) %>%
  ggplot(aes(x = source, y = prob_max, fill = source)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ interaction, scales = "free_y") +
  labs(x = "Sender (-> T cells)", y = "Max communication probability (per sample)",
       subtitle = "Baseline (PAPER)") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(rpath("diagnostics", "cellchat_fig_vs_other_celltypes_baseline.jpg"), p_vs_celltypes_baseline, width = 9, height = 4)

flip_marks <- call_compare %>% filter(call_flip) %>%
  distinct(sample, source, interaction) %>% mutate(flipped = TRUE)

bench_capped_marked <- res_capped$bench %>%
  filter(interaction %in% targets_main) %>%
  left_join(flip_marks, by = c("sample", "source", "interaction")) %>%
  mutate(flipped = coalesce(flipped, FALSE))

p_vs_celltypes_capped <- ggplot(bench_capped_marked, aes(x = source, y = prob_max, fill = source)) +
  geom_boxplot(outlier.size = 0.6) +
  geom_point(data = filter(bench_capped_marked, flipped), shape = 8, size = 2, colour = "firebrick") +
  facet_wrap(~ interaction, scales = "free_y") +
  labs(x = "Sender (-> T cells)", y = "Max communication probability (per sample)",
       subtitle = "Capped (SUPPLEMENTARY; red * = significance call flipped vs baseline)") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(rpath("diagnostics", "cellchat_fig_vs_other_celltypes_capped_supp.jpg"), p_vs_celltypes_capped, width = 9, height = 4)

## (ii) vs. matched negative controls, on MAGNITUDE 
build_neg_ctrl_mag <- function(per_sample_list, main_res_tbl) {
  ctrl_full_df <- map_dfr(per_sample_list, ~ mutate(.x$ctrl_full, sample = .x$sample))
  neg_ctrl_mag <- map_dfr(intersect(targets_main, lr_tab$interaction), function(tg) {
    ctl_pairs <- match_controls(tg, lr_tab)$interaction
    ctrl_full_df %>% filter(interaction %in% ctl_pairs) %>%
      transmute(sample, target = tg, prob_max, role = "Matched negative controls")
  })
  target_mag <- main_res_tbl %>%
    filter(source == cfg$ec_label, interaction %in% targets_main) %>%
    transmute(sample, target = interaction, prob_max = prob, role = "Target pathway (EC -> T)")
  bind_rows(target_mag, neg_ctrl_mag)
}

wilcox_by_target <- function(mag_tbl) {
  map_dfr(unique(mag_tbl$target), function(tg) {
    d  <- filter(mag_tbl, target == tg)
    tp <- d$prob_max[d$role == "Target pathway (EC -> T)"]
    nc <- d$prob_max[d$role == "Matched negative controls"]
    if (!length(tp) || !length(nc)) return(NULL)
    wt <- wilcox.test(tp, nc, alternative = "greater")
    tibble(target = tg, n_target_samples = length(tp), n_control_obs = length(nc),
           median_target = median(tp), median_control = median(nc), p_wilcox = wt$p.value)
  })
}

mag_compare_baseline <- build_neg_ctrl_mag(res_baseline$per_sample, res_baseline$main_res)
mag_compare_capped   <- build_neg_ctrl_mag(res_capped$per_sample,   res_capped$main_res)

p_neg_mag_baseline <- ggplot(mag_compare_baseline, aes(x = role, y = prob_max, fill = role)) +
  geom_boxplot(outlier.size = 0.6) + facet_wrap(~ target, scales = "free_y") +
  labs(x = NULL, y = "Communication probability (EC -> T, per sample)", subtitle = "Baseline (PAPER)") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(rpath("diagnostics", "cellchat_fig_vs_neg_controls_magnitude_baseline.jpg"), p_neg_mag_baseline, width = 9, height = 4)

p_neg_mag_capped <- ggplot(mag_compare_capped, aes(x = role, y = prob_max, fill = role)) +
  geom_boxplot(outlier.size = 0.6) + facet_wrap(~ target, scales = "free_y") +
  labs(x = NULL, y = "Communication probability (EC -> T, per sample)", subtitle = "Capped (SUPPLEMENTARY)") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(rpath("diagnostics", "cellchat_fig_vs_neg_controls_magnitude_capped_supp.jpg"), p_neg_mag_capped, width = 9, height = 4)

wilcox_results <- bind_rows(
  mutate(wilcox_by_target(mag_compare_baseline), cap_label = "baseline_uncapped"),
  mutate(wilcox_by_target(mag_compare_capped),   cap_label = "capped")
) %>% arrange(target, cap_label)
print(wilcox_results)

## ---------------------------------------------------------------------------
## 6. (B) T-CELL-SPLIT-RESOLVED BENCHMARK, NEGATIVE-CONTROL, AND DOT-PLOT TABLES
## ---------------------------------------------------------------------------

## (i) sender_types (EC, Myeloid, Malignant, Epithelial, Fibroblast) vs CD4-
## only or CD8-only T-cell states
build_bench_by_tcell <- function(main_res_tbl, target_filter) {
  main_res_tbl %>%
    filter(target_filter(target)) %>%
    group_by(cap_label, sample, dataset, histology, platform, source, interaction) %>%
    summarise(prob_max = max(prob), sig = any(pval < cfg$alpha & prob > 0), .groups = "drop")
}

bench_baseline_cd4 <- build_bench_by_tcell(res_baseline$main_res, is_cd4_subtype_group)
bench_baseline_cd8 <- build_bench_by_tcell(res_baseline$main_res, is_cd8_subtype_group)
bench_capped_cd4   <- build_bench_by_tcell(res_capped$main_res,   is_cd4_subtype_group)
bench_capped_cd8   <- build_bench_by_tcell(res_capped$main_res,   is_cd8_subtype_group)

build_bench_ec_aggregate_from_subtypes <- function(type_check_tbl, target_filter) {
  type_check_tbl %>%
    filter(is_ec_subtype_group(source), interaction %in% targets_main, target_filter(target)) %>%
    group_by(sample, dataset, histology, platform, interaction) %>%
    summarise(prob_max = max(prob), sig = any(pval < cfg$alpha & prob > 0), .groups = "drop") %>%
    mutate(source = cfg$ec_label, cap_label = "baseline_uncapped")
}

bench_ec_aggregate_cd4 <- build_bench_ec_aggregate_from_subtypes(type_check_baseline, is_cd4_subtype_group)
bench_ec_aggregate_cd8 <- build_bench_ec_aggregate_from_subtypes(type_check_baseline, is_cd8_subtype_group)

bench_baseline_cd4_ecagg <- bind_rows(filter(bench_baseline_cd4, source != cfg$ec_label), bench_ec_aggregate_cd4)
bench_baseline_cd8_ecagg <- bind_rows(filter(bench_baseline_cd8, source != cfg$ec_label), bench_ec_aggregate_cd8)

mark_flips_by_tcell <- function(bench_capped_tbl, bench_baseline_tbl) {
  compare <- bind_rows(bench_capped_tbl, bench_baseline_tbl) %>%
    select(cap_label, sample, source, interaction, sig, prob_max) %>%
    tidyr::pivot_wider(names_from = cap_label, values_from = c(sig, prob_max)) %>%
    filter(!is.na(sig_capped), !is.na(sig_baseline_uncapped)) %>%
    mutate(call_flip = sig_capped != sig_baseline_uncapped)
  flips <- compare %>% filter(call_flip) %>% distinct(sample, source, interaction) %>%
    mutate(flipped = TRUE)
  bench_capped_tbl %>% left_join(flips, by = c("sample", "source", "interaction")) %>%
    mutate(flipped = coalesce(flipped, FALSE))
}
bench_capped_cd4_marked <- mark_flips_by_tcell(bench_capped_cd4, bench_baseline_cd4)
bench_capped_cd8_marked <- mark_flips_by_tcell(bench_capped_cd8, bench_baseline_cd8)

## (ii) matched-negative-control MAGNITUDE, T-split-specific
build_neg_ctrl_mag_by_tcell <- function(per_sample_list, target_filter) {
  target_mag <- map_dfr(per_sample_list, function(ps) {
    ps$main %>% filter(source == cfg$ec_label, interaction %in% targets_main,
                       target_filter(target)) %>%
      transmute(sample = ps$sample, target = interaction, prob_max = prob,
                role = "Target pathway (EC -> T)")
  })
  neg_ctrl_mag <- map_dfr(intersect(targets_main, lr_tab$interaction), function(tg) {
    ctl_pairs <- match_controls(tg, lr_tab)$interaction
    map_dfr(per_sample_list, function(ps) {
      ps$ctrl_raw %>% filter(interaction %in% ctl_pairs, target_filter(target)) %>%
        transmute(sample = ps$sample, target = tg, prob_max = prob,
                  role = "Matched negative controls")
    })
  })
  bind_rows(target_mag, neg_ctrl_mag)
}

mag_compare_baseline_cd4 <- build_neg_ctrl_mag_by_tcell(res_baseline$per_sample, is_cd4_subtype_group)
mag_compare_baseline_cd8 <- build_neg_ctrl_mag_by_tcell(res_baseline$per_sample, is_cd8_subtype_group)
mag_compare_capped_cd4   <- build_neg_ctrl_mag_by_tcell(res_capped$per_sample,   is_cd4_subtype_group)
mag_compare_capped_cd8   <- build_neg_ctrl_mag_by_tcell(res_capped$per_sample,   is_cd8_subtype_group)

wilcox_results_cd4 <- bind_rows(
  mutate(wilcox_by_target(mag_compare_baseline_cd4), cap_label = "baseline_uncapped", t_split = "CD4"),
  mutate(wilcox_by_target(mag_compare_capped_cd4),   cap_label = "capped",            t_split = "CD4"))
wilcox_results_cd8 <- bind_rows(
  mutate(wilcox_by_target(mag_compare_baseline_cd8), cap_label = "baseline_uncapped", t_split = "CD8"),
  mutate(wilcox_by_target(mag_compare_capped_cd8),   cap_label = "capped",            t_split = "CD8"))
print(bind_rows(wilcox_results_cd4, wilcox_results_cd8) %>% arrange(target, t_split, cap_label), n = 100)

## (iii) EC-subtype x fine-T-subtype dot-plot table
build_dot_tab_fine <- function(tbl, target_filter, sample_meta_grp, interactions = targets_main) {
  tbl %>%
    filter(is_ec_subtype_group(source),
           interaction %in% interactions,
           target_filter(target)) %>%
    left_join(sample_meta_grp %>% select(sample, group), by = "sample") %>%
    filter(!is.na(group)) %>%
    group_by(source, target, interaction, group) %>%
    summarise(n_samples = n_distinct(sample),
              n_sig = n_distinct(sample[pval < cfg$alpha & prob > 0]),
              frac_sig = n_sig / n_samples,
              median_prob_sig = if (any(pval < cfg$alpha & prob > 0))
                median(prob[pval < cfg$alpha & prob > 0]) else NA_real_,
              .groups = "drop")
}

dot_tab_cd4_baseline <- build_dot_tab_fine(type_check_baseline, is_cd4_subtype_group, sample_meta_grp)
dot_tab_cd8_baseline <- build_dot_tab_fine(type_check_baseline, is_cd8_subtype_group, sample_meta_grp)

## Positive-control dot tables (figures 3C/4C/5C)
dot_tab_cd4_posctrl_baseline <- build_dot_tab_fine(type_check_posctrl_baseline, is_cd4_subtype_group,
                                                   sample_meta_grp, interactions = pos_ctrl_cd4)
dot_tab_cd8_posctrl_baseline <- build_dot_tab_fine(type_check_posctrl_baseline, is_cd8_subtype_group,
                                                   sample_meta_grp, interactions = pos_ctrl_cd8)

## (iv) matched negative-control MAGNITUDE, EC-SUBTYPE resolved 
ctrl_pairs_all <- unique(unlist(map(intersect(targets_main, lr_tab$interaction),
                                    function(tg) match_controls(tg, lr_tab)$interaction)))
tg_baseline_negctrl         <- run_type_grid(seu_full, "baseline_typecheck", primary_row, pairs = ctrl_pairs_all)
type_check_negctrl_baseline <- tg_baseline_negctrl$table

build_neg_ctrl_mag_by_tcell_ecsub <- function(target_tbl, negctrl_tbl, target_filter) {
  target_mag <- target_tbl %>%
    filter(is_ec_subtype_group(source), interaction %in% targets_main, target_filter(target)) %>%
    transmute(sample, source, target = interaction, prob_max = prob,
             role = "Target pathway (EC subtype -> T)")

  neg_ctrl_mag <- map_dfr(intersect(targets_main, lr_tab$interaction), function(tg) {
    ctl_pairs <- match_controls(tg, lr_tab)$interaction
    negctrl_tbl %>%
      filter(is_ec_subtype_group(source), interaction %in% ctl_pairs, target_filter(target)) %>%
      transmute(sample, source, target = tg, prob_max = prob,
               role = "Matched negative controls")
  })

  bind_rows(target_mag, neg_ctrl_mag)
}

mag_compare_baseline_cd4_ecsub <- build_neg_ctrl_mag_by_tcell_ecsub(type_check_baseline, type_check_negctrl_baseline, is_cd4_subtype_group)
mag_compare_baseline_cd8_ecsub <- build_neg_ctrl_mag_by_tcell_ecsub(type_check_baseline, type_check_negctrl_baseline, is_cd8_subtype_group)

## ---------------------------------------------------------------------------
## 9. SAVE OUTPUTS
## ---------------------------------------------------------------------------
fig_data <- list(
  cfg            = cfg,
  targets_main   = targets_main,
  t_split_names  = names(t_splits),

  prop_tabs_baseline      = prop_tabs_baseline,
  prop_tabs_capped_marked = prop_tabs_capped_marked,
  prop_compare            = prop_compare,
  sample_meta_grp         = sample_meta_grp,

  bench_baseline      = res_baseline$bench,
  bench_capped_marked = bench_capped_marked,
  flip_marks          = flip_marks,

  mag_compare_baseline = mag_compare_baseline,
  mag_compare_capped   = mag_compare_capped,
  wilcox_results       = wilcox_results,

  ctrl_out      = ctrl_out,
  lr_tab        = lr_tab,
  type_summary  = type_summary,
  bench_compare = bench_compare,
  call_compare  = call_compare,

  hp_models    = hp_models,
  perm_results = perm_results,

  bench_baseline_cd4       = bench_baseline_cd4,
  bench_baseline_cd8       = bench_baseline_cd8,
  bench_capped_cd4_marked  = bench_capped_cd4_marked,
  bench_capped_cd8_marked  = bench_capped_cd8_marked,
  mag_compare_baseline_cd4 = mag_compare_baseline_cd4,
  mag_compare_baseline_cd8 = mag_compare_baseline_cd8,
  mag_compare_capped_cd4   = mag_compare_capped_cd4,
  mag_compare_capped_cd8   = mag_compare_capped_cd8,
  wilcox_results_cd4       = wilcox_results_cd4,
  wilcox_results_cd8       = wilcox_results_cd8,
  dot_tab_cd4_baseline     = dot_tab_cd4_baseline,
  dot_tab_cd8_baseline     = dot_tab_cd8_baseline,

  pos_ctrl_cd4                 = pos_ctrl_cd4,
  pos_ctrl_cd8                 = pos_ctrl_cd8,
  dot_tab_cd4_posctrl_baseline = dot_tab_cd4_posctrl_baseline,
  dot_tab_cd8_posctrl_baseline = dot_tab_cd8_posctrl_baseline,

  mag_compare_baseline_cd4_ecsub = mag_compare_baseline_cd4_ecsub,
  mag_compare_baseline_cd8_ecsub = mag_compare_baseline_cd8_ecsub,

  bench_baseline_cd4_ecagg = bench_baseline_cd4_ecagg,
  bench_baseline_cd8_ecagg = bench_baseline_cd8_ecagg
)
saveRDS(fig_data, dpath("cellchat_figure_data.rds"))
message("Saved every figure-ready table/model to cellchat_figure_data.rds (",
        round(file.size(dpath("cellchat_figure_data.rds")) / 1e6, 1),
        " MB). figures/make_figures.R regenerates the paper figures from this file alone.")