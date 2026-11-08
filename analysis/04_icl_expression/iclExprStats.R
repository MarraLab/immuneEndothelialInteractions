# =============================================================================
# Statistical models for immune-checkpoint-ligand (ICL) expression
# (Figure 2, Supplementary Figure 2C, Supplementary Tables 2-4).
# =============================================================================

source("R/utils.R")
`%||%` <- function(a, b) if (length(a) == 0) b else a

suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans)
  library(edgeR); library(limma); library(speckle)
})
## ---------------------------------------------------------------------------
## 0. INPUT AND SETTINGS
## ---------------------------------------------------------------------------
# Merged Seurat object; "counts" layer must be raw counts. Required metadata columns:
#   sample, dataset, platform ("Flex"/"other"), histology ("LUAD"/"LUSC"),
#   tissue ("tumour"/"normal"), celltype_broad ("Endothelial", "Malignant", ...),
#   ec_subtype (ECs only), nCount_RNA

seu <- qs_read(dpath("nsclc_allCells_publicData_useThis_wNormal.qs2"))

ENDO <- "Endothelial"

cfg <- list(
  genes         = c(
    # "ACTB", "EPCAM", "PTPRC", "PECAM1"
    "CD274", "PDCD1LG2", "NECTIN2", "PVR"
    ),   
  min_ec_sample = 20,     # samples with fewer ECs are excluded from EC-level analyses
  min_cells_pb  = 20,     # minimum cells per pseudobulk profile
  cap_other     = 500,    # max cells per sample x cell type in the cell-type comparison (Fig 2A)
  cap_lusc_ec   = 200,    # max LUSC ECs per sample in each down-sampling draw
  n_iter        = 100,
  cores         = 5,    
  seed          = 1,
  alpha         = 0.05
)

## ---------------------------------------------------------------------------
## 1. HELPERS
## ---------------------------------------------------------------------------

# One row per cell: metadata + detection (det__GENE, count >= 1) and raw count (cnt__GENE).
build_cell_table <- function(seu, genes = cfg$genes, assay = "RNA") {
  miss <- setdiff(genes, rownames(seu))
  if (length(miss)) message("genes not in object (dropped): ", paste(miss, collapse = ", "))
  genes <- intersect(genes, rownames(seu))
  cnt <- GetAssayData(seu, assay = assay, layer = "counts")[genes, , drop = FALSE]
  cm  <- t(as.matrix(cnt))                                        # cells x genes
  md  <- seu@meta.data
  tab <- tibble(
    cell = rownames(md), sample = md$sample, dataset = md$dataset,
    platform  = factor(md$platform,  levels = c("Other", "FFPE")),
    histology = factor(md$histology, levels = c("LUAD", "LUSC")),
    tissue    = factor(md$tissue,    levels = c("normal", "tumour")),
    ctype = md$celltype_broad, ec_subtype = md$ec_subtype, nUMI = md$nCount_RNA)
  if ("is_tuhev" %in% names(md)) tab$is_tuhev <- md$is_tuhev
  tab$log_nUMI <- as.numeric(scale(log10(tab$nUMI)))              # detection depends on depth
  det <- as_tibble(cm > 0); names(det) <- paste0("det__", colnames(cm))
  cn  <- as_tibble(cm);     names(cn)  <- paste0("cnt__", colnames(cm))
  list(cells = bind_cols(tab, det, cn), genes = genes)
}
 
keep_samples_min_ec <- function(d, min_ec = cfg$min_ec_sample) {
  n_ec <- d %>% filter(ctype == ENDO) %>% count(sample, name = "n_ec")
  d %>% filter(sample %in% n_ec$sample[n_ec$n_ec >= min_ec])
}
 
cap_per_sample_type <- function(d, cap = cfg$cap_other, seed = cfg$seed) {
  set.seed(seed)
  d %>% group_by(sample, ctype) %>%
    filter(dplyr::row_number() %in% sample.int(dplyr::n(), min(cap, dplyr::n()))) %>%
    ungroup()
}
 
# glmmTMB with a fallback formula if the first fit fails or the Hessian is not positive definite
fit_glmm <- function(formula, data, family, fallback = NULL) {
  ok_fit <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess)
  formulas <- c(list(formula),
               if (!is.null(fallback)) (if (is.list(fallback)) fallback else list(fallback)))
  last_msg <- NULL
  bfgs_ctrl <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS", maxit = 1e4))
 
  try_one <- function(f, ctrl = glmmTMBControl()) {
    withCallingHandlers(
      tryCatch(glmmTMB(f, data = data, family = family, control = ctrl),
              error = function(e) { last_msg <<- conditionMessage(e); NULL }),
      warning = function(w) { last_msg <<- conditionMessage(w); invokeRestart("muffleWarning") })
  }
 
  for (i in seq_along(formulas)) {
    m <- try_one(formulas[[i]])
    if (ok_fit(m)) {
      attr(m, "fit_rung") <- paste0("formula ", i - 1, " (default optimizer)")
      attr(m, "used_fallback") <- i > 1
      return(m)
    }
    m2 <- try_one(formulas[[i]], bfgs_ctrl)
    if (ok_fit(m2)) {
      attr(m2, "fit_rung") <- paste0("formula ", i - 1, " (bfgs optimizer)")
      attr(m2, "used_fallback") <- i > 1
      return(m2)
    }
  }
  if (!is.null(last_msg)) {
    message("fit_glmm(): fitting failed across every formula/optimizer combination tried -- ", last_msg)
  } else {
    message("fit_glmm(): model(s) fit without error/warning but the Hessian was never ",
            "positive-definite across every formula/optimizer tried -- treat as a failed fit")
  }
  NULL
}
 
# exp(coefficient) with Wald 95% CI: an ODDS RATIO for a binomial model, a RATE RATIO for NB
coef_ratio <- function(m, term, effect = "OR") {
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients$cond
  if (!term %in% rownames(cf)) return(NULL)
  est <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]; z <- qnorm(0.975)
  tibble(effect = effect, term = term, estimate = est, se = se,
         ratio = exp(est), LCL = exp(est - z * se), UCL = exp(est + z * se),
         p = cf[term, "Pr(>|z|)"], used_fallback = isTRUE(attr(m, "used_fallback")),
         fit_rung = attr(m, "fit_rung") %||% "formula 0 (default optimizer)")
}
 
# One gene, one model. response "det" (binomial) or "cnt" (negative binomial)
fit_gene <- function(d, gene, rhs, term, fallback_rhs = NULL,
                     response = c("det", "cnt"), effect = "OR", family = binomial()) {
  response <- match.arg(response)
  f  <- as.formula(paste0(response, "__", gene, " ~ ", rhs))
  ff <- if (is.null(fallback_rhs)) NULL else
    lapply(fallback_rhs, function(r) as.formula(paste0(response, "__", gene, " ~ ", r)))
  r  <- coef_ratio(fit_glmm(f, d, family, ff), term, effect)
  if (is.null(r)) return(NULL)
  r$gene <- gene; r$n_cells <- nrow(d); r$n_samples <- n_distinct(d$sample)
  r
}
 
# Holm within `by` (default: across genes for one comparison)
add_padj <- function(df, by = character(0), method = "holm") {
  if (is.null(df) || !nrow(df) || !all(by %in% names(df))) return(df)
  df %>% group_by(across(all_of(by))) %>% mutate(p_adj = p.adjust(p, method)) %>% ungroup()
}
 
## ---------------------------------------------------------------------------
## 2. COMPOSITION (Fig 1F): speckle::propeller on per-sample EC-subtype proportions
## ---------------------------------------------------------------------------
comp_test <- function(ec, design_formula, plus, minus = NULL,
                      min_ec = cfg$min_ec_sample, transform = "asin") {
  ec <- ec %>% filter(ctype == ENDO, !is.na(ec_subtype_excl)) %>% keep_samples_min_ec(min_ec)
  props <- speckle::getTransformedProps(clusters = ec$ec_subtype_excl, sample = ec$sample,
                                        transform = transform)
  smp <- ec %>% distinct(sample, dataset, platform, histology, tissue) %>% as.data.frame()
  smp <- smp[match(colnames(props$TransformedProps), smp$sample), ]
  design <- model.matrix(design_formula, data = smp)
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    message("dropping non-estimable design columns: ", paste(ne, collapse = ", "))
    design <- design[, setdiff(colnames(design), ne), drop = FALSE]
  }
  contr <- setNames(rep(0, ncol(design)), colnames(design))
  stopifnot(plus %in% names(contr), is.null(minus) || minus %in% names(contr))
  contr[plus] <- 1; if (!is.null(minus)) contr[minus] <- -1
  out <- speckle::propeller.ttest(props, design = design, contrasts = contr,
                                  robust = TRUE, trend = FALSE, sort = TRUE)
  out %>% rownames_to_column("subtype") %>% mutate(n_samples = nrow(smp))
}
 
# LUAD vs LUSC (tumour only): histology adjusted for platform (Flex = LUSC only)
comp_hist <- function(cells)
  comp_test(filter(cells, tissue == "tumour"), ~ 0 + histology + platform,
            plus = "histologyLUSC", minus = "histologyLUAD")
comp_tissue <- function(cells)
  comp_test(cells, ~ 0 + tissue, plus = "tissuetumour", minus = "tissuenormal")
 
## ---------------------------------------------------------------------------
## 3. PREVALENCE (Fig 2A-C): binomial GLMMs -> ODDS RATIOS
## ---------------------------------------------------------------------------
# detected = raw count >= 1 in the cell. Every model includes scaled log10(nUMI) (detection
# rises with depth) and a random intercept for sample (cells in a sample are not independent).
 
# 3a. Fig 2A: each cell type vs endothelial cells (tumour). emmeans handles the contrasts;
#     Holm across cell types within gene. CIs are unadjusted. OR is reported as EC vs other.
prev_celltype <- function(cells, genes, cap = cfg$cap_other, seed = cfg$seed) {
  d <- cells %>% filter(tissue == "tumour", !is.na(ctype)) %>% cap_per_sample_type(cap, seed)
  d$ctype <- relevel(factor(d$ctype), ref = ENDO)
  map_dfr(genes, function(g) {
    m <- fit_glmm(as.formula(paste0("det__", g, " ~ ctype + log_nUMI + (1 | sample)")), d, binomial())
    if (is.null(m)) return(NULL)
    con <- contrast(emmeans(m, ~ ctype), method = "trt.vs.ctrl", ref = 1, adjust = "holm")
    s <- as_tibble(summary(con, type = "response", infer = TRUE))
    tibble(gene = g, contrast = as.character(s$contrast),           # "<other> / Endothelial"
           ratio = 1 / s$odds.ratio, LCL = 1 / s$asymp.UCL, UCL = 1 / s$asymp.LCL,
           p = s$p.value, p_adj = s$p.value, effect = "OR (EC vs other)")
  })
}
 
# 3b. Fig 2B/C: each EC subtype vs ALL OTHER ECs (the subtype is not part of its comparator)
prev_subtype <- function(cells, genes, tissue_sel = "tumour", min_cells_subtype = 15) {
  ec <- cells %>% filter(ctype == ENDO, tissue == tissue_sel, !is.na(ec_subtype_excl)) %>%
    keep_samples_min_ec()
  subs <- names(which(table(ec$ec_subtype_excl) >= min_cells_subtype))
  res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
    d <- ec %>% mutate(is_s = ec_subtype_excl == s)
    r <- fit_gene(d, g, "is_s + log_nUMI + (1 | sample)", "is_sTRUE")
    if (is.null(r)) return(NULL)
    r$subtype <- s; r
  }))
  add_padj(res, by = "gene")                                        # Holm across subtypes
}
 
# 3c. tumour vs normal ECs 
prev_tissue <- function(cells, genes) {
  ec <- cells %>% filter(ctype == ENDO) %>% keep_samples_min_ec()
  res <- map_dfr(genes, ~ fit_gene(ec, .x, "tissue + log_nUMI + (1 | dataset) + (1 | sample)",
                                   "tissuetumour",
                                   fallback_rhs = "tissue + log_nUMI + (1 | sample)"))
  add_padj(res)
}
 
# 3d. LUAD vs LUSC ECs (tumour). Three variants:
#   full         : histology adjusted for platform, dataset and sample
#   nonflex      : platform-matched (drops the Flex LUSC samples)
#   flex_effect  : Flex vs non-Flex within LUSC (size of the platform effect)
# prefilter = FALSE when `cells` is already the (possibly down-sampled) EC table.
prev_hist <- function(cells, genes, variant = c("full", "nonflex", "flex_effect"), prefilter = TRUE) {
  variant <- match.arg(variant)
  ec <- if (prefilter) cells %>% filter(ctype == ENDO, tissue == "tumour") %>% keep_samples_min_ec() else cells
  spec <- switch(variant,
    full        = list(d = ec,
                       rhs = "histology + platform + log_nUMI + (1 | dataset) + (1 | sample)",
                       fb  = "histology + platform + log_nUMI + (1 | sample)", term = "histologyLUSC"),
    nonflex     = list(d = filter(ec, platform == "other"),
                       rhs = "histology + log_nUMI + (1 | dataset) + (1 | sample)",
                       fb  = "histology + log_nUMI + (1 | sample)", term = "histologyLUSC"),
    flex_effect = list(d = filter(ec, histology == "LUSC"),
                       rhs = "platform + log_nUMI + (1 | sample)", fb = NULL, term = "platformFlex"))
  res <- map_dfr(genes, ~ fit_gene(spec$d, .x, spec$rhs, spec$term, fallback_rhs = spec$fb))
  if (is.null(res) || !nrow(res)) return(res)
  res$variant <- variant
  add_padj(res)                                                     # Holm across genes
}
 
## ---------------------------------------------------------------------------
## 3d. PREVALENCE BY GROUP (LUAD / LUSC / Normal)
## ---------------------------------------------------------------------------
 
# Pre-fit check for (quasi-)separation
check_separation_risk <- function(d, gene, min_events = 1) {
  det <- d[[paste0("det__", gene)]][d$is_s]
  n_pos <- sum(det); n_neg <- sum(!det)
  list(risky = n_pos < min_events || n_neg < min_events, n_pos = n_pos, n_neg = n_neg)
}
 
# Build the fixed/random-effects formula for ONE group's cells
build_rhs_ladder <- function(d, term, depth_term = "log_nUMI") {
  has_platform <- n_distinct(d$platform) > 1
  has_dataset  <- n_distinct(d$dataset) > 1
  rungs <- character(0)
  if (has_platform && has_dataset)
    rungs <- c(rungs, paste(term, "+ platform +", depth_term, "+ (1 | dataset) + (1 | sample)"))
  if (has_platform)
    rungs <- c(rungs, paste(term, "+ platform +", depth_term, "+ (1 | sample)"))
  if (has_dataset)
    rungs <- c(rungs, paste(term, "+", depth_term, "+ (1 | dataset) + (1 | sample)"))
  rungs <- c(rungs, paste(term, "+", depth_term, "+ (1 | sample)"))
  rungs <- c(rungs, paste(term, "+", depth_term))            # last resort: no random effects
  unique(rungs)
}
 
# Thin wrapper kept for the existing fit_gene_one_group() call site.
build_prev_rhs_ladder <- function(d) build_rhs_ladder(d, "is_s", "log_nUMI")
 
# One group (LUAD, LUSC, or Normal), one gene, one subtype.
fit_gene_one_group <- function(d, gene, term_main = "is_sTRUE", min_events = 1, label = "") {
  sep <- check_separation_risk(d, gene, min_events)
  if (sep$risky) {
    message(label, ": separation risk for ", gene, " (", sep$n_pos, " detected / ",
            sep$n_neg, " undetected among is_s==TRUE cells; need >= ", min_events,
            " of each) -- skipped, not fit")
    return(NULL)
  }
  ladder <- build_prev_rhs_ladder(d)
  fit_gene(d, gene, ladder[1], term_main, fallback_rhs = ladder[-1])
}
 
# Build the group column
add_ec_group <- function(cells) {
  cells %>% mutate(group = if_else(tissue == "normal", "Normal", as.character(histology)))
}

# Each EC subtype vs all other ECs, fit SEPARATELY within each group (LUAD,
# LUSC, Normal)
prev_subtype_by_group <- function(cells, genes, min_cells_subtype = 15,
                                  min_samples_subtype = 3, min_events = 1,
                                  group_levels = c("Normal", "LUAD", "LUSC")) {
  ec <- add_ec_group(cells) %>% filter(ctype == ENDO, !is.na(ec_subtype_excl), !is.na(group)) %>%
    keep_samples_min_ec()
  groups_present <- intersect(group_levels, unique(ec$group))
  message("Groups present: ", paste(groups_present, collapse = ", "), " (",
          paste(sapply(groups_present, function(g) paste0(g, ": ", n_distinct(ec$sample[ec$group == g]),
                                                           " samples, ", sum(ec$group == g), " ECs")),
                collapse = "; "), ")")
  all_subs <- unique(ec$ec_subtype_excl)
  flag_rows <- list()
 
  out <- map_dfr(groups_present, function(h) {
    ec_h <- filter(ec, group == h)
    n_cells_by_sub <- ec_h %>% count(ec_subtype_excl, name = "n_cells")
    subs <- n_cells_by_sub %>% filter(n_cells >= min_cells_subtype) %>% pull(ec_subtype_excl)
    dropped_cells <- setdiff(all_subs, subs)
    if (length(dropped_cells)) {
      message("[", h, "] subtypes dropped for < ", min_cells_subtype, " cells: ",
              paste(dropped_cells, collapse = ", "))
      flag_rows[[length(flag_rows) + 1]] <<- tibble(
        group = h, subtype = dropped_cells, gene = NA_character_, status = "dropped_cell_count",
        n_cells = coalesce(n_cells_by_sub$n_cells[match(dropped_cells, n_cells_by_sub$ec_subtype_excl)], 0L),
        n_samples = NA_integer_, n_pos = NA_integer_, n_neg = NA_integer_)
    }
 
    res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
      d <- ec_h %>% mutate(is_s = ec_subtype_excl == s)
      n_samp_pos <- n_distinct(d$sample[d$is_s])
      if (n_samp_pos < min_samples_subtype) {
        message("[", h, "] ", s, "/", g, ": only ", n_samp_pos, " samples carry this subtype, ",
                "skipped (< ", min_samples_subtype, ")")
        flag_rows[[length(flag_rows) + 1]] <<- tibble(
          group = h, subtype = s, gene = g, status = "dropped_sample_count",
          n_cells = sum(d$is_s), n_samples = n_samp_pos, n_pos = NA_integer_, n_neg = NA_integer_)
        return(NULL)
      }
      sep <- check_separation_risk(d, g, min_events)
      if (sep$risky) {
        flag_rows[[length(flag_rows) + 1]] <<- tibble(
          group = h, subtype = s, gene = g, status = "dropped_separation",
          n_cells = sum(d$is_s), n_samples = n_samp_pos, n_pos = sep$n_pos, n_neg = sep$n_neg)
      }
      r <- fit_gene_one_group(d, g, min_events = min_events, label = paste0("[", h, "] ", s, "/", g))
      if (is.null(r)) {
        if (!sep$risky)
          flag_rows[[length(flag_rows) + 1]] <<- tibble(
            group = h, subtype = s, gene = g, status = "dropped_fit_failed",
            n_cells = sum(d$is_s), n_samples = n_samp_pos, n_pos = sep$n_pos, n_neg = sep$n_neg)
        return(NULL)
      }
      flag_rows[[length(flag_rows) + 1]] <<- tibble(
        group = h, subtype = s, gene = g, status = "kept",
        n_cells = sum(d$is_s), n_samples = n_samp_pos, n_pos = sep$n_pos, n_neg = sep$n_neg)
      r$subtype <- s; r$label <- s; r$n_samples_pos <- n_samp_pos; r
    }))
    if (nrow(res)) res$group <- h
    res
  }) %>% mutate(group = factor(group, levels = group_levels)) %>%
    add_padj(by = c("group", "gene"))                      # Holm within EACH group's own comparison set
 
  attr(out, "flagged") <- bind_rows(flag_rows)
  out
}

get_flagged <- function(df) attr(df, "flagged")
 
# Formal test of whether the subtype's effect on detection DIFFERS between a
# PAIR of groups 
fit_binary_group_interaction <- function(d_all, gene, min_samples = 5, min_events = 1, label = "") {
  n_per_grp <- d_all %>% filter(is_flag) %>% distinct(sample, group) %>% count(group)
  qualifying <- n_per_grp$group[n_per_grp$n >= min_samples]
  dropped_grp <- setdiff(unique(d_all$group), qualifying)
  if (length(dropped_grp))
    message(label, ": group(s) excluded from the interaction model for < ", min_samples,
            " samples: ", paste(dropped_grp, collapse = ", "))
  if (length(qualifying) < 2) {
    message(label, ": fewer than 2 qualifying groups -- no pairwise comparison possible")
    return(NULL)
  }
  d <- d_all %>% filter(group %in% qualifying) %>% mutate(group = droplevels(factor(group)))
  sep_bad <- d %>% group_by(group) %>%
    summarise(n_pos = sum(.data[[paste0("det__", gene)]][is_flag]),
              n_neg = sum(!.data[[paste0("det__", gene)]][is_flag]), .groups = "drop") %>%
    filter(n_pos < min_events | n_neg < min_events)
  if (nrow(sep_bad)) {
    message(label, ": separation risk within group(s) ", paste(sep_bad$group, collapse = ", "),
            " -- interaction skipped, not fit")
    return(NULL)
  }

  ladder <- build_rhs_ladder(d, "is_flag * group", "log_nUMI")
  m <- fit_glmm(as.formula(paste0("det__", gene, " ~ ", ladder[1])), d, binomial(),
               fallback = lapply(ladder[-1], function(r) as.formula(paste0("det__", gene, " ~ ", r))))
  if (is.null(m)) return(NULL)
  emm <- tryCatch(emmeans(m, ~ is_flag * group), error = function(e) NULL)
  if (is.null(emm)) return(NULL)
  con <- tryCatch(contrast(emm, interaction = c("revpairwise", "revpairwise")), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  z <- qnorm(0.975)
  as_tibble(summary(con)) %>%
    transmute(pair = is_flag_revpairwise, group_pair = group_revpairwise,
             estimate, se = SE,
             ratio_of_ORs = exp(estimate), LCL = exp(estimate - z * se), UCL = exp(estimate + z * se),
             p = p.value) %>%
    mutate(n_groups_compared = length(qualifying), groups_compared = paste(qualifying, collapse = " vs "),
           fit_rung = attr(m, "fit_rung") %||% "formula 0 (default optimizer)",
           used_fallback = isTRUE(attr(m, "used_fallback")))
}
 
prev_subtype_group_interaction <- function(cells, genes, min_cells_subtype = 15,
                                           min_samples_subtype = 3, min_events = 1,
                                           group_levels = c("Normal", "LUAD", "LUSC")) {
  ec <- add_ec_group(cells) %>% filter(ctype == ENDO, !is.na(ec_subtype_excl), !is.na(group)) %>%
    keep_samples_min_ec()
  subs <- ec %>% count(ec_subtype_excl) %>% filter(n >= min_cells_subtype) %>% pull(ec_subtype_excl)
 
  res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
    d_all <- ec %>% mutate(is_flag = ec_subtype_excl == s)
    r <- fit_binary_group_interaction(d_all, g, min_samples_subtype, min_events,
                                      label = paste0(s, "/", g))
    if (is.null(r)) return(NULL)
    r %>% mutate(subtype = s, label = s, gene = g, effect = "interaction (ratio of odds ratios)")
  }))
  add_padj(res, by = "gene")   # Holm across subtype x group-pair tests within gene
}
 
 
## ---------------------------------------------------------------------------
## 3e. EC vs OTHER CELL TYPES (Fig 2A), SPLIT BY HISTOLOGY
## ---------------------------------------------------------------------------

check_separation_multilevel <- function(d, gene, level_col, min_events = 1) {
  det <- d[[paste0("det__", gene)]]
  d %>% mutate(.det = det) %>%
    group_by(.data[[level_col]]) %>%
    summarise(n_pos = sum(.det), n_neg = sum(!.det), .groups = "drop") %>%
    mutate(ok = n_pos >= min_events & n_neg >= min_events)
}
 
# Each OTHER cell type vs EC (the Fig 2A comparison), fit SEPARATELY within
# each group
prev_celltype_by_group <- function(cells, genes, cap = cfg$cap_other, seed = cfg$seed,
                                   min_cells_type = 15, min_samples_type = 3, min_events = 1,
                                   group_levels = c("Normal", "LUAD", "LUSC"),
                                   exclude_by_group = list(Normal = "Malignant")) {
  d0 <- add_ec_group(cells) %>% filter(!is.na(ctype), !is.na(group)) %>%
    cap_per_sample_type(cap, seed) %>% apply_group_exclusions(exclude_by_group)
  groups_present <- intersect(group_levels, unique(d0$group))
  message("Groups present: ", paste(groups_present, collapse = ", "))
  flag_rows <- list()
  for (grp in intersect(names(exclude_by_group), groups_present))
    flag_rows[[length(flag_rows) + 1]] <- tibble(
      group = grp, ctype = exclude_by_group[[grp]], gene = NA_character_,
      status = "excluded_by_design", n_cells = NA_integer_, n_samples = NA_integer_,
      n_pos = NA_integer_, n_neg = NA_integer_)
 
  out <- map_dfr(groups_present, function(h) {
    d_h <- filter(d0, group == h)
    type_n  <- d_h %>% count(ctype)
    type_ns <- d_h %>% distinct(sample, ctype) %>% count(ctype, name = "n_samples")
    type_info <- type_n %>% left_join(type_ns, by = "ctype")
    qualifying <- type_info %>% filter(n >= min_cells_type, n_samples >= min_samples_type) %>% pull(ctype)
    if (!ENDO %in% qualifying) {
      message("[", h, "] Endothelial reference does not meet the min_cells/min_samples bar ",
              "in this group -- skipping this group entirely")
      flag_rows[[length(flag_rows) + 1]] <<- tibble(
        group = h, ctype = as.character(type_info$ctype), gene = NA_character_,
        status = "dropped_reference_ec", n_cells = type_info$n, n_samples = type_info$n_samples,
        n_pos = NA_integer_, n_neg = NA_integer_)
      return(NULL)
    }
    dropped <- setdiff(unique(d_h$ctype), qualifying)
    if (length(dropped)) {
      message("[", h, "] cell types dropped for < ", min_cells_type, " cells or < ",
              min_samples_type, " samples: ", paste(dropped, collapse = ", "))
      di <- type_info %>% filter(ctype %in% dropped)
      flag_rows[[length(flag_rows) + 1]] <<- tibble(
        group = h, ctype = as.character(di$ctype), gene = NA_character_, status = "dropped_cell_or_sample_count",
        n_cells = di$n, n_samples = di$n_samples, n_pos = NA_integer_, n_neg = NA_integer_)
    }
    d_h <- d_h %>% filter(ctype %in% qualifying) %>% mutate(ctype = droplevels(factor(ctype)))
    d_h$ctype <- relevel(d_h$ctype, ref = ENDO)
 
    map_dfr(genes, function(g) {
      sep <- check_separation_multilevel(d_h, g, "ctype", min_events)
      bad <- sep$ctype[!sep$ok]
      d_g <- d_h
      if (length(bad)) {
        message("[", h, "] ", g, ": cell type(s) dropped for separation risk: ", paste(bad, collapse = ", "))
        sb <- sep %>% filter(ctype %in% bad)
        flag_rows[[length(flag_rows) + 1]] <<- tibble(
          group = h, ctype = as.character(sb$ctype), gene = g, status = "dropped_separation",
          n_cells = sb$n_pos + sb$n_neg, n_samples = NA_integer_, n_pos = sb$n_pos, n_neg = sb$n_neg)
        d_g <- d_g %>% filter(!ctype %in% bad) %>% mutate(ctype = droplevels(ctype))
        if (n_distinct(d_g$ctype) < 2 || !ENDO %in% levels(d_g$ctype)) return(NULL)
      }
      ladder <- build_rhs_ladder(d_g, "ctype", "log_nUMI")
      m <- fit_glmm(as.formula(paste0("det__", g, " ~ ", ladder[1])), d_g, binomial(),
                   fallback = lapply(ladder[-1], function(r) as.formula(paste0("det__", g, " ~ ", r))))
      if (is.null(m)) {
        flag_rows[[length(flag_rows) + 1]] <<- tibble(
          group = h, ctype = as.character(setdiff(levels(d_g$ctype), c(ENDO, bad))),
          gene = g, status = "dropped_fit_failed",
          n_cells = NA_integer_, n_samples = NA_integer_, n_pos = NA_integer_, n_neg = NA_integer_)
        return(NULL)
      }
      con <- contrast(emmeans(m, ~ ctype), method = "trt.vs.ctrl", ref = 1, adjust = "holm")
      s <- as_tibble(summary(con, type = "response", infer = TRUE))
      kept_types <- setdiff(levels(d_g$ctype), c(ENDO, bad))
      flag_rows[[length(flag_rows) + 1]] <<- tibble(
        group = h, ctype = kept_types, gene = g, status = "kept",
        n_cells = NA_integer_, n_samples = NA_integer_, n_pos = NA_integer_, n_neg = NA_integer_)
      tibble(group = h, gene = g, contrast = as.character(s$contrast),
             label = sub(" / .*$", "", as.character(s$contrast)),  
             ratio = 1 / s$odds.ratio, LCL = 1 / s$asymp.UCL, UCL = 1 / s$asymp.LCL,
             p = s$p.value, p_adj = s$p.value,                     
             effect = "OR (EC vs other)",                           
             fit_rung = attr(m, "fit_rung") %||% "formula 0 (default optimizer)",
             used_fallback = isTRUE(attr(m, "used_fallback")))
    })                                                              #   second add_padj() layer, matching the
  }) %>% mutate(group = factor(group, levels = group_levels))       #   original (pooled) prev_celltype()
 
  attr(out, "flagged") <- bind_rows(flag_rows)
  out
}
 
# Formal test: does a given (other cell type vs EC) comparison differ between
# groups? 
prev_celltype_group_interaction <- function(cells, genes, cap = cfg$cap_other, seed = cfg$seed,
                                            min_cells_type = 15, min_samples_type = 3, min_events = 1,
                                            group_levels = c("Normal", "LUAD", "LUSC"),
                                            exclude_by_group = list(Normal = "Malignant")) {
  d0 <- add_ec_group(cells) %>% filter(!is.na(ctype), !is.na(group)) %>%
    cap_per_sample_type(cap, seed) %>% apply_group_exclusions(exclude_by_group) %>%
    filter(group %in% group_levels)
  other_types <- setdiff(unique(d0$ctype[d0$ctype != ENDO]), NA)
  # a candidate type must clear min_cells/min_samples in at least one group 
  type_n  <- d0 %>% count(ctype)
  type_ns <- d0 %>% distinct(sample, ctype) %>% count(ctype, name = "n_samples")
  candidates <- type_n %>% left_join(type_ns, by = "ctype") %>%
    filter(ctype %in% other_types, n >= min_cells_type, n_samples >= min_samples_type) %>% pull(ctype)
 
  res <- map_dfr(genes, function(g) map_dfr(candidates, function(ct) {
    d_all <- d0 %>% filter(ctype %in% c(ENDO, ct)) %>% mutate(is_flag = ctype == ct)
    r <- fit_binary_group_interaction(d_all, g, min_samples_type, min_events, label = paste0(ct, "/", g))
    if (is.null(r)) return(NULL)
    r %>% mutate(ctype = ct, label = ct, gene = g, effect = "interaction (ratio of odds ratios)")
  }))
  add_padj(res, by = "gene")   # Holm across cell-type x group-pair tests within gene
}
 
## ---------------------------------------------------------------------------
## 4. MAGNITUDE
## ---------------------------------------------------------------------------
# 4a. Pseudobulk: sum counts per sample (x optional subtype)
pseudobulk <- function(counts, meta, by = "sample",
                       keep = c("dataset", "platform", "histology", "tissue"),
                       min_cells = cfg$min_cells_pb) {
  stopifnot(identical(colnames(counts), rownames(meta)))
  key <- do.call(paste, c(unname(as.list(meta[by])), list(sep = "|")))
  g   <- factor(key)
  M   <- Matrix::sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1,
                              dims = c(length(g), nlevels(g)),
                              dimnames = list(colnames(counts), levels(g)))
  pb   <- counts %*% M                                              # genes x profiles
  info <- meta %>% mutate(.key = key) %>% group_by(.key) %>%
    summarise(n_cells = dplyr::n(), across(all_of(c(by, keep)), ~ dplyr::first(.x)), .groups = "drop")
  info <- info[match(colnames(pb), info$.key), ]
  ok   <- info$n_cells >= min_cells
  info <- as.data.frame(info[ok, ]); rownames(info) <- info$.key
  list(counts = as.matrix(pb[, ok, drop = FALSE]), info = info)
}
 
# edgeR quasi-likelihood on pseudobulk profiles
run_edger <- function(pb, coef_name, design_formula = ~ dataset + histology,
                      panel = cfg$genes, min_flex_profiles = 3) {
  y <- DGEList(counts = pb$counts, samples = pb$info)
  design <- model.matrix(design_formula, data = y$samples)
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    message("dropping non-estimable design columns: ", paste(ne, collapse = ", "))
    design <- design[, setdiff(colnames(design), ne), drop = FALSE]
  }
  if (!coef_name %in% colnames(design))
    stop("'", coef_name, "' is not estimable with this design; try ~ platform + <group>")
  keep <- filterByExpr(y, design = design)
  flex <- y$samples$platform == "Flex"
  if (any(flex))
    keep <- keep & rowSums(y$counts[, flex, drop = FALSE] > 0) >= min(min_flex_profiles, sum(flex))
  force_keep <- rownames(y) %in% panel & rowSums(y$counts > 0) >= 3
  y <- y[keep | force_keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coef_name)
  tt  <- topTags(qlf, n = Inf)$table %>% rownames_to_column("gene") %>% as_tibble()
  inp <- tt$gene %in% panel
  tt$FDR_panel <- NA_real_; tt$FDR_panel[inp] <- p.adjust(tt$PValue[inp], "BH")
  tt$n_profiles <- ncol(y)
  tt
}
 
# LUAD vs LUSC ECs (tumour)
pb_hist <- function(cells, counts, panel = cfg$genes, prefilter = TRUE,
                    design_formula = ~ platform + histology) {
  ec <- if (prefilter) cells %>% filter(ctype == ENDO, tissue == "tumour") %>% keep_samples_min_ec() else cells
  meta <- as.data.frame(ec); rownames(meta) <- meta$cell
  pb <- pseudobulk(counts[, ec$cell, drop = FALSE], meta, by = "sample")
  run_edger(pb, "histologyLUSC", design_formula, panel)
}
 
# tumour vs normal ECs
pb_tissue <- function(cells, counts, panel = cfg$genes) {
  ec <- cells %>% filter(ctype == ENDO) %>% keep_samples_min_ec()
  meta <- as.data.frame(ec); rownames(meta) <- meta$cell
  pb <- pseudobulk(counts[, ec$cell, drop = FALSE], meta, by = "sample")
  run_edger(pb, "tissuetumour", ~ platform + tissue, panel)
}
 
# 4b. EC subtype vs other ECs
mag_subtype <- function(cells, genes, tissue_sel = "tumour", min_cells_subtype = 15) {
  ec <- cells %>% filter(ctype == ENDO, tissue == tissue_sel, !is.na(ec_subtype_excl)) %>%
    keep_samples_min_ec()
  subs <- names(which(table(ec$ec_subtype_excl) >= min_cells_subtype))
  res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
    d <- ec %>% mutate(is_s = ec_subtype_excl == s)
    r <- fit_gene(d, g, "is_s + offset(log(nUMI)) + (1 | sample)", "is_sTRUE",
                  response = "cnt", effect = "RR", family = nbinom2())
    if (is.null(r)) return(NULL)
    r$subtype <- s; r
  }))
  add_padj(res, by = "gene")
}
 
## ---------------------------------------------------------------------------
## 4c. MAGNITUDE BY GROUP
## ---------------------------------------------------------------------------

build_mag_rhs_ladder <- function(d) build_rhs_ladder(d, "is_s", "offset(log(nUMI))")
 
check_zero_risk <- function(d, gene, min_events = 1, expressors_only = TRUE) {
  cnt <- d[[paste0("cnt__", gene)]][d$is_s]
  n_cells <- length(cnt); n_nonzero <- sum(cnt > 0)
  risky <- n_cells < min_events || (!expressors_only && n_nonzero < min_events)
  list(risky = risky, n_cells = n_cells, n_nonzero = n_nonzero)
}
 
# One group, one gene, one subtype
fit_magnitude_one_group <- function(d, gene, term_main = "is_sTRUE", expressors_only = TRUE,
                                    min_events = 1, label = "") {
  d_fit <- if (expressors_only) filter(d, .data[[paste0("det__", gene)]]) else d
  chk <- check_zero_risk(d_fit, gene, min_events, expressors_only)
  if (chk$risky) {
    message(label, ": too few cells for ", gene, " (", chk$n_cells, " cells",
            if (!expressors_only) paste0(", ", chk$n_nonzero, " nonzero") else "",
            " among is_s==TRUE; need >= ", min_events, ") -- skipped, not fit")
    return(NULL)
  }
  ladder <- build_mag_rhs_ladder(d_fit)
  fit_gene(d_fit, gene, ladder[1], term_main, fallback_rhs = ladder[-1],
          response = "cnt", family = nbinom2(),
          effect = if (expressors_only) "RR (magnitude among expressors)" else "RR (magnitude)")
}
 
# Each EC subtype vs other ECs, MAGNITUDE
mag_subtype_by_group <- function(cells, genes, min_cells_subtype = 15, min_samples_subtype = 3,
                                 min_events = 1, expressors_only = TRUE,
                                 group_levels = c("Normal", "LUAD", "LUSC")) {
  ec <- add_ec_group(cells) %>% filter(ctype == ENDO, !is.na(ec_subtype_excl), !is.na(group)) %>%
    keep_samples_min_ec()
  groups_present <- intersect(group_levels, unique(ec$group))
  message("Groups present: ", paste(groups_present, collapse = ", "))
 
  map_dfr(groups_present, function(h) {
    ec_h <- filter(ec, group == h)
    subs <- ec_h %>% count(ec_subtype_excl) %>% filter(n >= min_cells_subtype) %>% pull(ec_subtype_excl)
    dropped <- setdiff(unique(ec$ec_subtype_excl), subs)
    if (length(dropped))
      message("[", h, "] subtypes dropped for < ", min_cells_subtype, " cells: ", paste(dropped, collapse = ", "))
 
    res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
      d <- ec_h %>% mutate(is_s = ec_subtype_excl == s)
      n_samp_pos <- n_distinct(d$sample[d$is_s])
      if (n_samp_pos < min_samples_subtype) {
        message("[", h, "] ", s, "/", g, ": only ", n_samp_pos, " samples carry this subtype, ",
                "skipped (< ", min_samples_subtype, ")")
        return(NULL)
      }
      r <- fit_magnitude_one_group(d, g, expressors_only = expressors_only, min_events = min_events,
                                   label = paste0("[", h, "] ", s, "/", g))
      if (is.null(r)) return(NULL)
      r$subtype <- s; r$label <- s; r$n_samples_pos <- n_samp_pos; r
    }))
    if (nrow(res)) res$group <- h
    res
  }) %>% mutate(group = factor(group, levels = group_levels)) %>%
    add_padj(by = c("group", "gene"))
}
 
# NB analog of fit_binary_group_interaction()
fit_magnitude_group_interaction <- function(d_all, gene, min_samples = 5, min_events = 1,
                                            expressors_only = TRUE, label = "") {
  n_per_grp <- d_all %>% filter(is_flag) %>% distinct(sample, group) %>% count(group)
  qualifying <- n_per_grp$group[n_per_grp$n >= min_samples]
  dropped_grp <- setdiff(unique(d_all$group), qualifying)
  if (length(dropped_grp))
    message(label, ": group(s) excluded from the interaction model for < ", min_samples,
            " samples: ", paste(dropped_grp, collapse = ", "))
  if (length(qualifying) < 2) {
    message(label, ": fewer than 2 qualifying groups -- no pairwise comparison possible")
    return(NULL)
  }
  d <- d_all %>% filter(group %in% qualifying) %>% mutate(group = droplevels(factor(group)))
  d_fit <- if (expressors_only) filter(d, .data[[paste0("det__", gene)]]) else d
 
  chk_bad <- d_fit %>% filter(is_flag) %>% group_by(group) %>%
    summarise(n = n(), n_nonzero = sum(.data[[paste0("cnt__", gene)]] > 0), .groups = "drop") %>%
    filter(n < min_events | (!expressors_only & n_nonzero < min_events))
  if (nrow(chk_bad)) {
    message(label, ": too few cells within group(s) ", paste(chk_bad$group, collapse = ", "),
            " -- interaction skipped, not fit")
    return(NULL)
  }
 
  ladder <- build_rhs_ladder(d_fit, "is_flag * group", "offset(log(nUMI))")
  m <- fit_glmm(as.formula(paste0("cnt__", gene, " ~ ", ladder[1])), d_fit, nbinom2(),
               fallback = lapply(ladder[-1], function(r) as.formula(paste0("cnt__", gene, " ~ ", r))))
  if (is.null(m)) return(NULL)
  emm <- tryCatch(emmeans(m, ~ is_flag * group), error = function(e) NULL)
  if (is.null(emm)) return(NULL)
  con <- tryCatch(contrast(emm, interaction = c("revpairwise", "revpairwise")), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  z <- qnorm(0.975)
  as_tibble(summary(con)) %>%
    transmute(pair = is_flag_revpairwise, group_pair = group_revpairwise,
             estimate, se = SE,
             ratio_of_RRs = exp(estimate), LCL = exp(estimate - z * se), UCL = exp(estimate + z * se),
             p = p.value) %>%
    mutate(n_groups_compared = length(qualifying), groups_compared = paste(qualifying, collapse = " vs "),
           fit_rung = attr(m, "fit_rung") %||% "formula 0 (default optimizer)",
           used_fallback = isTRUE(attr(m, "used_fallback")))
}
 
mag_subtype_group_interaction <- function(cells, genes, min_cells_subtype = 15,
                                          min_samples_subtype = 3, min_events = 1,
                                          expressors_only = TRUE,
                                          group_levels = c("Normal", "LUAD", "LUSC")) {
  ec <- add_ec_group(cells) %>% filter(ctype == ENDO, !is.na(ec_subtype_excl), !is.na(group)) %>%
    keep_samples_min_ec()
  subs <- ec %>% count(ec_subtype_excl) %>% filter(n >= min_cells_subtype) %>% pull(ec_subtype_excl)
 
  res <- map_dfr(genes, function(g) map_dfr(subs, function(s) {
    d_all <- ec %>% mutate(is_flag = ec_subtype_excl == s)
    r <- fit_magnitude_group_interaction(d_all, g, min_samples_subtype, min_events, expressors_only,
                                         label = paste0(s, "/", g))
    if (is.null(r)) return(NULL)
    r %>% mutate(subtype = s, label = s, gene = g, effect = "interaction (ratio of rate ratios)")
  }))
  add_padj(res, by = "gene")
}
 
# EC vs other cell types (Fig 2A analog), MAGNITUDE
mag_celltype_by_group <- function(cells, genes, cap = cfg$cap_other, seed = cfg$seed,
                                  min_cells_type = 15, min_samples_type = 3, min_events = 1,
                                  expressors_only = TRUE, group_levels = c("Normal", "LUAD", "LUSC"),
                                  exclude_by_group = list(Normal = "Malignant")) {
  d0 <- add_ec_group(cells) %>% filter(!is.na(ctype), !is.na(group)) %>%
    cap_per_sample_type(cap, seed) %>% apply_group_exclusions(exclude_by_group)
  groups_present <- intersect(group_levels, unique(d0$group))
  message("Groups present: ", paste(groups_present, collapse = ", "))
 
  map_dfr(groups_present, function(h) {
    d_h <- filter(d0, group == h)
    type_n  <- d_h %>% count(ctype)
    type_ns <- d_h %>% distinct(sample, ctype) %>% count(ctype, name = "n_samples")
    qualifying <- type_n %>% left_join(type_ns, by = "ctype") %>%
      filter(n >= min_cells_type, n_samples >= min_samples_type) %>% pull(ctype)
    if (!ENDO %in% qualifying) {
      message("[", h, "] Endothelial reference does not meet the min_cells/min_samples bar ",
              "in this group -- skipping this group entirely")
      return(NULL)
    }
    dropped <- setdiff(unique(d_h$ctype), qualifying)
    if (length(dropped))
      message("[", h, "] cell types dropped for < ", min_cells_type, " cells or < ",
              min_samples_type, " samples: ", paste(dropped, collapse = ", "))
    d_h <- d_h %>% filter(ctype %in% qualifying) %>% mutate(ctype = droplevels(factor(ctype)))
    d_h$ctype <- relevel(d_h$ctype, ref = ENDO)
 
    map_dfr(genes, function(g) {
      d_fit <- if (expressors_only) filter(d_h, .data[[paste0("det__", g)]]) else d_h
      chk <- d_fit %>% count(ctype)
      bad <- chk$ctype[chk$n < min_events]
      d_g <- d_fit
      if (length(bad)) {
        message("[", h, "] ", g, ": cell type(s) dropped for too few (expressing) cells: ",
                paste(bad, collapse = ", "))
        d_g <- d_g %>% filter(!ctype %in% bad) %>% mutate(ctype = droplevels(ctype))
        if (n_distinct(d_g$ctype) < 2 || !ENDO %in% levels(d_g$ctype)) return(NULL)
      }
      ladder <- build_rhs_ladder(d_g, "ctype", "offset(log(nUMI))")
      m <- fit_glmm(as.formula(paste0("cnt__", g, " ~ ", ladder[1])), d_g, nbinom2(),
                   fallback = lapply(ladder[-1], function(r) as.formula(paste0("cnt__", g, " ~ ", r))))
      if (is.null(m)) return(NULL)
      con <- contrast(emmeans(m, ~ ctype), method = "trt.vs.ctrl", ref = 1, adjust = "holm")
      s <- as_tibble(summary(con, type = "response", infer = TRUE))
      rcol <- intersect(c("ratio", "rate.ratio", "odds.ratio"), names(s))
      if (!length(rcol)) stop("emmeans summary has no recognized ratio column -- inspect names(s): ",
                              paste(names(s), collapse = ", "))
      tibble(group = h, gene = g, contrast = as.character(s$contrast),
             label = sub(" / .*$", "", as.character(s$contrast)),
             ratio = 1 / s[[rcol[1]]], LCL = 1 / s$asymp.UCL, UCL = 1 / s$asymp.LCL,
             p = s$p.value, p_adj = s$p.value,
             effect = if (expressors_only) "RR (magnitude among expressors, EC vs other)"
                      else "RR (magnitude, EC vs other)",
             fit_rung = attr(m, "fit_rung") %||% "formula 0 (default optimizer)",
             used_fallback = isTRUE(attr(m, "used_fallback")))
    })
  }) %>% mutate(group = factor(group, levels = group_levels))
}
 
mag_celltype_group_interaction <- function(cells, genes, cap = cfg$cap_other, seed = cfg$seed,
                                           min_cells_type = 15, min_samples_type = 3, min_events = 1,
                                           expressors_only = TRUE,
                                           group_levels = c("Normal", "LUAD", "LUSC"),
                                           exclude_by_group = list(Normal = "Malignant")) {
  d0 <- add_ec_group(cells) %>% filter(!is.na(ctype), !is.na(group)) %>%
    cap_per_sample_type(cap, seed) %>% apply_group_exclusions(exclude_by_group) %>%
    filter(group %in% group_levels)
  other_types <- setdiff(unique(d0$ctype[d0$ctype != ENDO]), NA)
  type_n  <- d0 %>% count(ctype)
  type_ns <- d0 %>% distinct(sample, ctype) %>% count(ctype, name = "n_samples")
  candidates <- type_n %>% left_join(type_ns, by = "ctype") %>%
    filter(ctype %in% other_types, n >= min_cells_type, n_samples >= min_samples_type) %>% pull(ctype)
 
  res <- map_dfr(genes, function(g) map_dfr(candidates, function(ct) {
    d_all <- d0 %>% filter(ctype %in% c(ENDO, ct)) %>% mutate(is_flag = ctype == ct)
    r <- fit_magnitude_group_interaction(d_all, g, min_samples_type, min_events, expressors_only,
                                         label = paste0(ct, "/", g))
    if (is.null(r)) return(NULL)
    r %>% mutate(ctype = ct, label = ct, gene = g, effect = "interaction (ratio of rate ratios)")
  }))
  add_padj(res, by = "gene")
}

## ---------------------------------------------------------------------------
## 5. LUSC DOWN-SAMPLING (sensitivity analyses against LUAD)
## ---------------------------------------------------------------------------
# LUSC ECs are far more numerous per sample than LUAD ECs (median ~384 vs ~52). Each draw
# keeps every LUAD EC and at most `cap_per_sample` randomly chosen ECs from each LUSC sample
draw_cells <- function(d, cap_per_sample = cfg$cap_lusc_ec, target_total = NULL, seed = 1) {
  set.seed(seed)
  keep    <- rep(TRUE, nrow(d))
  is_lusc <- d$histology == "LUSC"
  for (ix in split(which(is_lusc), d$sample[is_lusc])) {
    if (length(ix) > cap_per_sample) keep[setdiff(ix, sample(ix, cap_per_sample))] <- FALSE
  }
  k <- which(keep & is_lusc)
  if (!is.null(target_total) && length(k) > target_total)
    keep[setdiff(k, sample(k, target_total))] <- FALSE
  d[keep, ]
}
 
# Repeat any fit on n_iter independent draws. fit_fun(d) must return a data frame (or NULL).
run_downsample_iterations <- function(d, fit_fun, n_iter = cfg$n_iter, seed0 = 1000,
                                      cap_per_sample = cfg$cap_lusc_ec, target_total = NULL,
                                      cores = cfg$cores) {
  one <- function(i) {
    r <- tryCatch(fit_fun(draw_cells(d, cap_per_sample, target_total, seed = seed0 + i)),
                  error = function(e) { message("iteration ", i, ": ", conditionMessage(e)); NULL })
    if (is.null(r) || !nrow(r)) return(NULL)
    r$iter <- i; r
  }
  out <- if (cores > 1) parallel::mclapply(seq_len(n_iter), one, mc.cores = cores) else lapply(seq_len(n_iter), one)
  bind_rows(out)
}
 
# Distribution across draws. These are re-draws of the SAME data, so the spread describes
# sensitivity to which LUSC cells were kept
summarise_iterations <- function(res, by = "gene", est = "ratio", p = "p_adj",
                                 null = 1, alpha = cfg$alpha) {
  res %>% group_by(across(all_of(by))) %>%
    summarise(n_iter = n_distinct(iter),
              median_est = median(.data[[est]]),
              q025 = quantile(.data[[est]], 0.025), q975 = quantile(.data[[est]], 0.975),
              frac_sig = mean(.data[[p]] < alpha, na.rm = TRUE),
              frac_above_null = mean(.data[[est]] > null),
              median_p = median(.data[[p]], na.rm = TRUE), .groups = "drop")
}


## ---------------------------------------------------------------------------
## 7. DRIVER
## ---------------------------------------------------------------------------
control_genes <- c("ACTB", "GAPDH", "RPL13A", "B2M", "PPIA")   # housekeeping controls
tb    <- build_cell_table(seu, genes = c(cfg$genes, control_genes))
cells <- resolve_ec_labels(tb$cells); genes <- tb$genes
counts_all <- GetAssayData(seu, assay = "RNA", layer = "counts")

# ---- composition (speckle; sensitivity) ----
comp_h <- comp_hist(cells)
comp_t <- comp_tissue(cells)

# ---- prevalence: pooled ----
p2a   <- prev_celltype(cells, genes)
p2bc  <- prev_subtype(cells, genes)
p_tis <- prev_tissue(cells, genes)

# ---- prevalence: by disease state (Figure 2A/C, Supp Fig 2C, Supp Tables 2 and 4) ----
p2bc_grp <- prev_subtype_by_group(cells, genes)
p2bc_int <- prev_subtype_group_interaction(cells, genes)
p2a_grp  <- prev_celltype_by_group(cells, genes)
p2a_int  <- prev_celltype_group_interaction(cells, genes)

# ---- magnitude among expressors: by disease state (Figure 2B/D, Supp Tables 3 and 4) ----
m2bc_grp <- mag_subtype_by_group(cells, genes)
m2bc_int <- mag_subtype_group_interaction(cells, genes)
m2a_grp  <- mag_celltype_by_group(cells, genes)
m2a_int  <- mag_celltype_group_interaction(cells, genes)

icl        <- mag_celltype_by_group(cells, cfg$genes)
ctrl_panel <- mag_celltype_by_group(cells, control_genes)

# ---- LUAD vs LUSC variants and pseudobulk (sensitivity) ----
p_hist_full    <- prev_hist(cells, genes, "full")
p_hist_nonflex <- prev_hist(cells, genes, "nonflex")
p_hist_flex    <- prev_hist(cells, genes, "flex_effect")
m_hist   <- pb_hist(cells, counts_all)
m_tissue <- pb_tissue(cells, counts_all)
m_sub    <- mag_subtype(cells, genes)

# ---- LUSC down-sampling, 100 draws (sensitivity) ----
ec_t <- cells %>% filter(ctype == ENDO, tissue == "tumour") %>% keep_samples_min_ec()
it_prev <- run_downsample_iterations(ec_t, function(d) prev_hist(d, genes, "full", prefilter = FALSE))
it_pb   <- run_downsample_iterations(ec_t, function(d) pb_hist(d, counts_all, prefilter = FALSE) %>%
                                       filter(gene %in% cfg$genes))
it_comp <- run_downsample_iterations(ec_t, function(d) comp_test(d, ~ 0 + histology + platform,
                                       plus = "histologyLUSC", minus = "histologyLUAD"))
print(summarise_iterations(it_prev, by = "gene"))
print(summarise_iterations(it_pb, by = "gene", est = "logFC", p = "FDR_panel", null = 0))
print(summarise_iterations(it_comp, by = "subtype", est = "Tstatistic", p = "FDR", null = 0))

## ---------------------------------------------------------------------------
## 8. SAVE
## ---------------------------------------------------------------------------
all_results <- list(
  composition  = list(comp_h = comp_h, comp_t = comp_t),
  prevalence   = list(p2a = p2a, p2bc = p2bc, p_tis = p_tis,
                      p_hist_full = p_hist_full, p_hist_nonflex = p_hist_nonflex, p_hist_flex = p_hist_flex,
                      p2bc_grp = p2bc_grp, p2bc_int = p2bc_int, p2a_grp = p2a_grp, p2a_int = p2a_int),
  magnitude    = list(m_hist = m_hist, m_tissue = m_tissue, m_sub = m_sub,
                      m2bc_grp = m2bc_grp, m2bc_int = m2bc_int, m2a_grp = m2a_grp, m2a_int = m2a_int,
                      icl = icl, ctrl_panel = ctrl_panel),
  downsampling = list(it_prev = it_prev, it_pb = it_pb, it_comp = it_comp))
qs_save(all_results, dpath("senstivityStats_all_results.qs2"))

# Supplementary Tables 2-4 
bind_rows(mutate(p2a_grp, comparison = "EC vs other cell type"),
          mutate(p2bc_grp, comparison = "EC subtype vs other ECs")) %>%
  write_csv(rpath("SuppTable2_prevalence_odds_ratios.csv"))
bind_rows(mutate(m2a_grp, comparison = "EC vs other cell type"),
          mutate(m2bc_grp, comparison = "EC subtype vs other ECs")) %>%
  write_csv(rpath("SuppTable3_magnitude_rate_ratios.csv"))
bind_rows(mutate(p2a_int,  measure = "prevalence", comparison = "EC vs other cell type"),
          mutate(p2bc_int, measure = "prevalence", comparison = "EC subtype vs other ECs"),
          mutate(m2a_int,  measure = "magnitude",  comparison = "EC vs other cell type"),
          mutate(m2bc_int, measure = "magnitude",  comparison = "EC subtype vs other ECs")) %>%
  write_csv(rpath("SuppTable4_disease_state_interactions.csv"))