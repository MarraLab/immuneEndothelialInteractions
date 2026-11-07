# =============================================================================
# config/paths.R
#
# Every input/output location used by the analysis and figure scripts.
# All scripts are run from the repository root, e.g.
#   Rscript analysis/02_endothelial_cells/subsetECs.R
# =============================================================================

DATA_DIR       <- Sys.getenv("IEI_DATA_DIR",      "data")
# Per-sample processed Seurat objects for the internal (GSC/POG) cohort:
#   <INTERNAL_OBJ_DIR>/<sample>/rna_ProcessedSeurat.RDS
INTERNAL_OBJ_DIR <- Sys.getenv("IEI_INTERNAL_OBJ_DIR", file.path(DATA_DIR, "InternalSamples"))
# Per-donor Seurat objects built from the Salcher et al. (2022) CELLxGENE atlas
PUBLIC_OBJ_DIR <- Sys.getenv("IEI_PUBLIC_OBJ_DIR", file.path(DATA_DIR, "PublicData"))
# Per-sample inferCNV output folders
INFERCNV_DIR   <- Sys.getenv("IEI_INFERCNV_DIR",  file.path(DATA_DIR, "LungCohort", "InferCNV"))
# Raw CosMx SMI NSCLC FFPE bundle (Bruker): <COSMX_RAW_DIR>/<sample>/<sample>-Flat_files_and_images/
COSMX_RAW_DIR  <- Sys.getenv("IEI_COSMX_RAW_DIR", file.path(DATA_DIR, "CosMxRaw"))
# Reference files (cell-type annotation reference, gene-position file, biomaRt cache)
REF_DIR        <- Sys.getenv("IEI_REF_DIR",       file.path(DATA_DIR, "reference"))

RESULTS_DIR    <- Sys.getenv("IEI_RESULTS_DIR",   "results")
FIG_DIR        <- Sys.getenv("IEI_FIG_DIR",       file.path("figures", "output"))

# Shorthand constructors
dpath <- function(...) file.path(DATA_DIR, ...)
rpath <- function(...) file.path(RESULTS_DIR, ...)
fpath <- function(...) file.path(FIG_DIR, ...)

for (d in c(RESULTS_DIR, FIG_DIR, file.path(RESULTS_DIR, "diagnostics"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}