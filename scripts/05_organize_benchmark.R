#!/usr/bin/env Rscript
# Organize all database outputs into the benchmark directory tree.
#
# This script copies / symlinks DISCO RDS files and TabulaSapiens by_organ RDS
# files into organized_by_tissue/{Tissue}/ following the naming convention
# {Database}_{filename}.rds.
#
# HCL and MCA RDS files are already placed by script 04_convert_h5ad_to_seurat.R.
#
# Input:
#   01_DISCO/rds/*.rds                      (29 files)
#   06_TabulaSapiens/by_organ/TabulaSapiens_*.rds  (25 files)
#   organized_by_tissue/{Tissue}/HCL_HCL_*.rds     (already done by script 04)
#   organized_by_tissue/{Tissue}/MCA_MCA_*.rds     (already done by script 04)
#
# Output:
#   organized_by_tissue/  (112 files across 53 tissues)

library(Seurat)

BASE <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
setwd(BASE)

org_dir <- "organized_by_tissue"
dir.create(org_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(org_dir, "organization_summary.txt")
cat("", file = log_file)

log_msg <- function(...) {
  msg <- paste(Sys.time(), paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# ── DISCO filename → tissue mapping ──────────────────────────────────────────
DISCO_TISSUE_MAP <- list(
  "adipose_cell.rds"        = "Adipose",
  "adipose_nucleus.rds"     = "Adipose",
  "adrenal_gland.rds"       = "AdrenalGland",
  "basophil_mast_cell.rds"  = "ImmuneCell",
  "bladder.rds"             = "Bladder",
  "blood.rds"               = "Other_blood",
  "bone_marrow.rds"         = "BoneMarrow",
  "brain.rds"               = "Brain",
  "breast.rds"              = "MammaryGland",
  "breast_milk.rds"         = "MammaryGland",
  "dengue_blood.rds"        = "Blood_Disease",
  "eye.rds"                 = "Eye",
  "fibroblast.rds"          = "Fibroblast",
  "gingiva.rds"             = "Gingiva",
  "heart.rds"               = "Heart",
  "intestine.rds"           = "Intestine",
  "kidney.rds"              = "Kidney",
  "liver_cell.rds"          = "Liver",
  "liver_nucleus.rds"       = "Liver",
  "lung.rds"                = "Lung",
  "ovary.rds"               = "Gonad",
  "pancreas_cell.rds"       = "Pancreas",
  "placenta.rds"            = "Placenta",
  "skeletal_muscle.rds"     = "Muscle",
  "skin.rds"                = "Skin",
  "stomach.rds"             = "Stomach",
  "testis.rds"              = "Gonad",
  "thymus.rds"              = "Thymus",
  "tonsil.rds"              = "Tonsil"
)

# ── TabulaSapiens organ → tissue mapping ─────────────────────────────────────
TS_TISSUE_MAP <- list(
  "Bladder"         = "Bladder",
  "Blood"           = "Blood",
  "Bone_Marrow"     = "Bone_Marrow",
  "Eye"             = "Eye",
  "Fat"             = "Fat",
  "Heart"           = "Heart",
  "Intestine"       = "Intestine",
  "Kidney"          = "Kidney",
  "Large_Intestine" = "Large_Intestine",
  "Liver"           = "Liver",
  "Lung"            = "Lung",
  "Lymph_Node"      = "Lymph_Node",
  "Mammary"         = "Mammary",
  "Muscle"          = "Muscle",
  "Pancreas"        = "Pancreas",
  "Prostate"        = "Prostate",
  "Salivary_Gland"  = "Salivary_Gland",
  "Skin"            = "Skin",
  "Small_Intestine" = "Small_Intestine",
  "Spleen"          = "Spleen",
  "Thymus"          = "Thymus",
  "Tongue"          = "Tongue",
  "Trachea"         = "Trachea",
  "Uterus"          = "Uterus",
  "Vasculature"     = "Vasculature"
)

log_msg("=" * 60)
log_msg("Organizing benchmark dataset directory tree")
log_msg("=" * 60)

t0 <- Sys.time()
stats <- list(DISCO = 0, TabulaSapiens = 0, total = 0, total_size_mb = 0)

# ── DISCO ────────────────────────────────────────────────────────────────────
disco_src <- "01_DISCO/rds"
if (dir.exists(disco_src)) {
  disco_files <- list.files(disco_src, pattern = "\\.rds$", full.names = TRUE)
  for (f in disco_files) {
    fname <- basename(f)
    tissue <- DISCO_TISSUE_MAP[[fname]]
    if (is.null(tissue)) {
      log_msg(sprintf("  WARNING: unknown DISCO file mapping: %s", fname))
      next
    }
    dest_dir <- file.path(org_dir, tissue)
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(dest_dir, paste0("DISCO_", fname))
    file.copy(f, dest, overwrite = TRUE)
    stats$DISCO <- stats$DISCO + 1
    stats$total <- stats$total + 1
    stats$total_size_mb <- stats$total_size_mb + file.info(dest)$size / 1e6
    log_msg(sprintf("  DISCO  →  %s", dest))
  }
}

# ── TabulaSapiens ────────────────────────────────────────────────────────────
ts_dir <- "06_TabulaSapiens/by_organ"
if (dir.exists(ts_dir)) {
  ts_files <- list.files(ts_dir, pattern = "^TabulaSapiens_.*\\.rds$", full.names = TRUE)
  for (f in ts_files) {
    fname <- basename(f)
    organ <- sub("^TabulaSapiens_|\\.rds$", "", fname)
    tissue <- TS_TISSUE_MAP[[organ]]
    if (is.null(tissue)) {
      tissue <- organ  # use organ name directly if unmapped
    }
    dest_dir <- file.path(org_dir, tissue)
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(dest_dir, paste0("TabulaSapiens_", fname))
    file.copy(f, dest, overwrite = TRUE)
    stats$TabulaSapiens <- stats$TabulaSapiens + 1
    stats$total <- stats$total + 1
    stats$total_size_mb <- stats$total_size_mb + file.info(dest)$size / 1e6
    log_msg(sprintf("  TS     →  %s", dest))
  }
}

# ── Summary ──────────────────────────────────────────────────────────────────
total_elapsed <- difftime(Sys.time(), t0, units = "mins")
log_msg(sprintf("\nTotal files organized: %d", stats$total))
log_msg(sprintf("  DISCO: %d  |  TabulaSapiens: %d", stats$DISCO, stats$TabulaSapiens))
log_msg(sprintf("  HCL / MCA: already placed by script 04"))
log_msg(sprintf("Total size: %.1f GB", stats$total_size_mb / 1024))
log_msg(sprintf("Elapsed: %.1f min", total_elapsed))
log_msg(sprintf("\nOutput directory: %s", normalizePath(org_dir)))
sessionInfo()
