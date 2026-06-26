#!/usr/bin/env Rscript
# Convert HCL and MCA per-tissue h5ad files to Seurat RDS format.
#
# This is a batch conversion pipeline using easySCFr::easy.h5ad2seurat().
# Naming convention: {DB}_{original_name}.rds → organized_by_tissue/{Tissue}/DB_DB_Tissue.rds
#
# Input:
#   02_HCL/by_tissue/HCL_{tissue}.h5ad    (33 files)
#   02_HCL/MCA_by_tissue/MCA_{tissue}.h5ad  (25 files)
#
# Output:
#   organized_by_tissue/{Tissue}/HCL_HCL_{tissue}.rds
#   organized_by_tissue/{Tissue}/MCA_MCA_{tissue}.rds

library(easySCFr)
library(Seurat)

BASE <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
setwd(BASE)

log_file <- "conversion_log.txt"
cat("", file = log_file)

log_msg <- function(...) {
  msg <- paste(Sys.time(), paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

convert_and_save <- function(h5ad_path, rds_path, label) {
  log_msg(sprintf("[%s] Converting: %s", label, basename(h5ad_path)))
  t1 <- Sys.time()

  sobj <- easy.h5ad2seurat(h5ad_path)

  dir.create(dirname(rds_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(sobj, rds_path)

  elapsed <- difftime(Sys.time(), t1, units = "secs")
  size_mb <- file.info(rds_path)$size / 1e6
  log_msg(sprintf("  → %s  (%d cells, %.1f MB, %.1f s)",
                  basename(rds_path), ncol(sobj), size_mb, elapsed))
  return(TRUE)
}

log_msg("=" * 60)
log_msg("H5AD → Seurat RDS batch conversion")
log_msg("=" * 60)

# ── HCL files ────────────────────────────────────────────────────────────────
hcl_dir  <- "02_HCL/by_tissue"
hcl_files <- list.files(hcl_dir, pattern = "^HCL_.*\\.h5ad$", full.names = TRUE)
log_msg(sprintf("\nFound %d HCL h5ad files", length(hcl_files)))

for (f in hcl_files) {
  tissue <- sub("^HCL_|\\.h5ad$", "", basename(f))
  rds_path <- file.path("organized_by_tissue", tissue, sprintf("HCL_HCL_%s.rds", tissue))
  convert_and_save(f, rds_path, "HCL")
}

# ── MCA files ────────────────────────────────────────────────────────────────
mca_dir  <- "02_HCL/MCA_by_tissue"
mca_files <- list.files(mca_dir, pattern = "^MCA_.*\\.h5ad$", full.names = TRUE)
log_msg(sprintf("\nFound %d MCA h5ad files", length(mca_files)))

for (f in mca_files) {
  tissue <- sub("^MCA_|\\.h5ad$", "", basename(f))
  rds_path <- file.path("organized_by_tissue", tissue, sprintf("MCA_MCA_%s.rds", tissue))
  convert_and_save(f, rds_path, "MCA")
}

log_msg(sprintf("\nConversion complete. See %s for details.", log_file))
sessionInfo()
