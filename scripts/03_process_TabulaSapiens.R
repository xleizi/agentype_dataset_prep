#!/usr/bin/env Rscript
# Process Tabula Sapiens dataset and split by organ
#
# Input:
#   06_TabulaSapiens/exprMatrix.tsv  — sparse gene expression matrix (~50 GB)
#   06_TabulaSapiens/meta.tsv        — cell-level metadata (organ_tissue, cell_ontology_class, ...)
#
# Output:
#   06_TabulaSapiens/by_organ/TabulaSapiens_{organ}.rds  (24 organs)
#
# Reference: The Tabula Sapiens Consortium. Science (2022).

library(Seurat)
library(Matrix)
library(data.table)

BASE <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
setwd(BASE)

cat(sprintf("\n%s\n", paste(rep("=", 80), collapse="")))
cat("  Tabula Sapiens: Build Seurat Object & Split by Organ\n")
cat(sprintf("%s\n\n", paste(rep("=", 80), collapse="")))

# ── Step 1: load expression matrix ──────────────────────────────────────────
matrix_file <- "06_TabulaSapiens/exprMatrix.tsv"
meta_file   <- "06_TabulaSapiens/meta.tsv"
out_dir     <- "06_TabulaSapiens/by_organ"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

cat(sprintf("Step 1/3: Reading expression matrix (%s) ...\n", matrix_file))
cat("  Using data.table::fread() for speed (expect ~2-3 min) ...\n")
mat_dt <- fread(matrix_file, sep = "\t", header = TRUE, check.names = FALSE)
genes  <- mat_dt[[1]]
cells  <- colnames(mat_dt)[-1]
mat_dt <- mat_dt[, -1, with = FALSE]

cat(sprintf("  Converting to dgCMatrix (%d genes x %d cells) ...\n",
            nrow(mat_dt), ncol(mat_dt)))
counts <- as(as.matrix(mat_dt), "dgCMatrix")
rownames(counts) <- genes
colnames(counts) <- cells
rm(mat_dt)
gc()

# ── Step 2: build Seurat object ─────────────────────────────────────────────
cat(sprintf("\nStep 2/3: Loading metadata (%s) ...\n", meta_file))
meta <- fread(meta_file, sep = "\t", header = TRUE, data.table = FALSE)
rownames(meta) <- meta[, 1]

# Align cells between matrix and metadata
common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common]
meta   <- meta[common, ]

cat(sprintf("  Building Seurat object (%d cells) ...\n", length(common)))
sobj <- CreateSeuratObject(counts = counts, meta.data = meta, project = "TabulaSapiens")
rm(counts)
gc()

cat(sprintf("  Saving full Seurat object (%.1f M cells) ...\n", ncol(sobj) / 1e6))
saveRDS(sobj, "06_TabulaSapiens/TabulaSapiens_Seurat.rds")

# Parse cell ontology classes and add major cell types
source("06_TabulaSapiens/classify_cells.py", local = TRUE)

# ── Step 3: split by organ ──────────────────────────────────────────────────
cat(sprintf("\nStep 3/3: Splitting by organ_tissue ...\n"))
organ_col <- "organ_tissue"
organs <- sort(unique(sobj@meta.data[[organ_col]]))
cat(sprintf("  %d organs found\n", length(organs)))

summary_rows <- list()
for (i in seq_along(organs)) {
  organ <- organs[i]
  cat(sprintf("  [%2d/%d] %s ... ", i, length(organs), organ))
  t1 <- Sys.time()

  keep <- sobj@meta.data[[organ_col]] == organ
  sub <- sobj[, keep]
  out <- file.path(out_dir, sprintf("TabulaSapiens_%s.rds", gsub(" ", "_", organ)))
  saveRDS(sub, out)

  mb <- file.info(out)$size / 1e6
  elapsed <- difftime(Sys.time(), t1, units = "secs")
  cat(sprintf("%d cells  %.1f MB  (%.1f s)\n", ncol(sub), mb, elapsed))

  summary_rows[[organ]] <- data.frame(
    Organ = organ,
    N_Cells = ncol(sub),
    N_Genes = nrow(sub),
    File_MB = round(mb, 1),
    Output_File = basename(out),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(out_dir, "organ_split_summary.csv"), row.names = FALSE)

total <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\nCompleted in %.1f min — %d organ files written to %s\n",
            total, length(organs), out_dir))
cat("sessionInfo:\n")
sessionInfo()
