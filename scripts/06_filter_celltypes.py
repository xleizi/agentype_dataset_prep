#!/usr/bin/env python3
"""
Filter low-abundance cell types from benchmark Seurat RDS files.

Logic:
  Reads each RDS file, computes cell counts per cell type, removes types with
  fewer than MIN_CELLS cells, and saves the filtered RDS. Parallelized across
  N_CORES workers.

Input:
  organized_by_tissue/          (112 RDS files across 53 tissue directories)

Output:
  organized_by_tissue_filtered/ (106 RDS files; 6 files with only 1 remaining
                                 cell type after filtering are excluded)

Parameters (from filter_celltype_log.txt):
  MIN_CELLS = 100
  N_CORES   = 12
"""

import os
import subprocess
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE     = Path(__file__).resolve().parent.parent
SRC_DIR  = BASE / "organized_by_tissue"
OUT_DIR  = BASE / "organized_by_tissue_filtered"
MIN_CELLS = 100
N_CORES   = 12

R_FILTER_SCRIPT = """
suppressMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]
outfile <- args[2]
min_cells <- as.integer(args[3])

obj <- readRDS(infile)

# Detect the cell-type column name (celltype / cell_type / Cell_type / cell_ontology_class / free_annotation)
ct_cols <- intersect(c("celltype", "cell_type", "Cell_type", "cell_ontology_class", "free_annotation"),
                     colnames(obj@meta.data))
if (length(ct_cols) == 0) {
  stop(sprintf("No known cell-type column in %s. Columns: %s",
               infile, paste(colnames(obj@meta.data), collapse=", ")))
}
ct_col <- ct_cols[1]

ct_counts <- table(obj@meta.data[[ct_col]])
keep_types <- names(ct_counts[ct_counts >= min_cells])
n_before <- ncol(obj)
n_keep   <- sum(obj@meta.data[[ct_col]] %in% keep_types)
n_rm     <- n_before - n_keep
removed  <- names(ct_counts[ct_counts < min_cells])

if (length(keep_types) <= 1) {
  # If only one cell type remains, skip — not informative for benchmarking
  cat(sprintf("SKIP: only %d cell type(s) remaining after filter (%s)\\n",
              length(keep_types), paste(keep_types, collapse=", ")))
  quit(status = 2)
}

obj_filt <- subset(obj, subset = !!sym(ct_col) %in% keep_types)
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
saveRDS(obj_filt, outfile)

cat(sprintf("OK: %d → %d cells (%d removed, %d types kept)\\n",
            n_before, n_keep, n_rm, length(keep_types)))
if (length(removed) > 0) {
  cat("  Removed types:", paste(sprintf("%s (n=%d)", removed, ct_counts[removed]), collapse=", "), "\\n")
}
"""


def log(msg: str):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
    line = f"[{ts}] {msg}"
    print(line)


def scan_rds_files(root: Path) -> list[Path]:
    """Return sorted list of .rds files under root."""
    return sorted(root.rglob("*.rds"))


def filter_one(rds_path: Path, script_path: Path) -> dict:
    """
    Filter a single RDS file by launching an R subprocess.
    Returns a status dict.
    """
    rel    = rds_path.relative_to(SRC_DIR)
    out    = OUT_DIR / rel
    result = {"file": str(rel), "status": "unknown", "input_size_mb": rds_path.stat().st_size / 1e6}

    t1 = time.time()
    proc = subprocess.run(
        ["Rscript", str(script_path), str(rds_path), str(out), str(MIN_CELLS)],
        capture_output=True, text=True, timeout=600,
    )
    result["duration_s"] = round(time.time() - t1, 1)

    if proc.returncode == 2:
        result["status"] = "skipped"
        result["reason"] = proc.stdout.strip()
    elif proc.returncode == 0:
        result["status"] = "success"
        result["output_size_mb"] = out.stat().st_size / 1e6
    else:
        result["status"] = "error"
        result["error"] = proc.stderr.strip()[:200]

    # Parse log details
    for line in proc.stdout.splitlines():
        if "OK:" in line:
            result["details"] = line.strip()
        elif "Removed types:" in line:
            result["removed"] = line.strip()

    return result


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Write the R filter script to a temp file
    r_script = BASE / "scripts/_filter_one.R"
    r_script.write_text(R_FILTER_SCRIPT)

    rds_files = scan_rds_files(SRC_DIR)
    log(f"Found {len(rds_files)} RDS files in {SRC_DIR}")
    log(f"Output dir: {OUT_DIR}")
    log(f"Threshold: cell type < {MIN_CELLS} cells → removed")
    log(f"Parallel workers: {N_CORES}")

    t0 = time.time()
    results = []
    n_success = n_skipped = n_error = 0

    with ProcessPoolExecutor(max_workers=N_CORES) as ex:
        futures = {ex.submit(filter_one, f, r_script): f for f in rds_files}
        for future in as_completed(futures):
            r = future.result()
            results.append(r)
            status = r["status"]
            if status == "success":
                n_success += 1
            elif status == "skipped":
                n_skipped += 1
            else:
                n_error += 1
            log(f"  [{status:8s}] {r['file']}")

    total = time.time() - t0

    log(f"\n{'='*50}")
    log(f"Summary: {len(rds_files)} total  |  {n_success} success  |  "
        f"{n_skipped} skipped  |  {n_error} error")
    log(f"Output RDS files: {n_success}")
    log(f"Elapsed: {total/60:.1f} min")

    # Write summary JSON
    import json
    summary_path = OUT_DIR / "filter_summary.json"
    summary_path.write_text(json.dumps({
        "input_dir": str(SRC_DIR),
        "output_dir": str(OUT_DIR),
        "min_cells": MIN_CELLS,
        "n_cores": N_CORES,
        "total_input": len(rds_files),
        "n_success": n_success,
        "n_skipped": n_skipped,
        "n_error": n_error,
        "elapsed_min": round(total / 60, 1),
    }, indent=2))


if __name__ == "__main__":
    main()
