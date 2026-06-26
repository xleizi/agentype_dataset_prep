#!/usr/bin/env python3
"""
Pre-process HCL and MCA raw data: merge cell-type metadata from Excel/CSV into h5ad.

This is a prerequisite step before tissue splitting (scripts 01 and 02).
For HCL: reads HCL_Fig1_adata.h5ad + HCL_Fig1_cell_Info.xlsx
For MCA: reads MCA_BatchRemoved_Merge_dge.h5ad + MCA_CellAssignments.csv

Output:
  02_HCL/HCL_with_metadata.h5ad
  02_HCL/MCA_with_metadata.h5ad

Both output files contain obs columns: Tissue, celltype, batch, etc.
"""
import anndata as ad
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent


def process_hcl():
    """Merge HCL expression data with cell annotations from Excel."""
    h5ad_path = BASE / "02_HCL/HCL_Fig1_adata.h5ad"
    xlsx_path = BASE / "02_HCL/HCL_Fig1_cell_Info.xlsx"
    out_path  = BASE / "02_HCL/HCL_with_metadata.h5ad"

    print(f"[HCL] Loading {h5ad_path.name} ...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    print(f"[HCL] Merging {xlsx_path.name} ...")
    meta = pd.read_excel(xlsx_path, index_col=0)
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    for col in meta.columns:
        adata.obs[col] = meta.loc[common, col].values
    print(f"  Merged {len(common):,} cells with {len(meta.columns)} annotation columns")

    adata.write_h5ad(out_path)
    print(f"[HCL] Saved → {out_path.name}\n")


def process_mca():
    """Merge MCA expression data with cell assignments from CSV."""
    h5ad_path = BASE / "03_MCA/MCA_BatchRemoved_Merge_dge.h5ad"
    csv_path  = BASE / "03_MCA/MCA_CellAssignments.csv"
    out_path  = BASE / "02_HCL/MCA_with_metadata.h5ad"

    print(f"[MCA] Loading {h5ad_path.name} ...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    print(f"[MCA] Merging {csv_path.name} ...")
    meta = pd.read_csv(csv_path, index_col=0)
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    for col in meta.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = meta.loc[common, col]

    # MCA_CellAssignments.csv has columns: Batch, Tissue, Cell_type (or Annotation)
    # Rename "Cell_type" or "Annotation" → "celltype" for consistency
    for col in ["Cell_type", "Annotation"]:
        if col in adata.obs.columns:
            adata.obs["celltype"] = adata.obs[col]
            break

    print(f"  Merged {len(common):,} cells, columns: {adata.obs.columns.tolist()}")
    adata.write_h5ad(out_path)
    print(f"[MCA] Saved → {out_path.name}\n")


if __name__ == "__main__":
    process_hcl()
    process_mca()
    print("Done.")
