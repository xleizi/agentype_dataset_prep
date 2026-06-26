#!/usr/bin/env python3
"""Split MCA (Mouse Cell Atlas) by major tissue type.

Input:
  02_HCL/MCA_with_metadata.h5ad   - MCA expression data with cell-type metadata merged
  03_MCA/MCA_CellAssignments.csv  - cell type assignments (used if metadata not yet merged)

Output:
  02_HCL/MCA_by_tissue/MCA_{tissue}.h5ad   (25 major tissue classes)

The raw MCA contains 47 fine-grained tissue labels (e.g. "AdultKidney", "FetalKidney",
"NeonatalKidney", "AdultBoneMarrowcKit"). This script merges them into 25 anatomical
major classes and exports one h5ad per class.

Reference: Guo et al. "The Mouse Organogenesis Cell Atlas." Nature (2020).
"""
import anndata as ad
import pandas as pd
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

# 47 MCA sub-tissue labels → 25 anatomical major classes
MCA_TISSUE_MAP = {
    "AdultAdrenalGland": "AdrenalGland",
    "AdultBladder": "Bladder",
    "AdultBoneMarrow": "BoneMarrow",
    "AdultBoneMarrowcKit": "BoneMarrow",
    "AdultBrain": "Brain",
    "FetalBrain": "Brain",
    "NeontalBrain": "Brain",
    "NeonatalCalvaria": "Bone",
    "NeonatalRib": "Bone",
    "AdultPeripheralBlood": "Blood",
    "EmbryonicMesenchyme": "EmbryonicMesenchyme",
    "AdultOvary": "Gonad",
    "AdultTestis": "Gonad",
    "FetalFemaleGonad": "Gonad",
    "FetalMaleGonad": "Gonad",
    "FetalHeart": "Heart",
    "NeonatalHeart": "Heart",
    "AdultSmallIntestine": "Intestine",
    "AdultStomach": "Intestine",
    "FetalIntestine": "Intestine",
    "FetalStomach": "Intestine",
    "AdultKidney": "Kidney",
    "FetalKidney": "Kidney",
    "AdultLiver": "Liver",
    "FetalLiver": "Liver",
    "AdultLung": "Lung",
    "FetalLung": "Lung",
    "MammaryGland.Involution": "MammaryGland",
    "MammaryGland.Lactation": "MammaryGland",
    "MammaryGland.Pregnancy": "MammaryGland",
    "MammaryGland.Virgin": "MammaryGland",
    "AdultMuscle": "Muscle",
    "NeonatalMuscle": "Muscle",
    "AdultOmentum": "Omentum",
    "AdultPancreas": "Pancreas",
    "FetalPancreas": "Pancreas",
    "NeonatalPancreas": "Pancreas",
    "Placenta": "Placenta",
    "AdultPleura": "Pleura",
    "AdultProstate": "Prostate",
    "NeonatalSkin": "Skin",
    "AdultSpleen": "Spleen",
    "CulturedMesenchymalStemCells": "StemCells",
    "EmbryonicStemCells": "StemCells",
    "PrimaryMesenchymalStemCells": "StemCells",
    "AdultThymus": "Thymus",
    "AdultUterus": "Uterus",
}


def main():
    in_h5ad = BASE / "02_HCL/MCA_with_metadata.h5ad"
    out_dir = BASE / "02_HCL/MCA_by_tissue"
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"Loading {in_h5ad.name} ...")
    adata = ad.read_h5ad(in_h5ad)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  obs columns: {adata.obs.columns.tolist()}")

    # If metadata (Tissue, celltype columns) not yet in h5ad, merge from CSV
    csv_path = BASE / "03_MCA/MCA_CellAssignments.csv"
    if "Tissue" not in adata.obs.columns and csv_path.exists():
        print(f"  Merging cell annotations from {csv_path.name} ...")
        meta = pd.read_csv(csv_path, index_col=0)
        common = adata.obs_names.intersection(meta.index)
        adata = adata[common].copy()
        for col in meta.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = meta.loc[common, col]

    # Map sub-tissue → major tissue class
    tissue_col = "Tissue"
    if tissue_col not in adata.obs.columns:
        raise KeyError(f"'{tissue_col}' not found in adata.obs. Available: {adata.obs.columns.tolist()}")
    adata.obs["pre_tissue"] = adata.obs[tissue_col].map(MCA_TISSUE_MAP)
    unmapped = adata.obs.loc[adata.obs["pre_tissue"].isna(), tissue_col].unique()
    if len(unmapped):
        print(f"  WARNING: unmapped sub-tissues (excluded): {sorted(unmapped)}")

    tissues = sorted(adata.obs["pre_tissue"].dropna().unique())
    print(f"\nWriting {len(tissues)} tissue files → {out_dir}")
    rows = []
    for i, tissue in enumerate(tissues, 1):
        mask = adata.obs["pre_tissue"] == tissue
        out  = out_dir / f"MCA_{tissue}.h5ad"
        t1   = time.time()
        adata[mask].copy().write_h5ad(out)
        mb   = out.stat().st_size / 1e6
        print(f"  [{i:2d}/{len(tissues)}] {tissue:<22s} {int(mask.sum()):6,} cells  {mb:.0f} MB  ({time.time()-t1:.1f}s)")
        rows.append({"Tissue": tissue, "N_Cells": int(mask.sum()), "File_MB": round(mb, 1)})

    pd.DataFrame(rows).to_csv(out_dir / "split_summary.csv", index=False)
    print(f"\nCompleted in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
