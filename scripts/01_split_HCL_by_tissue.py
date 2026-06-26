#!/usr/bin/env python3
"""Split HCL (Human Cell Landscape) by major tissue type.

Input:
  02_HCL/HCL_Fig1_adata.h5ad      - expression data (599,926 cells x 27,341 genes)
  02_HCL/HCL_Fig1_cell_Info.xlsx  - cell annotations (Tissue, CellType columns)

Output:
  02_HCL/by_tissue/HCL_{tissue}.h5ad   (33 major tissue classes)

Reference: Han et al. Nature (2020). https://doi.org/10.1038/s41586-020-2157-4
"""
import anndata as ad
import pandas as pd
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

# 59 HCL sub-tissue labels (Excel "Tissue" column) → 33 anatomical major classes
TISSUE_MAP = {
    "AdultAdipose": "Adipose",
    "AdultAdrenalGland": "AdrenalGland", "FetalAdrenalGland": "AdrenalGland",
    "NeonatalAdrenalGland": "AdrenalGland",
    "AdultArtery": "Artery",
    "AdultBladder": "Bladder",
    "AdultPeripheralBlood": "Blood", "CordBlood": "Blood", "CordBloodCD34P": "Blood",
    "FetalCalvaria": "Bone", "FetalRib": "Bone",
    "AdultBoneMarrow": "BoneMarrow",
    "AdultCerebellum": "Brain", "AdultTemporalLobe": "Brain",
    "FetalBrain": "Brain", "FetalSpinalCord": "Brain",
    "AdultCervix": "Cervix",
    "AdultEsophagus": "Esophagus",
    "FetalEyes": "Eyes",
    "AdultFallopiantube": "Fallopiantube",
    "AdultGallbladder": "Gallbladder",
    "FetalFemaleGonad": "Gonad", "FetalMaleGonad": "Gonad",
    "AdultHeart": "Heart", "FetalHeart": "Heart",
    "hESC": "HESC",
    "AdultAscendingColon": "Intestine", "AdultDuodenum": "Intestine",
    "AdultEpityphlon": "Intestine", "AdultIleum": "Intestine",
    "AdultJejunum": "Intestine", "AdultRectum": "Intestine",
    "AdultSigmoidColon": "Intestine", "AdultTransverseColon": "Intestine",
    "FetalIntestine": "Intestine",
    "AdultKidney": "Kidney", "FetalKidney": "Kidney",
    "AdultLiver": "Liver", "FetalLiver": "Liver",
    "AdultLung": "Lung", "FetalLung": "Lung",
    "AdultMuscle": "Muscle", "FetalMuscle": "Muscle",
    "AdultOmentum": "Omentum",
    "AdultPancreas": "Pancreas", "FetalPancreas": "Pancreas", "NeonatalPancreas": "Pancreas",
    "Placenta": "Placenta", "ChorionicVillus": "Placenta",
    "AdultPleura": "Pleura",
    "AdultProstate": "Prostate",
    "AdultSkin": "Skin", "FetalSkin": "Skin",
    "AdultSpleen": "Spleen",
    "AdultThymus": "Thymus", "FetalThymus": "Thymus",
    "AdultThyroid": "Thyroid",
    "AdultTrachea": "Trachea",
    "AdultUreter": "Ureter",
    "AdultUterus": "Uterus",
}


def main():
    out_dir = BASE / "02_HCL/by_tissue"
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading HCL_Fig1_adata.h5ad ...")
    adata = ad.read_h5ad(BASE / "02_HCL/HCL_Fig1_adata.h5ad")
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    print("Merging HCL_Fig1_cell_Info.xlsx ...")
    meta = pd.read_excel(BASE / "02_HCL/HCL_Fig1_cell_Info.xlsx", index_col=0)
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    for col in meta.columns:
        adata.obs[col] = meta.loc[common, col].values

    adata.obs["pre_tissue"] = adata.obs["Tissue"].map(TISSUE_MAP)
    unmapped = adata.obs.loc[adata.obs["pre_tissue"].isna(), "Tissue"].unique()
    if len(unmapped):
        print(f"  WARNING: unmapped sub-tissues (excluded): {sorted(unmapped)}")

    tissues = sorted(adata.obs["pre_tissue"].dropna().unique())
    print(f"\nWriting {len(tissues)} tissue files → {out_dir}")
    rows = []
    for i, tissue in enumerate(tissues, 1):
        mask = adata.obs["pre_tissue"] == tissue
        out  = out_dir / f"HCL_{tissue}.h5ad"
        t1   = time.time()
        adata[mask].copy().write_h5ad(out)
        mb   = out.stat().st_size / 1e6
        print(f"  [{i:2d}/{len(tissues)}] {tissue:<22s} {int(mask.sum()):6,} cells  {mb:.0f} MB  ({time.time()-t1:.1f}s)")
        rows.append({"Tissue": tissue, "N_Cells": int(mask.sum()), "File_MB": round(mb, 1)})

    pd.DataFrame(rows).to_csv(out_dir / "split_summary.csv", index=False)
    print(f"\nCompleted in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
