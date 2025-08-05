from pathlib import Path
import argparse, numpy as np, pandas as pd

# Paths
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TIDY = DATA / "tidy" / "pp_tidy.csv"
SPLITS = DATA / "splits"

# Table
def show(name, *frames):
    def cnt(d):
        return (d["potency_bin"].value_counts().reindex(["<=5.0","5.0-6.0","6.0-7.5",">=7.5"], fill_value=0))
    c = list(map(cnt, frames))
    tags = ["train","val","test"] if len(frames)==3 \
           else ["warm","cold_cmpd","cold_str","cold_both"]
    print(f"\n[{name}]")
    hdr = "  set      rows   <=5.0   5.0-6.0   6.0-7.5   >=7.5"
    print(hdr, "\n  " + "-"*(len(hdr)-2))
    for tag, df_, ct in zip(tags, frames, c):
        print(f"  {tag:<10}{len(df_):7,} {ct['<=5.0']:5} {ct['5.0-6.0']:5} " f"{ct['6.0-7.5']:5} {ct['>=7.5']:5}")

# Main function
if __name__ == "__main__":
    seed = 0
    rng  = np.random.default_rng(seed)

    # DF and column helper
    df = pd.read_csv(TIDY, low_memory=False)
    df["potency_bin"] = pd.cut(df.pIC50, [-1,5,6,7.5,99], labels=["<=5.0","5.0-6.0","6.0-7.5",">=7.5"])

    # LOSO strains split
    print("\n Building LOSO splits (≥ 200 rows per strain)…")

    df["strain_norm"] = df["Strains"].str.split(r"[ /]").str[0]
    MIN_ROWS = 200
    sz = df["strain_norm"].value_counts()
    big = sz[sz >= MIN_ROWS].index
    print(f"Retained {len(big)} strains "f"({sz[big].sum():,} / {len(df):,} rows)")

    # Create LOSO splits
    for strain in big:
        mask_test = df["strain_norm"] == strain
        test_idx  = df[mask_test].index
        train_idx = df.index.difference(test_idx)

        # Random 10 % of train for validation
        val_idx = rng.choice(train_idx, size=int(0.10*len(train_idx)), replace=False)

        trn_df = df.loc[train_idx.difference(val_idx)]
        val_df = df.loc[val_idx]
        tst_df = df.loc[test_idx]

        # Create output directory and save splits
        out = SPLITS / f"loso_{strain}"
        out.mkdir(parents=True, exist_ok=True)
        trn_df.to_csv(out/"train.csv")
        val_df.to_csv(out/"val.csv")
        tst_df.to_csv(out/"test.csv")

        show(f"loso_{strain}", trn_df, val_df, tst_df)

    print("\n LOSO splits saved to data/splits/loso_*")

    # Cold-regime split 
    cold_cmpd = rng.choice(df["Molecule_ChEMBL_ID"].unique(), size=int(0.10*df["Molecule_ChEMBL_ID"].nunique()), replace=False)
    cold_str  = rng.choice(df["strain_norm"].unique(), size=int(0.10*df["strain_norm"].nunique()), replace=False)

    # Create cold-regime splits
    mask_cmpd = df["Molecule_ChEMBL_ID"].isin(cold_cmpd)
    mask_str = df["strain_norm"].isin(cold_str)

    # Create warm, cold_cmpd, cold_str, and cold_both splits
    warm = df[~mask_cmpd & ~mask_str]
    cold_cmpd_ = df[ mask_cmpd & ~mask_str]
    cold_str_ = df[~mask_cmpd &  mask_str]
    cold_both = df[ mask_cmpd &  mask_str]

    # Show distribution of cold-regime splits
    out = SPLITS / "cold_split"
    out.mkdir (parents=True, exist_ok=True)
    warm.to_csv (out/"warm.csv", index=False)
    cold_cmpd_.to_csv (out/"cold_cmpd.csv", index=False)
    cold_str_.to_csv (out/"cold_str.csv", index=False)
    cold_both.to_csv (out/"cold_both.csv", index=False)

    show("cold_split", warm, cold_cmpd_, cold_str_, cold_both)
    print("\n cold_split CSVs ", out)
