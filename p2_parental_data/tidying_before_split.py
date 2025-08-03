from pathlib import Path
import pandas as pd

# -------- configuration ----------------------------------------------------
INPUT_CSV  = Path("./pp.csv")
OUTPUT_CSV = Path("./pp_tidy.csv")
ESSENTIALS = ["Smiles", "Strains", "pIC50"]   # columns that must be present
# ---------------------------------------------------------------------------


def main() -> None:
    # 1 Load ---------------------------------------------------------------
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV.resolve()}")
    df = pd.read_csv(INPUT_CSV)

    # 2 Drop rows missing essentials --------------------------------------
    df.dropna(subset=ESSENTIALS, inplace=True)

    # 3 Potency bins (label-aligned) --------------------------------------
    df["potency_bin"] = pd.cut(
        df["pIC50"],
        bins=[-1, 6.0, 7.5, float("inf")],
        labels=["inactive", "active", "elite"],
        right=False              # 6.0 goes into "active"
    )

    # 4 Frequency counts ---------------------------------------------------
    mol_freq = df.groupby("Molecule_ChEMBL_ID").size()
    str_freq = df.groupby("Strains").size()
    df["mol_freq"]    = df["Molecule_ChEMBL_ID"].map(mol_freq)
    df["strain_freq"] = df["Strains"].map(str_freq)

    # 5 Duplicate-pair check ----------------------------------------------
    dup_mask = df.duplicated(subset=["Molecule_ChEMBL_ID", "Strains"])
    dup_count = dup_mask.sum()
    if dup_count:
        print(f"⚠️  Warning: {dup_count} duplicated molecule–strain rows kept.")
        # Optional: uncomment next line to drop duplicates, keeping first
        # df = df[~dup_mask]

    # 6 Save ---------------------------------------------------------------
    df.to_csv(OUTPUT_CSV, index=False)

    # 7 Report -------------------------------------------------------------
    print(f"[√] Saved tidy dataframe ➜ {OUTPUT_CSV.resolve()}")
    print(f"Rows after cleaning: {len(df):,}\n")

    print("Potency-bin distribution:")
    print(df["potency_bin"].value_counts(dropna=False), "\n")

    print("Top-5 molecules by frequency:")
    print(mol_freq.sort_values(ascending=False).head(), "\n")

    print("Top-5 strains by frequency:")
    print(str_freq.sort_values(ascending=False).head())


if __name__ == "__main__":
    main()