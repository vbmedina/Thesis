''' Description: This script generates violin and box plots to visualize the distribution of maximum Tanimoto similarities
between test molecules and their nearest neighbors in the training set across different data splitting methods and folds.

Preconditions:
1) "split_data" directory with subdirectories for each splitting method containing CSV files for folds 1-5 from step 2 
"Data splitting"
'''
# Conda activate molml 
# Imports
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator as FPG

RDLogger.DisableLog("rdApp.warning")

# Base directories
PROJECT   = Path("./p1_preprocessing/4 - Data splitting")
DATA_ROOT = PROJECT

BASE_DIRS = [DATA_ROOT / "2 - split_data"]
SPLITS    = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]

# Display names for nicer axis labels
DISPLAY_NAME = {
    "random": "Random (stratified)",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (k-means)",
    "umap_ward": "UMAP (ward)",
}

FOLDS     = [1, 2, 3, 4, 5]
SMI_COL   = "Smiles"

RADIUS = 2
N_BITS = 2048
USE_CHIRALITY = False

try:
    FP_GEN = FPG.GetMorganGenerator(
        radius=int(RADIUS),
        includeChirality=bool(USE_CHIRALITY),
        fpSize=int(N_BITS)
    )
    def ecfp4_bits(smi: str):
        m = Chem.MolFromSmiles(smi)
        return FP_GEN.GetFingerprint(m)
except Exception as e:
    print(f"[WARNING] MorganGenerator unavailable ({e}). Falling back to AllChem.GetMorganFingerprintAsBitVect.")
    def ecfp4_bits(smi: str):
        m = Chem.MolFromSmiles(smi)
        return AllChem.GetMorganFingerprintAsBitVect(m, int(RADIUS), nBits=int(N_BITS))

EXCLUDE_IDENTICAL = True  # drop exact SMILES matches when computing nearest train neighbor

# Output
FIG_DIR   = DATA_ROOT / "./3 - splits_data_visualization/tanimoto_violin_plots"
OUT_PNG   = FIG_DIR / f"max_tanimoto_violins_box_{N_BITS}.png"
STATS_CSV = FIG_DIR / f"max_tanimoto_stats_{N_BITS}.csv"

# Helpers - Return canonical SMILES, None if invalid
def canon(s: str) -> str | None:
    m = Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m, canonical=True) if m is not None else None

# Helpers - Find CSVs based on naming scheme
def find_csv(base_dirs, split: str, fold: int, which: str) -> Path | None:
    for root in base_dirs:
        d = root / split
        if not d.exists():
            continue
        for pat in (f"{split}_fold_{fold}_{which}.csv",
                    f"fold{fold}_{which}.csv",
                    f"{split}_fold{fold}_{which}.csv"):
            hits = list(d.glob(pat))
            if hits:
                return hits[0]
    return None

# Main loop
def main():
    rows = []
    fp_cache: dict[str, object] = {}

    # Build distributions on all splits/folds
    for sp in SPLITS:
        print(f"Processing split: {sp}")
        for fold in FOLDS:
            print(f"  Fold {fold}")
            # Find CSVs
            tr_path = find_csv(BASE_DIRS, sp, fold, "train")
            te_path = find_csv(BASE_DIRS, sp, fold, "test")
            if tr_path is None or te_path is None:
                print(f"[WARNING] Missing CSV for {sp} fold {fold}; skipping.")
                continue

            tr = pd.read_csv(tr_path)
            te = pd.read_csv(te_path)
            if SMI_COL not in tr or SMI_COL not in te:
                raise ValueError(f"CSV missing '{SMI_COL}': {tr_path}/{te_path}")

            train_smis = [canon(s) for s in tr[SMI_COL].astype(str)]
            test_smis  = [canon(s) for s in te[SMI_COL].astype(str)]
            train_smis = [s for s in train_smis if s is not None]
            test_smis  = [s for s in test_smis  if s is not None]

            # Cache train fps
            train_fps = []
            for s in train_smis:
                if s not in fp_cache:
                    fp_cache[s] = ecfp4_bits(s)
                train_fps.append(fp_cache[s])

            # Compute max similarity per test SMILES
            for s in test_smis:
                if s not in fp_cache:
                    fp_cache[s] = ecfp4_bits(s)
                sims = DataStructs.BulkTanimotoSimilarity(fp_cache[s], train_fps)
                if EXCLUDE_IDENTICAL:
                    sims = [val for val, tr_s in zip(sims, train_smis) if tr_s != s]
                if sims:
                    rows.append({"method": sp, "fold": fold, "max_tani": float(np.max(sims))})
    
    # If no values found
    if not rows:
        raise SystemExit("No values computed. Check BASE_DIRS and file names.")

    df = pd.DataFrame(rows)
    df["max_tani"] = df["max_tani"].clip(0, 1)

    # Stats for annotation
    stats = (df.groupby("method")["max_tani"]
             .agg(p05=lambda x: x.quantile(0.05),
                  q1 =lambda x: x.quantile(0.25),
                  median="median",
                  q3 =lambda x: x.quantile(0.75),
                  p95=lambda x: x.quantile(0.95),
                  n="count")
             .reset_index())

    # Print summary (of q1, q3, median, p05, and p95)
    print("\nPer-split summary (p05 / q1 / median / q3 / p95 / n):")
    print(stats[["method", "n", "p05", "q1", "median", "q3", "p95"]].to_string(
        index=False,
        float_format=lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
    ))

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    stats.to_csv(STATS_CSV, index=False)
    print(f"Saved stats to {STATS_CSV}")

    # Keep configured order
    order = [m for m in SPLITS if m in df["method"].unique()]

    # Reds palette
    reds_cmap = cm.get_cmap("Reds")
    colors = [reds_cmap(x) for x in np.linspace(0.35, 0.90, len(order))]
    PALETTE = dict(zip(order, colors))

    sns.set_style("whitegrid")
    sns.despine()

    # Create figure - figure size
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create figure - violin
    sns.violinplot(
        data=df, x="method", y="max_tani", order=order,
        inner=None, cut=0, bw_adjust=0.7, linewidth=1.0,
        palette=PALETTE, ax=ax
    )

    # Create figure - box plot
    sns.boxplot(
        data=df, x="method", y="max_tani", order=order,
        width=0.26, showcaps=True, showfliers=False, whis=(5, 95),
        boxprops=dict(facecolor="white", edgecolor="black", alpha=0.85, linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.5),
        ax=ax
    )

    # Use nicer method names on the x-axis
    ax.set_xticklabels([DISPLAY_NAME.get(m, m) for m in order], rotation=0, fontsize=13)

    # Labels
    ax.set_xlabel("Splitting Method", fontsize=15)
    ax.set_ylabel("Max Tanimoto Similarity to Training Set", fontsize=15)
    ax.set_ylim(0, 1.2)  # upper limit > 1 to leave space for text labels
    ax.set_title("Train-Test Tanimoto by Split", pad=12, fontsize = 18)

    # Annotate percentiles & median
    x_pos = {m: i for i, m in enumerate(order)}
    for _, r in stats.iterrows():
        m = r["method"]
        if m not in x_pos:
            continue
        x   = x_pos[m]
        p05 = float(r["p05"])
        q1  = float(r["q1"])
        med = float(r["median"])
        q3  = float(r["q3"])
        p95 = float(r["p95"])

        # Median
        ax.hlines(med, x-0.28, x-0.08, colors="black", linewidth=1.3)
        ax.text(x-0.35, med, f"{med:.2f}", ha="right", va="center", fontsize=10, fontweight="bold")

        # 5th/95th percentile
        ax.hlines(p95, x+0.08, x+0.28, colors="black", linewidth=1.3)
        ax.text(x+0.35, p95, f"{p95:.2f}", ha="left",  va="center", fontsize=10, fontweight="bold")
        ax.hlines(p05, x+0.08, x+0.28, colors="black", linewidth=1.3)
        ax.text(x+0.35, p05, f"{p05:.2f}", ha="left",  va="center", fontsize=10, fontweight="bold")

        # Q1 & Q3 labels
        ax.hlines(q1, x-0.08, x-0.28, colors="black", linewidth=1.3)
        ax.text(x-0.35, q1, f"{q1:.2f}", ha="right", va="center", fontsize=10, fontweight="bold")
        ax.hlines(q3, x-0.08, x-0.28, colors="black", linewidth=1.3)
        ax.text(x-0.35, q3, f"{q3:.2f}", ha="right", va="center",fontsize=10, fontweight="bold")

    # Out
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    print(f"\nSaved figure - {OUT_PNG}")

if __name__ == "__main__":
    main()