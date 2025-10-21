''' Description: This script generates violin and box plots to visualize the distribution of maximum Tanimoto similarities
between test molecules and their nearest neighbors in the training set across sexual and asexual data.

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

# Base directories for split/fold analysis between sexual and asexual data
PROJECT   = Path("./p1_preprocessing/4 - Data splitting")
DATA_ROOT = PROJECT
BASE_DIRS = [DATA_ROOT / "2 - split_data"]
SPLITS    = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]

SEXUALDATA = "./p0_all_csvs/postphase3_Asexual_Only.csv"

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
EXCLUDE_IDENTICAL = True  # drop exact SMILES matches when computing nearest neighbor

# Output for split/fold figure + stats
FIG_DIR   = DATA_ROOT / "./3 - splits_data_visualization/tanimoto_violin_plots"
OUT_PNG   = FIG_DIR / f"max_tanimoto_violins_box_{N_BITS}.png"
STATS_CSV = FIG_DIR / f"max_tanimoto_stats_{N_BITS}.csv"

# Explicit output paths you requested:
SEX_OUT_DIR = Path("/Users/victoriamedina/Thesis_Project/thesis/p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/tanimoto_violin_plots")
SEX_OUT_PNG = SEX_OUT_DIR / f"max_tanimoto_violins_box_{N_BITS}_sexual.png"
SEX_STATS_CSV  = SEX_OUT_DIR / f"max_tanimoto_stats_{N_BITS}_sexual.csv"
SEX_VALUES_CSV = SEX_OUT_DIR / f"max_tanimoto_values_{N_BITS}_sexual.csv"  # per-molecule output for traceability

# Fingerprints
try:
    FP_GEN = FPG.GetMorganGenerator(
        radius=int(RADIUS),
        includeChirality=bool(USE_CHIRALITY),
        fpSize=int(N_BITS)
    )
    def ecfp4_bits(smi: str):
        m = Chem.MolFromSmiles(smi)
        return FP_GEN.GetFingerprint(m) if m is not None else None
except Exception as e:
    print(f"[WARNING] MorganGenerator unavailable ({e}). Falling back to AllChem.GetMorganFingerprintAsBitVect.")
    def ecfp4_bits(smi: str):
        m = Chem.MolFromSmiles(smi)
        return AllChem.GetMorganFingerprintAsBitVect(m, int(RADIUS), nBits=int(N_BITS)) if m is not None else None

# ============ Helpers ============
def canon(s: str) -> str | None:
    m = Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m, canonical=True) if m is not None else None

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

# ============ ORIGINAL: Split/Fold main ============
def main_splits():
    rows = []
    fp_cache: dict[str, object] = {}

    for sp in SPLITS:
        print(f"Processing split: {sp}")
        for fold in FOLDS:
            print(f"  Fold {fold}")
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

            train_fps = []
            for s in train_smis:
                if s not in fp_cache:
                    fp_cache[s] = ecfp4_bits(s)
                train_fps.append(fp_cache[s])

            for s in test_smis:
                if s not in fp_cache:
                    fp_cache[s] = ecfp4_bits(s)
                sims = DataStructs.BulkTanimotoSimilarity(fp_cache[s], train_fps)
                if EXCLUDE_IDENTICAL:
                    sims = [val for val, tr_s in zip(sims, train_smis) if tr_s != s]
                if sims:
                    rows.append({"method": sp, "fold": fold, "max_tani": float(np.max(sims))})
    
    if not rows:
        print("No values computed for split/fold analysis. Check BASE_DIRS and file names.")
        return

    df = pd.DataFrame(rows)
    df["max_tani"] = df["max_tani"].clip(0, 1)

    stats = (df.groupby("method")["max_tani"]
             .agg(p05=lambda x: x.quantile(0.05),
                  q1 =lambda x: x.quantile(0.25),
                  median="median",
                  q3 =lambda x: x.quantile(0.75),
                  p95=lambda x: x.quantile(0.95),
                  n="count")
             .reset_index())

    print("\nPer-split summary (p05 / q1 / median / q3 / p95 / n):")
    print(stats[["method", "n", "p05", "q1", "median", "q3", "p95"]].to_string(
        index=False,
        float_format=lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
    ))

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    stats.to_csv(STATS_CSV, index=False)
    print(f"Saved stats to {STATS_CSV}")

    order = [m for m in SPLITS if m in df["method"].unique()]
    reds_cmap = cm.get_cmap("Reds")
    colors = [reds_cmap(x) for x in np.linspace(0.35, 0.90, len(order))]
    PALETTE = dict(zip(order, colors))

    sns.set_style("whitegrid")
    sns.despine()

    fig, ax = plt.subplots(figsize=(12, 6))

    sns.violinplot(
        data=df, x="method", y="max_tani", order=order,
        inner=None, cut=0, bw_adjust=0.7, linewidth=1.0,
        palette=PALETTE, ax=ax
    )

    sns.boxplot(
        data=df, x="method", y="max_tani", order=order,
        width=0.26, showcaps=True, showfliers=False, whis=(5, 95),
        boxprops=dict(facecolor="white", edgecolor="black", alpha=0.85, linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.5),
        ax=ax
    )

    ax.set_xticklabels([DISPLAY_NAME.get(m, m) for m in order], rotation=0, fontsize=13)
    ax.set_xlabel("Splitting Method", fontsize=15)
    ax.set_ylabel("Max Tanimoto Similarity to Training Set", fontsize=15)
    ax.set_ylim(0, 1.2)
    ax.set_title("Train-Test Tanimoto by Split", pad=12, fontsize=18)

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    print(f"\nSaved figure - {OUT_PNG}")

# ============ NEW: Sexual-only (no folds) ============
def main_sexual():
    # Load and canonicalize
    df = pd.read_csv(SEXUALDATA)
    if SMI_COL not in df.columns:
        raise ValueError(f"CSV '{SEXUALDATA}' missing '{SMI_COL}' column.")

    sexual_smis = [canon(s) for s in df[SMI_COL].astype(str)]
    sexual_smis = [s for s in sexual_smis if s is not None]

    # De-duplicate canonical SMILES within the set
    sexual_smis = list(dict.fromkeys(sexual_smis))
    print(f"Sexual unique SMILES: n={len(sexual_smis)}")

    # Build fingerprints once
    fp_cache: dict[str, object] = {}
    fps = []
    for s in sexual_smis:
        if s not in fp_cache:
            fp_cache[s] = ecfp4_bits(s)
        if fp_cache[s] is not None:
            fps.append(fp_cache[s])
        else:
            print(f"[WARN] Could not fingerprint: {s}")
    # Map index -> SMILES for exclusion logic
    idx2smi = {i: s for i, s in enumerate(sexual_smis)}

    # For each SMILES, compute max similarity to any *other* SMILES in the same set
    rows = []
    for i, s in enumerate(sexual_smis):
        fp_i = fp_cache.get(s)
        if fp_i is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp_i, fps)

        # Exclude self by comparing canonical string (and/or index alignment)
        # Here, filter by string to be robust if ordering changed
        sims_excl = []
        for j, (sim, smj) in enumerate(zip(sims, sexual_smis)):
            if EXCLUDE_IDENTICAL and smj == s:
                continue
            sims_excl.append(sim)

        if sims_excl:
            rows.append({"dataset": "Sexual", "smiles": s, "max_tani": float(np.max(sims_excl))})

    if not rows:
        raise SystemExit("No values computed for Sexual set. Check the CSV and SMILES parsing.")

    dvals = pd.DataFrame(rows)
    dvals["max_tani"] = dvals["max_tani"].clip(0, 1)

    # Stats
    stats = (dvals.groupby("dataset")["max_tani"]
            .agg(p05=lambda x: x.quantile(0.05),
                 q1 =lambda x: x.quantile(0.25),
                 median="median",
                 q3 =lambda x: x.quantile(0.75),
                 p95=lambda x: x.quantile(0.95),
                 n="count")
            .reset_index())

    # Save outputs
    SEX_OUT_DIR.mkdir(parents=True, exist_ok=True)
    dvals.to_csv(SEX_VALUES_CSV, index=False)
    stats.to_csv(SEX_STATS_CSV, index=False)
    print(f"Saved Sexual per-molecule values to: {SEX_VALUES_CSV}")
    print(f"Saved Sexual stats to: {SEX_STATS_CSV}")

    # Plot
    sns.set_style("whitegrid")
    sns.despine()

    order = ["Sexual"]
    reds_cmap = cm.get_cmap("Reds")
    PALETTE = {"Sexual": reds_cmap(0.65)}

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.violinplot(
        data=dvals, x="dataset", y="max_tani", order=order,
        inner=None, cut=0, bw_adjust=0.7, linewidth=1.0,
        palette=PALETTE, ax=ax
    )

    sns.boxplot(
        data=dvals, x="dataset", y="max_tani", order=order,
        width=0.32, showcaps=True, showfliers=False, whis=(5, 95),
        boxprops=dict(facecolor="white", edgecolor="black", alpha=0.9, linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.5),
        ax=ax
    )

    # Labeling
    ax.set_xlabel("")
    ax.set_ylabel("Tanimoto Similarity in Sexually Tested Molecules", fontsize=13)
    ax.set_ylim(0, 1.2)
    ax.set_title(f"Sexual Tested Molecules: Max In-Set Neighbor (ECFP4, {N_BITS} bits)", pad=10, fontsize=16)

    # Annotate percentiles & median
    r = stats.iloc[0]
    p05 = float(r["p05"]); q1 = float(r["q1"]); med = float(r["median"]); q3 = float(r["q3"]); p95 = float(r["p95"])
    x = 0
    ax.hlines(med, x-0.28, x-0.08, colors="black", linewidth=1.3); ax.text(x-0.35, med, f"{med:.2f}", ha="right", va="center", fontsize=10, fontweight="bold")
    ax.hlines(p95, x+0.08, x+0.28, colors="black", linewidth=1.3); ax.text(x+0.35, p95, f"{p95:.2f}", ha="left", va="center", fontsize=10, fontweight="bold")
    ax.hlines(p05, x+0.08, x+0.28, colors="black", linewidth=1.3); ax.text(x+0.35, p05, f"{p05:.2f}", ha="left", va="center", fontsize=10, fontweight="bold")
    ax.hlines(q1,  x-0.08, x-0.28, colors="black", linewidth=1.3); ax.text(x-0.35, q1,  f"{q1:.2f}",  ha="right", va="center", fontsize=10, fontweight="bold")
    ax.hlines(q3,  x-0.08, x-0.28, colors="black", linewidth=1.3); ax.text(x-0.35, q3,  f"{q3:.2f}",  ha="right", va="center", fontsize=10, fontweight="bold")

    fig.tight_layout()
    fig.savefig(SEX_OUT_PNG, dpi=300)
    print(f"Saved Sexual figure: {SEX_OUT_PNG}")

# ============ Entrypoint ============
if __name__ == "__main__":
    # 1) Run your original split/fold analysis (unchanged behavior)
    main_splits()
    # 2) Additionally produce the Sexual-only violin/statistics (no folds)
    main_sexual()
