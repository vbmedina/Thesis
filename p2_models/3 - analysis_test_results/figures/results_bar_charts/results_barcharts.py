"""
Make Sexual vs Asexual panels from a single wide CSV using seaborn Reds palette.

Inputs
- p2_models/3 - analysis_test_results/final_test_metrics.csv
  Expected wide format with duplicate columns (e.g., 'Split', 'Split.1', ...).

Outputs
- p2_models/analysis_outputs/figures/
    panel_hit_at_100.png
    panel_ef_at_100.png
    panel_rmse.png
    panel_spearman.png
    panel_delta_hit_at_100.png
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
INPUT_CSV = Path("p2_models/3 - analysis_test_results/final_test_metrics.csv")
OUT_DIR = Path("p2_models/3 - analysis_test_results/figures/results_bar_charts")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Visual style
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.linestyle": "--",
    "grid.alpha": 0.25,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# Requested brand colors
COLOR_ASEX = "#eca38c"   # asexual
COLOR_SEX  = "#b12d31"   # sexual
COLOR_GAP  = "#6e1a1d"   # gap/delta

# Split order & display labels (as requested)
SPLIT_ORDER = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
SPLIT_DISPLAY = {
    "random": "Random (stratified)",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (k-means)",
    "umap_ward": "UMAP (ward)",
}

EXACT_SPLIT_MAP = {
    "random (stratified)": "random",
    "random": "random",
    "scaffold": "scaffold",
    "butina": "butina",
    "umap (k-means)": "umap_kmeans",
    "umap k-means": "umap_kmeans",
    "umap (kmeans)": "umap_kmeans",
    "umap (ward)": "umap_ward",
}

def canon_split(val: str) -> str:
    if pd.isna(val):
        return None
    s = str(val).strip().lower()
    if s in EXACT_SPLIT_MAP:
        return EXACT_SPLIT_MAP[s]
    # fallback heuristics
    s2 = re.sub(r"[^a-z0-9 ]+", " ", s)
    if "random" in s2: return "random"
    if "scaffold" in s2: return "scaffold"
    if "butina" in s2: return "butina"
    if "kmeans" in s2 or "k means" in s2: return "umap_kmeans"
    if "ward" in s2: return "umap_ward"
    return None

def norm_dataset(v: str) -> str:
    s = str(v).lower()
    if "sex" in s or "external" in s or "ood" in s:
        return "sexual_data"
    return "fold_test"

# Load & reshape (stack each ".N" block)
if not INPUT_CSV.exists():
    raise FileNotFoundError(f"Cannot find: {INPUT_CSV}")

df_wide = pd.read_csv(INPUT_CSV)

base_names = ["Split", "Fold", "Dataset", "RMSE", "Spearman", "Hit Rate at 100", "EF at 100"]
suffixes = {""}
pat = re.compile(r"^(.*?)(\.\d+)$")
for c in df_wide.columns:
    m = pat.search(c)
    if m and m.group(1) in base_names:
        suffixes.add(m.group(2))

def extract_block(sfx: str):
    cols = {}
    for base in base_names:
        name = base + sfx if sfx else base
        if name in df_wide.columns:
            cols[base] = df_wide[name]

    if "Split" not in cols or "Dataset" not in cols:
        return None
    b = pd.DataFrame(cols)
    b = b.rename(columns={
        "Split": "split_raw",
        "Fold": "fold",
        "Dataset": "dataset_raw",
        "Hit Rate at 100": "hit_at_k",
        "EF at 100": "ef_at_k",
        "RMSE": "rmse",
        "Spearman": "spearman_rho",
    })
    # Normalize labels
    b["dataset"] = b["dataset_raw"].apply(norm_dataset)
    b["split_key"] = b["split_raw"].apply(canon_split)

    keep = ["split_key", "dataset", "fold", "hit_at_k", "ef_at_k", "rmse", "spearman_rho"]
    keep = [c for c in keep if c in b.columns]
    b = b[keep]

    metric_cols = [c for c in ["hit_at_k","ef_at_k","rmse","spearman_rho"] if c in b.columns]
    b = b[~b["split_key"].isna()]
    b = b[~b[metric_cols].isna().all(axis=1)]
    return b

# If not CSV found
blocks = []
for sfx in sorted(suffixes):
    blk = extract_block(sfx)
    if blk is not None and not blk.empty:
        blocks.append(blk)

if not blocks:
    raise SystemExit("No valid metric blocks found in the CSV.")

tidy = pd.concat(blocks, ignore_index=True)

# Aggregate mean and SD across folds/runs
def agg_mean_sd(tidy_df, value_col):
    return (tidy_df.groupby(["split_key", "dataset"], as_index=False)
                 .agg(val=(value_col, "mean"), sd=(value_col, "std")))

agg = {}
if "hit_at_k" in tidy.columns:      agg["hit"] = agg_mean_sd(tidy, "hit_at_k")
if "ef_at_k" in tidy.columns:       agg["ef"]  = agg_mean_sd(tidy, "ef_at_k")
if "rmse" in tidy.columns:          agg["rmse"] = agg_mean_sd(tidy, "rmse")
if "spearman_rho" in tidy.columns:  agg["rho"]  = agg_mean_sd(tidy, "spearman_rho")

# Ensure we have at least one plot to make
if not agg:
    raise SystemExit("No metrics found to plot (missing Hit/EF/RMSE/Spearman).")

# Ensure splits exist in order, but keep only those present
available = set()
for g in agg.values():
    available |= set(g["split_key"].unique())
SPLITS = [s for s in SPLIT_ORDER if s in available]
if not SPLITS:
    raise SystemExit("No recognized splits found in the data.")

# Plotting helpers (bars with error bars; legend on the right)
def grouped_bar_figure(metric_df, ylabel, title, out_path, ylim=None, zero_line=False):
    df_fold = metric_df[metric_df["dataset"]=="fold_test"].set_index("split_key").reindex(SPLITS)
    df_sex  = metric_df[metric_df["dataset"]=="sexual_data"].set_index("split_key").reindex(SPLITS)

    vals_a = np.nan_to_num(df_fold["val"].to_numpy(), nan=0.0)
    sds_a  = np.nan_to_num(df_fold["sd"].to_numpy(),  nan=0.0)
    vals_b = np.nan_to_num(df_sex["val"].to_numpy(),  nan=0.0)
    sds_b  = np.nan_to_num(df_sex["sd"].to_numpy(),   nan=0.0)

    x = np.arange(len(SPLITS))
    bar_w = 0.38
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.bar(x - bar_w/2, vals_a, width=bar_w, yerr=sds_a, capsize=3,
           color=COLOR_ASEX, label="Asexual")
    ax.bar(x + bar_w/2, vals_b, width=bar_w, yerr=sds_b, capsize=3,
           color=COLOR_SEX,  label="Sexual")

    ax.set_xticks(x)
    ax.set_xticklabels([SPLIT_DISPLAY[s] for s in SPLITS], rotation=15, ha="right")
    ax.set_xlabel("Split", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=16, pad=12)
    if ylim is not None:
        ax.set_ylim(ylim)
    if zero_line:
        ax.axhline(0.0, color="black", lw=1.2, zorder=3)

    # Legend to the right
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, ncol=1)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    print(f"Saved {out_path}")

def delta_bar_figure(df_hit_agg, ylabel, title, out_path):
    sex  = df_hit_agg[df_hit_agg["dataset"]=="sexual_data"].set_index("split_key").reindex(SPLITS)["val"]
    asex = df_hit_agg[df_hit_agg["dataset"]=="fold_test"].set_index("split_key").reindex(SPLITS)["val"]
    delta_vals = (sex - asex).fillna(0.0).to_numpy()

    x = np.arange(len(SPLITS))
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.bar(x, delta_vals, color=COLOR_GAP)
    ax.axhline(0.0, color="#666", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([SPLIT_DISPLAY[s] for s in SPLITS], rotation=15, ha="right")
    ax.set_xlabel("Split", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=16, pad=12)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    print(f"Saved {out_path}")

# Make panels

# 1) HitRate@100 (mean ± SD)
if "hit" in agg:
    grouped_bar_figure(
        metric_df=agg["hit"],
        ylabel="Hit Rate on Top 100 Molecules ",
        title="Cross-Stage Test Results: Hit Rate on the Top 100 Molecules",
        out_path=OUT_DIR / "panel_hit_at_100.png",
        ylim=(0,1)
    )

# 2) EF@100 (mean ± SD)
if "ef" in agg:
    grouped_bar_figure(
        metric_df=agg["ef"],
        ylabel="EF on Top 100 Molecules",
        title="Cross-Stage Test Results: Enrichment Factor on the Top 100 Molecules",
        out_path=OUT_DIR / "panel_ef_at_100.png",
        ylim=None
    )

# 3) RMSE (pIC50) (mean ± SD)
if "rmse" in agg:
    grouped_bar_figure(
        metric_df=agg["rmse"],
        ylabel="RMSE (pIC50)",
        title="Cross-Stage Test Results: Root Mean Squared Error",
        out_path=OUT_DIR / "panel_rmse.png",
    )

# 4) Spearman (rho) (mean ± SD), with zero reference line
if "rho" in agg:
    grouped_bar_figure(
        metric_df=agg["rho"],
        ylabel="Spearman Correlation",
        title="Cross-Stage Test Results: Rank Quality via Spearman Correlation",
        out_path=OUT_DIR / "panel_spearman.png",
        ylim=(-1, 1),
        zero_line=True
    )

# # 5) ΔHit@100 (Sexual - Asexual)
# if "hit" in agg:
#     delta_bar_figure(
#         df_hit_agg=agg["hit"],
#         ylabel="Change in Hit Rate at 100 Sexual between Asexual Tests",
#         title="Change in Hit Rate at 100 between Sexual and Asexual Test by Split",
#         out_path=OUT_DIR / "panel_delta_hit_at_100.png",
#     )

print("All done.")