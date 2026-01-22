'''
Description: This script creates data bins for each fold using each splitting method.

Requirements:
    1) Splitting methods completed from 2 - split_data
'''

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#  onfig
BASE = Path("./p1_preprocessing/4 - Data splitting/2 - split_data")
METHODS = ["random","scaffold","butina","umap_kmeans","umap_ward"]
FOLDS   = [1,2,3,4,5]
SUBSETS = ["train", "val", "test"]
TARGET  = "pIC50"
OUTDIR  = Path("./p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/activity_bins_barchat")
OUTDIR.mkdir(parents=True, exist_ok=True)

DISPLAY_NAME = {
    "random": "Random (stratified)",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (k-means)",
    "umap_ward": "UMAP (ward)",
}

# bins & tick labels (like the paper)
BINS = [4.0, 5.0, 6.0, 7.0, np.inf]
XLABEL = ["4–5", "5–6", "6–7", ">7"]

def find_csv(method: str, fold: int, subset: str) -> Path | None:
    """Try common filename patterns; return None if missing."""
    candidates = [
        f"{method}_fold_{fold}_{subset}.csv",
        f"{method}_fold{fold}_{subset}.csv",
        f"{method}_fold_{fold}-{subset}.csv",
        # just in case someone used 'validation'
        f"{method}_fold_{fold}_validation.csv" if subset == "val" else None,
        f"{method}_fold{fold}_validation.csv"  if subset == "val" else None,
    ]
    for name in filter(None, candidates):
        p = BASE / method / name
        if p.exists():
            return p
    return None

def bin_counts(values: pd.Series) -> np.ndarray:
    return np.histogram(values.astype(float).values, bins=BINS)[0]

def make_panel_for_subset(subset: str):
    counts_map, totals_map, ymax = {}, {}, 0
    for m in METHODS:
        for f in FOLDS:
            path = find_csv(m, f, subset)
            if path is None:
                continue
            df = pd.read_csv(path)
            cnts = bin_counts(df[TARGET])
            counts_map[(m, f)] = cnts
            totals_map[(m, f)] = int(len(df))
            ymax = max(ymax, cnts.max())

    n_rows, n_cols = len(METHODS), len(FOLDS)
    fig_w, fig_h = 3.0 * n_cols, 2.6 * n_rows
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), sharey=True)
    if n_rows == 1:
        axes = np.array([axes])

    bar_color = "#d14444"   # red
    edge_color = "#c73838"

    for r, method in enumerate(METHODS):
        for c, fold in enumerate(FOLDS, start=1):
            ax = axes[r, c-1]
            key = (method, fold)
            if key not in counts_map:
                ax.axis("off")
                continue

            cnts = counts_map[key]
            ax.bar(range(len(cnts)), cnts, width=0.8, color=bar_color, edgecolor=edge_color)

            ax.set_xticks(range(len(XLABEL)))
            ax.set_xticklabels(XLABEL, fontsize=9)
            ax.set_ylim(0, ymax * 1.1)
            ax.grid(axis="y", alpha=0.25, linestyle=":")

            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

            if r == 0:
                ax.set_title(f"Fold {fold}", fontsize=11, pad=6)
            if c == 1:
                ax.text(-0.55, 0.5, DISPLAY_NAME.get(method, method), rotation=90, va="center", ha="right",
                    transform=ax.transAxes, fontsize=11)

            ax.text(0.02, 0.92, f"Total: {totals_map[key]:,}", transform=ax.transAxes,
                    ha="left", va="top", fontsize=10)

    fig.suptitle(f"{TARGET} Distribution by Fold for Each Split ({subset} sets)", fontsize=16, y=0.975)
    fig.text(0.5, 0.02, f"{TARGET} bins", ha="center", fontsize=12)
    fig.text(0.05, 0.5, "Number of molecules", va="center", rotation=90, fontsize=12)

    fig.tight_layout(rect=[0.03, 0.05, 1, 0.97])
    outpng = OUTDIR / f"pIC50_distribution_by_split_{subset}.png"
    fig.savefig(outpng, dpi=300)
    plt.close(fig)
    print(f"Saved: {outpng}")

if __name__ == "__main__":
    for subset in SUBSETS:
        make_panel_for_subset(subset)

'''Sexual activity bin graph
comment out when not using'''

# # Sexual activity bins (comment out when using the top)
# from pathlib import Path
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# # ---- paths ----
# SEXUAL_CSV = Path("p2_models/0 - input_data/sexual_test.csv")
# OUT_DIR    = Path("p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/activity_bins_barchat/sexual_pIC50_distribution.png")
# OUT_DIR.mkdir(parents=True, exist_ok=True)

# # ---- config ----
# TARGET  = "pIC50"                     # case-insensitive match
# BINS    = [4.0, 5.0, 6.0, 7.0, np.inf]
# XLABELS = ["4–5", "5–6", "6–7", ">7"]
# BAR_COLOR  = "#e32f26"                # red
# EDGE_COLOR = "#a21e18"

# # ---- load ----
# df = pd.read_csv(SEXUAL_CSV)

# # case-insensitive column lookup for pIC50
# cols_lc = {c.lower(): c for c in df.columns}
# if TARGET.lower() not in cols_lc:
#     raise ValueError(f"Column '{TARGET}' not found in {SEXUAL_CSV}. Got: {list(df.columns)}")
# col = cols_lc[TARGET.lower()]

# y = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy()
# N = len(y)
# if N == 0:
#     raise ValueError(f"No numeric {TARGET} values in {SEXUAL_CSV}")

# # ---- bin counts ----
# counts, edges = np.histogram(y, bins=BINS)
# props = counts / counts.sum()

# # save counts table
# counts_df = pd.DataFrame({
#     "bin_label": XLABELS,
#     "lower": edges[:-1],
#     "upper": edges[1:],
#     "count": counts,
#     "proportion": props
# })
# counts_csv = OUT_DIR / "sexual_pIC50_bin_counts.csv"
# counts_df.to_csv(counts_csv, index=False)

# # ---- plot ----
# x = np.arange(len(counts))
# fig, ax = plt.subplots(figsize=(7.5, 4.75))
# ax.bar(x, counts, width=0.8, color=BAR_COLOR, edgecolor=EDGE_COLOR)

# ax.set_xticks(x)
# ax.set_xticklabels(XLABELS, fontsize=10)
# ax.set_ylabel("Number of molecules", fontsize=11)
# ax.set_xlabel("pIC50 bins", fontsize=11)
# ax.set_title(f"pIC50 Distribution — Sexual test set (N = {N:,})", fontsize=16, pad=12)
# ax.grid(axis="y", alpha=0.25, linestyle=":")

# # annotate counts above bars
# for xi, c in zip(x, counts):
#     ax.text(xi, c + max(counts)*0.02, f"{c:,}", ha="center", va="bottom", fontsize=10)

# for spine in ("top", "right"):
#     ax.spines[spine].set_visible(False)

# fig.tight_layout()
# out_png = OUT_DIR / "sexual_pIC50_distribution.png"
# fig.savefig(out_png, dpi=300)
# plt.close(fig)

# print("Saved figure:", out_png)
# print("Saved counts:", counts_csv)