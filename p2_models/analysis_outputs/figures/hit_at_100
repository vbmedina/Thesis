#!/usr/bin/env python3
# Save FIVE separate PNGs:
# 1) HitRate@100, 2) EF@100, 3) RMSE, 4) Spearman ρ, 5) ΔHitRate@100
# Legends are placed on the RIGHT of each figure.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- Paths --------------------
TOPK_CSV = "p2_models/analysis_outputs/topk_metrics_by_split_fold_rows.csv"
CORE_CSV = "p2_models/analysis_outputs/all_metrics_by_split_fold_core.csv"
OUT_DIR  = Path("p2_models/analysis_outputs/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------- Visuals --------------------
COLOR_ASEX = "#b82720"  # asexual (fold_test)
COLOR_SEX  = "#dc6d6d"  # sexual (external)
COLOR_GAP  = "#8d1115"  # delta bars
GRID_ALPHA = 0.2
plt.rcParams.update({
    "figure.dpi": 150,
    "axes.grid": True,
    "grid.alpha": GRID_ALPHA,
    "grid.linestyle": "--",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# Fixed left-to-right split order
SPLIT_ORDER  = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
SPLIT_LABELS = {
    "random": "random",
    "scaffold": "scaffold",
    "butina": "butina",
    "umap_kmeans": "umap kmeans",
    "umap_ward": "umap ward",
}

# -------------------- Helpers --------------------
def require_columns(df, cols, name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing columns: {missing}. Present: {list(df.columns)}")

def order_splits(available):
    return [s for s in SPLIT_ORDER if s in available]

def grouped_bar_figure(metric_df, value_col, sd_col, ylabel, title, splits, out_path, ylim=None, zero_line=False):
    """
    Draw grouped bars (Asexual vs Sexual) with error bars; legend on RIGHT; save to out_path.
    metric_df must have: split, dataset, value_col, sd_col
    """
    df_fold = metric_df[metric_df["dataset"]=="fold_test"].set_index("split").reindex(splits)
    df_sex  = metric_df[metric_df["dataset"]=="sexual_data"].set_index("split").reindex(splits)

    vals_a = np.nan_to_num(df_fold[value_col].to_numpy(), nan=0.0)
    sds_a  = np.nan_to_num(df_fold[sd_col].to_numpy(),   nan=0.0)
    vals_b = np.nan_to_num(df_sex[value_col].to_numpy(), nan=0.0)
    sds_b  = np.nan_to_num(df_sex[sd_col].to_numpy(),    nan=0.0)

    x = np.arange(len(splits))
    bar_w = 0.38
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.bar(x - bar_w/2, vals_a, width=bar_w, yerr=sds_a, capsize=3,
           color=COLOR_ASEX, label="Asexual (fold_test)")
    ax.bar(x + bar_w/2, vals_b, width=bar_w, yerr=sds_b, capsize=3,
           color=COLOR_SEX,  label="Sexual (external)")

    ax.set_xticks(x)
    ax.set_xticklabels([SPLIT_LABELS.get(s, s) for s in splits], rotation=15, ha="right")
    ax.set_xlabel("Split", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=16, pad=12)
    if ylim is not None:
        ax.set_ylim(ylim)

    # <<< add a black horizontal line at y=0 when requested >>>
    if zero_line:
        ax.axhline(0.0, color="black", lw=1.2, zorder=3)

    # Legend on the RIGHT (outside the axes)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, ncol=1)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    print("Saved", out_path)

def delta_bar_figure(delta_vals, splits, ylabel, title, out_path):
    x = np.arange(len(splits))
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.bar(x, delta_vals, color="#b11218")  # no label, no legend
    ax.axhline(0.0, color="#666", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([SPLIT_LABELS.get(s, s) for s in splits], rotation=15, ha="right")
    ax.set_xlabel("Split", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=16, pad=12)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    print("Saved", out_path)

    # Legend on the RIGHT (outside the axes)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=1)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    print("Saved", out_path)

# -------------------- Load --------------------
if not Path(TOPK_CSV).exists():
    raise FileNotFoundError(f"Cannot find {TOPK_CSV}")
if not Path(CORE_CSV).exists():
    raise FileNotFoundError(f"Cannot find {CORE_CSV}")

topk = pd.read_csv(TOPK_CSV)
core = pd.read_csv(CORE_CSV)

require_columns(topk, ["split","dataset","hit_at_k","ef_at_k"], "topk CSV")
require_columns(core, ["split","dataset","rmse","spearman_rho"], "core CSV")

# -------------------- Aggregate --------------------
agg_topk = (
    topk.groupby(["split","dataset"], as_index=False)
        .agg(hit_mean=("hit_at_k","mean"), hit_sd=("hit_at_k","std"),
             ef_mean=("ef_at_k","mean"),   ef_sd=("ef_at_k","std"))
)
agg_core = (
    core.groupby(["split","dataset"], as_index=False)
        .agg(rmse_mean=("rmse","mean"), rmse_sd=("rmse","std"),
             rho_mean=("spearman_rho","mean"), rho_sd=("spearman_rho","std"))
)

available = set(agg_topk["split"]).union(set(agg_core["split"]))
splits = order_splits(available)
if not splits:
    raise SystemExit("No recognized splits found in the CSVs.")

# Frames for plotting
df_hit = agg_topk.rename(columns={"hit_mean":"val","hit_sd":"sd"})
df_ef  = agg_topk.rename(columns={"ef_mean":"val","ef_sd":"sd"})
df_rmse= agg_core.rename(columns={"rmse_mean":"val","rmse_sd":"sd"})
df_rho = agg_core.rename(columns={"rho_mean":"val","rho_sd":"sd"})

# Delta (Sexual − Asexual) for Hit@100
hit_wide = agg_topk.pivot(index="split", columns="dataset", values="hit_mean").reindex(splits)
delta_hit = (hit_wide["sexual_data"] - hit_wide["fold_test"]).to_numpy()

# -------------------- Make separate PNGs --------------------
# 1) HitRate@100
grouped_bar_figure(
    metric_df=df_hit, value_col="val", sd_col="sd",
    ylabel="HitRate@100", title="Early Recognition — Hit Rate at 100 by Split",
    splits=splits,
    out_path=OUT_DIR / "panel_hit_at_100.png",
    ylim=(0,1)
)

# 2) EF@100
grouped_bar_figure(
    metric_df=df_ef, value_col="val", sd_col="sd",
    ylabel="EF@100", title="Early Recognition — Enrichment Factor at 100 by Split",
    splits=splits,
    out_path=OUT_DIR / "panel_ef_at_100.png",
    ylim=None
)

# 3) RMSE
grouped_bar_figure(
    metric_df=df_rmse, value_col="val", sd_col="sd",
    ylabel="RMSE (pIC50)", title="Regression Error — RMSE by Split",
    splits=splits,
    out_path=OUT_DIR / "panel_rmse.png",
    ylim=None
)

# 4) Spearman ρ  (full range -1..1) + black zero line
grouped_bar_figure(
    metric_df=df_rho, value_col="val", sd_col="sd",
    ylabel="Spearman ρ", title="Rank Quality — Spearman (rho) by Split",
    splits=splits,
    out_path=OUT_DIR / "panel_spearman.png",
    ylim=(-1, 1),
    zero_line=True
)

# 5) ΔHit@100 (Sexual − Asexual) — with right-side legend
delta_bar_figure(
    delta_vals=delta_hit, splits=splits,
    ylabel="Change in Hit Rate at 100 (Sexual - Asexual)",
    title="Generalization Gap — Difference in Hit Rate at 100",
    out_path=OUT_DIR / "panel_delta_hit_at_100.png"
)

print("All done.")
