# FIRST

#!/usr/bin/env python3
# Two-panel (Asexual vs Sexual) scatter for your chosen models only,
# with regression fit, right-side vertical legend, and metrics tables lower.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, mean_squared_error
from math import sqrt

# --------- CONFIG ---------
BASE_MODELS = Path("p2_models/models")
OUT_DIR     = Path("p2_models/analysis_outputs/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Only these split/fold pairs will be plotted:
SELECTED = [
    ("random", 1),
    ("scaffold", 1),
    ("butina", 2),
    ("umap_kmeans", 1),
    ("umap_ward", 5),
]

# Pretty names for the split label at the top
SPLIT_PRETTY = {
    "random": "Random",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (kmeans)",
    "umap_ward": "UMAP (ward)",
}

THRESHOLD = 6.0
TOP_K     = 100  # for HitRate@100

# Reds palette
PALETTE    = sns.color_palette("Reds", 6)
PT_COLOR   = PALETTE[4]            # points
IDEAL_CLR  = "#666"                # y=x & threshold guides
FIT_CLR    = "#000"                # regression fit
BG_TINT    = "#fde8ea"             # light red backdrop

# --------- I/O helpers ---------
def load_preds(split: str, fold: int, tag: str) -> pd.DataFrame:
    fn = "pred_vs_true_fold_test.csv" if tag == "fold_test" else "pred_vs_true_sexual_data.csv"
    p = BASE_MODELS / split / f"fold_{fold}" / fn
    if not p.exists():
        raise FileNotFoundError(f"Missing predictions: {p}")
    df = pd.read_csv(p)
    if not {"y_true","y_pred"}.issubset(df.columns):
        raise ValueError(f"{p} must contain y_true,y_pred columns")
    return df

# --------- metrics ---------
def confusion_counts(y_true, y_pred, thr):
    yb = (y_true >= thr).astype(int); yp = (y_pred >= thr).astype(int)
    tp = int(((yb==1)&(yp==1)).sum())
    fp = int(((yb==0)&(yp==1)).sum())
    tn = int(((yb==0)&(yp==0)).sum())
    fn = int(((yb==1)&(yp==0)).sum())
    return tp, fp, tn, fn

def mcc(tp, fp, tn, fn):
    den = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    return 0.0 if den == 0 else (tp*tn - fp*fn) / sqrt(den)

def hit_at_k(y_true, y_pred, thr, k):
    K = min(int(k), len(y_true))
    order = np.argsort(-y_pred)[:K]
    yb = (y_true >= thr).astype(int)
    return float(yb[order].mean())

def safe_roc_auc(y_true, y_pred, thr):
    yb = (y_true >= thr).astype(int)
    if len(np.unique(yb)) < 2:
        return np.nan
    return roc_auc_score(yb, y_pred)

def rmse(y_true, y_pred):
    return sqrt(mean_squared_error(y_true, y_pred))

# --------- fit helpers ---------
def linear_fit(y_true, y_pred):
    """Return slope b, intercept a, and R^2 for y_pred ≈ a + b * y_true."""
    b, a = np.polyfit(y_true, y_pred, 1)
    y_hat = a + b * y_true
    ss_res = float(np.sum((y_pred - y_hat)**2))
    ss_tot = float(np.sum((y_pred - np.mean(y_pred))**2))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return b, a, r2

# --------- annotated scatter panel ---------
def annotate_quadrants(ax, tp, fp, tn, fn):
    x0,x1 = ax.get_xlim(); y0,y1 = ax.get_ylim()
    offx = (x1-x0)*0.02; offy = (y1-y0)*0.02
    ax.text(x0+offx, y0+offy, f"TN = {tn:,}", fontsize=9,
            ha="left", va="bottom", bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.9))
    ax.text(x0+offx, y1-offy, f"FP = {fp:,}", fontsize=9,
            ha="left", va="top", bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.9))
    ax.text(x1-offx, y1-offy, f"TP = {tp:,}", fontsize=9,
            ha="right", va="top", bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.9))
    ax.text(x1-offx, y0+offy, f"FN = {fn:,}", fontsize=9,
            ha="right", va="bottom", bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.9))

def metrics_table_below(ax, hit, mcc_v, roc, rmse_v):
    """Place the table further *below* the axes so it doesn't crowd the plot."""
    cols = ["hit rate", "MCC", "ROC AUC", "RMSE"]
    vals = [f"{hit*100:.1f}%", f"{mcc_v:.3f}", f"{roc:.3f}" if not np.isnan(roc) else "—", f"{rmse_v:.3f}"]
    # bbox: [left, bottom, width, height] in axes coords; bottom < 0 puts it below.
    # moved lower (bottom=-0.62) and a bit taller (height=0.34)
    tbl = ax.table(cellText=[vals], colLabels=cols, cellLoc="center",
                   bbox=[0.0, -0.60, 1.0, 0.25])                            #HERE
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)

def draw_panel(ax, df, title):
    y_true = df["y_true"].to_numpy()
    y_pred = df["y_pred"].to_numpy()

    # square limits with a small margin
    lo = np.floor(min(y_true.min(), y_pred.min())*2)/2 - 0.2
    hi = np.ceil(max(y_true.max(), y_pred.max())*2)/2 + 0.2
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_facecolor(BG_TINT)

    # scatter
    ax.scatter(y_true, y_pred, s=6, alpha=0.35, color=PT_COLOR, edgecolor="none", label="Molecules")

    # ideal line y=x
    ax.plot([lo, hi], [lo, hi], ls="--", color=IDEAL_CLR, lw=1, label="Ideal: y = x")

    # regression fit y = a + b x
    b, a, r2 = linear_fit(y_true, y_pred)
    xline = np.array([lo, hi], dtype=float)
    ax.plot(xline, a + b*xline, color=FIT_CLR, lw=1.6,
            label=f"Fit: y = {b:.2f}x {a:+.2f}  (R²={r2:.3f})")

    # guides at threshold
    ax.axvline(THRESHOLD, color=IDEAL_CLR, lw=1)
    ax.axhline(THRESHOLD, color=IDEAL_CLR, lw=1)

    ax.set_xlabel("Measured pIC50"); ax.set_ylabel("Predicted pIC50")
    ax.set_title(title)

    # legend: vertical on the RIGHT, inside the axes
    ax.legend(loc="center right", bbox_to_anchor=(.74, -.22), ncol=1,
              frameon=True, facecolor="white", framealpha=0.9, fontsize=9)

    # metrics
    tp, fp, tn, fn = confusion_counts(y_true, y_pred, THRESHOLD)
    hit = hit_at_k(y_true, y_pred, THRESHOLD, TOP_K)
    mcc_v = mcc(tp, fp, tn, fn)
    roc = safe_roc_auc(y_true, y_pred, THRESHOLD)
    rmse_v = rmse(y_true, y_pred)

    annotate_quadrants(ax, tp, fp, tn, fn)
    metrics_table_below(ax, hit, mcc_v, roc, rmse_v)

# --------- run for the selected models ---------
for split, fold in SELECTED:
    try:
        df_asex = load_preds(split, fold, "fold_test")
        df_sex  = load_preds(split, fold, "sexual_data")
    except Exception as e:
        print(f"[skip] {split} fold {fold}: {e}")
        continue

    # a bit wider & taller to give room for tables below
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 7.0))
    # big label on top with split name + fold
    split_name = SPLIT_PRETTY.get(split, split.replace("_", " "))
    fig.suptitle(f"Split: {split_name} — fold {fold}", fontsize=16, y=0.95)

    draw_panel(axes[0], df_asex, "Asexual (fold test)")
    draw_panel(axes[1], df_sex,  "Sexual (external)")

    # extra space at bottom for the lower tables
    plt.subplots_adjust(top=0.86, bottom=0.34, wspace=0.28)

    out = OUT_DIR / f"scatter_true_vs_pred_{split}_fold{fold}.png"
    plt.tight_layout(pad=1.5)
    plt.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out)
