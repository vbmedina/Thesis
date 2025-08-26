#!/usr/bin/env python3
"""
02_choose_split_and_make_figs.py

After 01_score_models.py, choose the **showcase split** by a criterion and make the 3 figures:

Figures:
  1) scatter_<split>_bestfold.png  (two panels: asexual, sexual)
  2) pr_curves_<split>_bestfold.png (two panels)
  3) generalization_gap_bar.png     (ΔRMSE = Sexual − Asexual across splits)

Selection criterion (for the split and the fold within it):
  sexual_pr_auc   (default; maximize)
  sexual_mcc      (maximize)
  sexual_rmse     (minimize)
  delta_rmse      (minimize)

Usage:
  python 02_choose_split_and_make_figs.py \
    --base-root ~/Thesis/p2_models \
    --criterion sexual_pr_auc
"""

from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score
from scipy.stats import spearmanr

@dataclass
class Config:
    base_root: Path
    out_dir_name: str = "analysis_outputs"
    criterion: str = "sexual_pr_auc"   # sexual_pr_auc | sexual_mcc | sexual_rmse | delta_rmse
    threshold: float = 6.0

def parse_args() -> Config:
    p = argparse.ArgumentParser()
    p.add_argument("--base-root", type=str, required=True)
    p.add_argument("--out-dir-name", type=str, default="analysis_outputs")
    p.add_argument("--criterion", type=str, default="sexual_pr_auc",
                   choices=["sexual_pr_auc","sexual_mcc","sexual_rmse","delta_rmse"])
    p.add_argument("--threshold", type=float, default=6.0)
    a = p.parse_args()
    return Config(base_root=Path(a.base_root).expanduser(),
                  out_dir_name=a.out_dir_name,
                  criterion=a.criterion,
                  threshold=a.threshold)

def rmse(y_true, y_pred):
    from sklearn.metrics import mean_squared_error
    return float(mean_squared_error(y_true, y_pred, squared=False))

def spearman(y_true, y_pred):
    r, _ = spearmanr(y_true, y_pred)
    return float(r)

def main():
    cfg = parse_args()
    out_root = cfg.base_root / cfg.out_dir_name
    per_fold_path = out_root / "all_metrics_by_split_fold_core.csv"
    main_table    = out_root / "main_table_core_metrics.csv"

    if not per_fold_path.exists() or not main_table.exists():
        raise SystemExit("Run 01_score_models.py first to generate metrics CSVs.")

    per_fold_df = pd.read_csv(per_fold_path)
    main_df     = pd.read_csv(main_table)

    # Pick best split by criterion
    if cfg.criterion == "sexual_pr_auc":
        key = "sexual_data_pr_auc_mean"; maximize = True
    elif cfg.criterion == "sexual_mcc":
        key = "sexual_data_mcc_at_6_mean"; maximize = True
    elif cfg.criterion == "sexual_rmse":
        key = "sexual_data_rmse_mean"; maximize = False
    else:  # delta_rmse
        key = "delta_rmse_mean"; maximize = False

    best_idx = main_df[key].idxmax() if maximize else main_df[key].idxmin()
    best_split = str(main_df.loc[best_idx, "split"])
    print(f"Chosen showcase split: {best_split} (criterion={cfg.criterion}, value={main_df.loc[best_idx, key]:.4f})")

    # Best fold within split by same criterion on per-fold sexual_data metrics
    sub = per_fold_df[(per_fold_df["split"]==best_split) & (per_fold_df["dataset"]=="sexual_data")]
    if cfg.criterion == "sexual_pr_auc":
        best_fold = int(sub.sort_values("pr_auc", ascending=False).iloc[0]["fold"])
    elif cfg.criterion == "sexual_mcc":
        best_fold = int(sub.sort_values("mcc_at_6", ascending=False).iloc[0]["fold"])
    elif cfg.criterion == "sexual_rmse":
        best_fold = int(sub.sort_values("rmse", ascending=True).iloc[0]["fold"])
    else:
        # per-fold delta RMSE
        ft = per_fold_df[(per_fold_df["split"]==best_split) & (per_fold_df["dataset"]=="fold_test")][["fold","rmse"]].rename(columns={"rmse":"rmse_ft"})
        se = sub[["fold","rmse"]].rename(columns={"rmse":"rmse_se"})
        j = ft.merge(se, on="fold"); j["delta"] = j["rmse_se"] - j["rmse_ft"]
        best_fold = int(j.sort_values("delta", ascending=True).iloc[0]["fold"])

    print(f"Chosen fold: {best_fold}")

    # Figure 3: ΔRMSE bar
    plt.figure()
    x = np.arange(len(main_df))
    plt.bar(x, main_df["delta_rmse_mean"].values)
    plt.xticks(x, main_df["split"].values, rotation=15, ha="right")
    plt.ylabel("ΔRMSE (Sexual − Asexual)")
    plt.title("Generalization Gap by Split")
    plt.tight_layout()
    plt.savefig(out_root / "generalization_gap_bar.png", dpi=300)
    plt.close()

    # Load predictions for the chosen split/fold
    fold_dir = cfg.base_root / "models" / best_split / f"fold_{best_fold}"
    df_t = pd.read_csv(fold_dir / "pred_vs_true_fold_test.csv")
    df_s = pd.read_csv(fold_dir / "pred_vs_true_sexual_data.csv")

    # Figure 1: Scatter two-panel
    fig, axes = plt.subplots(1, 2, figsize=(9,4), dpi=150)
    for ax, df, ttl in zip(axes, [df_t, df_s], ["Asexual (fold test)", "Sexual (external)"]):
        yt, yp = df["y_true"].to_numpy(), df["y_pred"].to_numpy()
        lo, hi = float(min(yt.min(), yp.min())), float(max(yt.max(), yp.max()))
        ax.scatter(yt, yp, s=12, alpha=0.6, linewidths=0)
        ax.plot([lo,hi],[lo,hi], ls="--")
        if np.std(yt)>0:
            b = np.cov(yt, yp, ddof=0)[0,1] / np.var(yt, ddof=0)
            a = yp.mean() - b * yt.mean()
        else:
            a, b = 0.0, 1.0
        r = spearman(yt, yp)
        e = rmse(yt, yp)
        ax.set_title(ttl); ax.set_xlabel("True pIC50"); ax.set_ylabel("Predicted pIC50")
        ax.text(0.02, 0.98, f"RMSE={e:.3f}\nSpearman={r:.3f}\nfit: y={a:.2f}+{b:.2f}x",
                transform=ax.transAxes, va="top", ha="left")
    fig.suptitle(f"True vs Predicted — {best_split} (best fold={best_fold})", y=1.02)
    fig.tight_layout(); fig.savefig(out_root / f"scatter_{best_split}_bestfold.png", bbox_inches="tight"); plt.close(fig)

    # Figure 2: PR curves two-panel
    fig, axes = plt.subplots(1, 2, figsize=(9,4), dpi=150)
    for ax, df, ttl in zip(axes, [df_t, df_s], ["Asexual (fold test)", "Sexual (external)"]):
        yt, yp = df["y_true"].to_numpy(), df["y_pred"].to_numpy()
        yb = (yt >= cfg.threshold).astype(int)
        prec, rec, _ = precision_recall_curve(yb, yp)
        ap = average_precision_score(yb, yp)
        ax.plot(rec, prec)
        ax.set_xlabel("Recall"); ax.set_ylabel("Precision"); ax.set_title(f"{ttl} (AP={ap:.3f})")
    fig.suptitle(f"PR Curves — {best_split} (best fold={best_fold})", y=1.02)
    fig.tight_layout(); fig.savefig(out_root / f"pr_curves_{best_split}_bestfold.png", bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    main()
