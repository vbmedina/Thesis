"""
03_appendix_metrics.py

Appendix-only metrics computed from saved per-fold predictions:
- MAE, R²
- Calibration intercept & slope (Pred = a + b·True)
- ROC-AUC (scores), Balanced Accuracy @ 6.0
- Confusion matrix counts (TP, FP, TN, FN)
- EF@5% (sexual set only)

Requires that your per-fold prediction CSVs already exist:
  models/<split>/fold_<i>/pred_vs_true_fold_test.csv
  models/<split>/fold_<i>/pred_vs_true_sexual_data.csv

Outputs (in <base-root>/<out-dir-name>/):
  - appendix_metrics_by_split_fold.csv
  - appendix_summary_by_split.csv

Usage:
  python 03_appendix_metrics.py \
    --base-root ~/Thesis/p2_models
"""

from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import (
    mean_absolute_error, r2_score, balanced_accuracy_score,
    roc_auc_score, confusion_matrix
)

@dataclass
class Config:
    base_root: Path
    splits: list[str]
    folds: list[int]
    threshold: float = 6.0
    out_dir_name: str = "analysis_outputs"

def parse_args() -> Config:
    p = argparse.ArgumentParser()
    p.add_argument("--base-root", type=str, required=True)
    p.add_argument("--splits", type=str, default="random,scaffold,butina,umap_kmeans,umap_ward")
    p.add_argument("--folds", type=str, default="1,2,3,4,5")
    p.add_argument("--threshold", type=float, default=6.0)
    p.add_argument("--out-dir-name", type=str, default="analysis_outputs")
    a = p.parse_args()
    return Config(
        base_root=Path(a.base_root).expanduser(),
        splits=[s.strip() for s in a.splits.split(",") if s.strip()],
        folds=[int(x) for x in a.folds.split(",")],
        threshold=a.threshold,
        out_dir_name=a.out_dir_name,
    )

def load_preds(csv_path: Path):
    df = pd.read_csv(csv_path)
    return df["y_true"].to_numpy(), df["y_pred"].to_numpy()

def slope_intercept(y_true: np.ndarray, y_pred: np.ndarray):
    if np.var(y_true) <= 0:
        return float(y_pred.mean()), 0.0
    b = float(np.cov(y_true, y_pred, ddof=0)[0,1] / np.var(y_true, ddof=0))
    a = float(y_pred.mean() - b * y_true.mean())
    return a, b

def ef_at_k(y_true: np.ndarray, y_pred: np.ndarray, thr: float, frac: float) -> float:
    """Enrichment factor at top frac (e.g., 0.05). Active if y_true >= thr."""
    n = len(y_true)
    if n == 0:
        return float("nan")
    k = max(1, int(np.ceil(frac * n)))
    yb = (y_true >= thr).astype(int)
    order = np.argsort(-y_pred)  # descending by predicted score
    hits = int(yb[order[:k]].sum())
    prevalence = float(yb.mean())
    expected = prevalence * k
    return float(hits / expected) if expected > 0 else float("nan")

def main():
    cfg = parse_args()
    out_root = cfg.base_root / cfg.out_dir_name
    out_root.mkdir(parents=True, exist_ok=True)

    per_fold_rows = []

    for split in cfg.splits:
        for fold in cfg.folds:
            fold_dir = cfg.base_root / "models" / split / f"fold_{fold}"
            pt = fold_dir / "pred_vs_true_fold_test.csv"
            ps = fold_dir / "pred_vs_true_sexual_data.csv"
            if not (pt.exists() and ps.exists()):
                # silently skip missing folds
                continue

            # Compute appendix metrics for both datasets
            for tag, p in [("fold_test", pt), ("sexual_data", ps)]:
                yt, yp = load_preds(p)

                # Regression appendix metrics
                mae = float(mean_absolute_error(yt, yp))
                r2  = float(r2_score(yt, yp))
                a, b = slope_intercept(yt, yp)

                # Classification appendix metrics
                yb = (yt >= cfg.threshold).astype(int)
                roc = float(roc_auc_score(yb, yp)) if len(np.unique(yb)) > 1 else float("nan")
                yp_bin = (yp >= cfg.threshold).astype(int)
                bacc = float(balanced_accuracy_score(yb, yp_bin)) if len(np.unique(yb)) > 1 else float("nan")

                # Confusion matrix counts (labels [1,0] → rows: actual 1 then 0)
                cm = confusion_matrix(yb, yp_bin, labels=[1, 0])
                tp, fn, fp, tn = int(cm[0, 0]), int(cm[0, 1]), int(cm[1, 0]), int(cm[1, 1])

                # EF@5% only for sexual_data
                ef5 = ef_at_k(yt, yp, cfg.threshold, 0.05) if tag == "sexual_data" else float("nan")

                per_fold_rows.append({
                    "split": split, "fold": fold, "dataset": tag,
                    "mae": mae, "r2": r2,
                    "calib_intercept": a, "calib_slope": b,
                    "roc_auc": roc, "bacc_at_6": bacc,
                    "tp": tp, "fn": fn, "fp": fp, "tn": tn,
                    "ef5": ef5
                })

    per_fold_df = pd.DataFrame(per_fold_rows).sort_values(["split", "fold", "dataset"])
    per_fold_path = out_root / "appendix_metrics_by_split_fold.csv"
    per_fold_df.to_csv(per_fold_path, index=False)

    # Summary: mean ± sd across folds per split & dataset
    metrics = ["mae","r2","calib_intercept","calib_slope","roc_auc","bacc_at_6","ef5"]
    summary_rows = []
    for split in sorted(per_fold_df["split"].unique()):
        for dataset in ["fold_test", "sexual_data"]:
            sub = per_fold_df[(per_fold_df["split"] == split) & (per_fold_df["dataset"] == dataset)]
            if len(sub) == 0:
                continue
            row = {"split": split, "dataset": dataset}
            for m in metrics:
                vals = sub[m].dropna()
                row[f"{m}_mean"] = float(vals.mean()) if len(vals) > 0 else float("nan")
                row[f"{m}_sd"]   = float(vals.std(ddof=1)) if len(vals) > 1 else float("nan")
            summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows).sort_values(["split","dataset"])
    summary_path = out_root / "appendix_summary_by_split.csv"
    summary_df.to_csv(summary_path, index=False)

    print("Wrote:", per_fold_path)
    print("Wrote:", summary_path)

if __name__ == "__main__":
    main()
