'''
Compute early-recognition Top-K metrics (hit rate, enrichment factor, MCC) from per-fold prediction CSVs.

Requirements
- Per-fold prediction files produced by the evaluation step:
  ./p2_models/models/<split>/fold_<k>/
    * pred_vs_true_fold_test.csv
    * pred_vs_true_sexual_data.csv
- Required columns in each file: "y_true", "y_pred"

Splits & Folds
- Splits are discovered under ./p2_models/models/ unless restricted via --splits.
- Folds are discovered as fold_* unless restricted via --folds.
- Datasets evaluated per fold: ["fold_test", "sexual_data"] (skips missing files).

Metrics (computed on valid rows, keeping any replicates)
- Absolute Top-K (@K):
  * hit_at_k        = TP@K / K
  * ef_at_k         = hit_at_k / prevalence
  * mcc_at_k        = MCC when classifying only the Top-K as positive
- Fractional Top-K% (@K%):
  * hit_at_kpct, ef_at_kpct, mcc_at_kpct (same definitions using K% of rows)
- Additional fields:
  * prevalence, k_used, k_pct_used, n_rows
  * n_actives_total, n_inactives_total
  * tp_at_k, expected_tp_at_k (= prevalence * K)
  * tp_at_kpct, expected_tp_at_kpct
- Positives are defined by a pIC50 threshold (default 6.0): y_true ≥ threshold.

Procedure
1) Traverse ./p2_models/models/<split>/fold_<k>/ (optionally filtered by --splits/--folds).
2) For each dataset tag in {"fold_test","sexual_data"}:
   a. Load pred_vs_true_{tag}.csv.
   b. Sort rows by y_pred descending.
   c. Compute Top-K and Top-K% metrics.
3) Concatenate per-fold rows and write a single CSV.

Command-line Arguments
- --base-root     (default: ./p2_models)
- --out-dir-name  (default: analysis_outputs)
- --threshold     (float, default: 6.0)
- --k             (int absolute K, default: 100)
- --k-frac        (fractional K%, default: 0.05 for 5%; uses ceil, min 1 row)
- --splits        (comma-separated list to restrict splits, e.g., "random,scaffold")
- --folds         (comma-separated ints to restrict folds, e.g., "1,2,3")

Outputs
- ./p2_models/analysis_outputs/topk_metrics_by_split_fold_rows.csv
  Columns include: split, fold, dataset, hit_at_k/ef_at_k/mcc_at_k,
  hit_at_kpct/ef_at_kpct/mcc_at_kpct, prevalence, k_used, k_pct_used, n_rows, and counts.

Notes
- K is clamped to [1, N]; K% uses ceil(max(1, k_frac * N)).
- Files missing for a given (split, fold, tag) are skipped gracefully.
- Exits with a message if no prediction CSVs are found.
'''
from __future__ import annotations
import argparse, math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict
import numpy as np
import pandas as pd

@dataclass
class Config:
    base_root: Path
    out_dir_name: str
    threshold: float
    k_abs: int
    k_frac: float
    splits: List[str] | None
    folds: List[int] | None

def parse_args() -> Config:
    p = argparse.ArgumentParser()
    p.add_argument("--base-root", type=str, default="./p2_models")
    p.add_argument("--out-dir-name", type=str, default="analysis_outputs")
    p.add_argument("--threshold", type=float, default=6.0)
    p.add_argument("--k", type=int, default=100, help="Absolute K for top-K rows")
    p.add_argument("--k-frac", type=float, default=0.05, help="Fractional K (e.g., 0.05 for 5%)")
    p.add_argument("--splits", type=str, default=None, help="Comma list to restrict splits")
    p.add_argument("--folds", type=str, default=None, help="Comma list to restrict folds (ints)")
    a = p.parse_args()
    splits = [s.strip() for s in a.splits.split(",")] if a.splits else None
    folds = [int(x) for x in a.folds.split(",")] if a.folds else None
    return Config(
        base_root=Path(a.base_root).expanduser(),
        out_dir_name=a.out_dir_name,
        threshold=float(a.threshold),
        k_abs=int(a.k),
        k_frac=float(a.k_frac),
        splits=splits,
        folds=folds
    )

def load_pred_csv(base_root: Path, split: str, fold: int, tag: str) -> pd.DataFrame:
    """tag in {'fold_test','sexual_data'}"""
    fn = f"pred_vs_true_{'fold_test' if tag=='fold_test' else 'sexual_data'}.csv"
    p = base_root / "models" / split / f"fold_{fold}" / fn
    if not p.exists():
        raise FileNotFoundError(str(p))
    df = pd.read_csv(p)
    if not {"y_true","y_pred"}.issubset(df.columns):
        raise ValueError(f"{p} must contain columns: y_true, y_pred")
    return df

def mcc_from_confusion(tp, fp, tn, fn) -> float:
    den = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    if den <= 0: return 0.0
    return float((tp*tn - fp*fn) / (den ** 0.5))

def early_metrics_rows(df: pd.DataFrame, thr: float, k_abs: int, k_frac: float) -> Dict[str, float]:
    """Row-level early-recognition metrics (keep replicates)."""
    N = len(df)
    if N == 0:
        return {k: float("nan") for k in [
            "hit_at_k","ef_at_k","mcc_at_k",
            "hit_at_kpct","ef_at_kpct","mcc_at_kpct",
            "prevalence","k_used","k_pct_used","n_rows",
            "n_actives_total","n_inactives_total",
            "tp_at_k","expected_tp_at_k","tp_at_kpct","expected_tp_at_kpct"
        ]}

    # sort rows by predicted pIC50 (desc)
    order = np.argsort(-df["y_pred"].to_numpy())
    y_true = df["y_true"].to_numpy()
    yb = (y_true >= thr).astype(int)

    pos = int(yb.sum())
    prevalence = pos / N if N>0 else float("nan")

    # Absolute K
    K = max(1, min(int(k_abs), N))
    top_idx = order[:K]
    tp_k = int(yb[top_idx].sum())
    hit_k = tp_k / K
    ef_k = (hit_k / prevalence) if prevalence > 0 else float("nan")

    yp = np.zeros(N, dtype=int); yp[top_idx] = 1
    tp = int(((yb==1) & (yp==1)).sum())
    tn = int(((yb==0) & (yp==0)).sum())
    fp = int(((yb==0) & (yp==1)).sum())
    fn = int(((yb==1) & (yp==0)).sum())
    mcc_k = mcc_from_confusion(tp, fp, tn, fn)

    # Fractional K%
    Kp = max(1, int(math.ceil(k_frac * N)))
    top_idx_p = order[:Kp]
    tp_kp = int(yb[top_idx_p].sum())
    hit_kp = tp_kp / Kp
    ef_kp = (hit_kp / prevalence) if prevalence > 0 else float("nan")

    yp_p = np.zeros(N, dtype=int); yp_p[top_idx_p] = 1
    tp2 = int(((yb==1) & (yp_p==1)).sum())
    tn2 = int(((yb==0) & (yp_p==0)).sum())
    fp2 = int(((yb==0) & (yp_p==1)).sum())
    fn2 = int(((yb==1) & (yp_p==0)).sum())
    mcc_kp = mcc_from_confusion(tp2, fp2, tn2, fn2)

    return {
        "hit_at_k": hit_k, "ef_at_k": ef_k, "mcc_at_k": mcc_k,
        "hit_at_kpct": hit_kp, "ef_at_kpct": ef_kp, "mcc_at_kpct": mcc_kp,
        "prevalence": prevalence, "k_used": K, "k_pct_used": Kp, "n_rows": N,
        "n_actives_total": pos, "n_inactives_total": N - pos,
        "tp_at_k": tp_k, "expected_tp_at_k": prevalence * K,
        "tp_at_kpct": tp_kp, "expected_tp_at_kpct": prevalence * Kp
    }

def main():
    cfg = parse_args()
    out_root = cfg.base_root / cfg.out_dir_name
    out_root.mkdir(parents=True, exist_ok=True)

    rows = []
    models_root = cfg.base_root / "models"
    if not models_root.exists():
        raise SystemExit(f"Models directory not found: {models_root}")

    for split_dir in sorted(models_root.iterdir()):
        if not split_dir.is_dir(): continue
        split = split_dir.name
        if cfg.splits and split not in cfg.splits: continue

        for fold_dir in sorted(split_dir.glob("fold_*")):
            try:
                fold = int(fold_dir.name.split("_")[-1])
            except Exception:
                continue
            if cfg.folds and fold not in cfg.folds: continue

            for tag in ["fold_test","sexual_data"]:
                try:
                    df = load_pred_csv(cfg.base_root, split, fold, tag)
                except FileNotFoundError:
                    continue
                mets = early_metrics_rows(df, cfg.threshold, cfg.k_abs, cfg.k_frac)
                rows.append({"split": split, "fold": fold, "dataset": tag, **mets})

    if not rows:
        raise SystemExit("No prediction CSVs found. Expected pred_vs_true_{fold_test,sexual_data}.csv under models/<split>/fold_*/")

    per_df = pd.DataFrame(rows).sort_values(["split","fold","dataset"])
    per_path = out_root / "topk_metrics_by_split_fold_rows.csv"
    per_df.to_csv(per_path, index=False)
    print("Wrote:", per_path)

if __name__ == "__main__":
    main()