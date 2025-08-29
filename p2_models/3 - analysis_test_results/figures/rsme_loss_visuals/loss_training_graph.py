'''
This script plots per-epoch RMSE curves (train & val) for each split/fold from training logs, and save one PNG per run.

Requirements
- One CSV per run at: <base>/\<split>/fold_<k>/metrics_per_epoch.csv
  * Must contain: "epoch"
  * Optional (used if present): "train/rmse", "val/rmse"

Splits & Folds
- Fixed split order: ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
- Folds are provided via --folds (default "1,2,3,4,5")
- Plot titles use pretty names: Random, Scaffold, Butina, UMAP (kmeans), UMAP (ward)

What it does
1) Loads each metrics_per_epoch.csv, keeps the last record per epoch (handles repeated epochs),
   and casts "epoch" to int.
2) Draws RMSE vs epoch lines using seaborn/matplotlib:
   - Train curve: semi-transparent line (if "train/rmse" exists)
   - Validation curve: solid line (if "val/rmse" exists)
   - Best validation epoch (min "val/rmse"): vertical dashed line + point marker
3) Uses a consistent Reds palette—one distinct shade per split; no legend.
4) Saves figures to the output folder with informative filenames.

Command-line Arguments
- --base   Path to the models root containing split/fold dirs
           (default: /Users/victoriamedina/Thesis_Project/thesis/p2_models/models)
- --out    Output directory for PNGs
           (default: /Users/victoriamedina/Thesis_Project/thesis/p2_models/analysis_outputs/figures)
- --folds  Comma-separated list of fold indices to plot (default: 1..5)

Outputs
- One PNG per (split, fold) that has a valid CSV:
  <out>/loss_rmse_<split>_fold_<k>.png

Notes
- Runs are skipped if the CSV is missing or lacks both "train/rmse" and "val/rmse".
- Axes use a white-grid theme; titles are "RMSE per Epoch: <Split> Fold <k>".
- The script prints which files were saved and a final count of generated figures.
'''

#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, gzip
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

SP_ORDER = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
PRETTY = {
    "random": "Random",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (kmeans)",
    "umap_ward": "UMAP (ward)",
}

FILE_CANDIDATES = ["metrics_per_epoch.csv"]

def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")

def _sep(path: Path) -> str:
    suf = path.suffix.lower()
    inner = Path(path.stem).suffix.lower() if suf == ".gz" else suf
    return "\t" if inner == ".tsv" else ","

def find_metrics_file(run_dir: Path) -> Path | None:
    # try exact first, then globs
    checks = []
    for name in FILE_CANDIDATES:
        if "*" in name:
            checks += sorted(run_dir.glob(name))
        else:
            p = run_dir / name
            if p.exists(): return p
    return checks[0] if checks else None

def read_per_epoch(path: Path) -> pd.DataFrame:
    with _open(path) as f:
        df = pd.read_csv(f, sep=_sep(path))
    if "epoch" not in df.columns:
        raise ValueError(f'No "epoch" col in {path}')
    df = df[df["epoch"].notna()].copy()
    df["epoch"] = df["epoch"].astype(int)
    # keep last row for each epoch
    df = df.groupby("epoch", as_index=False).last().sort_values("epoch")
    return df

def pick_val_rmse(df: pd.DataFrame) -> pd.Series | None:
    # preferred columns in order
    if "val/rmse" in df.columns and df["val/rmse"].notna().any():
        return df["val/rmse"]
    if "val/mse" in df.columns and df["val/mse"].notna().any():
        return np.sqrt(df["val/mse"])
    # cautious fallback: if only val_loss exists, try sqrt(val_loss)
    if "val_loss" in df.columns and df["val_loss"].notna().any():
        return np.sqrt(df["val_loss"])
    # final fallback: training RMSE if validation is missing
    if "train/rmse" in df.columns and df["train/rmse"].notna().any():
        return df["train/rmse"]
    if "train/mse" in df.columns and df["train/mse"].notna().any():
        return np.sqrt(df["train/mse"])
    return None

def rmse_to_learning_rate(rmse: pd.Series, smooth_window: int | None = 3) -> pd.Series:
    # ΔRMSE per epoch: previous − current (positive = improvement)
    rate = rmse.shift(1) - rmse
    # optional simple smoothing to reduce noise
    if smooth_window and smooth_window > 1:
        rate = rate.rolling(window=smooth_window, min_periods=1, center=False).mean()
    return rate

def plot_combo(base_dir: Path, out_dir: Path, runs: list[tuple[str, int]],
               smooth_window: int = 3) -> None:
    sns.set_theme(style="whitegrid")
    reds = sns.color_palette("Reds", n_colors=len(SP_ORDER))
    color_map = {sp: reds[i] for i, sp in enumerate(SP_ORDER)}

    # 1) Derived learning-rate (from RMSE deltas)
    plt.figure(figsize=(8.5, 5))
    drew = 0
    for split, fold in runs:
        run_dir = base_dir / split / f"fold_{fold}"
        mpath = find_metrics_file(run_dir)
        if not mpath:
            print(f"[skip] metrics file missing in {run_dir}")
            continue
        df = read_per_epoch(mpath)
        y = pick_val_rmse(df)
        if y is None:
            print(f"[skip] no usable RMSE in {mpath}")
            continue
        rate = rmse_to_learning_rate(y, smooth_window=smooth_window)
        label = f"{PRETTY.get(split, split)} (fold {fold})"
        plt.plot(df["epoch"], rate, label=label, linewidth=2.2,
                 alpha=0.95, color=color_map.get(split))
        drew += 1
    if drew == 0:
        print("No curves plotted. Ensure your files have val/rmse or val/mse.")
    plt.axhline(0.0, linestyle="--", linewidth=1.0, alpha=0.5)
    plt.title("Rate of Learning from RMSE (ΔRMSE per epoch)", fontsize=15)
    plt.xlabel("Epoch"); plt.ylabel("RMSE improvement per epoch  (↑ is better)")
    plt.legend(title="Run", fontsize=9)
    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_lr = out_dir / "learning_rate_from_rmse_selected.png"
    plt.savefig(out_lr, dpi=220); plt.close()
    print(f"saved {out_lr}")

    # 2) Reference: combined RMSE curves
    plt.figure(figsize=(8.5, 5))
    drew = 0
    for split, fold in runs:
        run_dir = base_dir / split / f"fold_{fold}"
        mpath = find_metrics_file(run_dir)
        if not mpath:
            continue
        df = read_per_epoch(mpath)
        y = pick_val_rmse(df)
        if y is None:
            continue
        # mark best epoch (min RMSE)
        best_idx = int(np.nanargmin(y.values))
        best_epoch = int(df.iloc[best_idx]["epoch"])
        best_val = float(y.iloc[best_idx])
        label = f"{PRETTY.get(split, split)} (fold {fold})"
        plt.plot(df["epoch"], y, label=label, linewidth=2.2,
                 alpha=0.95, color=color_map.get(split))
        plt.axvline(best_epoch, linestyle="--", alpha=0.4, color=color_map.get(split))
        plt.scatter([best_epoch], [best_val], s=30, zorder=5, color=color_map.get(split))
        drew += 1
    if drew == 0:
        print("No RMSE curves plotted.")
    plt.title("Validation RMSE per Epoch: Selected Runs", fontsize=15)
    plt.xlabel("Epoch"); plt.ylabel("RMSE")
    plt.legend(title="Run", fontsize=9)
    plt.tight_layout()
    out_rmse = out_dir / "rmse_selected.png"
    plt.savefig(out_rmse, dpi=220); plt.close()
    print(f"saved {out_rmse}")

def parse_runs(s: str) -> list[tuple[str, int]]:
    runs = []
    for tok in [t for t in s.split(",") if t.strip()]:
        t = tok.strip().replace("=", ":").replace(" ", ":")
        parts = [p for p in t.split(":") if p]
        if len(parts) >= 2:
            runs.append((parts[0], int(parts[1])))
    return runs

def main():
    parser = argparse.ArgumentParser(description="Plot ΔRMSE-per-epoch (rate of learning) and RMSE for selected runs.")
    parser.add_argument("--base", type=Path,
        default=Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/2 - train_val_test_results"),
        help="Root containing <split>/fold_<k>/metrics_per_epoch*")
    parser.add_argument("--out", type=Path,
        default=Path("./p2_models/3 - analysis_test_results/figures/rsme_loss_visuals"),
        help="Where to save the PNGs")
    parser.add_argument("--runs", type=str,
        default="random:1,scaffold:1,butina:2,umap_kmeans:1,umap_ward:5",
        help="Comma-separated split:fold pairs (spaces/= ok)")
    parser.add_argument("--smooth", type=int, default=3,
        help="Rolling window (epochs) to smooth ΔRMSE; 1 disables smoothing")
    args = parser.parse_args()

    plot_combo(args.base, args.out, parse_runs(args.runs), smooth_window=args.smooth)

if __name__ == "__main__":
    main()