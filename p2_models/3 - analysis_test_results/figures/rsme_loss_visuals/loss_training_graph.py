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

from __future__ import annotations
from pathlib import Path
import argparse
import gzip
import re
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

# Anything matching these is considered a learning-rate column
LR_REGEXES = [
    re.compile(r"(^|[\/\-\_:])lr($|[\/\-\_:]\d+?$)"),   # lr, lr-Adam, lr-0, foo/lr
    re.compile(r"learning[_/ \-]?rate"),               # learning_rate, learning-rate
]

METRIC_CANDIDATES = [
    "metrics_per_epoch.csv",
    "metrics_per_epoch.tsv",
    "metrics.csv",
    "metrics_per_step.csv",
    "metrics*.csv",
    "metrics*.csv.gz",
    "metrics*.tsv",
    "metrics*.tsv.gz",
    "metrics_per_epoch",  # no extension
]

def pretty_split(s: str) -> str:
    return PRETTY.get(s, s.replace("_", " ").title())

def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

def _sep_for(path: Path) -> str:
    # tsv if final suffix (or inner suffix for .gz) is .tsv
    suf = path.suffix.lower()
    inner = Path(path.stem).suffix.lower() if suf == ".gz" else suf
    return "\t" if inner == ".tsv" else ","

def _find_lr_cols(columns) -> list[str]:
    cols = []
    for c in columns:
        lc = str(c).lower()
        for rx in LR_REGEXES:
            if rx.search(lc):
                cols.append(c)
                break
    return cols

def _read_df(path: Path) -> pd.DataFrame:
    with _open_text(path) as f:
        return pd.read_csv(f, sep=_sep_for(path))

def _group_last_per_epoch(df: pd.DataFrame) -> pd.DataFrame:
    if "epoch" not in df.columns:
        raise ValueError('No "epoch" column found.')
    d = df[df["epoch"].notna()].copy()
    d["epoch"] = d["epoch"].astype(int)
    # keep last record per epoch
    return d.groupby("epoch", as_index=False).last()

def find_lr_source_df(run_dir: Path) -> tuple[pd.DataFrame, str] | tuple[None, None]:
    """
    Look for a file in run_dir that has both 'epoch' and at least one LR column.
    Returns (df, lr_col_name_to_use) or (None, None)
    """
    # First try specific filenames, then globs
    files = []
    for name in METRIC_CANDIDATES:
        if "*" in name:
            files += sorted(run_dir.glob(name))
        else:
            p = run_dir / name
            if p.exists():
                files.append(p)

    seen = set()
    for p in files:
        if p in seen: 
            continue
        seen.add(p)
        try:
            df = _read_df(p)
        except Exception:
            continue

        if "epoch" not in df.columns:
            continue

        lr_cols = _find_lr_cols(df.columns)
        if not lr_cols:
            # maybe it's per-epoch metrics with no LR; skip
            continue

        # Prefer a single stable column: choose the first sorted by name
        lr_cols_sorted = sorted(lr_cols, key=lambda x: str(x))
        lr_col = lr_cols_sorted[0]
        # collapse to last row per epoch (handles per-step logging)
        per_epoch = _group_last_per_epoch(df[["epoch", lr_col]].copy())
        # guard: drop rows with all-NaN LR
        if per_epoch[lr_col].notna().any():
            return per_epoch, lr_col

    return None, None

def plot_lr_combo(base_dir: Path, out_dir: Path, runs: list[tuple[str, int]], logy: bool = True) -> bool:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 5))
    reds = sns.color_palette("Reds", n_colors=len(SP_ORDER))
    color_map = {sp: reds[i] for i, sp in enumerate(SP_ORDER)}

    drew_any = False
    for split, fold in runs:
        run_dir = base_dir / split / f"fold_{fold}"
        df, lr_col = find_lr_source_df(run_dir)
        if df is None:
            print(f"[skip] no LR data found under {run_dir}")
            continue

        label = f"{pretty_split(split)} (fold {fold})"
        plt.plot(df["epoch"], df[lr_col], linewidth=2.2, alpha=0.95, label=label, color=color_map.get(split))
        drew_any = True

    if not drew_any:
        print("No learning-rate curves were plotted. Make sure LR is being logged.")
        plt.close()
        re
