# plot_25_rmse_sns_reds.py
from __future__ import annotations
from pathlib import Path
import argparse, re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

SP_ORDER = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
REQUIRED_FILE = "metrics_per_epoch.csv"

def pretty_split_name(s: str) -> str:
    # "umap_kmeans" -> "Umap kmeans"
    return s.replace("_", " ").capitalize()

def read_per_epoch(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df[df["epoch"].notna()].copy()
    df["epoch"] = df["epoch"].astype(int)
    return df.groupby("epoch", as_index=False).last()

def plot_one(base_models_dir: Path, out_dir: Path, split: str, fold: int, color):
    fold_dir = Path(base_models_dir) / split / f"fold_{fold}"
    csv_path = fold_dir / REQUIRED_FILE
    if not csv_path.exists():
        print(f"[skip] missing {csv_path}")
        return False

    per_epoch = read_per_epoch(csv_path)

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(6, 4))

    did_any = False

    # Train line (lighter / transparent)
    if "train/rmse" in per_epoch.columns:
        sns.lineplot(x=per_epoch["epoch"], y=per_epoch["train/rmse"],
                     color=color, alpha=0.55, linewidth=2)
        did_any = True

    # Val line (same hue, solid). Keep marker/line for best, but no text label.
    if "val/rmse" in per_epoch.columns:
        yv = per_epoch["val/rmse"]
        sns.lineplot(x=per_epoch["epoch"], y=yv,
                     color=color, linewidth=2.2)
        if yv.notna().any():
            best_idx = int(np.nanargmin(yv.values))
            best_epoch = int(per_epoch.iloc[best_idx]["epoch"])
            best_val   = float(yv.iloc[best_idx])
            plt.axvline(best_epoch, linestyle="--", color=color, alpha=0.7)
            plt.scatter([best_epoch], [best_val], s=35, zorder=5, color=color)
        did_any = True

    if not did_any:
        print(f"[skip] no RMSE columns in {csv_path}")
        plt.close()
        return False

    split_title = pretty_split_name(split)
    plt.title(f"{split_title} fold {fold} — RMSE per epoch", size=15)
    plt.tight_layout()
    plt.xlabel("Epoch", size=11)
    plt.ylabel("RMSE", size=11)
    # No legend
    plt.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"loss_rmse_{split}_fold_{fold}.png"
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"saved {out_png}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Create 25 RMSE figures (5 splits × 5 folds) with Reds palette, no legend, capitalized titles.")
    parser.add_argument("--base", type=Path,
        default=Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/models"),
        help="Path to .../p2_models/models")
    parser.add_argument("--out", type=Path,
        default=Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/analysis_outputs/figures"),
        help="Where to save the 25 PNGs")
    parser.add_argument("--folds", type=str, default="1,2,3,4,5",
        help="Comma-separated folds to plot (default 1..5)")
    args = parser.parse_args()

    folds = [int(x) for x in args.folds.split(",") if x.strip()]
    # Five distinct shades of red in the requested split order
    reds = sns.color_palette("Reds", n_colors=len(SP_ORDER))
    split_to_color = {sp: reds[i] for i, sp in enumerate(SP_ORDER)}

    made = 0
    for split in SP_ORDER:
        color = split_to_color[split]
        for f in folds:
            ok = plot_one(args.base, args.out, split, f, color)
            if ok:
                made += 1
    print(f"Completed: {made} figure(s).")

if __name__ == "__main__":
    main()
