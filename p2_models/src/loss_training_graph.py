# plot_training_curves.py
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_curves(base_models_dir: Path, split: str, fold: int, metrics=("rmse","mae","mse")):
    """
    Make per-epoch plots for the chosen model: train & val curves, and mark best val.
    Saves PNGs next to the checkpoint.
    """
    fold_dir = Path(base_models_dir) / split / f"fold_{fold}"
    log_csv  = fold_dir / "metrics_per_epoch.csv"
    if not log_csv.exists():
        raise FileNotFoundError(f"Missing: {log_csv}")

    df = pd.read_csv(log_csv)

    # Some loggers write multiple rows per epoch (steps + epoch summary).
    # Take the last row per epoch to represent the epoch summary.
    df = df[df["epoch"].notna()].copy()
    df["epoch"] = df["epoch"].astype(int)
    per_epoch = df.groupby("epoch", as_index=False).last()

    def smooth(s, w=1):
        # rolling average for nicer curves; set w=1 to disable smoothing
        return s.rolling(w, min_periods=1, center=True).mean()

    for metric in metrics:
        tcol = f"train/{metric}"
        vcol = f"val/{metric}"
        have_any = (tcol in per_epoch.columns) or (vcol in per_epoch.columns)
        if not have_any:
            continue

        plt.figure()
        x = per_epoch["epoch"]

        if tcol in per_epoch.columns:
            plt.plot(x, smooth(per_epoch[tcol], w=1), label="train")

        if vcol in per_epoch.columns:
            yv = per_epoch[vcol]
            plt.plot(x, smooth(yv, w=1), label="val")
            # mark best val (minimum)
            best_idx = int(np.nanargmin(yv.values))
            best_epoch = int(per_epoch.iloc[best_idx]["epoch"])
            best_val   = float(yv.iloc[best_idx])
            plt.axvline(best_epoch, linestyle="--", alpha=0.7)
            plt.scatter([best_epoch], [best_val], s=30, zorder=5)
            plt.text(best_epoch, best_val, f"  best {metric.upper()}={best_val:.3f} @ {best_epoch}", va="bottom")

        plt.xlabel("Epoch")
        plt.ylabel(metric.upper())
        plt.title(f"{split} · fold {fold} — {metric.upper()} per epoch")
        plt.legend()
        plt.tight_layout()

        out_png = fold_dir / f"curve_{metric}.png"
        plt.savefig(out_png, dpi=200)
        plt.close()
        print(f"saved {out_png}")

    # Optional: if your CSV has explicit 'train/loss' and 'val/loss', plot those too.
    if {"train/loss","val/loss"}.intersection(per_epoch.columns):
        plt.figure()
        if "train/loss" in per_epoch.columns:
            plt.plot(per_epoch["epoch"], per_epoch["train/loss"], label="train/loss")
        if "val/loss" in per_epoch.columns:
            yv = per_epoch["val/loss"]
            plt.plot(per_epoch["epoch"], yv, label="val/loss")
            best_idx = int(np.nanargmin(yv.values))
            be = int(per_epoch.iloc[best_idx]["epoch"])
            bv = float(yv.iloc[best_idx])
            plt.axvline(be, linestyle="--", alpha=0.7)
            plt.scatter([be], [bv], s=30, zorder=5)
            plt.text(be, bv, f"  best loss={bv:.3f} @ {be}", va="bottom")
        plt.xlabel("Epoch"); plt.ylabel("LOSS"); plt.title(f"{split} · fold {fold} — Loss per epoch")
        plt.legend(); plt.tight_layout()
        out_png = Path(base_models_dir) / split / f"fold_{fold}" / "curve_loss.png"
        plt.savefig(out_png, dpi=200); plt.close()
        print(f"saved {out_png}")

if __name__ == "__main__":
    # Example usage:
    base = Path.home() / "Thesis" / "p2_models" / "models"
    plot_curves(base_models_dir=base, split="butina", fold=1, metrics=("rmse","mae"))
