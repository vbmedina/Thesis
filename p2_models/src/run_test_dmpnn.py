#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import re
import numpy as np
import pandas as pd
import torch
import lightning.pytorch as pl
from chemprop import data, models
from sklearn import metrics
from scipy.stats import spearmanr
from rdkit import Chem 

# CPU-only settings & paths
DEVICE = torch.device("cpu") 
BASE_ROOT = Path("./p2_models").expanduser()
EXTERNAL_CSV = BASE_ROOT / "input_data" / "sexual_test.csv"
OUT_DIR_NAME = "analysis_outputs"
SMILES_COL = "Smiles"
TARGET_COL = "pIC50"
THRESHOLD = 6.0
NUM_WORKERS = 0
SEED = 1337

TARGET_MODELS = {
    "random":      1,
    "scaffold":    1,
    "butina":      2,
    "umap_kmeans": 1,
    "umap_ward":   5,
}

# Config
@dataclass
class Config:
    base_root: Path = BASE_ROOT
    external_csv: Path = EXTERNAL_CSV
    out_dir_name: str = OUT_DIR_NAME
    smiles_col: str = SMILES_COL
    target_col: str = TARGET_COL
    threshold: float = THRESHOLD

# Helpers for finding paths to best-v1.ckpt
def discover_splits_and_folds(base_root: Path) -> Dict[str, List[int]]:
    """Find all <split>/fold_<n>/best-v1.ckpt under <base_root>/models."""
    models_root = base_root / "models"
    if not models_root.exists():
        raise FileNotFoundError(f"Missing models dir: {models_root}")

    splits: Dict[str, List[int]] = {}
    for split_dir in sorted(p for p in models_root.iterdir() if p.is_dir()):
        split = split_dir.name
        folds: List[int] = []
        for fd in sorted(split_dir.glob(f"fold_{TARGET_MODELS[split]}")):
            m = re.fullmatch(r"fold_(\d+)", fd.name)
            if not m:
                continue
            fold_num = int(m.group(1))
            if (fd / "best-v1.ckpt").exists():
                folds.append(fold_num)
        if folds:
            splits[split] = sorted(folds)

    if not splits:
        raise FileNotFoundError(f"No split/fold checkpoints found under {models_root}")
    return splits

# Helpers for finding paths to test
def find_test_csv(base_root: Path, split: str, fold: int) -> Path:
    split_dir = base_root / "input_data" / "split_data" / split
    for name in (f"{split}_fold_{fold}_test.csv", f"{split}_fold{fold}_test.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Missing test CSV for split={split} fold={fold} in {split_dir}")

# Data and Test Model
def smiles_to_mols_and_mask(df: pd.DataFrame, smiles_col: str) -> Tuple[List, np.ndarray]:
    # Convert SMILES into RDKit Mol
    smi = df[smiles_col].astype(str).tolist()
    mols: List = []
    mask = np.zeros(len(smi), dtype=bool)
    bad_ix = []
    for i, s in enumerate(smi):
        m = Chem.MolFromSmiles(s)
        if m is None:
            bad_ix.append(i)
        else:
            mask[i] = True
            mols.append(m)
    if bad_ix:
        print(f" [warn] Dropping {len(bad_ix)} invalid SMILES (first few idx: {bad_ix[:5]})")
    return mols, mask

def build_dataset_from_mols(df: pd.DataFrame, smiles_col: str, target_col: str):
    # Build a Chemprop MoleculeDataset from RDKit Mol
    from chemprop.data import MoleculeDataset, MoleculeDatapoint

    mols, mask = smiles_to_mols_and_mask(df, smiles_col)
    y_all = df[target_col].to_numpy()
    y = y_all[mask]

    dpoints = []
    for mol, yv in zip(mols, y):
        tgt = [float(yv)] if pd.notna(yv) else [None]
        dpoints.append(MoleculeDatapoint(mol, tgt))
    return MoleculeDataset(dpoints), mask

def make_dataloader(dset, shuffle=False, num_workers=NUM_WORKERS):
    if hasattr(data, "build_dataloader"):
        return data.build_dataloader(dset, shuffle=shuffle, num_workers=num_workers)
    raise RuntimeError("chemprop.data.build_dataloader is not available, and this build lacks MoleculeDataLoader." 
                       "Install a Chemprop version that provides one.")

# Run Lightning predict
def predict(trainer: pl.Trainer, model: models.MPNN, loader) -> np.ndarray:
    def _to_flat_cpu_numpy(x):
        if isinstance(x, dict):
            # try common keys first, else take the first value
            for k in ("y_hat", "preds", "y_pred", "pred", "logits", "output", "out"):
                if k in x:
                    x = x[k]
                    break
            else:
                x = next(iter(x.values()))
        if isinstance(x, (list, tuple)):
            x = x[0]
        x = torch.as_tensor(x).detach().cpu().numpy()
        return np.ravel(x)

    with torch.no_grad():
        outs = trainer.predict(model=model, dataloaders=loader, ckpt_path=None)

    flats = [_to_flat_cpu_numpy(b) for b in outs]
    return np.concatenate(flats, axis=0)

# Metrics 

def compute_core_metrics(y_true: np.ndarray, y_pred: np.ndarray, thr: float):
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()

    # RMSE
    mse  = float(metrics.mean_squared_error(y_true, y_pred))
    rmse = float(np.sqrt(mse))

    # Spearman
    rho, _ = spearmanr(y_true, y_pred)
    spearman_rho = float(rho)

    # PR-AUC 
    yb = (y_true >= thr).astype(int)
    pr_auc = float(metrics.average_precision_score(yb, y_pred))

    # MCC @ 6 pIC50 threshold
    yp = (y_pred >= thr).astype(int)
    tp = int(((yb == 1) & (yp == 1)).sum())
    tn = int(((yb == 0) & (yp == 0)).sum())
    fp = int(((yb == 0) & (yp == 1)).sum())
    fn = int(((yb == 1) & (yp == 0)).sum())
    num = tp * tn - fp * fn
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc_at_6 = float(num / den) if den > 0 else 0.0

    return {"rmse": rmse, "spearman_rho": spearman_rho, "pr_auc": pr_auc, "mcc_at_6": mcc_at_6}

# Main loop

def main():
    cfg = Config()

    if not cfg.external_csv.exists():
        raise FileNotFoundError(f"External CSV not found: {cfg.external_csv}")
    split_folds = discover_splits_and_folds(cfg.base_root)

    out_root = cfg.base_root / cfg.out_dir_name
    out_root.mkdir(parents=True, exist_ok=True)

    pl.seed_everything(SEED, workers=True)

    # Force CPU Testing
    trainer = pl.Trainer(
        accelerator="cpu",
        devices=1,
        logger=False,
        enable_progress_bar=False,
        precision=32,
    )

    per_fold_rows = []

    for split, folds in split_folds.items():
        print(f"\n=== SPLIT: {split} ===")
        models_dir = cfg.base_root / "models" / split
        for fold in folds:
            fold_dir = models_dir / f"fold_{fold}"
            ckpt = fold_dir / "best-v1.ckpt"
            print(f"  - fold {fold}: {ckpt}")

            # Data
            test_csv = find_test_csv(cfg.base_root, split, fold)
            df_test = pd.read_csv(test_csv)
            df_ext  = pd.read_csv(cfg.external_csv)

            # Build datasets from RDKit Mols
            dtest, mask_t = build_dataset_from_mols(df_test, cfg.smiles_col, cfg.target_col)
            dext,  mask_e = build_dataset_from_mols(df_ext,  cfg.smiles_col, cfg.target_col)

            test_loader = make_dataloader(dtest, shuffle=False, num_workers=NUM_WORKERS)
            ext_loader  = make_dataloader(dext,  shuffle=False, num_workers=NUM_WORKERS)

            # model & predict (CPU)
            model = models.MPNN.load_from_checkpoint(str(ckpt), map_location=DEVICE)
            model.to(DEVICE)
            model.eval()

            # Filter truths to valid rows before computing metrics
            y_true_t = df_test[cfg.target_col].to_numpy()[mask_t]
            y_true_e = df_ext[cfg.target_col].to_numpy()[mask_e]

            y_pred_t = predict(trainer, model, test_loader)
            y_pred_e = predict(trainer, model, ext_loader)

            # Sanity check
            assert len(y_true_t) == len(y_pred_t), "Mismatch test y_true vs y_pred after SMILES filtering."
            assert len(y_true_e) == len(y_pred_e), "Mismatch external y_true vs y_pred after SMILES filtering."

            # save per-fold predictions next to checkpoint (only valid rows)
            df_test_valid = df_test.loc[mask_t, [cfg.smiles_col]].copy()
            df_test_valid["y_true"] = y_true_t
            df_test_valid["y_pred"] = y_pred_t
            df_test_valid.to_csv(fold_dir / "pred_vs_true_fold_test.csv", index=False)

            df_ext_valid = df_ext.loc[mask_e, [cfg.smiles_col]].copy()
            df_ext_valid["y_true"] = y_true_e
            df_ext_valid["y_pred"] = y_pred_e
            df_ext_valid.to_csv(fold_dir / "pred_vs_true_sexual_data.csv", index=False)

            # Record invalid SMILES rows
            inv_test = df_test.loc[~mask_t, [cfg.smiles_col]]
            inv_ext  = df_ext.loc[~mask_e,  [cfg.smiles_col]]
            if len(inv_test):
                inv_test.to_csv(fold_dir / "invalid_smiles_fold_test.csv", index=False)
            if len(inv_ext):
                inv_ext.to_csv(fold_dir / "invalid_smiles_sexual_data.csv", index=False)

            # Metrics for both datasets
            for tag, yt, yp in (("fold_test", y_true_t, y_pred_t), ("sexual_data", y_true_e, y_pred_e)):
                row = {"split": split, "fold": fold, "dataset": tag}
                row.update(compute_core_metrics(yt, yp, cfg.threshold))
                per_fold_rows.append(row)

    if not per_fold_rows:
        raise SystemExit("No (split, fold) pairs produced results. Check your folders and filenames.")

    # Aggregate and save tables 
    per_fold_df = pd.DataFrame(per_fold_rows).sort_values(["split", "fold", "dataset"])
    (out_root / "all_metrics_by_split_fold_core.csv").write_text(per_fold_df.to_csv(index=False))
    print(f"\nSaved to {out_root / 'all_metrics_by_split_fold_core.csv'}")

    agg = (
        per_fold_df
        .groupby(["split", "dataset"])
        .agg(rmse_mean=("rmse", "mean"), rmse_sd=("rmse", "std"), spearman_rho_mean=("spearman_rho", "mean"), 
            spearman_rho_sd=("spearman_rho", "std"), pr_auc_mean=("pr_auc", "mean"), pr_auc_sd=("pr_auc", "std"),
            mcc_at_6_mean=("mcc_at_6", "mean"),      mcc_at_6_sd=("mcc_at_6", "std"))
        .reset_index())

    wide = agg.pivot(index="split", columns="dataset")
    delta_rmse_mean = wide["rmse_mean"]["sexual_data"] - wide["rmse_mean"]["fold_test"]

    def name(metric_stat: str, dataset: str) -> str:
        return f"{dataset}_{metric_stat}"

    main_df = pd.DataFrame({"split": wide.index})
    for metric_stat in ["rmse_mean", "rmse_sd",
                        "spearman_rho_mean", "spearman_rho_sd",
                        "pr_auc_mean", "pr_auc_sd",
                        "mcc_at_6_mean", "mcc_at_6_sd"]:
        for dataset in ["fold_test", "sexual_data"]:
            main_df[name(metric_stat, dataset)] = wide[metric_stat][dataset].values

    main_df["delta_rmse_mean"] = delta_rmse_mean.values

    ordered = ["split",
               "fold_test_rmse_mean", "fold_test_rmse_sd",
               "fold_test_spearman_rho_mean", "fold_test_spearman_rho_sd",
               "fold_test_pr_auc_mean", "fold_test_pr_auc_sd",
               "fold_test_mcc_at_6_mean", "fold_test_mcc_at_6_sd",
               "sexual_data_rmse_mean", "sexual_data_rmse_sd",
               "sexual_data_spearman_rho_mean", "sexual_data_spearman_rho_sd",
               "sexual_data_pr_auc_mean", "sexual_data_pr_auc_sd",
               "sexual_data_mcc_at_6_mean", "sexual_data_mcc_at_6_sd",
               "delta_rmse_mean"]
    main_df = main_df[ordered]
    (out_root / "main_table_core_metrics.csv").write_text(main_df.to_csv(index=False))
    print(f"Saved to {out_root / 'main_table_core_metrics.csv'}")

    def fmt(m, s):
        return f"{m:.3f} ± {s:.3f}" if np.isfinite(m) and np.isfinite(s) else "nan ± nan"

    metrics_rows = []
    for _, row in main_df.iterrows():
        metrics_rows.append({
            "split": row["split"],
            "fold_test_rmse": fmt(row["fold_test_rmse_mean"], row["fold_test_rmse_sd"]),
            "fold_test_spearman_rho": fmt(row["fold_test_spearman_rho_mean"], row["fold_test_spearman_rho_sd"]),
            "fold_test_pr_auc": fmt(row["fold_test_pr_auc_mean"], row["fold_test_pr_auc_sd"]),
            "fold_test_mcc_at_6": fmt(row["fold_test_mcc_at_6_mean"], row["fold_test_mcc_at_6_sd"]),
            "sexual_data_rmse": fmt(row["sexual_data_rmse_mean"], row["sexual_data_rmse_sd"]),
            "sexual_data_spearman_rho": fmt(row["sexual_data_spearman_rho_mean"], row["sexual_data_spearman_rho_sd"]),
            "sexual_data_pr_auc": fmt(row["sexual_data_pr_auc_mean"], row["sexual_data_pr_auc_sd"]),
            "sexual_data_mcc_at_6": fmt(row["sexual_data_mcc_at_6_mean"], row["sexual_data_mcc_at_6_sd"]),
            "delta_rmse_mean": row["delta_rmse_mean"],
        })
    metrics_df = pd.DataFrame(metrics_rows)
    (out_root / "main_table_core_metrics.csv").write_text(metrics_df.to_csv(index=False))
    print(f"Saved to {out_root / 'main_table_core_metrics.csv'}")


if __name__ == "__main__":
    main()
