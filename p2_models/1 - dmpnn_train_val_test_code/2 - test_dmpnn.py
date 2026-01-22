'''
Description: Evaluates Chemprop D-MPNN best checkpoints on fold test sets and the
external sexual test set, and exports predictions + core metrics.

Core test metrics (imbalanced-aware)
- RMSE
- Spearman rho
- PR-AUC (Average Precision; threshold-free)
- HR@100 (top-100 hit rate at pIC50 ≥ threshold)
- EF@100 (enrichment factor at top-100 relative to test-set prevalence)
- Prevalence (fraction actives in each evaluated set)

Requirements
- Model checkpoints per split/fold at: ./p2_models/models/<split>/fold_<k>/
  - Accepted names (checked in this order): best-v1.ckpt / best.ckpt / last-v1.ckpt / last.ckpt
- Fold test CSVs from the data-splitting step: ./p2_models/input_data/split_data/{split}/
  - Filenames accepted: {split}_fold_{k}_test.csv or {split}_fold{k}_test.csv
- External evaluation CSV: ./p2_models/input_data/sexual_test.csv
- Required columns: "Smiles", "pIC50"

Outputs
- Per split/fold directory: ./p2_models/models/<split>/fold_<k>/
  * pred_vs_true_fold_test.csv
  * pred_vs_true_sexual_data.csv
- Aggregated table (per-fold rows) in: ./p2_models/3 - analysis_outputs/
  * Overall final test metrics
    (columns: split, fold, dataset, rmse, spearman_rho, pr_auc, hr_at_100, ef_at_100, prevalence)

Notes
- Evaluation is CPU-only and uses a fixed seed for determinism of non-GPU components.
'''

# Imports
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

# CPU only dir
DEVICE = torch.device("cpu") 
BASE_ROOT = Path("./p2_models").expanduser()
EXTERNAL_CSV = BASE_ROOT / "0 - input_data" / "sexual_test.csv"
OUT_DIR_NAME = "./3 - analysis_test_results"
SMILES_COL = "Smiles"
TARGET_COL = "pIC50"

THRESHOLD = 6.0 
K_HIT_RATE = 100 

NUM_WORKERS = 0
SEED = 2858808528

HR_AT_K_COL = f"hr_at_{K_HIT_RATE}"
EF_AT_K_COL = f"ef_at_{K_HIT_RATE}"

TARGET_MODELS = {
    "random": 1,
    "scaffold": 1,
    "butina": 2,
    "umap_kmeans": 1,
    "umap_ward": 5,
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

# Helper for data and heckpoints
def discover_splits_and_folds(base_root: Path) -> Dict[str, List[int]]:
    models_root = base_root / "2 - train_val_test_results"
    if not models_root.exists():
        raise FileNotFoundError(f"Missing models dir: {models_root}")

    splits: Dict[str, List[int]] = {}
    ckpt_names = ("best-v1.ckpt", "best.ckpt", "last-v1.ckpt", "last.ckpt")

    for split_dir in sorted(p for p in models_root.iterdir() if p.is_dir()):
        split = split_dir.name
        folds: List[int] = []
        for fd in sorted(split_dir.glob(f"fold_{TARGET_MODELS[split]}")):  # kept as-is (minimal change)
            m = re.fullmatch(r"fold_(\d+)", fd.name)
            if not m:
                continue
            fold_num = int(m.group(1))
            if any((fd / n).exists() for n in ckpt_names):
                folds.append(fold_num)
        if folds:
            splits[split] = sorted(folds)

    if not splits:
        raise FileNotFoundError(f"No split/fold checkpoints found under {models_root}")
    return splits

def find_test_csv(base_root: Path, split: str, fold: int) -> Path:
    split_dir = base_root / "0 - input_data" / "split_data" / split
    for name in (f"{split}_fold_{fold}_test.csv", f"{split}_fold{fold}_test.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Missing test CSV for split={split} fold={fold} in {split_dir}")

def pick_checkpoint(fold_dir: Path) -> Path:
    for name in ("best-v1.ckpt", "best.ckpt", "last-v1.ckpt", "last.ckpt"):
        p = fold_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No checkpoint found in {fold_dir} (checked best-v1.ckpt, best.ckpt, last-v1.ckpt, last.ckpt).")

# Build dataset
def smiles_to_mols_and_mask(df: pd.DataFrame, smiles_col: str) -> Tuple[List, np.ndarray]:
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
    raise RuntimeError("chemprop.data.build_dataloader is not available. Install a Chemprop build that provides one.")

# Predict
def predict(trainer: pl.Trainer, model: models.MPNN, loader) -> np.ndarray:
    def _to_flat_cpu_numpy(x):
        if isinstance(x, dict):
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
def compute_core_metrics(y_true: np.ndarray, y_pred: np.ndarray, thr: float, k: int):
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    n = len(y_true)

    # RMSE
    mse  = float(metrics.mean_squared_error(y_true, y_pred))
    rmse = float(np.sqrt(mse))

    # Spearman
    rho, _ = spearmanr(y_true, y_pred)
    spearman_rho = float(rho)

    # Binary labels at threshold
    yb = (y_true >= thr).astype(int)
    pos = int(yb.sum())
    prevalence = float(pos / n) if n > 0 else float("nan")

    # PR-AUC (Average Precision); handle degenerate single-class cases
    if pos == 0 or pos == n:
        pr_auc = float("nan")
    else:
        pr_auc = float(metrics.average_precision_score(yb, y_pred))

    # HR@K (top-K by predicted score)
    order = np.argsort(-y_pred)  # descending
    k_eff = int(min(k, n))
    hits_at_k = int(yb[order[:k_eff]].sum()) if k_eff > 0 else 0
    hr_at_k = float(hits_at_k / k_eff) if k_eff > 0 else float("nan")

    # EF@K = (HR@K / prevalence)
    ef_at_k = float(hr_at_k / prevalence) if prevalence > 0 else float("nan")

    return {
        "rmse": rmse,
        "spearman_rho": spearman_rho,
        "pr_auc": float(pr_auc),
        HR_AT_K_COL: hr_at_k,
        EF_AT_K_COL: ef_at_k,
        "prevalence": prevalence,
    }

# Main def
def main():
    cfg = Config()

    if not cfg.external_csv.exists():
        raise FileNotFoundError(f"External CSV not found: {cfg.external_csv}")
    split_folds = discover_splits_and_folds(cfg.base_root)

    out_root = cfg.base_root / cfg.out_dir_name
    out_root.mkdir(parents=True, exist_ok=True)

    pl.seed_everything(SEED, workers=True)

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
        models_dir = cfg.base_root / "2 - train_val_test_results" / split
        for fold in folds:
            fold_dir = models_dir / f"fold_{fold}"
            ckpt = pick_checkpoint(fold_dir)
            print(f"  - fold {fold}: using {ckpt}")

            # Data
            test_csv = find_test_csv(cfg.base_root, split, fold)
            df_test = pd.read_csv(test_csv)
            df_ext  = pd.read_csv(cfg.external_csv)

            # Build datasets
            dtest, mask_t = build_dataset_from_mols(df_test, cfg.smiles_col, cfg.target_col)
            dext,  mask_e = build_dataset_from_mols(df_ext,  cfg.smiles_col, cfg.target_col)

            test_loader = make_dataloader(dtest, shuffle=False, num_workers=NUM_WORKERS)
            ext_loader  = make_dataloader(dext,  shuffle=False, num_workers=NUM_WORKERS)

            # model & predict (CPU)
            model = models.MPNN.load_from_checkpoint(str(ckpt), map_location=DEVICE)
            model.to(DEVICE)
            model.eval()

            # Filter truths to valid rows
            y_true_t = df_test[cfg.target_col].to_numpy()[mask_t]
            y_true_e = df_ext[cfg.target_col].to_numpy()[mask_e]

            y_pred_t = predict(trainer, model, test_loader)
            y_pred_e = predict(trainer, model, ext_loader)

            assert len(y_true_t) == len(y_pred_t), "Mismatch test y_true vs y_pred after SMILES filtering."
            assert len(y_true_e) == len(y_pred_e), "Mismatch external y_true vs y_pred after SMILES filtering."

            # save per-fold predictions (only valid rows)
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
                row.update(compute_core_metrics(yt, yp, cfg.threshold, K_HIT_RATE))
                per_fold_rows.append(row)

    if not per_fold_rows:
        raise SystemExit("No (split, fold) pairs produced results. Check your folders and filenames.")

    # Save per-fold metrics table
    per_fold_df = pd.DataFrame(per_fold_rows).sort_values(["split", "fold", "dataset"])
    (out_root / "final_test_metrics.csv").write_text(per_fold_df.to_csv(index=False))
    print(f"\nSaved to {out_root / 'final_test_metrics.csv'}")

    # Mean and SD
    agg = (
        per_fold_df
        .groupby(["split", "dataset"])
        .agg(
            rmse_mean=("rmse", "mean"), rmse_sd=("rmse", "std"),
            spearman_rho_mean=("spearman_rho", "mean"), spearman_rho_sd=("spearman_rho", "std"),
            pr_auc_mean=("pr_auc", "mean"), pr_auc_sd=("pr_auc", "std"),
            **{f"{HR_AT_K_COL}_mean": (HR_AT_K_COL, "mean"),
               f"{HR_AT_K_COL}_sd":   (HR_AT_K_COL, "std"),
               f"{EF_AT_K_COL}_mean": (EF_AT_K_COL, "mean"),
               f"{EF_AT_K_COL}_sd":   (EF_AT_K_COL, "std"),
               "prevalence_mean": ("prevalence", "mean"),
               "prevalence_sd":   ("prevalence", "std")}
        )
        .reset_index()
    )

if __name__ == "__main__":
    main()