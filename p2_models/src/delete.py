'''
Description: This script trains and validate Chemprop DMPNN model for pIC50 across multiple data-split strategies while saving 
best ad last checkpoints.

Requirements:
    * Split CSVs from Data Splitting preprocessing step: "~/p2_models/input_data/split_data/{split}" with required columns:
    "Smiles", "pIC50"

Procedure:
1) For each split in ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
   and each fold k=1-5
   a. Load train/val CSVs.
   b. Featurize molecules with Chemprop's SimpleMoleculeMolGraphFeaturizer.
   c. Normalize targets using the TRAIN set's scaler; apply the same scaler to val.
      (Predictions are unscaled back to the original pIC50 space via an UnscaleTransform.)
   d. Build a Chemprop DMPNN: BondMessagePassing to MeanAggregation to RegressionFFN,
      with metrics RMSE/MSE/MAE.
   e. Train with PyTorch Lightning:
        - Early stopping on "val/rmse" (patience=15)
        - Checkpoint the best model (min val/rmse)

OUTPUTS (per split and fold)
Directory: ~/p2_models/models/<split>/fold_<i>/
    * best.ckpt (Lightning checkpoint) and last.ckpt
    * TensorBoard logs (view with: tensorboard --logdir ~/p2_models/models/<split>)
'''
from pathlib import Path
import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
from lightning import pytorch as pl
from lightning.pytorch.callbacks import ModelCheckpoint, EarlyStopping
from lightning.pytorch.loggers import TensorBoardLogger, CSVLogger
from chemprop import data, featurizers, models, nn
from sklearn import metrics
from scipy.stats import spearmanr
import shutil

# Config
ALL_SPLITS = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
p2models_dir = Path.home() / "Thesis" / "p2_models"
BASE_SPLITS_DIR = p2models_dir / "input_data" / "split_data"
BASE_OUT_DIR = p2models_dir / "input_data"

# Inputs
SMILES_COL = "Smiles"
TARGET_COLS = ["pIC50"]

# Settings
NUM_WORKERS = 0
MAX_EPOCHS = 1
PATIENCE = 15
use_gpu = torch.cuda.is_available()

# Helpers to file and molecules/targets
def find_file(split_dir: Path, split: str, i: int, kind: str) -> Path:
    for name in (f"{split}_fold_{i}_{kind}.csv", f"{split}_fold{i}_{kind}.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {kind} file for fold {i} in {split_dir}")

def to_points(df):
    return [data.MoleculeDatapoint.from_smi(s, y)
            for s, y in zip(df[SMILES_COL].values, df[TARGET_COLS].values)]

# Metrics helpers
def _predict_flat(trainer: pl.Trainer, model, loader, ckpt_path=None) -> np.ndarray:
    with torch.no_grad():
        outs = trainer.predict(model=model, dataloaders=loader, ckpt_path=ckpt_path)
    flats = []
    for o in outs:
        x = o
        if isinstance(x, dict):
            for k in ("y_hat", "preds", "y_pred", "pred", "logits", "output", "out"):
                if k in x:
                    x = x[k]; break
            else:
                x = next(iter(x.values()))
        if isinstance(x, (list, tuple)):
            x = x[0]
        x = torch.as_tensor(x).detach().cpu().numpy()
        flats.append(np.ravel(x))
    return np.concatenate(flats, axis=0)

# Regression metrics helper
def _regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    y_true = np.asarray(y_true, float).ravel()
    y_pred = np.asarray(y_pred, float).ravel()
    mse = float(metrics.mean_squared_error(y_true, y_pred))
    rmse = float(np.sqrt(mse))
    mae = float(metrics.mean_absolute_error(y_true, y_pred))
    r2  = float(metrics.r2_score(y_true, y_pred))
    rho, _ = spearmanr(y_true, y_pred)
    return {"rmse": rmse, "mae": mae, "mse": mse, "r2": r2, "spearman_rho": float(rho), "n": int(len(y_true))}

# Helper to evaluate and save
def evaluate_and_save(trainer: pl.Trainer, model, train_set, val_set, df_train: pd.DataFrame, df_val: pd.DataFrame, ckpt_path: str, out_csv_path: Path):
    
    train_loader_eval = data.build_dataloader(train_set, num_workers=NUM_WORKERS, shuffle=False)
    val_loader_eval   = data.build_dataloader(val_set,   num_workers=NUM_WORKERS, shuffle=False)

    y_pred_train = _predict_flat(trainer, model, train_loader_eval, ckpt_path=ckpt_path)
    y_pred_val   = _predict_flat(trainer, model, val_loader_eval,   ckpt_path=ckpt_path)

    y_true_train = df_train[TARGET_COLS].to_numpy().ravel()
    y_true_val   = df_val[TARGET_COLS].to_numpy().ravel()

    m_train = _regression_metrics(y_true_train, y_pred_train)
    m_val   = _regression_metrics(y_true_val,   y_pred_val)

    rows = []
    rows.append({"dataset": "train", **m_train})
    rows.append({"dataset": "val",   **m_val})
    pd.DataFrame(rows).to_csv(out_csv_path, index=False)

# Training loop
for SPLIT in ALL_SPLITS:
    print(f"\n SPLIT: {SPLIT.upper()}", flush=True)

    SPLIT_DIR = BASE_SPLITS_DIR / SPLIT
    OUTPUTDIR = (BASE_OUT_DIR / SPLIT); OUTPUTDIR.mkdir(parents=True, exist_ok=True)

    for i in range(1, 6):
        print(f"\n===== Fold {i} ({SPLIT}) =====", flush=True)

        train_path = find_file(SPLIT_DIR, SPLIT, i, "train")
        val_path   = find_file(SPLIT_DIR, SPLIT, i, "val")

        df_train = pd.read_csv(train_path)
        df_val   = pd.read_csv(val_path)
        print(f"Loaded | train={len(df_train):,}  val={len(df_val):,}", flush=True)

        # Prepare data
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        train_set = data.MoleculeDataset(to_points(df_train), featurizer)
        scaler = train_set.normalize_targets()

        val_set  = data.MoleculeDataset(to_points(df_val), featurizer)
        val_set.normalize_targets(scaler)

        train_loader = data.build_dataloader(train_set, num_workers=NUM_WORKERS)
        val_loader   = data.build_dataloader(val_set,   num_workers=NUM_WORKERS, shuffle=False)

        # Model
        mp = nn.BondMessagePassing()
        agg = nn.MeanAggregation()
        out_tf = nn.UnscaleTransform.from_standard_scaler(scaler)
        ffn = nn.RegressionFFN(output_transform=out_tf)
        metric_list = [nn.metrics.RMSE(), nn.metrics.MSE(), nn.metrics.MAE()]
        model = models.MPNN(mp, agg, ffn, batch_norm=True, metrics=metric_list)

        fold_out = OUTPUTDIR / f"fold_{i}"
        fold_out.mkdir(parents=True, exist_ok=True)

        # Loggers: keep TB + add CSV
        tb_logger  = TensorBoardLogger(save_dir=str(OUTPUTDIR), name=f"fold_{i}")
        csv_logger = CSVLogger(save_dir=str(OUTPUTDIR),     name=f"fold_{i}")

        ckpt = ModelCheckpoint(
            dirpath=fold_out,
            filename="best",
            monitor="val/rmse",
            mode="min",
            save_top_k=1,
            save_last=True)
        es = EarlyStopping(monitor="val/rmse", mode="min", patience=PATIENCE)

        trainer = pl.Trainer(
            accelerator="gpu" if use_gpu else "auto",
            devices=1,
            precision="16-mixed" if use_gpu else "32-true",
            max_epochs=MAX_EPOCHS,
            callbacks=[ckpt, es],
            logger=[tb_logger, csv_logger],
            enable_checkpointing=True,
            enable_progress_bar=True)

        print(f"Training (max_epochs={MAX_EPOCHS}, patience={PATIENCE})...", flush=True)
        trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
        print(f"Best checkpoint: {ckpt.best_model_path}", flush=True)

        # Evaluate best ckpt on train and val: save summary CSV 
        summary_csv = fold_out / "final_eval_metrics.csv"
        evaluate_and_save(trainer, model, train_set, val_set, df_train, df_val, ckpt.best_model_path, summary_csv)
        print(f"Saved metrics summary to {summary_csv}", flush=True)

        # Copy per-epoch CSV logs next to checkpoint for convenience
        per_epoch_src = Path(csv_logger.log_dir) / "metrics.csv"
        per_epoch_dst = fold_out / "metrics_per_epoch.csv"
        if per_epoch_src.exists():
            shutil.copy(per_epoch_src, per_epoch_dst)
            print(f"Saved per-epoch metrics -> {per_epoch_dst}", flush=True)
        else:
            print("WARN: CSV per-epoch log not found; check logger path.", flush=True)

print(f"\nDone. Per-split logs & ckpts under: {BASE_OUT_DIR}/<split>/fold_<i>/")