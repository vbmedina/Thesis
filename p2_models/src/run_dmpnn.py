'''
Description: Train and validate Chemprop DMPNN model for pIC50 across multiple data-split strategies,
saving checkpoints (best and last).

Requirements:
    * Split CSVs from Data Splitting preprocessing step: "~/p2_models/input_data/split_data/{split}" with required columns:
    "Smiles", "pIC50"

Procedure:
1) For each split in ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
   and each fold i=1..5:
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
from lightning.pytorch.loggers import TensorBoardLogger
from chemprop import data, featurizers, models, nn

# Config
ALL_SPLITS = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
p2models_dir = Path.home() / "Thesis" / "p2_models"
BASE_SPLITS_DIR = p2models_dir / "input_data" / "split_data"
BASE_OUT_DIR = p2models_dir / "models"

# Inputs
SMILES_COL = "Smiles"
TARGET_COLS = ["pIC50"]

# Settings
NUM_WORKERS = 0
MAX_EPOCHS = 200
PATIENCE = 15
use_gpu = torch.cuda.is_available()

# Helpers: Find files
def find_file(split_dir: Path, split: str, i: int, kind: str) -> Path:
    for name in (f"{split}_fold_{i}_{kind}.csv", f"{split}_fold{i}_{kind}.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {kind} file for fold {i} in {split_dir}")

def to_points(df):
    return [data.MoleculeDatapoint.from_smi(s, y)
            for s, y in zip(df[SMILES_COL].values, df[TARGET_COLS].values)]

# Main loop for splits
for SPLIT in ALL_SPLITS:
    print(f"\n SPLIT: {SPLIT.upper()}", flush=True)

    SPLIT_DIR = BASE_SPLITS_DIR / SPLIT
    OUTPUTDIR = (BASE_OUT_DIR / SPLIT); OUTPUTDIR.mkdir(parents=True, exist_ok=True)

    # Find CSVs
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

        val_set  = data.MoleculeDataset(to_points(df_val), featurizer);  val_set.normalize_targets(scaler)

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
        logger = TensorBoardLogger(save_dir=str(OUTPUTDIR), name=f"fold_{i}")

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
            logger=logger,
            enable_checkpointing=True,
            enable_progress_bar=True)

        print(f"Training (max_epochs={MAX_EPOCHS}, patience={PATIENCE})...", flush=True)
        trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
        print(f"Best checkpoint: {ckpt.best_model_path}", flush=True)

print(f"\nDone. Per-split logs & ckpts under: {BASE_OUT_DIR}/<split>/fold_<i>/")