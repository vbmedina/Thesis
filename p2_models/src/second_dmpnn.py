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
p2models_dir = Path.home() / "p2_models"
BASE_SPLITS_DIR = p2models_dir / "input_data" / "split_data"
BASE_OUT_DIR = p2models_dir / "models" 

# Sexual Test Data
EXTERNAL_TEST_CSV = p2models_dir / "input_data" / "sexual_test.csv"

# Inputs
SMILES_COL = "Smiles"
TARGET_COLS = ["pIC50"]

# Settings
NUM_WORKERS = 0
MAX_EPOCHS = 200
PATIENCE = 15
use_gpu = torch.cuda.is_available()

# Helpers
# Finding file paths
def find_file(split_dir: Path, split: str, i: int, kind: str) -> Path:
    for name in (f"{split}_fold_{i}_{kind}.csv", f"{split}_fold{i}_{kind}.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {kind} file for fold {i} in {split_dir}")

def to_points(df):
    return [data.MoleculeDatapoint.from_smi(s, y)
            for s, y in zip(df[SMILES_COL].values, df[TARGET_COLS].values)]

def pull_rmse(d):
    for k, v in d.items():
        if "rmse" in k.lower():
            return float(v)
    return np.nan

# Ploting True vs Predicted Values per fold
def save_pred_plot(trainer, model, loader, df_true, out_png, split_name: str, fold: int, dataset_label: str):
    with torch.no_grad():
        preds_batches = trainer.predict(model=model, dataloaders=loader, ckpt_path="best")
    y_hat = np.concatenate([p.detach().cpu().numpy().ravel() for p in preds_batches])
    y_true = df_true["pIC50"].to_numpy().ravel()

    # fit + stats
    m, b = np.polyfit(y_true, y_hat, 1)
    r = np.corrcoef(y_true, y_hat)[0, 1]
    r2 = float(r * r)
    rmse = float(np.sqrt(np.mean((y_hat - y_true) ** 2)))

    lo = float(min(y_true.min(), y_hat.min()))
    hi = float(max(y_true.max(), y_hat.max()))
    xline = np.linspace(lo, hi, 100)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=170)
    ax.scatter(y_true, y_hat, s=14, alpha=0.6, color="#d92525", marker="o",
               linewidths=0, edgecolors="none",
               label=f"Data (n={len(y_true):,}, RMSE={rmse:.3f})")
    ax.plot(xline, m * xline + b, color="black", lw=2,
            label=f"Fit: y = {m:.2f}x + {b:.2f}  (R²={r2:.3f})")
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=1.25, ls="--", label="Ideal: y = x")

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("True pIC50"); ax.set_ylabel("Predicted pIC50")
    ax.set_title(f"True vs Predicted Values on {dataset_label}: {split_name.upper()} Fold {fold}", pad=15, fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.25)
    ax.legend(loc="best", frameon=True, facecolor="white", framealpha=0.95)
    fig.tight_layout(); fig.savefig(out_png, dpi=300); plt.close(fig)

    pd.DataFrame({"Smiles": df_true[SMILES_COL], "y_true": y_true, "y_pred": y_hat}) \
      .to_csv(Path(out_png).with_suffix(".csv"), index=False)

# Main loop for splits
df_external = pd.read_csv(EXTERNAL_TEST_CSV)
all_rows = []

for SPLIT in ALL_SPLITS:
    print(f"\n SPLIT: {SPLIT.upper()}", flush=True)

    SPLIT_DIR = BASE_SPLITS_DIR / SPLIT
    OUTPUTDIR = (BASE_OUT_DIR / SPLIT); OUTPUTDIR.mkdir(parents=True, exist_ok=True)

    per_split_rows = []

    for i in range(1, 6):
        print(f"\n===== Fold {i} ({SPLIT}) =====", flush=True)

        train_path = find_file(SPLIT_DIR, SPLIT, i, "train")
        val_path   = find_file(SPLIT_DIR, SPLIT, i, "val")
        test_path  = find_file(SPLIT_DIR, SPLIT, i, "test")

        df_train = pd.read_csv(train_path)
        df_val   = pd.read_csv(val_path)
        df_test  = pd.read_csv(test_path)
        print(f"Loaded | train={len(df_train):,}  val={len(df_val):,}  "
              f"test={len(df_test):,}  external={len(df_external):,}", flush=True)

        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

        train_set = data.MoleculeDataset(to_points(df_train), featurizer)
        scaler = train_set.normalize_targets()

        val_set  = data.MoleculeDataset(to_points(df_val), featurizer);  val_set.normalize_targets(scaler)
        test_set = data.MoleculeDataset(to_points(df_test), featurizer); test_set.normalize_targets(scaler)
        ext_set  = data.MoleculeDataset(to_points(df_external), featurizer); ext_set.normalize_targets(scaler)

        train_loader = data.build_dataloader(train_set, num_workers=NUM_WORKERS)
        val_loader   = data.build_dataloader(val_set,   num_workers=NUM_WORKERS, shuffle=False)
        test_loader  = data.build_dataloader(test_set,  num_workers=NUM_WORKERS, shuffle=False)
        ext_loader   = data.build_dataloader(ext_set,   num_workers=NUM_WORKERS, shuffle=False)

        mp = nn.BondMessagePassing()
        agg = nn.MeanAggregation()
        out_tf = nn.UnscaleTransform.from_standard_scaler(scaler)
        ffn = nn.RegressionFFN(output_transform=out_tf)
        metric_list = [nn.metrics.RMSE(), nn.metrics.MSE(), nn.metrics.MAE()]  # logs 'val/rmse' etc.
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
            save_last=True,
        )
        es = EarlyStopping(monitor="val/rmse", mode="min", patience=PATIENCE)

        trainer = pl.Trainer(
            accelerator="gpu" if use_gpu else "auto",
            devices=1,
            precision="16-mixed" if use_gpu else "32-true",
            max_epochs=MAX_EPOCHS,
            callbacks=[ckpt, es],
            logger=logger,
            enable_checkpointing=True,
            enable_progress_bar=True,
        )
        
        # Print statement for training and checkpoints comparing test outputs
        print(f"Training (max_epochs={MAX_EPOCHS}, patience={PATIENCE})...", flush=True)
        trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
        print(f"Best checkpoint: {ckpt.best_model_path}", flush=True)

        test_metrics = trainer.test(dataloaders=test_loader, ckpt_path="best")[0]
        ext_metrics  = trainer.test(dataloaders=ext_loader,  ckpt_path="best")[0]
        test_rmse = pull_rmse(test_metrics)
        ext_rmse  = pull_rmse(ext_metrics)
        print(f"Fold test RMSE={test_rmse:.4f} | External RMSE={ext_rmse:.4f}", flush=True)

        # Visuals
        save_pred_plot(trainer, model, test_loader, df_test,
                       fold_out / "pred_vs_true_fold_test.png", SPLIT, i, "Fold Test Data")
        save_pred_plot(trainer, model, ext_loader,  df_external,
                       fold_out / "pred_vs_true_sexual_data.png", SPLIT, i, "Sexual Data")

        row = {"split": SPLIT, "fold": i, "fold_test_rmse": test_rmse, "external_rmse": ext_rmse}
        per_split_rows.append(row)
        all_rows.append(row)

    # Save per-split summary
    per_split_df = pd.DataFrame(per_split_rows).set_index(["split", "fold"])
    per_split_df.to_csv(OUTPUTDIR / "cv_results_with_external.csv")
    print(f"\n>>> {SPLIT} summary:\n", per_split_df.groupby(level=0).mean())

# Save combined summary across all splits
all_df = pd.DataFrame(all_rows).set_index(["split", "fold"]).sort_index()
all_df.to_csv(BASE_OUT_DIR / "cv_results_with_external_ALL.csv")
print("\n====================  ALL SPLITS — MEANS  ====================")
print(all_df.groupby(level=0)[["fold_test_rmse", "external_rmse"]].mean())
print(f"\nCombined summary saved to: {BASE_OUT_DIR / 'cv_results_with_external_ALL.csv'}")
print(f"Per-split logs & ckpts under: {BASE_OUT_DIR}/<split>/fold_<i>/")
