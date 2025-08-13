# Activate molml env
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
SPLIT = "umap"  # or: "random", "scaffold", "butina", "umap_hdb"

# Directory to splits
p2models_dir = Path.home() / "Thesis/p2_models"
SPLIT_DIR = p2models_dir / "data" / "splits" / SPLIT
OUTPUTDIR = p2models_dir / "models" / "checkpoints" / SPLIT
OUTPUTDIR.mkdir(parents=True, exist_ok=True)

# Directory to sexual data test set
EXTERNAL_TEST_CSV = p2models_dir / "data" / "splits" / "sexual_test.csv"

# Find targets
SMILES_COL = "Smiles"
TARGET_COLS = ["pIC50"]

# Settings for time and compute
NUM_WORKERS = 0
MAX_EPOCHS = 1
PATIENCE = 15
use_gpu = torch.cuda.is_available()

# Path helper for csv's
def find_file(split_dir: Path, split: str, i: int, kind: str) -> Path:
    for name in (f"{split}_fold_{i}_{kind}.csv", f"{split}_fold{i}_{kind}.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {kind} file for fold {i} in {split_dir}")

# Load CSV for sexual test set
df_ext_master = pd.read_csv(EXTERNAL_TEST_CSV)

# Results storage for later
rows = []

# Training and validation over folds
for i in range(1, 6):
    print(f"\n===== Fold {i} ({SPLIT}) =====", flush=True)

    train_path = find_file(SPLIT_DIR, SPLIT, i, "train")
    val_path   = find_file(SPLIT_DIR, SPLIT, i, "val")
    test_path  = find_file(SPLIT_DIR, SPLIT, i, "test")

    df_train = pd.read_csv(train_path)
    df_val   = pd.read_csv(val_path)
    df_test  = pd.read_csv(test_path)

    print(f"Loaded | train={len(df_train):,}  val={len(df_val):,}  "
          f"test={len(df_test):,}  external={len(df_ext_master):,}", flush=True)

    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

    def to_points(df):
        return [data.MoleculeDatapoint.from_smi(s, y)
                for s, y in zip(df[SMILES_COL].values, df[TARGET_COLS].values)]

    # Datasets - train, validation, test, and sexual test sets
    train_set = data.MoleculeDataset(to_points(df_train), featurizer)
    scaler = train_set.normalize_targets()

    val_set  = data.MoleculeDataset(to_points(df_val), featurizer)
    val_set.normalize_targets(scaler)

    test_set = data.MoleculeDataset(to_points(df_test), featurizer)
    test_set.normalize_targets(scaler)

    ext_set  = data.MoleculeDataset(to_points(df_ext_master), featurizer)
    ext_set.normalize_targets(scaler)

    # Dataloaders
    train_loader = data.build_dataloader(train_set, num_workers=NUM_WORKERS)
    val_loader   = data.build_dataloader(val_set,   num_workers=NUM_WORKERS, shuffle=False)
    test_loader  = data.build_dataloader(test_set,  num_workers=NUM_WORKERS, shuffle=False)
    ext_loader   = data.build_dataloader(ext_set,   num_workers=NUM_WORKERS, shuffle=False)

    # Model
    mp = nn.BondMessagePassing()
    agg = nn.MeanAggregation()  # change to SumAggregation/NormAggregation to compare
    out_tf = nn.UnscaleTransform.from_standard_scaler(scaler)
    ffn = nn.RegressionFFN(output_transform=out_tf)
    metric_list = [nn.metrics.RMSE(), nn.metrics.MSE(), nn.metrics.MAE()]  # logs keys like 'val/rmse'
    model = models.MPNN(mp, agg, ffn, batch_norm=True, metrics=metric_list)

    # Per-fold output and logger
    fold_out = OUTPUTDIR / f"fold_{i}"
    fold_out.mkdir(parents=True, exist_ok=True)

    # Try TensorBoard; if not available, fall back to CSVLogger
    try:
        from lightning.pytorch.loggers import TensorBoardLogger
        logger = TensorBoardLogger(save_dir=str(OUTPUTDIR), name=f"fold_{i}")
        print("∙ Using TensorBoard logger", flush=True)
    except Exception:
        from lightning.pytorch.loggers import CSVLogger
        logger = CSVLogger(save_dir=str(OUTPUTDIR), name=f"fold_{i}")
        print("∙ TensorBoard not available — using CSVLogger", flush=True)

    # Callbacks for modelcheck point and save top result
    ckpt = ModelCheckpoint(
        dirpath=fold_out,
        filename="best",
        monitor="val/rmse",
        mode="min",
        save_top_k=1,
        save_last=True,
    )
    es = EarlyStopping(monitor="val/rmse", mode="min", patience=PATIENCE)
#
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

    print(f"Training (max_epochs={MAX_EPOCHS}, patience={PATIENCE})...", flush=True)
    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
    print(f"Best checkpoint: {ckpt.best_model_path}", flush=True)

    # Evaluate on fold test
    test_metrics = trainer.test(dataloaders=test_loader, ckpt_path="best")[0]
    print(f"Fold test RMSE: " +
          ", ".join([f"{k}={v:.4f}" for k, v in test_metrics.items() if 'rmse' in k.lower()]), flush=True)

    # Evaluate on sexual test set
    ext_metrics = trainer.test(dataloaders=ext_loader, ckpt_path="best")[0]
    print(f"External test RMSE: " +
          ", ".join([f"{k}={v:.4f}" for k, v in ext_metrics.items() if 'rmse' in k.lower()]), flush=True)

    # Visuals for predicted vs true for both test sets
    def save_pred_plot(loader, df_true, out_png, split_name: str, fold: int, dataset_label: str):
        with torch.no_grad():
            preds_batches = trainer.predict(model=model, dataloaders=loader, ckpt_path="best")

        y_hat = np.concatenate([p.detach().cpu().numpy().ravel() for p in preds_batches])
        y_true = df_true["pIC50"].to_numpy().ravel()

        # fit + summary stats
        m, b = np.polyfit(y_true, y_hat, 1)
        r = np.corrcoef(y_true, y_hat)[0, 1]
        r2 = float(r * r)
        rmse = float(np.sqrt(np.mean((y_hat - y_true) ** 2)))

        lo = float(min(y_true.min(), y_hat.min()))
        hi = float(max(y_true.max(), y_hat.max()))
        xline = np.linspace(lo, hi, 100)

        fig, ax = plt.subplots(figsize=(5, 5), dpi=170)

        # RED dots (each point is a circle)
        ax.scatter(
            y_true, y_hat,
            s=14, alpha=0.6, color="#d92525",
            marker="o", linewidths=0, edgecolors="none",
            label=f"Data (n={len(y_true):,}, RMSE={rmse:.3f})"
        )

        # BLACK regression line
        ax.plot(xline, m * xline + b, color="black", lw=2,
                label=f"Fit: y = {m:.2f}x + {b:.2f}  (R²={r2:.3f})")

        # GRAY y=x reference
        ax.plot([lo, hi], [lo, hi], color="0.6", lw=1.25, ls="--", label="Ideal: y = x")

        # cosmetics
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("True pIC50")
        ax.set_ylabel("Predicted pIC50")
        ax.set_title(f"True vs Predicted Values on {dataset_label}: {split_name.upper()} Fold {fold}")
        ax.grid(True, linestyle=":", alpha=0.25)
        ax.legend(loc="best", frameon=True, facecolor="white", framealpha=0.95)

        fig.tight_layout()
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

        # Save predictions to CSV next to the plot
        pd.DataFrame(
            {"Smiles": df_true[SMILES_COL], "y_true": y_true, "y_pred": y_hat}
        ).to_csv(Path(out_png).with_suffix(".csv"), index=False)

    save_pred_plot(test_loader, df_test, fold_out / "pred_vs_true_fold_test.png")
    save_pred_plot(ext_loader,  df_ext_master, fold_out / "pred_vs_true_sexual_data.png")

    # Keep metrics for summary table
    def pull_rmse(d):
        for k, v in d.items():
            if "rmse" in k.lower():
                return float(v)
        return np.nan

    rows.append({
        "fold": i,
        "fold_test_rmse": pull_rmse(test_metrics),
        "external_rmse":  pull_rmse(ext_metrics),
        **{f"fold_{k}": v for k, v in test_metrics.items()},
        **{f"external_{k}": v for k, v in ext_metrics.items()},
    })

# Write summary
summary = pd.DataFrame(rows).set_index("fold")
summary.to_csv(OUTPUTDIR / "cv_results_with_external.csv")
print("\nPer-fold summary (RMSE columns):\n", summary[["fold_test_rmse","external_rmse"]])
print("\nMeans:\n", summary[["fold_test_rmse","external_rmse"]].mean())
print(f"\nLogs & checkpoints at: {OUTPUTDIR}")
print(f"(If you installed TensorBoard: tensorboard --logdir {OUTPUTDIR})")