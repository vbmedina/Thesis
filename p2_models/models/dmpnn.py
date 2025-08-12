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

# Configeration
SPLIT = "umap"  # or: "random", "scaffold", "butina", "umap_hdb"
p2models_dir = Path.cwd().parent
SPLIT_DIR = p2models_dir / "data" / "splits" / SPLIT
OUTPUTDIR = p2models_dir / "models" / "checkpoints" / SPLIT
OUTPUTDIR.mkdir(parents=True, exist_ok=True)

# External test file
EXTERNAL_TEST_CSV = p2models_dir / "data" / "splits" / "sexual_test.csv"

SMILES_COL = "Smiles"
TARGET_COLS = ["pIC50"]
NUM_WORKERS = 0
MAX_EPOCHS = 200
PATIENCE = 15
use_gpu = torch.cuda.is_available()

# File paths
def find_file(split_dir: Path, split: str, i: int, kind: str) -> Path:
    for name in (f"{split}_fold_{i}_{kind}.csv", f"{split}_fold{i}_{kind}.csv"):
        p = split_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {kind} file for fold {i} in {split_dir}")

# Load the external test dataframe
df_ext_master = pd.read_csv(EXTERNAL_TEST_CSV)

# For table at the end
rows = []

# Print statements for splits
for i in range(1, 6):
    print(f"\n===== Fold {i} ({SPLIT}) =====")

    train_path = find_file(SPLIT_DIR, SPLIT, i, "train")
    val_path   = find_file(SPLIT_DIR, SPLIT, i, "val")
    test_path  = find_file(SPLIT_DIR, SPLIT, i, "test")

    df_train = pd.read_csv(train_path)
    df_val   = pd.read_csv(val_path)
    df_test  = pd.read_csv(test_path)

    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

    def to_points(df):
        return [data.MoleculeDatapoint.from_smi(s, y)
                for s, y in zip(df[SMILES_COL].values, df[TARGET_COLS].values)]

    # Train, validatio and test sets
    train_set = data.MoleculeDataset(to_points(df_train), featurizer)
    scaler = train_set.normalize_targets()

    val_set  = data.MoleculeDataset(to_points(df_val), featurizer)
    val_set.normalize_targets(scaler)

    test_set = data.MoleculeDataset(to_points(df_test), featurizer)
    test_set.normalize_targets(scaler)

    # Sexual test set
    ext_set = data.MoleculeDataset(to_points(df_ext_master), featurizer)
    ext_set.normalize_targets(scaler)

    train_loader = data.build_dataloader(train_set, num_workers=NUM_WORKERS)
    val_loader   = data.build_dataloader(val_set,   num_workers=NUM_WORKERS, shuffle=False)
    test_loader  = data.build_dataloader(test_set,  num_workers=NUM_WORKERS, shuffle=False)
    ext_loader   = data.build_dataloader(ext_set,   num_workers=NUM_WORKERS, shuffle=False)

    # Model
    mp = nn.BondMessagePassing()
    agg = nn.MeanAggregation()  # try SumAggregation()/NormAggregation() if you want to compare
    out_tf = nn.UnscaleTransform.from_standard_scaler(scaler)
    ffn = nn.RegressionFFN(output_transform=out_tf)
    metric_list = [nn.metrics.RMSE(), nn.metrics.MSE(), nn.metrics.MAE()]  # monitors val_rmse
    model = models.MPNN(mp, agg, ffn, batch_norm=True, metric_list=metric_list)

    # Logging and callbacks
    fold_out = OUTPUTDIR / f"fold_{i}"
    fold_out.mkdir(parents=True, exist_ok=True)
    logger = TensorBoardLogger(save_dir=OUTPUTDIR, name=f"fold_{i}")

    ckpt = ModelCheckpoint(
        dirpath=fold_out,
        filename="best-{epoch}-{val_rmse:.4f}",
        monitor="val_rmse", mode="min", save_top_k=1, save_last=True,
    )
    es = EarlyStopping(monitor="val_rmse", mode="min", patience=PATIENCE)

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

    # Train
    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
    print("Best checkpoint:", ckpt.best_model_path)

    # Evaluate on fold test
    test_metrics = trainer.test(dataloaders=test_loader, ckpt_path="best")[0]
    print("Fold test metrics:", test_metrics)

    # Evaluate on external test
    ext_metrics = trainer.test(dataloaders=ext_loader, ckpt_path="best")[0]
    print("External test metrics:", ext_metrics)

    # Visuals: predicted vs true for both test sets
    def save_pred_plot(loader, df_true, out_png):
        with torch.no_grad():
            preds_batches = trainer.predict(model=model, dataloaders=loader, ckpt_path="best")
        y_hat = np.concatenate([p.detach().cpu().numpy().ravel() for p in preds_batches])
        y_true = df_true["pIC50"].to_numpy().ravel()
        plt.figure(figsize=(4,4))
        plt.scatter(y_true, y_hat, s=10, alpha=0.5)
        lo = float(min(y_true.min(), y_hat.min())); hi = float(max(y_true.max(), y_hat.max()))
        plt.plot([lo, hi], [lo, hi], lw=2)
        plt.xlabel("True pIC50"); plt.ylabel("Predicted pIC50"); plt.tight_layout()
        plt.savefig(out_png, dpi=200); plt.close()
        # also save predictions
        pd.DataFrame({"Smiles": df_true[SMILES_COL], "y_true": y_true, "y_pred": y_hat}).to_csv(
            Path(out_png).with_suffix(".csv"), index=False
        )

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
print("\nPer-fold summary:\n", summary[["fold_test_rmse","external_rmse"]])
print("\nMeans:\n", summary[["fold_test_rmse","external_rmse"]].mean())
