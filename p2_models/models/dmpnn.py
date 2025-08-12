from pathlib import Path
from lightning import pytorch as pl
from lightning.pytorch.callbacks import ModelCheckpoint
import pandas as pd
from chemprop import data, featurizers, models, nn

OUTPUTDIR = "/rds/general/user/vbm24/home/Thesis/p2_models/models/checkpoints"

p2models_dir = Path.cwd().parent
train_path = p2models_dir / "data" / "splits" / "umap" / "umap_fold_1_train.csv"
val_path = p2models_dir / "data" / "splits" / "umap" / "umap_fold1_val.csv"
test_path = p2models_dir / "data" / "splits" / "umap" / "umap_fold1_test.csv"
num_workers = 0 # number of workers for dataloader.
smiles_column = 'Smiles' # name of the column containing SMILES strings
target_columns = ['pIC50'] # list of names of the columns containing targets


df_train = pd.read_csv(train_path)
df_val = pd.read_csv(val_path)
df_test = pd.read_csv(test_path)


smis_train = df_train.loc[:, smiles_column].values
ys_train = df_train.loc[:, target_columns].values

smis_val = df_val.loc[:, smiles_column].values
ys_val = df_val.loc[:, target_columns].values

smis_test = df_test.loc[:, smiles_column].values
ys_test = df_test.loc[:, target_columns].values


train_data = [data.MoleculeDatapoint.from_smi(smi, y) for smi, y in zip(smis_train, ys_train)]
val_data = [data.MoleculeDatapoint.from_smi(smi, y) for smi, y in zip(smis_val, ys_val)]
test_data = [data.MoleculeDatapoint.from_smi(smi, y) for smi, y in zip(smis_test, ys_test)]

featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

train_dset = data.MoleculeDataset(train_data, featurizer)
scaler = train_dset.normalize_targets()

val_dset = data.MoleculeDataset(val_data, featurizer)
val_dset.normalize_targets(scaler)

test_dset = data.MoleculeDataset(test_data, featurizer)

train_loader = data.build_dataloader(train_dset, num_workers=num_workers)
val_loader = data.build_dataloader(val_dset, num_workers=num_workers, shuffle=False)
test_loader = data.build_dataloader(test_dset, num_workers=num_workers, shuffle=False)

mp = nn.BondMessagePassing()

# Type of aggregation
agg = nn.MeanAggregation()

output_transform = nn.UnscaleTransform.from_standard_scaler(scaler)

ffn = nn.RegressionFFN(output_transform=output_transform)

batch_norm = True


metric_list = [nn.metrics.RMSE(), nn.metrics.MSE(), nn.metrics.MAE()] # Only the first metric is used for training and early stopping

mpnn = models.MPNN(mp, agg, ffn, batch_norm, metric_list)

checkpointing = ModelCheckpoint(
    OUTPUTDIR,  # Directory where model checkpoints will be saved
    "best-{epoch}-{val_loss:.2f}",  # Filename format for checkpoints, including epoch and validation loss
    "val_loss",  # Metric used to select the best checkpoint (based on validation loss)
    mode="min",  # Save the checkpoint with the lowest validation loss (minimization objective)
    save_last=True,  # Always save the most recent checkpoint, even if it's not the best
)

trainer = pl.Trainer(
    logger=True,
    enable_checkpointing=True, # Use `True` if you want to save model checkpoints. The checkpoints will be saved in the `checkpoints` folder.
    enable_progress_bar=True,
    accelerator="auto",
    devices=1,
    max_epochs=20, # number of epochs to train for
    callbacks=[checkpointing], # Use the configured checkpoint callback
)

trainer.fit(mpnn, train_loader, val_loader)

results = trainer.test(dataloaders=test_loader)