from pathlib import Path
import pandas as pd
import sys

# Paths
ROOT   = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "data" / "splits"

# Help message
print("\nRunning split-integrity assertions…")

# Cold‑regime: no ID overlap warm cold_cmpd
cold_dir = SPLITS / "cold_regime"
if cold_dir.exists():
    warm       = pd.read_csv(cold_dir / "warm.csv")
    cold_cmpd_ = pd.read_csv(cold_dir / "cold_cmpd.csv")
    overlap = set(warm["Molecule_ChEMBL_ID"]) & set(cold_cmpd_["Molecule_ChEMBL_ID"])
    assert not overlap, f"Molecule leak between warm and cold_cmpd: {len(overlap)} IDs"

# LOSO: held‑out strain is truly unseen during training
for loso_dir in SPLITS.glob("loso_*"):
    strain = loso_dir.name.replace("loso_", "")
    trn = pd.read_csv(loso_dir / "train.csv")
    tst = pd.read_csv(loso_dir / "test.csv")
    assert strain not in trn["strain_norm"].unique(), f"Strain {strain} leaked into train"
    assert strain in  tst["strain_norm"].unique(),  f"Strain {strain} missing from test"

# Chemistry folds: SMILES are unique across train, val, test
for split in ["random", "scaffold", "umap", "butina"]:
    base = SPLITS / split
    if not base.exists():
        continue
    for fold_dir in sorted(base.glob("fold*")):
        train_smiles = pd.read_csv(fold_dir / "train.csv")["Smiles"]
        val_smiles   = pd.read_csv(fold_dir / "val.csv")["Smiles"]
        test_smiles  = pd.read_csv(fold_dir / "test.csv")["Smiles"]

        assert train_smiles.isin(test_smiles).sum() == 0, f"SMILES leak train↔test in {fold_dir}"
        assert train_smiles.isin(val_smiles).sum()  == 0, f"SMILES leak train↔val  in {fold_dir}"
        assert val_smiles.isin(test_smiles).sum()   == 0, f"SMILES leak val ↔test in {fold_dir}"

print("All split-integrity checks passed.\n")
