import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import SaltRemover
from pathlib import Path

# Asexual Data
# in_csv = Path("./do_not_touch/postphase3_Asexual_Only.csv")
# out_fir = Path("./do_not_touch/postphase4_validSmiles.csv")
# out_inv = Path("./do_not_touch/postphase4_invalidSmiles.csv")

# Sexual Data
in_csv = Path("./p0_do_not_touch/postphase3_Sexual_Only.csv")
out_fir = Path("./p2_models/sexual_test.csv")
out_inv = Path("./p0_do_not_touch/invalidSmiles_sexual_data.csv")
               
df = pd.read_csv(in_csv)
before_removal = len(df)

df = df[df['Smiles'].apply(lambda x: isinstance(x, str) and len(x) > 0)]

# Initialize the salt remover with default salts
remover = SaltRemover.SaltRemover()
# Remove salts from SMILES
df['Smiles'] = df['Smiles'].apply(lambda x: Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(x))))

# Delete invalid SMILES
def is_valid_smiles(smi):
    """Returns True if SMILES is valid, False otherwise."""
    if not isinstance(smi, str) or len(smi.strip()) == 0:
        return False
    mol = Chem.MolFromSmiles(smi)
    return mol is not None

# Apply the validation function
valid_mask = df['Smiles'].apply(is_valid_smiles)

# Split into valid and invalid DataFrames
valid_df = df[valid_mask].copy()
invalid_df = df[~valid_mask].copy()

# Print summary
print(f"Valid SMILES: {len(valid_df)}")
print(f"Invalid SMILES removed: {len(invalid_df)}")
print(f"Empty rows removed: {before_removal - len(valid_df)}")

# Save the valid and invalid DataFrames
valid_df.to_csv(out_fir, index=False)
if len(invalid_df) > 0:
    invalid_df.to_csv(out_inv, index=False)