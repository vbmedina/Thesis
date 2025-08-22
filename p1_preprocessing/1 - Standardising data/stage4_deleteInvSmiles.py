''' Description: This script is the fourth stage to the pipeline (SMILES). It works on the SMILES.
This script:
1) removes rows with invalid SMILES.

Preconditions:
1) "postphase3_Asexual_Only.csv" - generated from stage3 sexual filter
2) "postphase3_Sexual_Only.csv" - generated from stage3 sexual filter
'''

# Imports
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover
from pathlib import Path

for sexType in ["sexual", "asexual"]:

    if sexType == "asexual":
        # Asexual Data
        in_csv = Path("./p0_all_csvs/postphase3_Asexual_Only.csv")
        out_fir = Path("./p0_all_csvs/postphase4_validSmiles.csv")
        out_inv = Path("./p0_all_csvs/postphase4_invalidSmiles.csv")

    else:
        # Sexual Data
        in_csv = Path("./p0_all_csvs/postphase3_Sexual_Only.csv")
        out_fir = Path("./p0_all_csvs/sexual_test.csv")
        out_inv = Path("./p0_all_csvs/invalidSmiles_sexual_data.csv")

    # Load DataFrame 
    df = pd.read_csv(in_csv)
    before_removal = len(df)

    # Remove 'Smiles' column if it isnt an empty string
    df = df[df['Smiles'].apply(lambda x: isinstance(x, str) and len(x) > 0)]

    # Initialize the salt remover with default salts
    remover = SaltRemover.SaltRemover()

    # Remove salts from SMILES
    df['Smiles'] = df['Smiles'].apply(lambda x: Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(x))))

    # Returns True if the SMILES is valid (is a string and has a mol), False otherwise
    def is_valid_smiles(smi):
        if not isinstance(smi, str) or len(smi.strip()) == 0:
            return False
        mol = Chem.MolFromSmiles(smi)
        return mol is not None

    # Split into valid and invalid smiles
    valid_mask = df['Smiles'].apply(is_valid_smiles)
    valid_df = df[valid_mask].copy()
    invalid_df = df[~valid_mask].copy()

    # Print summary
    print(f"Valid SMILES: {len(valid_df)}")
    print(f"Invalid SMILES removed: {len(invalid_df)}")
    print(f"Empty rows removed: {before_removal - len(valid_df)}")

    # Save the valid and invalid DFs to CSV files
    valid_df.to_csv(out_fir, index=False)
    if len(invalid_df) > 0:
        invalid_df.to_csv(out_inv, index=False)