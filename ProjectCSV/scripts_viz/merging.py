import pandas as pd

# 1) Load chembl.csv with only the ID and SMILES columns
chembl = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv", low_memory=False, usecols=["Molecule_ChEMBL_ID", "Smiles"])

# 2) Load the pairings file
pair_df = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/pairings_strains_IC50_edits.csv")
