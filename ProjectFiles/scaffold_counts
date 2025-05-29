import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV
df = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl_scaf.csv", low_memory=False)

# Generate scaffold column if not present
def smiles_to_scaffold(smiles):
    if not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) if mol else None

df['Scaffold'] = df['Smiles'].apply(smiles_to_scaffold)

# Count molecules per scaffold
scaffold_counts = df['Scaffold'].value_counts()
counts = scaffold_counts.values

# Save scaffold counts to a new CSV
scaffold_counts.to_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/scaffold_counts_after"
".csv", header=True)