import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV
df = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv", low_memory=False)

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

sns.set_theme(style="whitegrid", palette="viridis")

# Log-log histogram
bins = np.logspace(np.log10(1), np.log10(counts.max()), num=50)
plt.figure()
plt.hist(counts, bins=bins)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Molecules per scaffold (log scale)')
plt.ylabel('Number of scaffolds (log scale)')
plt.title('Distribution of Scaffold Frequency')
plt.tight_layout()

# Save the histogram
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/scaffold_frequency_histogram.png", dpi=300)