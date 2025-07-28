import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib.pyplot as plt
import seaborn as sns

# Load CSV
df = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl_scaf.csv", low_memory=False)

# Generate scaffolds
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

# Log-log histogram with inverted axes
bins = np.logspace(np.log10(1), np.log10(counts.max()), num=50)
plt.figure()
plt.hist(counts, bins=bins, orientation='horizontal')  # Invert axes by setting orientation to horizontal
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Molecules per scaffold (log scale)')
plt.xlabel('Frequency of scaffolds (log scale)')
plt.title('Distribution of Scaffold Frequency')
plt.tight_layout()

plt.savefig("./scaffold_freq_hist_scaff_play.png", dpi=300)