import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.manifold import TSNE
from pathlib import Path

# Load data
path= "./pp.csv"
df = pd.read_csv(path)
print(f"Loaded {len(df):,} rows from {path}")

# Filter to valid SMILES and pChEMBL values
df = df.dropna(subset=["Smiles", "pChEMBL Value"]).reset_index(drop=True)
print(f"{len(df):,} rows with valid SMILES and pChEMBL")

# Sample subset for viz
df = df.sample(n=40324)

# Convert SMILES to Morgan fingerprints
def smiles_to_fp(Smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(Smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

# Extract scaffold from SMILES
def extract_scaffold(Smiles):
    mol = Chem.MolFromSmiles(Smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)

# Generate fingerprints and scaffolds
fps = []
valid_indices = []
scaffolds = []

# Iterate over SMILES to generate fingerprints and scaffolds
for i, smi in enumerate(df["Smiles"]):
    fp = smiles_to_fp(smi)
    scaffold = extract_scaffold(smi)
    if fp:
        fps.append(np.array(fp))
        valid_indices.append(i)
        scaffolds.append(scaffold)

# Convert to numpy array
fps = np.array(fps)
df = df.iloc[valid_indices].copy()
df["Scaffold"] = scaffolds
print(f"Generated fingerprints and scaffolds for {len(df):,} molecules")

# # Bucket pChEMBL into potency categories
# def bucket_pchembl(p):
#     if p >= 8:
#         return "High"
#     elif p >= 6:
#         return "Moderate"
#     else:
#         return "Low"

# df["Potency_Bucket"] = df["pChEMBL Value"].apply(bucket_pchembl)

# Run t-SNE
print("Running t-SNE...")
tsne = TSNE(n_components=2, perplexity=30, n_iter=1000)
embedding = tsne.fit_transform(fps)
print("t-SNE completed.")

print(embedding.shape)

# Plot the results
plt.figure(figsize=(10, 6))
colors = {"High": "#C43032", "Moderate": "#e8d5cb", "Low": "#455DCE"}

for category, color in colors.items():
    mask = df["Potency_Bucket"] == category
    plt.scatter(embedding[mask, 0], embedding[mask, 1],
                c=color, label=category, alpha=0.6, s=10)

# Customize plot
plt.title("t-SNE of Morgan Fingerprints Colored by Potency")
plt.xlabel("t-SNE-1")
plt.ylabel("t-SNE-2")
plt.legend(title="Potency")
plt.grid(True)
plt.tight_layout()

# Save
output_path = Path("./tsne_chemspace_scaff_viz.png")
plt.savefig(output_path, dpi=300)
print(f"Plot saved to {output_path}")
