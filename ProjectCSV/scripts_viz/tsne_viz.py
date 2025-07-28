import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from pathlib import Path

# Load data
DATA_PATH = "./pp.csv"
df = pd.read_csv(DATA_PATH)
print(f"Loaded {len(df):,} rows from {DATA_PATH}")

# Filter to valid SMILES and pChEMBL values
df = df.dropna(subset=["Smiles", "pChEMBL Value"]).reset_index(drop=True)
print(f"{len(df):,} rows with valid SMILES and pChEMBL")


# Sample a subset for viz
df = df.sample(n=39907)

# Convert SMILES to Morgan fingerprints
def smiles_to_fp(Smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(Smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

fps = []
valid_indices = []

for i, smi in enumerate(df["Smiles"]):
    fp = smiles_to_fp(smi)
    if fp:
        fps.append(np.array(fp))
        valid_indices.append(i)

fps = np.array(fps)
df = df.iloc[valid_indices].copy()
print(f"Generated fingerprints for {len(df):,} molecules")

# Bucket pChEMBL into potency categories
def bucket_pchembl(p):
    if p >= 8:
        return "High"
    elif p >= 6:
        return "Moderate"
    else:
        return "Low"

df["Potency_Bucket"] = df["pChEMBL Value"].apply(bucket_pchembl)

# Run t-SNE
print(" Running t-SNE...")
tsne = TSNE(n_components=2, perplexity=30, n_iter=1000)
embedding = tsne.fit_transform(fps)
print("t-SNE completed.")

# Plot the results
plt.figure(figsize=(10, 6))
colors = {"High": "#C43032", "Moderate": "#e8d5cb", "Low": "#455DCE"}

for category, color in colors.items():
    mask = df["Potency_Bucket"] == category
    plt.scatter(embedding[mask, 0], embedding[mask, 1],
                c=color, label=category, alpha=0.6, s=10)

# Add a legend and labels
plt.title("t-SNE of Morgan Fingerprints Colored by Potency")
plt.xlabel("t-SNE-1")
plt.ylabel("t-SNE-2")
plt.legend(title="Potency")
plt.grid(True)
plt.tight_layout()

# Save and show
output_path = Path("./tsne_chemspace_viz2.png")
plt.savefig(output_path, dpi=300)
print(f"Plot saved to {output_path}")
