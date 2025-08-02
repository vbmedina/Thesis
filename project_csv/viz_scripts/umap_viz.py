import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.preprocessing import StandardScaler
import umap

# Load data
data = ".pp.csv"
df = pd.read_csv(data)
print(f"Loaded {len(df):,} rows from {data}")

# Filter to valid SMILES and pChEMBL values
df = df.dropna(subset=["Smiles", "pChEMBL Value"]).reset_index(drop=True)
print(f"{len(df):,} rows with valid SMILES and pChEMBL")

# 3Sample a subset for viz
df = df.sample(n=40324)

# Convert SMILES to Morgan fingerprints
def smiles_to_fp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

# Extract scaffold from SMILES
fps = []
valid_indices = []

# Iterate over SMILES to generate fingerprints
for i, smi in enumerate(df["Smiles"]):
    fp = smiles_to_fp(smi)
    if fp:
        fps.append(np.array(fp))
        valid_indices.append(i)

# Convert to numpy array
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

# Run UMAP
print("Running UMAP...")
scaled_fps = StandardScaler().fit_transform(fps)
reducer = umap.UMAP()
embedding = reducer.fit_transform(scaled_fps)
print("UMAP completed.")

# Plot the results
plt.figure(figsize=(10, 6))
colors = {"High": "#FDE725FF", "Moderate": "#7AD151FF", "Low": "#22A884FF"}

# Scatter plot with colors based on potency
for category, color in colors.items():
    mask = df["Potency_Bucket"] == category
    plt.scatter(embedding[mask, 0], embedding[mask, 1],
                c=color, label=category, alpha=0.6, s=10)

# Add a legend and labels
plt.title("UMAP of Morgan Fingerprints Colored by Potency")
plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.legend(title="Potency")
plt.grid(True)
plt.tight_layout()

# 8. Save
output_path = Path("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/umap_chemspace_viz.png")
plt.savefig(output_path, dpi=300)
print(f"Plot saved to {output_path}")
plt.show()
