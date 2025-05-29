import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Scatterplot of Molecular Weight (MW) vs AlogP ---------------------------------
# Load the data
data_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl_scaf.csv"
data = pd.read_csv(data_path)

# Count number of assays per molecule
molecule_counts = data["Molecule_ChEMBL_ID"].value_counts().rename("Number of Molecules")
data = data.merge(molecule_counts, left_on="Molecule_ChEMBL_ID", right_index=True)

#Count how many fall in drug-like zone
druglike_mask = (data["Molecular Weight"] <= 500) & (data["AlogP"] >= 1) & (data["AlogP"] <= 4)
num_in_zone = druglike_mask.sum()
total = len(data)
percent_in_zone = 100 * num_in_zone / total

print(f"{num_in_zone:,} out of {total:,} molecules fall within the drug-like zone")
print(f"That’s {percent_in_zone:.2f}% of the dataset")

# Plot setup
plt.figure(figsize=(12, 8))

# Add shaded "drug-like" zone: AlogP 1 to 4, MW ≤ 500
plt.axvspan(1, 4, ymin=0, ymax=400/1000, color='green', alpha=0.3, label="Drug-like zone")

# Scatterplot
sns.scatterplot(
    data=data,
    x="AlogP",
    y="Molecular Weight",
    size="Number of Molecules",
    hue="pChEMBL Value",
    palette="viridis",
    sizes=(10, 300),
    alpha=0.7
)

# Customize
plt.xlabel('AlogP (lipophilicity)', fontsize=12)
plt.ylabel('Molecular Weight (Da)', fontsize=12)
plt.title('MW vs AlogP, Colored by pIC50 and Sized by # of Molecules', fontsize=14)
plt.legend(title='pIC50 and Assays', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()

# Save + show
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/mw_vs_alogp_scat_scaff.png", dpi=300)
plt.show()