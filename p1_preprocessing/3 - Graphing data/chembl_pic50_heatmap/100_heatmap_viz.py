''' Description: Using the data from Step 2: pIC50 binning and average, this script uses a copy of the file 
"pairings_center.csv" to generate a heatmap of the top 100 strain-molecule pairs in the whole ChEMBL dataset
using a Representative pIC50. (reference: 2- pIC50 binning and average.py for more information on binning and averaging).'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# File Path
df_heatmap = pd.read_csv("./p1_preprocessing/3 - Graphing data/chembl_pic50_heatmap/pairings_center_copy.csv")

# Drop rows with missing Representative_pIC50
df_heatmap_clean = df_heatmap.dropna(subset=["Representative_pIC50"])

# Sort by pIC50 and select top 100 unique pairs
top_100 = df_heatmap_clean.sort_values(by="Representative_pIC50", ascending=False).head(100)
print(top_100["Representative_pIC50"])

# Pivot the table for heatmap
heatmap_data = top_100.pivot_table(index="Molecule", columns="Strain", values="Representative_pIC50")

# Plot heatmap 
plt.figure(figsize=(11, 7))
sns.heatmap(
    heatmap_data,
    cmap="Reds",
    annot=True,
    fmt=".2f",
    cbar_kws={'label': 'Representative pIC50'},
    annot_kws={"size": 6}
)

# Customize
plt.title("Top 100 Strain-Molecule Pairs by Representative pIC50")
plt.xlabel("Strain", labelpad=1)
plt.ylabel("Molecule", labelpad=10)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()


plt.savefig("./p1_preprocessing/3 - Graphing data/chembl_pic50_heatmap/top_100_pIC_hm.png", dpi=300)