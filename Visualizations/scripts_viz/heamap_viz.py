import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the newly uploaded file
df_heatmap = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/pairings_strains_IC50_scaff.csv")

# Drop rows with missing Representative_pIC50
df_heatmap_clean = df_heatmap.dropna(subset=["Representative_pIC50"])

# Sort by Representative_pIC50 (descending) and select top 100 unique strain-chemical pairs
top_100 = df_heatmap_clean.sort_values(by="Representative_pIC50", ascending=False).head(100)
print(top_100["Representative_pIC50"])
# Pivot the table for heatmap: rows as chemicals, columns as strains
heatmap_data = top_100.pivot_table(index="Chemical", columns="Strain", values="Representative_pIC50")

# Plot the heatmap 
plt.figure(figsize=(11, 7))
sns.heatmap(
    heatmap_data,
    cmap="coolwarm",
    annot=True,
    fmt=".2f",
    cbar_kws={'label': 'Representative pIC50'},
    annot_kws={"size": 6}
)

plt.title("Top 100 Strain-Molecule Pairs by Representative pIC50")
plt.xlabel("Strain", labelpad=1)
plt.ylabel("Molecule", labelpad=10)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()

# Save the heatmap
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/top_100_pairs_hm_scaff.png", dpi=300)