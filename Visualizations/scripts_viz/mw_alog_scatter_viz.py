# Scatterplot of Molecular Weight (MW) vs AlogP ---------------------------------
# Load the data
data_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv"
data = pd.read_csv(data_path)

# Count number of assays per molecule
assay_counts = data["Molecule_ChEMBL_ID"].value_counts().rename("Number of Assays")
data = data.merge(assay_counts, left_on="Molecule_ChEMBL_ID", right_index=True)

# Plot setup
plt.figure(figsize=(12, 8))

# Add shaded "drug-like" zone: AlogP –1 to 5, MW ≤ 500
plt.axvspan(1, 4, ymin=0, ymax=500/1000, color='green', alpha=0.3, label="Drug-like zone")

# Scatterplot
sns.scatterplot(
    data=data,
    x="AlogP",
    y="Molecular Weight",
    size="Number of Assays",
    hue="pChEMBL Value",
    palette="viridis",
    sizes=(10, 300),
    alpha=0.7
)

# Customize
plt.xlabel('AlogP (lipophilicity)', fontsize=12)
plt.ylabel('Molecular Weight (Da)', fontsize=12)
plt.title('MW vs AlogP, Colored by pIC50 and Sized by # Assays', fontsize=14)
plt.legend(title='pIC50 and Assays', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()

# Save + show
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/mw_vs_alogp_scatterplot.png", dpi=300)
plt.show()