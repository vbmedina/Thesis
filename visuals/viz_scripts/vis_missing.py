import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load datasets
matrix_df = pd.read_csv("./pairings_table.csv", index_col=0)
center_df = pd.read_csv("./pairings_center.csv", index_col=0)

# Drop "Grand Total"
matrix_df = matrix_df.drop(index="Grand Total", errors="ignore")
matrix_df = matrix_df.drop(columns=["Grand Total"], errors="ignore")

# Calculate completeness and rank molecules (rows) and strains (columns)
row_completeness = matrix_df.notna().sum(axis=1).sort_values(ascending=False)
col_completeness = matrix_df.notna().sum(axis=0).sort_values(ascending=False)

# DF as top 50 most complete rows and columns
top_rows = row_completeness.head(50).index
top_cols = col_completeness.head(50).index

# Filter matrix and sort by completeness
filtered_df = matrix_df.loc[top_rows, top_cols]
binary_matrix = filtered_df.notna().astype(int)

# Sort rows and columns based on completeness
sorted_rows = binary_matrix.sum(axis=1).sort_values(ascending=False).index
sorted_cols = binary_matrix.sum(axis=0).sort_values(ascending=False).index
binary_matrix = binary_matrix.loc[sorted_rows, sorted_cols]

# Generate pIC50 heatmap matrix using layout
pic50_pivot = center_df.pivot_table(index="Molecule", columns=center_df.index, values="Representative_pIC50")
ordered_pic50 = pic50_pivot.loc[sorted_rows, sorted_cols]

# Compute completeness
def compute_matrix_completeness(df):
    total_values = df.size
    missing = df.isna().sum().sum()
    return round(100 * (1 - missing / total_values), 2)

completeness = compute_matrix_completeness(ordered_pic50)

# Plot
plt.figure(figsize=(14, 10))
hp = sns.heatmap(
    ordered_pic50,
    cmap="Reds",
    cbar=True,
    cbar_kws={"label": "pIC50"},
    linewidths=0.1,
    linecolor="lightgray"
)

# Rotate colorbar 
colorbar = hp.collections[0].colorbar
colorbar.ax.set_ylabel("pIC50", rotation=-90, labelpad=15)

# Customizations
plt.xticks(fontsize=6, rotation=90)
plt.yticks(fontsize=6)
plt.title(f"Top 50 Molecules Strain Pairings: pIC50 Heatmap\nHeatmap Completeness - {completeness}%", fontsize=14)
plt.xlabel("Strains")
plt.ylabel("Molecules")


plt.tight_layout()
plt.show()