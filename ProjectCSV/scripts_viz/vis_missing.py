# def compute_matrix_completeness(df):
#     """
#     Compute matrix completeness percentage. Either provide a DataFrame, or provide the components directly:
#     - unique_cells: list/set or int
#     - unique_compounds: list/set or int
#     - data_points: number of measurements
#     :return: Completeness percentage (0-100)
#     """
#     #number of na cells
#     total_missing = df.isna().values.sum()
#     return 100 - (total_missing / (df.shape[0] * df.shape[1]) * 100)

# import pandas as pd
# import missingno as msno
# import matplotlib.pyplot as plt

# # Load your dataset
# df = pd.read_csv("./pairings_table.csv")

# # Create two subsets: first 100 and last 100 rows
# df_first_100 = df.head(100)
# df_last_100 = df.tail(100)

# fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# # Plot first 100 rows
# msno.matrix(df_first_100, ax=axes[0])
# axes[0].set_title("First 100 Rows")

# # Plot last 100 rows
# msno.matrix(df_last_100, ax=axes[1])
# axes[1].set_title("Last 100 Rows")

# plt.tight_layout()
# plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def compute_matrix_completeness(df):
    total_values = df.shape[0] * df.shape[1]
    total_missing = df.isna().sum().sum()
    completeness = 100 - (total_missing / total_values) * 100
    return round(completeness, 2)

# Load your dataset
df = pd.read_csv("./pairings_table.csv", index_col=0)

# Drop "Grand Total" row and column if they exist
df = df.drop(index="Grand Total", errors="ignore")
df = df.drop(columns=["Grand Total"], errors="ignore")

# Rank rows and columns by how many values they have
row_completeness = df.notna().sum(axis=1).sort_values(ascending=False)
col_completeness = df.notna().sum(axis=0).sort_values(ascending=False)

# Choose top N most complete
top_rows = row_completeness.head(50).index
top_cols = col_completeness.head(50).index
filtered_df = df.loc[top_rows, top_cols]

# Binary mask (1 = present, 0 = missing)
binary_matrix = filtered_df.notna().astype(int)

# Sort again by completeness
binary_matrix = binary_matrix.loc[
    binary_matrix.sum(axis=1).sort_values(ascending=False).index,
    binary_matrix.sum(axis=0).sort_values(ascending=False).index
]

# Calculate % completeness
completeness = compute_matrix_completeness(filtered_df)

# Plot
plt.figure(figsize=(14, 10))
sns.heatmap(
    binary_matrix,
    cmap="Greys",
    cbar=False,
    linewidths=0.1,
    linecolor="lightgray"
)

# Rotate axis labels and set font size
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)


plt.title(f"Top 50 Molecules x Strains — Data Completeness: {completeness}%", fontsize=14)
plt.xlabel("Strains")
plt.ylabel("Molecules")
plt.tight_layout()
plt.show()

# #--------------------------------------
# import pandas as pd
# import matplotlib.pyplot as plt
# import missingno as msno

# def compute_matrix_completeness(df):
#     total = df.shape[0] * df.shape[1]
#     missing = df.isna().sum().sum()
#     return round(100 - (missing / total) * 100, 2)

# # Load and clean
# df = pd.read_csv("./Thesis/ProjectCSV/pairings_table.csv", index_col=0)
# df = df.drop(index="Grand Total", errors="ignore")
# df = df.drop(columns=["Grand Total"], errors="ignore")

# # Transpose: strains become rows
# df_transposed = df.T.dropna(how="all")

# # Sort strains by completeness
# sorted_strains = df_transposed.notna().sum(axis=1).sort_values(ascending=False)

# # Most complete strains
# top_strains = sorted_strains.head(50).index
# bottom_strains = sorted_strains.tail(50).index

# df_top = df_transposed.loc[top_strains]
# df_bottom = df_transposed.loc[bottom_strains]

# # Compute completeness
# comp_top = compute_matrix_completeness(df_top)
# comp_bottom = compute_matrix_completeness(df_bottom)

# # Plot
# fig, axes = plt.subplots(1, 2, figsize=(24, 12))

# msno.matrix(df_top, ax=axes[0], sparkline=False)
# axes[0].set_title(f"Top 50 Most Complete Strains — Completeness: {comp_top}%")
# axes[0].set_xlabel("Molecules")
# axes[0].set_ylabel("Strains")

# msno.matrix(df_bottom, ax=axes[1], sparkline=False)
# axes[1].set_title(f"Bottom 50 Least Complete Strains — Completeness: {comp_bottom}%")
# axes[1].set_xlabel("Molecules")
# axes[1].set_ylabel("")

# plt.tight_layout()
# plt.show()

