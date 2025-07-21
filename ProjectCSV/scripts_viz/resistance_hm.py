import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load and clean the data
df = pd.read_csv('./pp.csv')
df_clean = df.dropna(subset=['Representative_pIC50'])

# # First
# #----------------------------------
# # Filter for "CQR" in Flag Groups
# df_cqr = df_clean[df_clean['Flag Groups'].str.contains('CQR', na=False)]

# # Select top 100 by Representative_pIC50
# top_molecules = (
#     df_cqr
#     .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
#     .max()
#     .sort_values(ascending=False)
#     .head(100)
#     .index
# )

# # Subset the original CQR data to only these top molecules
# df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# # Pivot to matrix for heatmap
# heatmap_data = df_top100_mols.pivot_table(
#     index='Molecule_ChEMBL_ID',
#     columns='Strains',
#     values='Representative_pIC50'
# )

# # Plot heatmap
# plt.figure(figsize=(14, 20))


# # Define vmin/vmax for color scaling
# vmin = heatmap_data.min().min()
# vmax = heatmap_data.max().max()

# ax = sns.heatmap(
#     heatmap_data,
#     cmap='Reds',
#     vmin=vmin,
#     vmax=vmax,
#     annot=True,
#     fmt='.2f',
#     cbar_kws={'label': 'Representative pIC50'},
#     annot_kws={'size': 6},
#     linewidths=0.25,
#     linecolor='gray'
# )

# # Rotate axis labels
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=10)
# ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4)

# # 8. Adjust layout and title
# plt.title("Drug Potency of Top 100 Molecules Tested Against Chloroquine Resistant Strains", pad=20)
# plt.xlabel('Strain')
# plt.ylabel('Molecule ChEMBL ID')
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.25)
# plt.show()


# #  Second
# #----------------------------------
# # Filter for "MEFR" in Flag Groups
# df_cqr = df_clean[df_clean['Flag Groups'].str.contains('MEFR', na=False)]

# # Select top 100 by Representative_pIC50
# top_molecules = (
#     df_cqr
#     .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
#     .max()
#     .sort_values(ascending=False)
#     .head(100)
#     .index
# )

# # Subset the original CQR data to only these top molecules
# df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# # Pivot to matrix for heatmap
# heatmap_data = df_top100_mols.pivot_table(
#     index='Molecule_ChEMBL_ID',
#     columns='Strains',
#     values='Representative_pIC50'
# )

# # Plot heatmap
# plt.figure(figsize=(14, 20))


# # Define vmin/vmax for color scaling
# vmin = heatmap_data.min().min()
# vmax = heatmap_data.max().max()

# ax = sns.heatmap(
#     heatmap_data,
#     cmap='Reds',
#     vmin=vmin,
#     vmax=vmax,
#     annot=True,
#     fmt='.2f',
#     cbar_kws={'label': 'Representative pIC50'},
#     annot_kws={'size': 6},
#     linewidths=0.25,
#     linecolor='gray'
# )

# # Rotate axis labels
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=10)
# ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4)

# # 8. Adjust layout and title
# plt.title("Drug Potency of Top 100 Molecules Tested Against Mefloquine Resistant Strains", pad=20)
# plt.xlabel('Strain')
# plt.ylabel('Molecule ChEMBL ID')
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.25) 
# plt.show()


# # Third
# #----------------------------------
# # Filter for "PYRR" in Flag Groups
# df_cqr = df_clean[df_clean['Flag Groups'].str.contains('PYRR', na=False)]

# # Select top 100 by Representative_pIC50
# top_molecules = (
#     df_cqr
#     .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
#     .max()
#     .sort_values(ascending=False)
#     .head(100)
#     .index
# )

# # Subset the original CQR data to only these top molecules
# df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# # Pivot to matrix for heatmap
# heatmap_data = df_top100_mols.pivot_table(
#     index='Molecule_ChEMBL_ID',
#     columns='Strains',
#     values='Representative_pIC50'
# )

# # Plot heatmap
# plt.figure(figsize=(14, 20))


# # Define vmin/vmax for color scaling
# vmin = heatmap_data.min().min()
# vmax = heatmap_data.max().max()

# ax = sns.heatmap(
#     heatmap_data,
#     cmap='Reds',
#     vmin=vmin,
#     vmax=vmax,
#     annot=True,
#     fmt='.2f',
#     cbar_kws={'label': 'Representative pIC50'},
#     annot_kws={'size': 6},
#     linewidths=0.25,
#     linecolor='gray'
# )

# # Rotate axis labels
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=10)
# ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4)

# # 8. Adjust layout and title
# plt.title("Drug Potency of Top 100 Molecules Tested Against Pyrimethamine Resistant Strains", pad=20)
# plt.xlabel('Strain')
# plt.ylabel('Molecule ChEMBL ID')
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.25) 
# plt.show()


# # Fourth
# #----------------------------------
# # Filter for "CQS" in Flag Groups
# df_cqr = df_clean[df_clean['Flag Groups'].str.contains('CQS', na=False)]

# # Select top 100 by Representative_pIC50
# top_molecules = (
#     df_cqr
#     .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
#     .max()
#     .sort_values(ascending=False)
#     .head(100)
#     .index
# )

# # Subset the original CQR data to only these top molecules
# df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# # Pivot to matrix for heatmap
# heatmap_data = df_top100_mols.pivot_table(
#     index='Molecule_ChEMBL_ID',
#     columns='Strains',
#     values='Representative_pIC50'
# )

# # Plot heatmap
# plt.figure(figsize=(14, 20))


# # Define vmin/vmax for color scaling
# vmin = heatmap_data.min().min()
# vmax = heatmap_data.max().max()

# ax = sns.heatmap(
#     heatmap_data,
#     cmap='Reds',
#     vmin=vmin,
#     vmax=vmax,
#     annot=True,
#     fmt='.2f',
#     cbar_kws={'label': 'Representative pIC50'},
#     annot_kws={'size': 6},
#     linewidths=0.25,
#     linecolor='gray'
# )

# # Rotate axis labels
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=10)
# ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4)

# # 8. Adjust layout and title
# plt.title("Drug Potency of Top 100 Molecules Tested Against Chloroquine Sensitive Strains", pad=20)
# plt.xlabel('Strain')
# plt.ylabel('Molecule ChEMBL ID')
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.25) 
# plt.show()


# Fifth
#----------------------------------
# Filter for "MEFS" in Flag Groups
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('MEFS', na=False)]

# Select top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset the original CQR data to only these top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot heatmap
plt.figure(figsize=(14, 20))


# Define vmin/vmax for color scaling
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

ax = sns.heatmap(
    heatmap_data,
    cmap='Reds',
    vmin=vmin,
    vmax=vmax,
    annot=True,
    fmt='.2f',
    cbar_kws={'label': 'Representative pIC50'},
    annot_kws={'size': 6},
    linewidths=0.25,
    linecolor='gray'
)

# Rotate axis labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=10)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4)

# 8. Adjust layout and title
plt.title("Drug Potency of Top 100 Molecules Tested Against Mefloquine Sensitive Strains", pad=20)
plt.xlabel('Strain')
plt.ylabel('Molecule ChEMBL ID')
plt.tight_layout()
plt.subplots_adjust(bottom=0.25) 
plt.show()

