''' Description: This script uses the CSV file generated from the Step 1. 1_mapping_drug_resistance.py
to create heatmaps for different drug resistance flags. It visualizes the drug potency of the top 100 molecules
tested against various strains, categorized by their resistance flags.'''

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load and clean the data
df = pd.read_csv('./p1_preprocessing/3 - Graphing data/drug_resistance/csv_with_resistance_flags.csv')
df_clean = df.dropna(subset=['Representative_pIC50'])

# 1. CQR HM ---------------------------------------------------------------
# Filter for "CQR"
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('CQR', na=False)]

# Top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset original CQR data to only top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot
plt.figure(figsize=(14, 20))

# Define vmin/vmax for color bar
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

# Create heatmap
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

# Layout and title
plt.title("Drug Potency of Top 100 Molecules Tested Against Chloroquine Resistant Strains", pad=20)
plt.xlabel('Strain')
plt.ylabel('Molecule ChEMBL ID')
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.show()

plt.savefig('./p1_preprocessing/3 - Graphing data/drug_resistance/cqr_hm.png', dpi=300, bbox_inches='tight')

# 2. MRFR HM ---------------------------------------------------------------
# Filter for "MEFR"
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('MEFR', na=False)]

# Top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset original CQR data to only top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]

# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot heatmap
plt.figure(figsize=(14, 20))

# Define vmin/vmax for color bar
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

# Create heatmap
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

# Adjust layout and title
plt.title("Drug Potency of Top 100 Molecules Tested Against Mefloquine Resistant Strains", pad=20)
plt.xlabel('Strain')
plt.ylabel('Molecule ChEMBL ID')
plt.tight_layout()
plt.subplots_adjust(bottom=0.25) 
plt.show()

plt.savefig('./p1_preprocessing/3 - Graphing data/drug_resistance/mefr_hm.png', dpi=300, bbox_inches='tight')

# 3. PYRR HM ---------------------------------------------------------------
# Filter for "PYRR"
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('PYRR', na=False)]

# Top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset original CQR data to only top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]

# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot heatmap
plt.figure(figsize=(14, 20))

# Define vmin/vmax for color bar
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

# Create heatmap
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

# Adjust layout and title
plt.title("Drug Potency of Top 100 Molecules Tested Against Pyrimethamine Resistant Strains", pad=20)
plt.xlabel('Strain')
plt.ylabel('Molecule ChEMBL ID')
plt.tight_layout()
plt.subplots_adjust(bottom=0.25) 
plt.show()

plt.savefig('./p1_preprocessing/3 - Graphing data/drug_resistance/pyrr_hm.png', dpi=300, bbox_inches='tight')

# 4. CQS HM ---------------------------------------------------------------
# Filter for "CQS"
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('CQS', na=False)]

# Top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset the original CQR data to only top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]

# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot heatmap
plt.figure(figsize=(14, 20))

# Define vmin/vmax for color bar
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

# Create heatmap
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

# Adjust layout and title
plt.title("Drug Potency of Top 100 Molecules Tested Against Chloroquine Sensitive Strains", pad=20)
plt.xlabel('Strain')
plt.ylabel('Molecule ChEMBL ID')
plt.tight_layout()
plt.subplots_adjust(bottom=0.25) 
plt.show()

plt.savefig('./p1_preprocessing/3 - Graphing data/drug_resistance/cqs_hm.png', dpi=300, bbox_inches='tight')

# 5. MEFS HM ---------------------------------------------------------------
# Filter for "MEFS"
df_cqr = df_clean[df_clean['Flag Groups'].str.contains('MEFS', na=False)]

# Top 100 by Representative_pIC50
top_molecules = (
    df_cqr
    .groupby('Molecule_ChEMBL_ID')['Representative_pIC50']
    .max()
    .sort_values(ascending=False)
    .head(100)
    .index
)

# Subset the original CQR data to only top molecules
df_top100_mols = df_cqr[df_cqr['Molecule_ChEMBL_ID'].isin(top_molecules)]


# Pivot to matrix for heatmap
heatmap_data = df_top100_mols.pivot_table(
    index='Molecule_ChEMBL_ID',
    columns='Strains',
    values='Representative_pIC50'
)

# Plot heatmap
plt.figure(figsize=(14, 20))


# Define vmin/vmax for color bar
vmin = heatmap_data.min().min()
vmax = heatmap_data.max().max()

# Create heatmap
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
plt.savefig('./p1_preprocessing/3 - Graphing data/drug_resistance/mefs_hm.png', dpi=300, bbox_inches='tight')