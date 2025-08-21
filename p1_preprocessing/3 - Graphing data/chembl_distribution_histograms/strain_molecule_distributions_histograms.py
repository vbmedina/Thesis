import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Top 50 Strains by Molecule Pairings (Non-Unique) ---------------------------------------------------------------
# Load data
strain_data = pd.read_csv("./str_counts.csv")

# Sort top 50 strains
top_strains = strain_data.nlargest(50, 'Count').sort_values(by='Count', ascending=False)

# Create color palette: Reds inverted
strain_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Make barplot
plt.figure(figsize=(14, 8))
sns.barplot(data=top_strains, x='Strain', y='Count', palette=strain_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.xlabel('Strain Identifier', fontsize=12)
plt.ylabel('Number of Times Tested', fontsize=12)
plt.title('Top 50 Most Frequently Tested Strains', fontsize=14)
plt.tight_layout()

# Save and display
plt.savefig("./top_50_strains_hm.png", dpi=300)
plt.show()

# Top 50 Molecules by Strain Pairings (Non-Unique) ---------------------------------------------------------------
# Load data
molecule_data = pd.read_csv("./mol_counts.csv")

# Sort top 50 molecules
top_molecules = molecule_data.nlargest(50, 'Count').sort_values(by='Count', ascending=False)

# Reversed red palette
molecule_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Generate molecule barplot
plt.figure(figsize=(14, 8))
sns.barplot(data=top_molecules, x='Molecule', y='Count', palette=molecule_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=7)
plt.xlabel('ChEMBL Molecule Identifier', fontsize=12)
plt.ylabel('Number of Times Tested', fontsize=12)
plt.title('Top 50 Most Frequently Tested Molecules', fontsize=14)
plt.tight_layout()

# Save and display
plt.savefig("./top_50_molecules_hm.png", dpi=300)
plt.show()

# Top 50 Molecules by Unique Strain Coverage ---------------------------------------------------------------
# Load data
pairings_data = pd.read_csv("./pairings_center.csv")

# Calculate unique strains tested per molecule
unique_strain_counts = (pairings_data.groupby('Molecule')['Strain'].nunique()
    .nlargest(50)
    .reset_index()
    .sort_values(by='Strain', ascending=False))

unique_strain_counts.columns = ['Chemical', 'UniqueStrains']

# Color palette
diversity_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Generate barplot 
plt.figure(figsize=(14, 8))
sns.barplot(data=unique_strain_counts, x='Chemical', y='UniqueStrains', palette=diversity_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=7)
plt.xlabel('ChEMBL Molecule Identifier', fontsize=12)
plt.ylabel('Number of Unique Strains Tested', fontsize=12)
plt.title('Top 50 Molecules by Most Unique Strain Coverage', fontsize=14)
plt.tight_layout()

# Save and display
plt.savefig("./molecules_unique_strain_coverage.png", dpi=300)
plt.show()