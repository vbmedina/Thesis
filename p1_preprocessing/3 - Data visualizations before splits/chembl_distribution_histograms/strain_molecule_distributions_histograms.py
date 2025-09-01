''' Description: This script generates a set of histograms visualizing the distribution of strain and molecule counts.
This script:
1) Generates histograms for the top 50 strains with most amount of pIC50 scores (non-unique molecules).
2) Generates a histogram for the top 50 molecules with most amount of pIC50 score (non-unique strain coverage)
2) Creates a histogram for the top 50 molecules based on unique strain coverage.

Preconditions:
1) "strain_counts.csv" - generated from strain_counts.py
2) "molecule_counts.csv" - generated from molecule_counts.py
3) "pairings_center_copy.csv" - generated from pairings_center.py made in step 2 of the pipeline
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Top 50 Strains by Molecule Pairings (Non-Unique) ---------------------------------------------------------------
# Load data
strain_data = pd.read_csv("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/strain_counts.csv")

# Sort top 50 strains
top_strains = strain_data.nlargest(50, 'Count').sort_values(by='Count', ascending=False)

# Create color palette: Reds inverted
strain_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Make barplot
plt.figure(figsize=(14, 8))
sns.barplot(data=top_strains, x='Strain', y='Count', palette=strain_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=20)
plt.xlabel('Strain Identifier', fontsize=25)
plt.ylabel('Number of pIC50 Scores', fontsize=25)
plt.title('Top 50 Strains with the Greatest Amount of pIC50 Scores\nin CHEMBL364 Dataset (non-unique molecules)', fontsize=28)
plt.tight_layout()

# Save and display
plt.savefig("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/top_50_strains_histogram.png", dpi=300)
plt.show()

# Top 50 Molecules by Strain Pairings (Non-Unique) ---------------------------------------------------------------
# Load data
molecule_data = pd.read_csv("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/molecule_counts.csv")

# Sort top 50 molecules
top_molecules = molecule_data.nlargest(50, 'Count').sort_values(by='Count', ascending=False)

# Reversed red palette
molecule_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Generate molecule barplot
plt.figure(figsize=(14, 8))
sns.barplot(data=top_molecules, x='Molecule', y='Count', palette=molecule_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=20)
plt.xlabel('ChEMBL Molecule Identifier', fontsize=25)
plt.ylabel('Number of IC50 Scores', fontsize=25)
plt.title('Top 50 Molecules with the Greatest Amount of pIC50 Scores\nin CHEMBL364 Dataset (non-unique strain coverage)', fontsize=28)
plt.tight_layout()

# Save and display
plt.savefig("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/top_50_molecules_histogram.png", dpi=300)
plt.show()

# Top 50 Molecules by Unique Strain Coverage ---------------------------------------------------------------
# Load data
pairings_data = pd.read_csv("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/pairings_center_copy.csv")

# Calculate unique strains tested per molecule
unique_strain_counts = (pairings_data.groupby('Molecule')['Strain'].nunique()
    .nlargest(50)
    .reset_index()
    .sort_values(by='Strain', ascending=False))

unique_strain_counts.columns = ['Molecule', 'UniqueStrains']

# Color palette
diversity_palette = sns.color_palette("Reds", n_colors=50)[::-1]

# Generate barplot 
plt.figure(figsize=(14, 8))
sns.barplot(data=unique_strain_counts, x='Molecule', y='UniqueStrains', palette=diversity_palette)

# Customizations
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=20)
plt.xlabel('ChEMBL Molecule Identifier', fontsize=25)
plt.ylabel('Number of Unique Strains Tested', fontsize=25)
plt.title('Top 50 Molecules by Most Unique Strain Coverage\nin CHEMBL364 Dataset', fontsize=28)
plt.tight_layout()

# Save and display
plt.savefig("./p1_preprocessing/3 - Data visualizations before splits/chembl_distribution_histograms/molecules_unique_strain_coverage.png", dpi=300)
plt.show()