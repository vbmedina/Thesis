import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

#Histogram of the top 20 strains on molecules ---------------------------------
# Load the data
data_path = "./str_counts_final_scaff.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'Count')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='Strain', y='Count', palette='coolwarm')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Target Strain', fontsize=12)
plt.ylabel('Number of Tested Molecules', fontsize=12)
plt.title('Top 20 Strains Tested on the Most Molecules', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("./top_20_strains_histogram_scaff.png", dpi=300)
plt.show()

# Histogram of the top 20 molecules on strains ---------------------------------
# Load the data
data_path = "./chem_counts_final_scaff.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'Count')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='Chemical', y='Count', palette='coolwarm')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Molecule ID', fontsize=12)
plt.ylabel('Number of Tested Strains', fontsize=12)
plt.title('Top 20 Molecules Tested on the Most Strains', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("./top_20_mols_histogram.png", dpi=300)
plt.show()

# Histogram of the top 20  ---------------------------------# 
# Load the data
data_path = "./pairings_center.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'uniq_str')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='mol_id', y='uniq_str', palette='coolwarm')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Molecule ID', fontsize=12)
plt.ylabel('Number of Unique Strains Tested On', fontsize=12)
plt.title('Top 20 Molecules Tested on Unique Strains', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("./molecules_unique_strains_histogram_scaff.png", dpi=300)
plt.show()