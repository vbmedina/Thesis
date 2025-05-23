import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

#Histogram of the top 20 strains ---------------------------------
# Load the data
data_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/str_counts_final.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'Count')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='Strain', y='Count', palette='viridis')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Target Strain', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title('Frequency of Top 20 Strains', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/top_20_strains_histogram.png", dpi=300)
plt.show()

# Histogram of the top 20 chemicals ---------------------------------
# Load the data
data_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chem_counts_final.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'Count')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='Chemical', y='Count', palette='viridis')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Molecule ID', fontsize=12)
plt.ylabel('Number of Molecules', fontsize=12)
plt.title('Frequency of Top 20 Molecules', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/top_20_chemicals_histogram.png", dpi=300)
plt.show()

# Histogram of the top 20 chemicals ---------------------------------# 
# Load the data
data_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/pairings_strains_IC50.csv"
data = pd.read_csv(data_path)

# Sort and filter the top 20 strains
top_strains = data.nlargest(20, 'Unique_Strains_Tested_On')

# Create the histogram
plt.figure(figsize=(12, 8))
sns.barplot(data=top_strains, x='Chemical_ID', y='Unique_Strains_Tested_On', palette='viridis')

# Customize the plot
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel('Molecule ID', fontsize=12)
plt.ylabel('Number of Unique Strains Tested On', fontsize=12)
plt.title('Molecules and Unique Chemicals Tested', fontsize=14)
plt.tight_layout()

# Save and show the plot
plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/molecules_unique_strains_histogram.png", dpi=300)
plt.show()