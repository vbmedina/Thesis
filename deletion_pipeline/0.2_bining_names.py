import re
import itertools
from pathlib import Path
from difflib import SequenceMatcher
import pandas as pd
from collections import Counter

# File path 
path = "/Users/victoriamedina/Thesis_Project/Thesis/TESTING.csv"

data = pd.read_csv(path)

print(data.columns)

# Verify the column "stand_strain" exists
if "Standardized Strain" in data.columns:
    # Count occurrences of each item in the "stand_strain column
    strain_counts = Counter(data["Standardized Strain"].dropna())
    
    # Print the counts
    for strain, count in strain_counts.items():
        print(f"{strain}: {count}")
else:
    print("The column 'Standardized Strain' does not exist in the CSV file.")
    strain_counts = Counter()  

# Create DF from the strain counts
strain_counts_df = pd.DataFrame(strain_counts.items(), columns=["Strain", "Count"])


# Save DF to new CSV
strain_counts_df.to_csv(path, index=False)

print(f"Strain counts have been saved to {path}")
 