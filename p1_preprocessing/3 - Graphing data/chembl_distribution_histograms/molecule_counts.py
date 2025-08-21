''' Description: This script generates a CSV file containing the counts of unique molecules and their occurrences
in a dataset, specifically focusing on the "Molecule_ChEMBL_ID" column.
Preconditions:
1) "final_data_copy.csv"
'''

import pandas as pd
from collections import Counter

# Path
path = "./p1_preprocessing/3 - Graphing data/chembl_distribution_histograms/final_data_copy.csv"

data = pd.read_csv(path)

print(data.columns)

# Verify the column "stand_strain" exists
if "Molecule_ChEMBL_ID" in data.columns:
    # Count item in the "stand_strain" column
    molecule_counts = Counter(data["Molecule_ChEMBL_ID"].dropna())
    
    # Print counts
    for molecule, count in molecule_counts.items():
        print(f"{molecule}: {count}")
else:
    print("The column 'Molecule_ChEMBL_ID' does not exist in the CSV file.")
    molecule_counts = Counter() 

# Create a DF from the strain counts
molecule_counts_df = pd.DataFrame(molecule_counts.items(), columns=["Molecule", "Count"])

# File path for the new CSV
output_csv_path = "./p1_preprocessing/3 - Graphing data/chembl_distribution_histograms/molecule_counts.csv"

# Save the DF to a new CSV file
molecule_counts_df.to_csv(output_csv_path, index=False)

print(f"Strain counts have been saved to {output_csv_path}")
