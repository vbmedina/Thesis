import pandas as pd
from collections import Counter

# File path to the CSV
csv_file_path = "/Users/victoriamedina/Thesis_Project/Thesis/ProjectFiles/chembl_final.csv"

data = pd.read_csv(csv_file_path)

print(data.columns)

# Verify the column "stand_strain" exists
if "Standardized_Strain" in data.columns:
    # Count occurrences of each item in the "stand_strain column
    strain_counts = Counter(data["Standardized_Strain"].dropna())
    
    # Print the counts
    for strain, count in strain_counts.items():
        print(f"{strain}: {count}")
else:
    print("The column 'Standardized_Strain' does not exist in the CSV file.")
    strain_counts = Counter()  # Initialize an empty Counter if the column doesn't exist

# Create a DataFrame from the strain counts
strain_counts_df = pd.DataFrame(strain_counts.items(), columns=["Strain", "Count"])

# File path for the new CSV
output_csv_path = "/Users/victoriamedina/Thesis_Project/Thesis/str_counts_final.csv"

# Save the DataFrame to a new CSV file
strain_counts_df.to_csv(output_csv_path, index=False)

print(f"Strain counts have been saved to {output_csv_path}")
