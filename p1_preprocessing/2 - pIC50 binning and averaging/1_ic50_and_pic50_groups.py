''' Description: This script is the 1st part in finding a measure of center for molecules found in the cleaned 
ChEMBL dataset. Many molecules in the dataset have multiple pIC50 values for the same strain. This section aims 
to aggregate these strain-molecule pairs with all of their IC50 and pIC50 values. The output is a CSV file that 
contains the number of each pairing, and their respective IC50 and pIC50 scores. 
'''

# Imports
import pandas as pd
import numpy as np

# Load DF
df = pd.read_csv("./p1_preprocessing/2 - pIC50 binning and averaging/final_data_copy.csv", low_memory=False)

#Dicts
ic50_dict = {}
pic50_dict = {}

for i, row in df.iterrows():
    strain = row['standardized_strain']
    molecule = row['Molecule_ChEMBL_ID']
    ic50 = row['Standard_Value']
    pic50 = row['pIC50']
    key = strain + "_" + molecule
    
    if key not in ic50_dict:
        ic50_dict[key] = []
        pic50_dict[key] = []

    ic50_dict[key].append(ic50)
    pic50_dict[key].append(pic50)

new_data = pd.DataFrame(columns=["Strain", "Molecule", "Frequency", "IC50 Scores", "pIC50 Scores"])

# Populate new DataFrame with collected data
for key in ic50_dict:
    strain, molecule = key.split("_")
    frequency = len(ic50_dict[key])
    ic50_scores = ic50_dict[key]
    pic50_scores = pic50_dict[key]
    
    ic50_str = " : ".join([str(x) for x in ic50_scores])
    pic50_str = " : ".join([str(x) for x in pic50_scores])
    
    new_data.loc[len(new_data)] = {
    "Strain": strain,
    "Molecule": molecule,
    "Frequency": frequency,
    "IC50 Scores": ic50_str,
    "pIC50 Scores": pic50_str
    }

# Save the new DataFrame to a CSV file
new_data.to_csv("./p1_preprocessing/2 - pIC50 binning and averaging/pairing_IC50_pIC50.csv", index=False)