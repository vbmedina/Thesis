''' Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
This script:
1) Converts IC50 ("Standard_Value") values to pIC50 values - new column "pIC50"


Preconditions:
1) "postphase2_deleteEquiv.csv" - generated from stage2 pt3'''

# Imports
import pandas as pd
import numpy as np

# Path
IN_CSV = "./p0_all_csvs/postphase2_3.csv"
OUT_CSV = "./p0_all_csvs/postphase2_4.csv"

# Load the data
df = pd.read_csv(IN_CSV)

# Sanity check: Numeric values should remain in "Standard Value". Print any non-numeric values
for i, value in df["Standard_Value"].items():
    if not isinstance(value, (int, float)):
        print(f"Row {i} has non-numeric Standard Value: {value}")

# Convertion tool to convert IC50 values in nM to pIC50 values
def convert(ic50_nM):
    if ic50_nM <= 0 or pd.isna(ic50_nM):
        return np.nan
    ic50_M = ic50_nM * 1e-9
    return -np.log10(ic50_M)  

# Apply the tool on 'Standard_Value' to create 'pIC50'. Round to 2 decimal places
df['pIC50'] = df['Standard_Value'].apply(convert).round(2)

# Sanity check for pIC50 values
print(df[['Standard_Value','pIC50']].head())

# Save the dataframe to a new CSV file
df.to_csv(OUT_CSV, index=False)