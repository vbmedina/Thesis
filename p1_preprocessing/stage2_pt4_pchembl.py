# Description: For stage 2 of preprocessing (IC50), this script converts IC50 ("Standard_Value") values to pIC50 values 
# and adds them to a new column in the CSV file a new column called "pIC50".

# Imports
import pandas as pd
import numpy as np

# Path
IN_CSV = "./do_not_touch/postphase2_deleteEquiv.csv"
OUT_CSV = "./do_not_touch/postphase2_CorrectpIC50.csv"

# Load the data
df = pd.read_csv(IN_CSV)

# Sanity check: Numeric values should remain in "Standard Value"
for i, value in df["Standard_Value"].items():
    if not isinstance(value, (int, float)):
        print(f"Row {i} has non-numeric Standard Value: {value}")

# Convertion tool to convert IC50 values in nM to pIC50 values
def convert(ic50_nM):
    if ic50_nM <= 0 or pd.isna(ic50_nM):
        return np.nan
    ic50_M = ic50_nM * 1e-9
    return -np.log10(ic50_M)  

# Apply the tool on 'Standard_Value' to create 'pIC50'
df['pIC50'] = df['Standard_Value'].apply(convert).round(2)

# Sanity check for pIC50 values
print(df[['Standard_Value','pIC50']].head())

# Save the dataframe to a new CSV file
df.to_csv(OUT_CSV, index=False)