import pandas as pd
import numpy as np

#Path
IN_CSV = "./Do Not Touch/postphase5.csv"

# Load the data
df = pd.read_csv(IN_CSV)

# Filter convert IC50 to pIC50
def convert(ic50_nM):
    if ic50_nM <= 0 or pd.isna(ic50_nM):
        return np.nan
    ic50_M = ic50_nM * 1e-9
    return -np.log10(ic50_M)  # Convert IC50 in nanomolar to pIC50 in M

# Apply the conversion to new column
df['pIC50'] = df['Standard_Value'].apply(convert).round(2)

# Sanity check for pIC50 values
print(df[['Standard_Value','pIC50']].head())

# Save the dataframe to a new CSV file
df.to_csv("./ProjectCSV/ppchembl.csv", index=False)