''' Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
This script: 
1) deletes rows where "Standard Relation" column contains "~", ">", "≈", "NA", and "0".

Preconditions:
1) "postphase2_deleteDuplicates.csv" - generated from stage2 pt2'''

# Imports
from pathlib import Path
import pandas as pd

# Paths
csv_in = Path("./p0_all_csvs/postphase2_2.csv")
csv_out = Path("./p0_all_csvs/postphase2_3.csv")

# Load DF
df = pd.read_csv(csv_in)

# If the "Standard Relation" row contains a ~, >, ≈, NA or 0, remove it
before = len(df)
df = df.dropna(subset=["Standard Relation"])
df = df[~df["Standard Relation"].str.contains("~", na=False)]
df = df[~df["Standard Relation"].str.contains(">", na=False)]
df = df[~df["Standard Relation"].str.contains("≈", na=False)]
df = df[~df["Standard Relation"].str.contains("NA", na=False)]
df = df[~df["Standard Relation"].str.contains("0", na=False)]
after = len(df)

print(df["Standard Relation"].unique())

# Save CSV 
df.to_csv(csv_out, index=False)

# Print results
print(f"Overwrote {csv_out}")
print(f" Removed {before - after:,} rows with incorrect 'Standard Relation' values"
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")

# Print all unique values in the "Standard Relation" column
unique_values = df["Standard Relation"].unique()
print("Unique values in 'Standard Relation':", unique_values)