''' 
Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
1) This script: Adds activity boundaries and defines active and inactive.

Preconditions:
1) "postphase2_CorrectpIC50.csv" - generated from stage2 pt4
'''

# Import
import pandas as pd
import numpy as np

# Paths
in_path = "./p0_all_csvs/postphase2_4.csv"      
out_path = "./p0_all_csvs/postphase2_5.csv"               

# Load the data
df = pd.read_csv(in_path)
before = len(df)                              


lower_bound = 4.0
upper_bound = 11.0

# Drop pIC50 values outside of window
df = df[df["pIC50"].between(lower_bound, upper_bound)].copy()
after = len(df)

# Label the remaining rows as inactive, active, or high potency
df["potency_label"] = "inactive"
df.loc[df["pIC50"] >= 6.0, "potency_label"]  = "active"
df.loc[df["pIC50"] >= 7.5, "potency_label"]  = "high_potency"

# Save the table
df.to_csv(out_path, index=False)

# Output results
print(f"Overwrote {out_path}")
print(f"Removed {before - after:,} rows with implausible pIC50 values "
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")