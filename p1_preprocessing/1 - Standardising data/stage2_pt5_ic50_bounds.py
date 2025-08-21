''' Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
This script:
1) deletes rows with implausible pIC50 values (above 11 or below 4)
2) label rows as inactive, active, or high potency.

Preconditions:
1) "postphase2_CorrectpIC50.csv" - generated from stage2 pt4'''

# Import
import pandas as pd
import numpy as np

# Paths
in_path   = "./p0_all_csvs/postphase2_4.csv"      
out_path  = "./p0_all_csvs/postphase2_5.csv"               

# Load the data
df = pd.read_csv(in_path)
before = len(df)                              

# Define chemically plausible window for pIC50
lower_bound = 4.0      #  pIC50  < 4  ⇒  IC50 > 100 uM  = biologically meaningless
upper_bound = 11.0     #  pIC50  > 11 ⇒  IC50 < 0.01 nM = almost always an artefact

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