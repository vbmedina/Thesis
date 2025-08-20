# Description: For stage 2 of preprocessing (IC50), delete any rows with implausible pIC50 values (above 11 or below 4)
# and label the remaining rows as inactive, active, or high potency.

# Import
import pandas as pd
import numpy as np

# Paths
in_path   = "./do_not_touch/postphase2_CorrectpIC50.csv"      
out_path  = "./do_not_touch/postphase2_cleaned.csv"               

# Load the data
df = pd.read_csv(in_path)
before = len(df)                              

# Define chemically plausible window for pIC50
lower_bound = 4.0      #  pIC50  < 4  ⇒  IC50 > 100 uM  = biologically meaningless
upper_bound = 11.0     #  pIC50  > 11 ⇒  IC50 < 0.01 nM = almost always an artefact

# Drop rows outside of window
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