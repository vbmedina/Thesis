import pandas as pd
import numpy as np

# Paths
in_path   = "./do_not_touch/postphase2_CorrectpIC50.csv"      # raw file
out_path  = "./do_not_touch/postphase2_cleaned.csv"               # cleaned file

# DF
df = pd.read_csv(in_path)
before = len(df)                               # initial row count

# Define chemically plausible window for pIC50
lower_bound = 4.0      #  pIC50  < 4  ⇒  IC50 > 100 µM  = biologically meaningless
upper_bound = 11.0     #  pIC50  > 11 ⇒  IC50 < 0.01 nM = almost always an artefact

# Drop rows outside of window
df = df[df["pIC50"].between(lower_bound, upper_bound)].copy()
after = len(df)

# Label the remaining rows as inactive, active, or high potency
df["potency_label"] = "inactive"
df.loc[df["pIC50"] >= 6.0, "potency_label"]  = "active"
df.loc[df["pIC50"] >= 7.5, "potency_label"]  = "high_potency"

# Save the cleaned table
df.to_csv(out_path, index=False)

# Output results
print(f"Overwrote {out_path}")
print(f"Removed {before - after:,} rows with implausible pIC50 values "
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")