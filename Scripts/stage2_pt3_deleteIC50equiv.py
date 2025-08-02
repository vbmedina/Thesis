from pathlib import Path
import pandas as pd

# Paths
csv_in = Path("./Do Not Touch/postphase2_deleteDuplicates.csv")
csv_out = Path("./Do Not Touch/postphase2_deleteEquiv.csv")

# Load DF
df = pd.read_csv(csv_in)

# Filter on "Standard Relation" 
before = len(df)
df = df.dropna(subset=["Standard Relation"])
df = df[~df["Standard Relation"].str.contains("~", na=False)]
df = df[~df["Standard Relation"].str.contains(">", na=False)]
df = df[~df["Standard Relation"].str.contains("≈", na=False)]
df = df[~df["Standard Relation"].str.contains("NA", na=False)]
df = df[~df["Standard Relation"].str.contains("0", na=False)]
after = len(df)

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