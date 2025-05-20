from pathlib import Path
import pandas as pd

# ── 1.  FILE LOCATION  ──────────────────────────────────────────
CSV_PATH = Path("/Users/victoriamedina/Thesis_Project/Thesis/CHEMBL_Feb_5_Final_Cleaned_copy.csv")

# ── 2.  LOAD  ──────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# ── 3.  FILTER ROWS BASED ON "Standard Relation"  ──────────────
# Keep rows where "Standard Relation" does not contain "~".
before = len(df)
df = df[~df["Standard Relation"].str.contains("~", na=False)]
after = len(df)

# ── 4.  SAVE (OVERWRITE)  ──────────────────────────────────────
df.to_csv(CSV_PATH, index=False)

print(f"Overwrote {CSV_PATH}")
print(f" Removed {before - after:,} rows without a numeric IC50 "
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")

#Print all unique values in the "Standard Relation" column
unique_values = df["Standard Relation"].unique()
print("Unique values in 'Standard Relation':", unique_values)