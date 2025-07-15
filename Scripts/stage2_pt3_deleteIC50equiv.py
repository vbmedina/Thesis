from pathlib import Path
import pandas as pd

# ── 1.  FILE LOCATION  ──────────────────────────────────────────
CSV_PATH = Path("./Do Not Touch/postphase2_deleteDuplicates.csv")
CSV_PATH_OUT = Path("./Do Not Touch/postphase2_deleteEquiv.csv")

# ── 2.  LOAD  ──────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# ── 3.  FILTER ROWS BASED ON "Standard Relation"  ──────────────
# Keep rows where "Standard Relation" does not contain "~".
before = len(df)
df = df.dropna(subset=["Standard Relation"])
df = df[~df["Standard Relation"].str.contains("~", na=False)]
df = df[~df["Standard Relation"].str.contains(">", na=False)]
df = df[~df["Standard Relation"].str.contains("≈", na=False)]
df = df[~df["Standard Relation"].str.contains("NA", na=False)]
df = df[~df["Standard Relation"].str.contains("0", na=False)]
after = len(df)

# ── 4.  SAVE (OVERWRITE)  ──────────────────────────────────────
df.to_csv(CSV_PATH_OUT, index=False)

print(f"Overwrote {CSV_PATH_OUT}")
print(f" Removed {before - after:,} rows with incorrect 'Standard Relation' values"
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")

#Print all unique values in the "Standard Relation" column
unique_values = df["Standard Relation"].unique()
print("Unique values in 'Standard Relation':", unique_values)