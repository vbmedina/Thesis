from pathlib import Path
import pandas as pd

# ── 1.  FILE LOCATION  ──────────────────────────────────────────
CSV_PATH = Path("/Users/victoriamedina/Thesis_Project/Thesis/"
                "CHEMBL_Feb_5_Asexual_Only.csv")

# ── 2.  LOAD  ──────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# ── 3.  COERCE → NUMERIC & DROP NaNs  ──────────────────────────
# Any blank, "n/a", or non-numeric text becomes NaN, then is dropped.
df["Standard_Value"] = pd.to_numeric(df["Standard_Value"], errors="coerce")

before = len(df)
df     = df.dropna(subset=["Standard_Value"])
after  = len(df)

# ── 4.  SAVE (OVERWRITE)  ──────────────────────────────────────
df.to_csv(CSV_PATH, index=False)

print(f"Overwrote {CSV_PATH}")
print(f" Removed {before - after:,} rows without a numeric IC50 "
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")