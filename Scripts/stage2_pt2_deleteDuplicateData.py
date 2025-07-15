from pathlib import Path
import pandas as pd

# ── 1.  FILE LOCATION  ──────────────────────────────────────────
CSV_PATH = Path("./Do Not Touch/postphase2_convertedUnits.csv")
CSV_PATH_OUT = Path("./Do Not Touch/postphase2_deleteDuplicates.csv")

# ── 2.  LOAD  ──────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# ── 3.  delete duplicate data  ──────────────────────────
# Delete duplicate rows
before = len(df)
df = df.drop_duplicates()
after  = len(df)

# ── 4.  SAVE (OVERWRITE)  ──────────────────────────────────────
df.to_csv(CSV_PATH_OUT, index=False)

print(f"Wrote {CSV_PATH_OUT}")
print(f" Removed {before - after:,} rows that were duplicates"
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")