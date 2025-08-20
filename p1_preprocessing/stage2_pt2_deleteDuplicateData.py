# Description: This script removes exact duplicate rows from a CSV file.

from pathlib import Path
import pandas as pd

# Paths
csv_in = Path("./Do Not Touch/postphase2_convertedUnits.csv")
csv_out = Path("./Do Not Touch/postphase2_deleteDuplicates.csv")

# Load DF
df = pd.read_csv(csv_in)

# Delete duplicates
before = len(df)
df = df.drop_duplicates()
after  = len(df)

# Save
df.to_csv(csv_out, index=False)

print(f"Wrote {csv_out}")
print(f" Removed {before - after:,} rows that were duplicates"
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")