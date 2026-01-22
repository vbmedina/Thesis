''' 
Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
This script:
1) Removes exact duplicate rows from a CSV file.

Preconditions:
1) "postphase2_convertedUnits.csv" - generated from stage2 pt1
'''

from pathlib import Path
import pandas as pd

# Paths
csv_in = Path("./p0_all_csvs/postphase2_1.csv")
csv_out = Path("./p0_all_csvs/postphase2_2.csv")

# Load DF
df = pd.read_csv(csv_in)

# Delete duplicate rows
# Note: This will keep the first occurrence of each duplicate row and remove subsequent duplicates.
before = len(df)
df = df.drop_duplicates()
after = len(df)

# Save
df.to_csv(csv_out, index=False)

print(f"Wrote {csv_out}")
print(f" Removed {before - after:,} rows that were duplicates"
      f"({(before - after) / before:.2%})")
print(f"Remaining rows: {after:,}")