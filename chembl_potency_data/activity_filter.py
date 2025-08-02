from pathlib import Path
import pandas as pd

# Path to CSV file
csv = Path("./potency_str.csv")

# Load CSV
df = pd.read_csv(csv)

# col = 'Comment'
# s = df[col]

# # Function to determine the type of comment
# type_counts = s.map(type).value_counts(dropna=False)
# print(type_counts)

# Delete inconclusive data
df = df[~df[comment].str.contains('inconclusive', case=False, na=False)]

# ── 4.  SAVE (OVERWRITE)  ──────────────────────────────────────
df.to_csv(csv, index=False)
