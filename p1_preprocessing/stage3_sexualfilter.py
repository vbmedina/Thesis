# Description: For stage 3 of preprocessing (Filtering Sexual Assays), this script filters out rows from a CSV file based
# on specific keywords in assay descriptions to separate asexual and sexual assays.

# Import
import re
from pathlib import Path
import pandas as pd

# Paths
csv_in  = Path("./do_not_touch/postphase2_cleaned.csv")
out_fir = Path("./do_not_touch/postphase3_Asexual_Only.csv")
out_sec = Path("./do_not_touch/postphase3_Sexual_Only.csv")

# Load DataFrame
df = pd.read_csv(csv_in)

# Keywords for identifying sexual stage assays
stage_pat = re.compile(r"""
    (gamet|gametocyte|stage\s?[IVV]+|transmission|block|zygot|ook|
     pfs16|pfs25|pfs230|luc|luciferase|gfp|cbg99?|promoter|reporter|attb|elo1|
     mosq|oocyst|midgut|sporo|sporozoite|
     liver|hep\s?[0-9]*|hepatocyte|frg|huh|hc04|invasion)
""", re.IGNORECASE | re.VERBOSE)

assay_cols = df.columns.values

# Flag and remove sexual assay rows 
is_sec = df[assay_cols].fillna("").apply(
    lambda row: any(stage_pat.search(str(cell)) for cell in row), axis=1)

# Helper to count and separate rows
n_total     = len(df)
n_sec       = int(is_sec.sum())
n_asec      = n_total - n_sec

# Output results with print statements
print(f"• Scanned {n_total:,} rows")
print(f"• Removing {n_sec:,} non-asexual rows "
      f"({n_sec/n_total:.2%})")

df_asec = df[~is_sec].copy()
df_sec = df[is_sec].copy()

# Save the filtered DataFrames
df_sec.to_csv(out_sec, index=False)
print(f"Saved sexual dataset: {out_sec} "
      f"({n_sec:,} rows)")

df_asec.to_csv(out_fir, index=False)
print(f"Saved strictly-asexual dataset: {out_fir} "
      f"({n_asec:,} rows)")