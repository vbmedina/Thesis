''' Description: This script is the third stage to the pipeline (Filtering Sexual Assays).
This script:
1) Clasifies rows as asexual or sexual assays. - clasifies by specific keywords in assay descriptions

Preconditions:
1) "postphase2_cleaned.csv" - generated from stage2 pt5'''

# Import
import re
from pathlib import Path
import pandas as pd

# Paths
csv_in  = Path("./p0_all_csvs/postphase2_5.csv")
out_fir = Path("./p0_all_csvs/postphase3_Asexual_Only.csv")
out_sec = Path("./p0_all_csvs/postphase3_Sexual_Only.csv")

# Load DataFrame
df = pd.read_csv(csv_in)

# Keywords for identifying sexual stage assays. All these keywords are related to sexual stages of the malaria parasite.
stage_pat = re.compile(r"""
    (gamet|gametocyte|stage\s?[IVV]+|transmission|block|zygot|ook|
     pfs16|pfs25|pfs230|luc|luciferase|gfp|cbg99?|promoter|reporter|attb|elo1|
     mosq|oocyst|midgut|sporo|sporozoite|
     liver|hep\s?[0-9]*|hepatocyte|frg|huh|hc04|invasion)
""", re.IGNORECASE | re.VERBOSE)

# Creates a bool per row. True indicates a keyword was found in any of the row's cells
assay_cols = df.columns.values
is_sex = df[assay_cols].fillna("").apply(
    lambda row: any(stage_pat.search(str(cell)) for cell in row), axis=1)

# Split the DataFrame into sexual and asexual
df_asec = df[~is_sex].copy()
df_sec = df[is_sex].copy()

# Count number of sexual and asexual rows
n_total     = len(df)
n_sex       = int(is_sex.sum())
n_asec      = n_total - n_sex

# Print results
print(f"• Scanned {n_total:,} rows")
print(f"• Removing {n_sex:,} non-asexual rows "
      f"({n_sex/n_total:.2%})")

# Save the filtered DataFrames
df_sec.to_csv(out_sec, index=False)
print(f"Saved sexual dataset: {out_sec} "
      f"({n_sex:,} rows)")

df_asec.to_csv(out_fir, index=False)
print(f"Saved strictly-asexual dataset: {out_fir} "
      f"({n_asec:,} rows)")