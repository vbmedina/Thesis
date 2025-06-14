#!/usr/bin/env python3
"""
apply_strain_mapping.py
───────────────────────
Adds a 'Standardized Strain' column to CHEMBL_Feb_5_Phase_3.csv
based on strain_mapping_v1.csv, then writes a new copy.

Run:
    python apply_strain_mapping.py
"""

from pathlib import Path
import pandas as pd
import re, datetime as dt

# ── 1.  FILE LOCATIONS ──────────────────────────────────────────
IN_CSV     = Path("/Users/victoriamedina/Thesis_Project/Thesis/"
                  "CHEMBL_Feb_5_Phase_3.csv")
MAPPING_CSV = Path("/Users/victoriamedina/Thesis_Project/Thesis/"
                   "strain_mapping_v1.csv")
OUT_CSV    = IN_CSV.with_name("CHEMBL_Feb_5_Phase_3_with_standard_strain.csv")

# ── 2.  LOAD DATA ───────────────────────────────────────────────
df       = pd.read_csv(IN_CSV)
mapping  = pd.read_csv(MAPPING_CSV)

# build dict: every raw variant → strain_clean
map_dict = {}
for _, row in mapping.iterrows():
    clean   = row["strain_clean"]
    for raw in str(row["raw_variants"]).split(";"):
        map_dict[raw.strip()] = clean

# ── 3.  APPLY MAPPING ───────────────────────────────────────────
def standardize(raw):
    return map_dict.get(raw, raw)      # default: leave as-is if not found

df["Standardized Strain"] = df["Original_Strain"].apply(standardize)

# ── 4.  SAVE NEW COPY ───────────────────────────────────────────
df.to_csv(OUT_CSV, index=False)
print(f"✅  Wrote file with standardized strains → {OUT_CSV}")
print("    Rows:", len(df), "| unique strains:",
      df['Standardized Strain'].nunique())



# In every stage of the pipeline, we should keep track of the number of deleted rows