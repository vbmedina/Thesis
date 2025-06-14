import re
from pathlib import Path
import pandas as pd

# ── 1. FILE LOCATIONS ────────────────────────────────────────────
IN_FILE  = Path("/Users/victoriamedina/Thesis_Project/Thesis/CHEMBL_Feb_5_Sexual_and_Asexual.csv")
OUT_FILE = Path("/Users/victoriamedina/Thesis_Project/Thesis/"
                "CHEMBL_Feb_5_Asexual_Only.csv")

# ── 2. LOAD DATA ────────────────────────────────────────────────
df = pd.read_csv(IN_FILE)

# ── 3. KEYWORD PATTERN FOR NON-ASEXUAL ASSAYS ──────────────────
STAGE_PAT = re.compile(r"""
    (gamet|gametocyte|stage\s?[IVV]+|transmission|block|zygot|ook|
     pfs16|pfs25|pfs230|luc|luciferase|gfp|cbg99?|promoter|reporter|attb|elo1|
     mosq|oocyst|midgut|sporo|sporozoite|
     liver|hep\s?[0-9]*|hepatocyte|frg|huh|hc04|invasion)
""", re.IGNORECASE | re.VERBOSE)

# ── 4. COLUMNS THAT CONTAIN ASSAY TEXT ─────────────────────────
ASSAY_COLS = [c for c in ["Assay Description",
                          "Assay Name",
                          "Description",
                          "Assay_Title"]
              if c in df.columns]

if not ASSAY_COLS:
    raise ValueError(
        "No assay-description column found."
        "Add its name to ASSAY_COLS before running."
    )

# ── 5. FLAG & REMOVE NON-ASEXUAL ROWS ───────────────────────────
is_non_asec = df[ASSAY_COLS].fillna("").apply(
    lambda row: any(STAGE_PAT.search(str(cell)) for cell in row), axis=1
)

n_total     = len(df)
n_non_asec  = int(is_non_asec.sum())
n_asec      = n_total - n_non_asec

print(f"• Scanned {n_total:,} rows")
print(f"• Removing {n_non_asec:,} non-asexual rows "
      f"({n_non_asec/n_total:.2%})")

df_asec = df[~is_non_asec].copy()

# ── 6. SAVE CLEAN DATA ─────────────────────────────────────────
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
df_asec.to_csv(OUT_FILE, index=False)

print(f"Saved strictly-asexual dataset: {OUT_FILE} "
      f"({n_asec:,} rows)")

#CHANGES: add column to define sexual type. split into two files, one for sexual and one for asexual.