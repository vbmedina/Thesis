import re
from pathlib import Path
import pandas as pd

# ── 1. FILE LOCATIONS ────────────────────────────────────────────
IN_FILE  = Path("./Do Not Touch/postphase2_deleteEquiv.csv")
OUT_FILE_ASEC = Path("./Do Not Touch/postphase3_Asexual_Only.csv")
OUT_FILE_SEC = Path("./Do Not Touch/postphase3_Sexual_Only.csv")

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
ASSAY_COLS = df.columns.values

# ── 5. FLAG & REMOVE NON-ASEXUAL ROWS ───────────────────────────
is_sec = df[ASSAY_COLS].fillna("").apply(
    lambda row: any(STAGE_PAT.search(str(cell)) for cell in row), axis=1
)

n_total     = len(df)
n_sec       = int(is_sec.sum())
n_asec      = n_total - n_sec

print(f"• Scanned {n_total:,} rows")
print(f"• Removing {n_sec:,} non-asexual rows "
      f"({n_sec/n_total:.2%})")

df_asec = df[~is_sec].copy()
df_sec = df[is_sec].copy()

# ── 6. SAVE CLEAN DATA ─────────────────────────────────────────
df_sec.to_csv(OUT_FILE_SEC, index=False)

print(f"Saved sexual dataset: {OUT_FILE_SEC} "
      f"({n_sec:,} rows)")

df_asec.to_csv(OUT_FILE_ASEC, index=False)

print(f"Saved strictly-asexual dataset: {OUT_FILE_ASEC} "
      f"({n_asec:,} rows)")