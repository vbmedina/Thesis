from pathlib import Path
import re, shutil, datetime as dt
import pandas as pd

# ── 1.  PATHS ───────────────────────────────────────────────────
CSV = Path("./Do Not Touch/postphase1_strainsMerged.csv")
CSVOUT = Path("./Do Not Touch/postphase2_convertedUnits.csv")

# ── 2.  LOAD ────────────────────────────────────────────────────
df = pd.read_csv(CSV, low_memory=False)   # low_memory=False suppresses dtype warning

unit_col = "Standard Units"               # change if your header differs
val_col  = "Standard_Value"

before = len(df)
df = df.dropna(subset=[val_col])
after = len(df)
totalDropped = before - after
# ── 2.1  DELETE NA ────────────────────────────────────────
print(f"   Empty IC50   : {totalDropped:,}  ({(totalDropped)/before:.2%})")
print(f"   Rows remaining : {after:,}")

# ── 3.  HELPER: canonicalise unit strings ───────────────────────
def canon(u: str) -> str:
    """
    Upper-case, strip non-alphanumerics
    e.g. "uM well-1" → "UMWELL1"
    """
    if pd.isna(u):
        return ""
    u = u.upper()
    return re.sub(r"[^A-Z0-9]", "", u)

# direct factors (canon_string → multiplier for value → nM)
FACTOR = {
    "NM":        1,
    "NANOMOLAR": 1,
    "UM":        1_000,         # µM
    "UMWELL1":   1_000,         # µM well-1
    "MM":        1_000_000,     # mM
    "PM":        0.001,         # pM  (rare)
}

# regex for strings like "10^-2microM", "10^-5 uM"
exp_pat = re.compile(r"10\^?-?(-?\d+)\s*(?:U|MICRO)?M", re.I)
# regex for strings like "10'-5 g/L"
exp_pat_concentration = re.compile(r"(?:10\^?'?-?(-?\d+)\s*)?(U?)(?:G)/(M?)L", re.I)

# ── 4.  CONVERT row-by-row ──────────────────────────────────────
unknown_units = {}     # unit_string → count
convertible   = 0
dropped       = 0

def convert_row(row):
    global convertible, dropped
    raw_unit = str(row[unit_col]).strip()

    # Case A – direct mapping
    cu = canon(raw_unit)
    if cu in FACTOR:
        convertible += 1
        return row[val_col] * FACTOR[cu]

    # Case B – "10^-x µM"
    m = exp_pat.fullmatch(raw_unit.replace(" ", ""))
    if m:
        exponent = int(m.group(1))                 # e.g. −2
        factor   = (10 ** exponent) * 1_000        # µM → nM
        convertible += 1
        return row[val_col] * factor

    # Case C – "10^-x g/L"
    m = exp_pat_concentration.fullmatch(raw_unit.replace(" ", ""))
    if m:
        exponent = int(m.group(1))                 # e.g. −2
        is_ug = m.group(2).upper() == "U"          # e.g. "U" for µ
        is_mL = m.group(3).upper() == "M"          # e.g. "M" for mL
        factor = (10 ** exponent) * 1e9 * (0.000001 if is_ug else 1) * (1000 if is_mL else 1)      # ug/mL → g/L or g/L → g/L
        divisor = row["Molecular Weight"]          # g/L → nM
        convertible += 1
        return row[val_col] * factor / (divisor if divisor else 1)  # avoid division by zero
    
    # Case D – unknown / non-convertible
    unknown_units[raw_unit] = unknown_units.get(raw_unit, 0) + 1
    dropped += 1
    return None

df[val_col] = df.apply(convert_row, axis=1)

# ── 5.  DROP non-convertible rows & fix unit column ─────────────
before = len(df)
df = df.dropna(subset=[val_col])
after = len(df)
df[unit_col] = "nM"

# ── 6.  SAVE (overwrite) ────────────────────────────────────────
df.to_csv(CSVOUT, index=False)
print(f"\n  Wrote {CSVOUT}")
print(f"   Converted rows : {convertible:,}")
print(f"   Dropped rows   : {dropped:,}  ({dropped/before:.2%})")
print(f"   Rows remaining : {before:,}")

if unknown_units:
    print("\n  Units NOT converted (dropped):")
    for u, c in sorted(unknown_units.items(), key=lambda x: (-x[1], x[0])):
        print(f"  • {u:<12}  {c:>6} rows")
print()
print(f"Total {dropped + totalDropped:,} removed.")
#needs to try to salvage units that have ug/mL, 10^15 g/L, 10^-4 g/l, 10^-3 g/L should check if thse are all the units (Standard Units) that are used using Standard_Value and Molecular Weight column
#can delete ug/well, and nM/g, and uMxhr