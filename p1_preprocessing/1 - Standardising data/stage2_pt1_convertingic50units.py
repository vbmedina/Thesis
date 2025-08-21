''' Description: This script is the second stage to the pipeline (IC50). It works on the IC50 scores.
This script:
1) converts the IC50 values to nM
2) removes rows with an empty, negative, or 0 IC50 value.

Preconditions:
1) "postphase1_strainsmerged.csv" - generated from stage1'''


# Imports
from pathlib import Path
import re, shutil, datetime as dt
import pandas as pd

# Path
csv = Path("./p0_all_csvs/postphase1_strainsmerged.csv")
out = Path("./p0_all_csvs/postphase2_1.csv")

# Load
df = pd.read_csv(csv, low_memory=False) 

unit_col = "Standard Units" 
val_col  = "Standard_Value"

# Check columns and if "Standard Value" rows are empty if not, drop them. Also drop rows with values <= 0
before = len(df)
df = df.dropna(subset=[val_col])
df = df[df[val_col] > 0]
after = len(df)
totalDropped = before - after

# Print statement for deleted rows with empty IC50 values
print(f"   Empty IC50   : {totalDropped:,}  ({(totalDropped)/before:.2%})")
print(f"   Rows remaining : {after:,}")

# Canonicalize unit strings
def canon(u: str) -> str:
    if pd.isna(u):
        return ""
    u = u.upper()
    return re.sub(r"[^A-Z0-9]", "", u)

# Factors to multiply by to turn units to nM
FACTOR = {
    "NM":        1,
    "NANOMOLAR": 1,
    "UM":        1_000,         
    "MM":        1_000_000,    
    "PM":        0.001,
}


# Regex for strings like "10^-2microM", "10^-5 uM"
exp_pat = re.compile(r"10\^?-?(-?\d+)\s*(?:U|MICRO)?M", re.I)
# Regex for strings like "10'-5 g/L"
exp_pat_concentration = re.compile(r"(?:10\^?'?-?(-?\d+)\s*)?(U?)(?:G)/(M?)L", re.I)

# Convert
unknown_units = {}    
convertible   = 0
dropped       = 0

# Function to convert units to nM. There are three cases:
# A. Direct mapping uM or MM to nM
# B. Change "10^-2microM", "10^-x µM" to nM
# C. Change "10^-x g/L" to nM
# D. Remove unknown units
def convert_row(row):
    global convertible, dropped
    raw_unit = str(row[unit_col]).strip()

    # Case A – direct mapping
    cu = canon(raw_unit)
    if cu in FACTOR:
        convertible += 1
        return row[val_col] * FACTOR[cu]
    
# Dropping this and need to make it ug/mL,10^-5, g/l, 10$^-4g/L and 10^-3g/L and make it usuable by multiplying by MW in csv
    # Case B – "10^-x µM"
    m = exp_pat.fullmatch(raw_unit.replace(" ", ""))
    if m:
        exponent = int(m.group(1))  
        factor   = (10 ** exponent) * 1_000        # µM to nM
        convertible += 1
        return row[val_col] * factor

    # Case C – "10^-x g/L"
    m = exp_pat_concentration.fullmatch(raw_unit.replace(" ", ""))
    if m:
        exponent = int(m.group(1))
        is_ug = m.group(2).upper() == "U"
        is_mL = m.group(3).upper() == "M"
        factor = (10 ** exponent) * 1e9 * (0.000001 if is_ug else 1) * (1000 if is_mL else 1)
        convertible += 1
        return row[val_col] * factor
    
    # Case D – unknown / non-convertible
    unknown_units[raw_unit] = unknown_units.get(raw_unit, 0) + 1
    dropped += 1
    return None

# Run the conversion
df[val_col] = df.apply(convert_row, axis=1)

# DROP rows that couldn't be converted & overwrite unit column to be "nM"
before = len(df)
df = df.dropna(subset=[val_col])
after = len(df)
df[unit_col] = "nM"

# SAVE
df.to_csv(out, index=False)
print(f"\n  Wrote {out}")
print(f"   Converted rows : {convertible:,}")
print(f"   Dropped rows   : {dropped:,}  ({dropped/before:.2%})")
print(f"   Rows remaining : {before:,}")

# Print unknown units
if unknown_units:
    print("\n  Units NOT converted (dropped):")
    for u, c in sorted(unknown_units.items(), key=lambda x: (-x[1], x[0])):
        print(f"  • {u:<12}  {c:>6} rows")
print()
print(f"Total {dropped + totalDropped:,} removed.")
