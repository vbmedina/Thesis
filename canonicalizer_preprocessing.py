import pandas as pd
import numpy as np
import re
from difflib import SequenceMatcher
import itertools

counts = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/strain_counts_3.csv")     # Strain, Count
chembl = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/CHEMBL_Feb_5_Canon.csv")        # your cleaned big file

# tell the script which column holds the raw strain label
STRAIN_COL = "stand_strain"          # ← CHANGE HERE if your file differs

if STRAIN_COL not in chembl.columns:
    raise KeyError(
        f"‘{STRAIN_COL}’ not found. Available columns:\n{chembl.columns.tolist()}"
    )

# keep naming tidy
chembl = chembl.rename(columns={STRAIN_COL: "raw_strain"})

# ───────────────────────────────────────────────────────────────
# 2.  Canonicalise every label
# ───────────────────────────────────────────────────────────────
def canonical(name: str) -> str:
    name = name.strip().upper()
    return re.sub(r"[^A-Z0-9]", "", name)

counts["canon"]       = counts["Strain"].apply(canonical)
chembl["canon_raw"]   = chembl["raw_strain"].apply(canonical)

# ───────────────────────────────────────────────────────────────
# 3.  Exact-synonym merge (spelling / punctuation only)
# ───────────────────────────────────────────────────────────────
syn_groups = (
    counts.groupby("canon")["Strain"]
    .apply(list)
    .to_dict()
)
rep = {c: sorted(names, key=len)[0] for c, names in syn_groups.items()}

# ───────────────────────────────────────────────────────────────
# 4.  Flag “possible derivative” pairs (HB3-leuR1 vs HB3, etc.)
#     – uses the rules table we discussed last time
# ───────────────────────────────────────────────────────────────
SIM_THRESHOLD = 0.95
EXTRA_CHARS   = 2
suspects = []

for c1, c2 in itertools.combinations(rep.keys(), 2):
    ratio  = SequenceMatcher(None, c1, c2).ratio()
    ldiff  = abs(len(c1) - len(c2))
    if ratio > SIM_THRESHOLD and ldiff <= EXTRA_CHARS:
        suspects.append((rep[c1], rep[c2]))

pd.DataFrame(suspects, columns=["parent", "possible_derivative"]) \
  .to_csv("strain_derivative_review.csv", index=False)

print("▶ Review", len(suspects), "questionable pairs in strain_derivative_review.csv")

# ───────────────────────────────────────────────────────────────
# 5.  Build mapping  (edit `rep` after manual review if needed)
# ───────────────────────────────────────────────────────────────
mapping = {orig: rep[canon]
           for canon, names in syn_groups.items()
           for orig in names}

# ───────────────────────────────────────────────────────────────
# 6.  Apply mapping and save
# ───────────────────────────────────────────────────────────────
chembl["Canon_Strain"] = chembl["Standard_Strain"].map(mapping) \
                                              .fillna(chembl["Before_Strain"])

chembl.to_csv("CHEMBL_Feb_5_with_clean_strains.csv", index=False)
print("Cleaned file written: CHEMBL_Feb_5_with_clean_strains.csv")