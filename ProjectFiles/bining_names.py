import re
import itertools
from pathlib import Path
from difflib import SequenceMatcher
import pandas as pd

# # ── 1. CONFIGURE PATHS ─────────────────────────────────────────
# IN_CSV  = Path("/Users/victoriamedina/Thesis_Project/Thesis/Do Not Touch/CHEMBL_Feb_5_Before_Deleting_NAN.csvcsv")
# OUT_CSV = Path("/Users/victoriamedina/Thesis_Project/Thesis/TESTING.csv")

# # ── 2. LOAD DATA ────────────────────────────────────────────────
# df = pd.read_csv(IN_CSV, low_memory=False)

# # ── 3. CANONICALISE FUNCTION ────────────────────────────────────
# def canonical(s: str) -> str:
#     """Upper-case, replace Greek µ→U, strip all non-alphanumeric."""
#     if pd.isna(s):
#         return ""
#     u = s.strip().upper().replace("Μ","U").replace("μ","U").replace("µ","U")
#     return re.sub(r"[^A-Z0-9]", "", u)

# # compute canonical keys
# df["canon"] = df["Original_Strain"].apply(canonical)

# # ── 4. SYNONYM GROUPING (Rules 1 & 3) ─────────────────────────────
# # group all raw labels by their canonical key
# syn_groups = df.groupby("canon")["Original_Strain"].unique().to_dict()
# # pick a representative for each canonical key (shortest spelling)
# rep = {c: sorted(v, key=len)[0] for c, v in syn_groups.items()}

# # count datapoints per canon
# canon_counts = df["canon"].value_counts().to_dict()

# # ── 5. FLAG SPECIAL (Rule 4) ─────────────────────────────────────
# # any label containing these tokens stays separate
# RES_TOKENS = [
#     "MEF","CQ","PYR","LUM","R1","R2","RES","LEU","MDR","KO",
#     "PFS","GFP","LUC","CBG","ATTB","PFCRT"
# ]
# res_pat = re.compile("|".join(RES_TOKENS), re.I)
# special_canons = {c for c,r in rep.items() if res_pat.search(r)}

# # ── 6. INITIAL MAPPING (Synonyms merged) ─────────────────────────
# mapping = {raw: rep[c] for c, raws in syn_groups.items() for raw in raws}

# # ── 7. MINOR SUFFIX MERGING (Rule 2) ──────────────────────────────
# THRESHOLD = 0.10  # 10%
# # for each canonical group that looks like parent+one-char
# for c in list(syn_groups):
#     if c in special_canons:
#         continue
#     # skip pure synonyms (len(canon)==len(rep) automatically)
#     # find any child canon of length parent+1
#     if len(c) <= 1:
#         continue
#     parent = c[:-1]
#     if parent in syn_groups:
#         # do we have few datapoints?
#         if canon_counts.get(c, 0) < THRESHOLD * canon_counts.get(parent, 1):
#             for raw in syn_groups[c]:
#                 mapping[raw] = rep[parent]

# # ── 8. TYPOGRAPHICAL SIMILARITY (Rule 5) ───────────────────────────
# SIM_THRESHOLD = 0.95
# # consider only non-special canons
# c_non_special = [c for c in syn_groups if c not in special_canons]
# for a, b in itertools.combinations(c_non_special, 2):
#     if SequenceMatcher(None, a, b).ratio() > SIM_THRESHOLD \
#        and abs(len(a) - len(b)) <= 2:
#         # merge the rarer into the more common (or shorter rep)
#         rep_a, rep_b = rep[a], rep[b]
#         # choose which representative to keep
#         keep = rep_a if canon_counts.get(a,0) >= canon_counts.get(b,0) else rep_b
#         for raw in syn_groups[a] + syn_groups[b]:
#             mapping[raw] = keep

# # ── 9. FINAL ASSIGNMENT (Rules 6-8) ───────────────────────────────
# # any raw not in mapping keeps itself
# def standardize(raw: str) -> str:
#     return mapping.get(raw, raw)

# df["Standardized Strain"] = df["Original_Strain"].apply(standardize)

# # ── 10. SAVE OUTPUT ──────────────────────────────────────────────
# df.to_csv(OUT_CSV, index=False)
# print("✅ Wrote standardized file:", OUT_CSV)
# print("   Unique strains:", df["Standardized Strain"].nunique())

import pandas as pd
from collections import Counter

# File path to the CSV
csv_file_path = "/Users/victoriamedina/Thesis_Project/Thesis/TESTING.csv"

data = pd.read_csv(csv_file_path)

print(data.columns)

# Verify the column "stand_strain" exists
if "Standardized Strain" in data.columns:
    # Count occurrences of each item in the "stand_strain column
    strain_counts = Counter(data["Standardized Strain"].dropna())
    
    # Print the counts
    for strain, count in strain_counts.items():
        print(f"{strain}: {count}")
else:
    print("The column 'Standardized Strain' does not exist in the CSV file.")
    strain_counts = Counter()  # Initialize an empty Counter if the column doesn't exist

# Create a DataFrame from the strain counts
strain_counts_df = pd.DataFrame(strain_counts.items(), columns=["Strain", "Count"])


# Save the DataFrame to a new CSV file
strain_counts_df.to_csv(csv_file_path, index=False)

print(f"Strain counts have been saved to {csv_file_path}")
