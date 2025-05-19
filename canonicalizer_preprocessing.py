import pandas as pd
import numpy as np
import re
from difflib import SequenceMatcher

counts = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/strain_counts_3.csv")     # Strain, Count
chembl = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/CHEMBL_Feb_5_Canon.csv")        # your cleaned big file

# ------------------------------------------------------------------
# 2.  Canonicalise every strain name
# ------------------------------------------------------------------
def canonical(name: str) -> str:
    """
    Upper-case, strip leading/trailing blanks, then
    delete every character that is NOT A–Z or 0–9.
    """
    name = name.strip().upper()
    return re.sub(r"[^A-Z0-9]", "", name)

counts["canon"] = counts["Strain"].apply(canonical)

# ------------------------------------------------------------------
# 3.  First-pass merge: exact canonical matches
# ------------------------------------------------------------------
syn_groups = (
    counts.groupby("canon")["Strain"]
    .apply(list)
    .to_dict()
)

# Keep a representative label for each group
# (here: the shortest original spelling, or pick the most frequent one)
rep = {
    canon: sorted(names, key=len)[0]            # customise if you prefer
    for canon, names in syn_groups.items()
}

# ------------------------------------------------------------------
# 4.  Second-pass: look for *near* matches that differ by
#     a couple of extra chars (“derivatives”) and let you decide.
# ------------------------------------------------------------------
SIM_THRESHOLD = 0.95          # SequenceMatcher ratio
EXTRA_CHARS   = 2             # how many extra alphanumeric chars
suspects = []                 # (parent, possible_derivative)

for canon1, canon2 in itertools.combinations(rep.keys(), 2):
    s1, s2 = rep[canon1], rep[canon2]
    ratio   = SequenceMatcher(None, canon1, canon2).ratio()
    length_diff = abs(len(canon1) - len(canon2))
    if ratio > SIM_THRESHOLD and length_diff <= EXTRA_CHARS:
        # Same canonical root plus ≤2 new chars (e.g. HB3 vs HB3R1)
        suspects.append((s1, s2))

suspect_df = pd.DataFrame(suspects, columns=["parent", "possible_derivative"])
suspect_df.to_csv("strain_derivative_review.csv", index=False)
print("▶ Review", len(suspect_df), "questionable pairs in strain_derivative_review.csv")

# ------------------------------------------------------------------
# 5.  Build the final mapping dict
#     (After manual review, extend/modify `rep` as you see fit)
# ------------------------------------------------------------------
mapping = {original: rep[canon] for canon, names in syn_groups.items() for original in names}

# ------------------------------------------------------------------
# 6.  Apply to the full ChEMBL dataset
# ------------------------------------------------------------------
chembl["strain_clean"] = chembl["stran_strain"].map(mapping).fillna(chembl["stran_strain"])

#                         ──────►  ready for model training
chembl.to_csv("CHEMBL_Feb_5_with_clean_strains.csv", index=False)