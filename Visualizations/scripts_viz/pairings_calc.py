#!/usr/bin/env python3

import numpy as np
import pandas as pd
from itertools import combinations

# 1) File paths
in_path  = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/pairings_strains_IC50.csv"
out_path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/pairings_strains_IC50_edits.csv"

# 2) Helpers
def parse_scores(score_str):
    """Convert a colon-separated pIC50 string into a list of floats."""
    return [float(x) for x in score_str.split(':') if x.strip()]

def compute_rmse(scores):
    """Compute RMSE of all pairwise differences in a list of floats."""
    if len(scores) < 2:
        return np.nan
    sq_diffs = [(a - b)**2 for a, b in combinations(scores, 2)]
    a = float(np.sqrt(np.mean(sq_diffs)))
    return np.round(a, 6)  

def choose_representative_pic50(scores, is_reliable):
    """Mean if reliable, else median."""
    if not scores:
        return np.nan
    a = float(np.mean(scores)) if is_reliable else float(np.median(scores))
    return np.round(a, 6) 

# 3) Load the data
df = pd.read_csv(in_path)

# 4) Compute new columns
rmse_vals        = []
reliable_flags   = []
representatives  = []

for txt in df['pIC50 Scores']:
    scores = parse_scores(txt)
    rmse    = compute_rmse(scores)
    is_rel  = (not np.isnan(rmse)) and (rmse < 1)
    
    rmse_vals.append(rmse)
    reliable_flags.append(is_rel)
    representatives.append(choose_representative_pic50(scores, is_rel))

# 5) Assign back into the DataFrame as floats
df['delta_pIC50_RMSE']      = pd.Series(rmse_vals)
df['Reliable_pIC50']        = pd.Series(reliable_flags)
df['Representative_pIC50']  = pd.Series(representatives)


# 6) Save to a new CSV
df.to_csv(out_path, index=False)
print(f"Written processed CSV to: {out_path}")
