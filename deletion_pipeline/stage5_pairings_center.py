# Run first then comment out the code below
# ------------------------------------------------------------------------------------------

# import numpy as np
# import pandas as pd
# from itertools import combinations

# # File paths
# in_path  = "./pairing_IC50_pIC50.csv"
# out_path = "./pairings_center.csv"

# # Helpers
# def parse_scores(score_str):
#     """Convert a colon-separated pIC50 string into a list of floats."""
#     return [float(x) for x in score_str.split(':') if x.strip()] if isinstance(score_str, str) else [score_str]

# def compute_rmse(scores):
#     """Compute RMSE of all pairwise differences in a list of floats."""
#     if len(scores) < 2:
#         return np.nan
#     sq_diffs = [(a - b)**2 for a, b in combinations(scores, 2)]
#     a = float(np.sqrt(np.mean(sq_diffs)))
#     return np.round(a, 6)  

# def choose_representative_pic50(scores, is_reliable):
#     """Mean if reliable, else median."""
#     if not scores:
#         return np.nan
#     a = float(np.mean(scores)) if is_reliable else float(np.median(scores))
#     return np.round(a, 6) 

# # Load the data
# df = pd.read_csv(in_path)

# # Compute new columns
# rmse_vals        = []
# reliable_flags   = []
# representatives  = []

# for txt in df['pIC50 Scores']:
#     scores = parse_scores(txt)
#     rmse    = compute_rmse(scores)
#     is_rel  = (not np.isnan(rmse)) and (rmse < 1)
    
#     rmse_vals.append(rmse)
#     reliable_flags.append(is_rel)
#     representatives.append(choose_representative_pic50(scores, is_rel))

# # Assign back into the DataFrame as floats
# df['delta_pIC50_RMSE']      = pd.Series(rmse_vals)
# df['Reliable_pIC50']        = pd.Series(reliable_flags)
# df['Representative_pIC50']  = pd.Series(representatives)


# # Save to a new CSV
# df.to_csv(out_path, index=False)
# print(f"Written processed CSV to: {out_path}")

# ------------------------------------------------------------------------------------------
# Adding run second and comment top out... this is to add pIC50 to pp csv

import pandas as pd

pp_path       = 'pp.csv'
pairings_path = 'pairings_center.csv'

iter_pp = pd.read_csv(pp_path)
pairings_df = pd.read_csv(pairings_path)

def match_pIC50(row):
    #return the pIC50 value from parings_center.csv if the strain and chemical match
    if int(row.name) % 1000 == 0:
        print(row.name)
    strain = row['Strains']
    chemical = row['Molecule_ChEMBL_ID']
    match = pairings_df[(pairings_df['Strain'] == strain) & (pairings_df['Chemical'] == chemical)]
    if not match.empty:
        return match['Representative_pIC50'].values[0]

iter_pp['Representative_pIC50'] = iter_pp.apply(match_pIC50, axis=1)

iter_pp.to_csv("test.csv", index=False)
