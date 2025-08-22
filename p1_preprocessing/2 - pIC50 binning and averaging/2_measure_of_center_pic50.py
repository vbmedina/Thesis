''' Description: This script is split into 2 parts. It aims to calculate the best measure of center for the pIC50 scores of 
each molecule-strain pairing. Using the CSV file generated in the previous step, "pairing_IC50_pIC50.csv", it computes the 
RMSE of the pIC50 scores for each pairing. In this case, RMSE (Root Mean Square Error) quantifies the typical deviation of 
pIC50 values from their mean and serves as an indicator of data variability and outlier presence. Since the mean is highly 
sensitive to outliers and can be skewed by extreme values, this script uses RMSE as a diagnostic tool: when RMSE > 1 
(indicating high variability likely due to outliers), it resorts to the median as a robust measure of center that is resistant 
to extreme values. When RMSE <= 1 (indicating low variability with minimal outliers), it uses the mean as the most precise 
estimate of central tendency using all available data points. 

Preconditions:
1) "pairing_IC50_pIC50.csv" - generated from stage 2 of the pipeline.
2) "pairings_center.csv" - gnerated from stage 2 of the pipeline.
'''



import numpy as np
import pandas as pd
from itertools import combinations

'''---------------------------------------------------------------------------------------------------------------------'''
''' Part 1: This script is the 1st part of the script, which goes through the "pairing_IC50_pIC50.csv" file and 
computes a representative pIC50 for each strain-chemical pair based on the RMSE of their pIC50 scores. If the
RMSE is less than 1, the mean of the pIC50 scores is used; otherwise, the median is used.'''

# File paths
in_path  = "./p1_preprocessing/2 - pIC50 binning and averaging/pairing_IC50_pIC50.csv"
out_path = "./p1_preprocessing/2 - pIC50 binning and averaging/pairings_center.csv"

# Helpers
def parse_scores(score_str):
    # Convert a colon-separated pIC50 string into a list of floats
    return [float(x) for x in score_str.split(':') if x.strip()] if isinstance(score_str, str) else [score_str]

def compute_rmse(scores):
    # Compute RMSE of all pairwise differences in a list of floats
    if len(scores) < 2:
        return np.nan
    sq_diffs = [(a - b)**2 for a, b in combinations(scores, 2)]
    a = float(np.sqrt(np.mean(sq_diffs)))
    return np.round(a, 6)  

def choose_representative_pic50(scores, is_reliable):
    # Mean if reliable, else median
    if not scores:
        return np.nan
    a = float(np.mean(scores)) if is_reliable else float(np.median(scores))
    return np.round(a, 6) 

# Load the data
df = pd.read_csv(in_path)

# Compute new columns
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

# Assign back into the DataFrame as floats
df['delta_pIC50_RMSE']      = pd.Series(rmse_vals)
df['Reliable_pIC50']        = pd.Series(reliable_flags)
df['Representative_pIC50']  = pd.Series(representatives)

# Save to a new CSV
df.to_csv(out_path, index=False)
print(f"Written processed CSV to: {out_path}")

'''---------------------------------------------------------------------------------------------------------------------'''
''' Part 2: This 2nd part of the script saves the "Representative_pIC50" result to final_data_copy.csv. '''

final_data = './p1_preprocessing/2 - pIC50 binning and averaging/final_data_copy.csv'

iter_fd = pd.read_csv(final_data)
pairings_df = pd.read_csv(out_path)

def match_pIC50(row):
    #return the pIC50 value from parings_center.csv if the strain and chemical match
    if int(row.name) % 1000 == 0:
        print(row.name)
    strain = row['standardized_strain']
    chemical = row['Molecule_ChEMBL_ID']
    match = pairings_df[(pairings_df['Strain'] == strain) & (pairings_df['Molecule'] == chemical)]
    if not match.empty:
        return match['Representative_pIC50'].values[0]

iter_fd['Representative_pIC50'] = iter_fd.apply(match_pIC50, axis=1)

iter_fd.to_csv("./p1_preprocessing/2 - pIC50 binning and averaging/final_data_copy.csv", index=False)
