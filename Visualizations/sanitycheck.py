import sys
import re
import numpy as np
import pandas as pd
from pathlib import Path

chembl = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv")
couples = pd.read_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/IC50_freq_with_pIC50.csv")

ic50_scores = []
for i, row in enumerate(couples['IC50 Scores']):
    scores = row.split(":")
    if(couples['Frequency'][i] != len(scores)):
        print(i)
        print(scores)
        print(couples['IC50 Scores'][i])
    for score in scores:
        ic50_scores.append(float(score))
ic50_scores = np.array(ic50_scores)

pIC50_scores = []
for row in couples['pIC50 Scores']:
    scores = row.split(":")
    for score in scores:
        pIC50_scores.append(float(score))
pIC50_scores = np.array(pIC50_scores)

# print (f"IC50 scores: {ic50_scores}")
# print (f"pIC50 scores: {pIC50_scores}")
# print(len(ic50_scores))
# print(len(pIC50_scores))

print("IC50 scores below 0: ", np.where(ic50_scores < 0)[0])
print("pIC50 scores above 15: ", np.where(pIC50_scores > 15)[0])

already_checked = {}

copy_rows = []

for i, row in chembl.iterrows():
    strain = row['Standardized_Strain']
    chemical = row['Molecule_ChEMBL_ID']
    ic50 = row['Standard_Value']
    doc = row['Document ChEMBL ID']

    key = strain + chemical + str(ic50) + doc
    if key not in already_checked:
        already_checked[key] = []

    already_checked[key].append(i)

temp = [v for key, v in already_checked.items() if len(v) > 1]
temp = np.concatenate(temp)
print(temp)
copy_rows = [[value+2, chembl.iloc[value].to_dict()] for value in temp]
print(copy_rows)
# Save the duplicate rows to a new CSV file
copy_df = pd.DataFrame(copy_rows)
copy_df.to_csv("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/duplicates.csv", index=False)

# # Print summary of duplicates
# for key, count in copy_counts.items():
#     print(f"{key}: {count} copies")