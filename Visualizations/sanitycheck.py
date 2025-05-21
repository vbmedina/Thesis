import sys
import re
import numpy as np
import pandas as pd
from pathlib import Path

chembl = pd.read_csv("/Users/victoriamedina/Thesis_Project/Visualizations/chembl.csv")
couples = pd.read_csv("/Users/victoriamedina/Thesis_Project/Visualizations/IC50_freq_with_pIC50.csv")

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