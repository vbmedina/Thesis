
import pandas as pd
import numpy as np

# ---------- 1. Load & clean raw data ----------
df = pd.read_csv("./Do Not Touch/postphase5.csv")

ic50_dict = {}
pic50_dict = {}

for i, row in df.iterrows():
    strain = row['stand_strain_bin']
    chemical = row['Molecule_ChEMBL_ID']
    ic50 = row['Standard_Value']
    pic50 = row['pChEMBL Value']
    key = strain + "_" + chemical
    
    if key not in ic50_dict:
        ic50_dict[key] = []
        pic50_dict[key] = []

    ic50_dict[key].append(ic50)
    pic50_dict[key].append(pic50)

new_data = pd.DataFrame(columns=["Strain", "Chemical", "Frequency", "IC50 Scores", "pIC50 Scores"])

for key in ic50_dict:
    strain, chemical = key.split("_")
    frequency = len(ic50_dict[key])
    ic50_scores = ic50_dict[key]
    pic50_scores = pic50_dict[key]
    
    ic50_str = " : ".join([str(x) for x in ic50_scores])
    pic50_str = " : ".join([str(x) for x in pic50_scores])
    
    new_data.loc[len(new_data)] = {
    "Strain": strain,
    "Chemical": chemical,
    "Frequency": frequency,
    "IC50 Scores": ic50_str,
    "pIC50 Scores": pic50_str
    }

new_data.to_csv("./Visualizations/pairings_strains_IC50_scaff.csv", index=False)