
import pandas as pd
import numpy as np

# ---------- 1. Load & clean raw data ----------
df = pd.read_csv("/Users/victoriamedina/Thesis_Project/Visualizations/chembl.csv")

ic50_dict = {}
pic50_dict = {}

for i, row in df.iterrows():
    strain = row['Standardized_Strain']
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
    
    new_data = new_data.append({
        "Strain": strain,
        "Chemical": chemical,
        "Frequency": frequency,
        "IC50 Scores": ic50_str,
        "pIC50 Scores": pic50_str
    }, ignore_index=True)

pd.save_csv(new_data, "/Users/victoriamedina/Thesis_Project/Visualizations/IC50_freq_with_pIC50.csv", index=False)

# keep only positive numeric IC50s
# df = df[pd.to_numeric(df["Standard_Value"], errors="coerce") > 0].copy()
# df["Standard_Value"] = pd.to_numeric(df["Standard_Value"])

# # if Standard_Type exists, restrict to true IC50 rows
# if "Standard_Type" in df.columns:
#     df = df[df["Standard_Type"].str.upper() == "IC50"]

# # ---------- 2. Helper functions ----------
# def join_values(values, fmt="ic50"):
#     """
#     Join an iterable of floats into a colon-separated string.
#     fmt = 'ic50'  -> raw values
#     fmt = 'pic50' -> −log10 scaled values
#     """
#     if fmt == "pic50":
#         values = [9 - np.log10(v) for v in values if v > 0]
#     return " : ".join(f"{x:.3g}".rstrip(".") for x in sorted(values))

# # ---------- 3. Build the pair-level table ----------
# summary = (
#     df.groupby(["Standardized_Strain", "Molecule_ChEMBL_ID"])
#       .agg(
#            Frequency        = ("Standard_Value", "size"),
#            **{"IC50 Scores" : ("Standard_Value", lambda v: join_values(v, "ic50")),
#               "pIC50 Scores": ("Standard_Value", lambda v: join_values(v, "pic50"))}
#       )
#       .reset_index()
#       .rename(columns={"Standardized_Strain": "Strain",
#                        "Molecule_ChEMBL_ID": "Chemical"})
# )

# # sort so the most frequently measured pairs come first
# summary = summary.sort_values(["Frequency", "Strain"],
#                               ascending=[False, True])

# # ---------- 4. Save to CSV ----------
# out_file = "strain_chemical_frequency_ic50scores_with_pIC50.csv"
# summary.to_csv(out_file, index=False)

# print(f"✓ Saved: {out_file}")
