import pandas as pd
import numpy as np

csv_path = "/Users/victoriamedina/Thesis_Project/IC50_Plasmodium.csv"

known_strains = ["W2", "NF54", "D6", "HB3", "K1", "D10", "3D7", "Dd2", "7G8", "FCR3", "FJB D9", "FcB1R", "FCB", "FSH 14", "D2", "Haiti 135", "W-2", "F08B40", "F09A41", "F09N1", "TM91C235", "D-8", "IMT Guy", "IMT 16332", "IMT K2", "SB1", "Bre1", "TM90C2B", "GB4", "CsI-2", "V1/S", "GB4", "F-32 Tanzania", "FcBIR", "wtTM4/8.2", "FCA 20/Ghana", "F32", "FcB1", "FCK2", "DD2", "F32-Tanzania", "K1CB1", "FCA 20/ Ghana", "VS/1", "Indochina I", "NHP1337", "Thai", "Dd2-luc", "subclone 3C20", "PF18", "F09A21", "FCC1/HN", "Ghana", "IMT K14", "IMT 10500", "IMT K4", "IMT Bres", "TM1C235", "NF54-pfs16-GFP", "GCO3", "FcB1/Columbia", "W-2 Indochina III", "GCO3 SND", "3BA6 SND", "GCO3", "FcB1R/Columbia", "IMT K4", "W2/Indochina", "FCA 20/ Ghana", " Plasmodium yoelii N67", "Plasmodium falciparum FcM29-Cameroon", "FCM29", "VS1", "MRC-02", "C235", "3D7A", "NF 54", "FCR-3", "GHA", "KI", "K14", "IMT A4", "IMT 31", "IMT 8425", "IMT 10336", "IMT Vol", "FCR-3/A2", "C235", "FcM29", "3D7A"] 

def find_strains(description, known_strains):
    if pd.isna(description):
        return "ERROR"
    found = [strain for strain in known_strains if strain in description]
    found = np.array(found)
    found = found[np.argsort([len(x) for x in found])]

    return found[-1] if len(found)>0 else "ERROR"

df = pd.read_csv(csv_path, low_memory=False)

assay = [col for col in df.columns if "Assay Description" in col]

if assay:
    df['Assay'] = df[assay[0]].apply(lambda x: find_strains(x, known_strains))
else:
    df['Assay'] = "ERROR"

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_path, index=False)

print(f"CSV file updated and saved successfully: {csv_path}")
