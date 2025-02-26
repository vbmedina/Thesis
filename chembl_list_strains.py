import pandas as pd
import numpy as np
import re
from english_words import get_english_words_set

english_words = get_english_words_set(['web2'], lower="True")

custom_words = ["Antiplasmodial", "multidrug-resistant", "cells", "drug-resistant", "mefloquine-susceptible", "3", "Antimicrobial", "harboring", "CDC", "Leone", "3-200 nM", "3-200", "nM", "1", "uL", "LDH", "mefloquine-resistant", "37C", "3H-hypoxanthine", "IC50", "A-positive","3-200 nM1", "required", "50", "RBC", "parasites", "SYBR", "dose-response", "nanoGlo", "72h", "chloroquine-resistant", "microdilution", "chloroquine-sensitive", "ELISA", "vitro", "48", "incubated", "CQ-resistance", "24", "3H", "hrs", "hr", "erythrocytes", "choloroquine-sensitive", "HRP2", "72", "DAPI"]

csv_path = "/Users/victoriamedina/Thesis_Project/IC50_Plasmodium_practice.csv"

known_strains = ["W2", "NF54", "D6", "HB3", "K1", "D10", "3D7", "Dd2", "7G8", "FCR3", "FJB D9", "FcB1R", "FCB", "FSH 14", "D2", "Haiti 135", "W-2", "F08B40", "F09A41", "F09N1", "TM91C235", "D-8", "IMT Guy", "IMT 16332", "IMT K2", "SB1", "Bre1", "TM90C2B", "GB4", "CsI-2", "V1/S", "GB4", "F-32 Tanzania", "FcBIR", "wtTM4/8.2", "FCA 20/Ghana", "F32", "FcB1", "FCK2", "DD2", "F32-Tanzania", "K1CB1", "FCA 20/ Ghana", "VS/1", "Indochina I", "NHP1337", "Thai", "Dd2-luc", "subclone 3C20", "PF18", "F09A21", "FCC1/HN", "Ghana", "IMT K14", "IMT 10500", "IMT K4", "IMT Bres", "TM1C235", "NF54-pfs16-GFP", "GCO3", "FcB1/Columbia", "W-2 Indochina III", "GCO3 SND", "3BA6 SND", "GCO3", "FcB1R/Columbia", "IMT K4", "W2/Indochina", "FCA 20/ Ghana", " Plasmodium yoelii N67", "Plasmodium falciparum FcM29-Cameroon", "FCM29", "VS1", "MRC-02", "C235", "3D7A", "NF 54", "FCR-3", "GHA", "KI", "K14", "IMT A4", "IMT 31", "IMT 8425", "IMT 10336", "IMT Vol", "FCR-3/A2", "C235", "FcM29", "3D7A"] 

def extract_terms(description):
    if not isinstance(description, str):
        return description
    
    tokens = re.findall(r"\b\w[\w-]*\b", description)

    filtered_words = [word for word in tokens if word.lower() not in english_words and word not in custom_words]

    found_strains = [word for word in filtered_words if word in known_strains]
    
    # If any known strain is found, return only those strains
    if found_strains:
        return " ".join(found_strains)
    else:
        return " ".join(filtered_words)

df = pd.read_csv(csv_path, low_memory=False)

df["cleaned_description"] = df["Assay Description"].apply(extract_terms)

print(df["cleaned_description"])

# assay = [col for col in df.columns if "Assay Description" in col]

# if assay:
#     df['Assay'] = df[assay[0]].apply(extract_terms)
# else:
#     df['Assay'] = "ERROR"

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_path, index=False)

print(f"CSV file updated and saved successfully: {csv_path}")
