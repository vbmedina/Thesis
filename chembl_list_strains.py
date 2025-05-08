import pandas as pd
import numpy as np
import re
from english_words import get_english_words_set

english_words = get_english_words_set(['web2'], lower="True")

custom_words = ["Antiplasmodial", "DNA", "PicoGreen", "SYBR-Green-1", "gametocytes", "microfluorimetry", "containing", "pfmdr1-containing", "parasitaemia", "56", "10", "5", "amplicon", "using", "giemsa",  "FM3A", "pretreated", "parasitaemia", "anti-plasmodial" "antiplasmodial", "8", "DNDI", "EC50", "2", "18", "20", "Nluc", "followed", "35S", "Cys", "N51I", "C59R", "I164L", "S108N", "DHFR", "multidrug-resistant", "cells", "drug-resistant", " G-3H", "proguanil-resistant", "cytometry", "pyrimethamine", "mefloquine-susceptible", "3", "Antimicrobial", "harboring", "CDC", "Leone", "3-200 nM", "3-200", "nM", "1", "uL", "LDH", "mefloquine-resistant", "37C", "3H-hypoxanthine", "IC50", "A-positive","3-200 nM1", "required", "50", "RBC", "parasites", "SYBR", "dose-response", "nanoGlo", "72h", "chloroquine-resistant", "microdilution", "chloroquine-sensitive", "ELISA", "vitro", "48", "incubated", "CQ-resistance", "24", "3H", "hrs", "hr", "erythrocytes", "propidium","choloroquine-sensitive", "HRP2", "72", "DAPI"]
csv_path = "/Users/victoriamedina/Thesis_Project/IC50_Plasmodium_practice.csv"

known_strains = ["W2", "CDC1", "GCO3", "RKL-9 30", "FCR-3/Gambia", "FCA 20 Ghana", "CAMP", "Cam3.II ", "MRA1239", "D-6 Sierra Leone", "Kenyan", "Tm90C2b", "FCR-3", "PA", "TM6", "FCM29", "FCR-3", "3D7r_MMV848", "Cam 3.1 R539", "IMT L1", "NF54-Mal8p1.16-GFP-Luc", "FCA10/GHA", "RCS", "wild-type TM4/8.2", "F07-3", "F09A9", "NF54", "W2mef", "D6 (CQS)Sierra Leone", "T9-96", "F07-27", "D6", "HB3", "TM90-C2A", "F09A54", "V1S", "Indochina I", "Nigerian", "subclone 3A26", "HB", "K1", "VS/1", "D10", "FCA", "3D7", "Dd2", "106/1'76I", "FCB1-R/Colombia", "7G8", "FCR3", "FCA 20/ Ghana", "F09A61", "F09N22", "F08B63", "FJB D9", "HB3-leuR1", "FCM29", "FCB1-R/Colombia", "FcB1R", "FCB", "FSH 14", "D2", "Haiti 135", "F09N18", "W-2", "F08B40", "F09A41", "F09N1", "TM91C235", "IMT 10354", "DAN", "CDC/I/HB-3", "F07-7", "FDL-HD","D-8", "IMT Guy", "IMT 16332", "IMT K2", "SB1", "Bre1", "TM90C2B", "GB4", "CsI-2", "V1/S", "GB4", "TM4", "CSL-2", "F-32 Tanzania","FCR-3", "T9/94","FcBIR", "wtTM4/8.2", "FCA 20", "FLD-B", "F32", "FcB1", "FCK2", "DD2", "F32-Tanzania", "K1CB1", "VS/1", "Indochina I", "NHP1337", "Thai", "Dd2-luc", "AZ10011008", "subclone 3C20", "3D7elo1-pfs16-CBG99", "PF18", "F09A21", "FCCI-HN", "FCC1/HN", "FVO", "Ghana", "Palo Alto/Ouganda", "IMT K14", "IMT 10500", "IMT K4", "IMT Bres", "TM1C235", "NF54-pfs16-GFP", "GCO3", "FcB1/Columbia", "W-2 Indochina III", "GCO3 SND", "3BA6 SND", "GCO3", "FCA 20/Ghana", "Cam 3.1 R539 ", "FcB1R/Columbia", "IMT K4", "NF54-pfs16-GPF", "SRIV35", "RKL09", "FCB125", "W2/Indochina", "RKL9", "SKF58", "Plasmodium yoelii N67", "Plasmodium falciparum FcM29-Cameroon", "FCM29", "VS1", "MRC-02", "C235", "3D7A", "NF 54", "FCR-3", "GHA", "KI", "PH1263-C", "F07-11", "K14", "IMT A4", "T9.94", "NF54-R", "IMT 31", "IMT 8425", "CP286", "IMT 10336", "IMT Vol", "FCR-3/A2", "C235", "FcM29", "3D7A"] 

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
