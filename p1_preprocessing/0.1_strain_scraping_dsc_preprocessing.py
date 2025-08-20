# Description: Need to finish description... udpate paths


# Imports
import pandas as pd
import numpy as np

# Path
df = pd.read_csv("./potency.csv", low_memory=False)

# Dictionary of strains
lst_strains = ["3D7","NF54 3D7","3D7 NF54", "IPC4912","3d7","D7","W2","W-2","K1","KI","Kl","NF54","NF 54","NFS4","NF S4","NF-54","NF5","NF54-Mal8p1.16-GFP-Luc","Dd2","DD2","dD2","Dd2-luc","D6","D-6","D10","D-10","D10 CRT","HB3","HB","HB-3","CDC/I/HB-3","7G8","7GB","FcB1","FCB1","FCB125","Fcb1","GB4","PC49","TM91-C235", "FCR3F86", "TM91C235","TM91-C23","TM91 c235","3D7elo1-pfs16-CBG99", "TM90C2b", "TM91 C235","C235","TM1C235","FcR3","FCR-3","FCR3","FRC-3","FRC3","3D7A","F32","F-32","FcB1-R","FcB1R","FCB1-R","FcBIR","TM4","wtTM4","WTTM4","WT TM4","wt TM4","INDO","Indo","GHA","Ghana", "ghana","(Ghana)","FcM29","FCM29","K1CB1","K14","FCB","GCO3","GC03","GCO3 C2","FCK2","RKL9","RKL-9","RKL09","CsL-2","CsI-2","CSL-2","D2","D-2","FCA20","FCA 20","FCA","T9/96","T9-96","T9.96","T9 96","T996","T9/94","T9-94","T9.94","T9 94","T994","V1/S","VS/1","V1S","VS1","Indochina I","TM6","Haiti 135","W2mef","W2MEF","W2Mef","W2-Mef","W2-MEF","NF54-R","Bre1","IMT 10354","IMT10354","MRC-2","MRC-02","MRC2","FSH14","FSH 14","FJB-D9","FJB D9","THAI","Thai","THAI/Thailand", "ITG2", "FcBR","FCBR","PA","SRIV35","106/1","PH3","PH 3","D8","D-8","IMT 16332","IMT16332","TM267TR","IMT A4","FC27","IMT Bres","IMT 10500","MRC 20","MRC20","NHP1337","IMT 9996","IMT 31","IMT K4","IMT 10336","IMT 8425","IMT 9881","IMT L1","K2","K-2","FDL-B","K39","K-39","FLD-B","FLDB","FLD-NG","C2A","MP-14","MP14","RKL-2","RKL2","RKL02","FAC8","FDL-NG","FDL-NG","FDL-HD","FCS","Cam 3.II","Cam3.II","SB1","IMT031","IMT 31","Mad20","MAD20","ItG","6320","IMT Guy","IMT 10354","IMT10354","PFB","CAM","Cam 3.1^R539T","Cam 3.1 R539","6218","FcR3F86","FcR 3F86","F32-TEM","F32A","CS2","L1","IMT Vol","IMT Vo1","Palo Alto","PaloAlto","NF-54HGL","NF54 HGL","FCC1","GP1","VNS","CDC1","IMT 10500","IMT10500","RSA11","IMT L1","IMTL1","VOLL","Voll","BRES","Bres","8425","SGE2","J164","3BA6","TM90-C2B","C2B","Tm90-C2B","TM90 c2B","TM90 C2B","Tm90","Tm90C2b","TM90C2B","BHz 26/86","BHz26/86","IMT K14","TM93-C1088","SB1-A6","Cam 3.I","Cam3.1","Cam 3.l","FVO","CAMP","Camp","E8B","FCC-1/HN","FCC1/HN","FCCI-HN","FCC1-HN","FCC1-HN","FCCI/HN","Dd 2048","Dd2048", "DD2048","DD 2048","Fab9","Fab 9","2087","Cam 3.II^rev","Cam3.II^rev","Itg2","Itg-2","Itg 2","MRA1239","FcB2","FcB-2","FCB2","FCB-2","KT1","KT-1","KT 1","FG4","FG 4","FG-4","FG3","FG1","FG2","FCM 17","FCM17","FCM6","FCM 6","SGE2","SGE 2","FcR3TC","FRC 3TC","FCR3TC","D3","D 3","D-3","Bres1","Bres 1","Pow","PoW","poW","IMT 10336","IMT10336","FC27","FC 27","IMT 6311","IMT6311","K1Mef","K1MEF","K1mef","K1-mef","F32","FC32","PH1263-C","IMT K2","IMTK2","KT1","KT3","IMT 16116","IMT16116","TD7","Nigerian","Nigeria","RF12","RKL9","S20","CP286","FCA20/GHA","FCA10/GHA","FCA 10/GHA","FCA 20/GHA","SKF58","NF54-HGL","NF54 HGL","TM4/8.2","TM4-8.2","TM 4/8.2","wtTM4/8.2","T9/94 RC17","T9/94RC17","TM90-C2A","C2A","Tm90-C2A","TM90 c2A","TM90 C2A","Tm90C2a","TM90C2A"] 

# Dictionary of isolates
lst_ic = ["WC4","XE7", "3A3", "F07-46", "AK022D-0", "JC3","JF11","JH6","3C20","XD8","KH7","JON","KC2","F09A21","XB3","K1AM","K1-AM","K1 AM", "3A3", "HB3-leuR1","AUD","LA10","NF10","NIC","TF1","KA6","QF5","XG10","DAN","DEV","JC9","JB8","KC5","K1HF","K1Hf","KB8","KT3","LC12","Dd2-B2","DD2-B2","FcR3-A","FCR3-A","WE2","C2GC03","C2","MRA1240","MRA 1240","WF12","106/1^76T","106/1'76T","AL2","34-1/E","106/1^76I","106/1'76I","Smith","JF6","FcR-3/A","FCR-3/A","FCR-3/A","FCR3/A","FcR3/A","FcR-3/A2","FCR-3/A2","FCR3/A2","GA3","106/1^76N","106/1'76N","CH12","JB12","AZ10011003","K1H6/2","R1H6/2","F09N66","F09N40","F09N33","TM3036","F07-59","F08B41","AK144","F09A55","F07-56","KC98","F07-6","F09N68","F07-13","F08B53","F07-10","F09N35","F09N6","F07-47","F08B32","F09N58","AZ10011008","F07-23","F07-29","SB","KC98 358","KC98358","3A7","F07-11","F09A41","F08B40","F09N1","PF18","F08B63","F09A61","F09N22","F07-7","F09N18","F09A54","F09A9","F07-3","F07-27","F07-28","Cam3.II^C580Y","3A26","F09N72","F07-8","F07-35","F07-40","F09N44","F09N64","SM","F09A10","Feng","MB","3C18","3A1","Kil-164","F08B27","F07-42","AK127","AK033 D-49","NHP1337","W2AL80","8968","6816","8973","208432","3A8","3A12","3C1","AK158","F07-58","AK182","XF12","AK183","AK022 D-28","F07-31","208432","3A12","AK062","F08B72","8885","9067","AK033 D-0","AZ10011017","AK121","AK152 D-56","9070","8977","F08B9","AK222","AK249", "8948","3C10","F07-9","3A17","F08B2","F08B61","F09N29","F09N78","F07-34","F08B7","F07-25","AK150 D-42","F08B60","F07-50","AK227","AK167","AK018","8966","TM90C6B","AK150 D-0","3C20","F09A21", "AZ10011022","CDA-1553"]

# Extracting terms with nomenclature rules
def extract_terms(description):

    if not isinstance(description, str):
        return "ERROR"
    
    description = " " + description + " "

    # Rules for extracting strains and isolates
    strains = [strain for strain in lst_strains if (" " + strain + " " in description) or (" " + strain + "-" in description) or (" " + strain + "/" in description) or (" " + strain + ":" in description) or (" " + strain + ";" in description) or (" " + strain + "." in description) or (" " + strain + "(" in description)]
    ics = [ic for ic in lst_ic if (" " + ic + " " in description) or (" " + ic + "-" in description) or (" " + ic + "/" in description) or (" " + ic + ":" in description) or (" " + ic + ";" in description) or (" " + ic + "." in description) or (" " + ic + "(" in description)]

    strains = np.array(strains)
    ics = np.array(ics)

    compounds = np.concatenate((strains, ics))

    unique_compounds = []
    for compound1 in compounds:
        #Find smallest compound
        smallest = compound1
        for compound2 in compounds:
            if compound2 in compound1:
                if len(compound2) < len(smallest):
                    smallest = compound2

        #Find largest compound
        similar = []
        for compound2 in compounds:
            if smallest in compound2:
                similar.append(compound2)

        largest = max(similar, key=len)

        if largest not in unique_compounds:
            unique_compounds.append(largest)


# Specific rules for strain names
    if "NF54 3D7" in unique_compounds:
        unique_compounds.remove("NF54 3D7")
        unique_compounds.append("3D7")

    if "3D7 NF54" in unique_compounds:
        unique_compounds.remove("3D7 NF54")
        unique_compounds.append("3D7")
        
    if len(unique_compounds) > 1:
        if "Ghana" in unique_compounds:
            unique_compounds.remove("Ghana")
        if "Thai" in unique_compounds:
            unique_compounds.remove("Thai")
        if "THAI" in unique_compounds:
            unique_compounds.remove("THAI")
        if "Thailand" in unique_compounds:
            unique_compounds.remove("Thailand")
        if "Nigerian" in unique_compounds:
            unique_compounds.remove("Nigerian")
        if "Nigeria" in unique_compounds:
            unique_compounds.remove("Nigeria")

    unique_compounds = np.array(unique_compounds)
    
    return " ".join(unique_compounds)

# Create new column with extracted terms
df["Strains"] = df["Assay Description"].apply(extract_terms)

# Save path
save_path = "./potency_str.csv"
df.to_csv(save_path, index=False)
