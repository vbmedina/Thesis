# Script to verify that every key and alias appearing in dictionaries A and B
# exists as a key in dictionary C.  Missing entries are printed.

from pprint import pprint

# ------------------------
# Dictionaries A, B, and C
# ------------------------
lst_strains = {
    "3D7":  ["NF54 3D7", "3D7 NF54", "3d7", "D7"],
    "W2": ["W-2"],
    "K1":  ["KI", "Kl"],
    "NF54": ["NF 54","NFS4","NF S4","NF-54","NF5", "NF54-Mal8p1.16-GFP-Luc"],
    "DD2": ["Dd2", "dD2", "Dd2-luc"],
    "D6": ["D-6"],
    "D10": ["D-10", "D10 CRT"],
    "HB3": ["HB", "HB-3", "CDC/I/HB-3"],
    "7G8": ["7GB"],
    "FcB1": ["FCB1", "FCB125", "Fcb1"],
    "GB4": [],
    "PC49": [],
    "TM91-C235": ["TM91C235", "TM91-C23", "TM91 c235", "TM91 C235", "C235", "TM1C235"],
    "FcR3": ["FCR-3", "FCR3", "FRC-3", "FRC3"],
    "3D7A": [],
    "F32": ["F-32", "FC32"],
    "FcB1-R": ["FcB1R", "FCB1-R", "FcBIR"],
    "TM4": ["wtTM4", "WTTM4", "WT TM4", "wt TM4"],
    "INDO": ["Indo"],
    "GHA": ["Ghana", "ghana", "(Ghana)"],
    "FcM29": ["FCM29"],
    "K1CB1": [],
    "K14": [],
    "FCB": [],
    "GCO3": ["GC03", "GCO3 C2"],
    "FCK2": [],
    "RKL9": ["RKL-9", "RKL09"],
    "CsL-2": ["CsI-2", "CSL-2"],
    "D2": ["D-2"],
    "FCA20": ["FCA 20", "FCA"],
    "T9/96": ["T9-96", "T9.96", "T9 96", "T996"],
    "T9/94": ["T9-94", "T9.94", "T9 94", "T994"],
    "V1/S": ["VS/1", "V1S", "VS1"],
    "Indochina I": [],
    "TM6": [],
    "Haiti 135": [],
    "W2mef": ["W2MEF", "W2Mef", "W2-Mef", "W2-MEF"],
    "NF54-R": [],
    "Bre1": [],
    "MRC-2": ["MRC-02", "MRC2"],
    "FSH14": ["FSH 14"],
    "FJB-D9": ["FJB D9"],
    "THAI": ["Thai", "THAI/Thailand"],
    "FcBR": ["FCBR"],
    "PA": [],
    "SRIV35": [],
    "106/1": [],
    "PH3": ["PH 3"],
    "D8": ["D-8"],
    "IMT 16332": ["IMT16332"],
    "TM267TR": [],
    "IMT A4": [],
    "IMT Bres": [],
    "MRC 20": ["MRC20"],
    "NHP1337": [],
    "IMT 9996": [],
    "IMT K4": [],
    "IMT 10336": [],
    "IMT 8425": [],
    "IMT 9881": [],
    "IMT L1": [],
    "K2": ["K-2"],
    "FDL-B": [],
    "K39": ["K-39"],
    "FLD-B": ["FLDB"],
    "FLD-NG": [],
    "MP-14": ["MP14"],
    "RKL-2": ["RKL2", "RKL02"],
    "FAC8": [],
    "FDL-NG": [],
    "FDL-HD": [],
    "RCS": [],
    "Cam 3.II": ["Cam3.II"],
    "SB1": [],
    "IMT 31": ["IMT031"],
    "Mad20": ["MAD20"],
    "ItG": [],
    "6320": [],
    "IMT Guy": [],
    "IMT 10354": ["IMT10354"],
    "PFB": [],
    "CAM": [],
    "Cam 3.1^R539T": ["Cam 3.1 R539"],
    "6218": [],
    "FcR3F86": ["FcR 3F86"],
    "F32-TEM": [],
    "F32A": [],
    "CS2": [],
    "L1": [],
    "IMT Vol": ["IMTVo1"],
    "Palo Alto": ["PaloAlto"],
    "NF54-HGL": ["NF54 HGL"],
    "FCC1": [],
    "GP1": [],
    "VNS": [],
    "CDC1": [],
    "IMT 10500": ["IMT10500"],
    "RSA11": [],
    "IMT L1": ["IMTL1"],
    "VOLL": ["Voll"],
    "BRES": ["Bres"],
    "8425": [],
    "J164": [],
    "3BA6": [],
    "TM90-C2B": ["C2B", "Tm90-C2B", "TM90 c2B", "TM90 C2B", "Tm90", "Tm90C2b", "TM90C2B"],
    "TM90-C2A": ["C2A", "Tm90-C2A", "TM90 c2A", "TM90 C2A", "Tm90C2a", "TM90C2A"],
    "BHz 26/86": ["BHz26/86"],
    "IMT K14": [],
    "TM93-C1088": [],
    "SB1-A6": [],
    "Cam 3.I": ["Cam3.1", "Cam 3.l"],
    "FVO": [],
    "CAMP": ["Camp"],
    "E8B": [],
    "FCC-1/HN": ["FCC1/HN", "FCCI-HN", "FCC1-HN", "FCCI/HN"],
    "Dd 2048": ["Dd2048", "DD2048", "DD 2048"],
    "Fab9": ["Fab 9"],
    "2087": [],
    "Cam 3.II^rev": ["Cam3.II^rev"],
    "Itg2": ["Itg-2", "Itg 2"],
    "MRA1239": [],
    "FcB2": ["FcB-2", "FCB2", "FCB-2"],
    "KT1": ["KT-1", "KT 1"],
    "FG4": ["FG 4", "FG-4"],
    "FG3": [],
    "FG1": [],
    "FG2": [],
    "FCM 17": ["FCM17"],
    "FCM6": ["FCM 6"],
    "SGE2": ["SGE 2"],
    "FcR3TC": ["FRC 3TC", "FCR3TC"],
    "D3": ["D 3", "D-3"],
    "Bres1": ["Bres 1"],
    "Pow": ["PoW", "poW"],
    "IMT 10336": ["IMT10336"],
    "FC27": ["FC 27"],
    "IMT 6311": ["IMT6311"],
    "K1Mef": ["K1MEF", "K1mef", "K1-mef"],
    "PH1263-C": [],
    "IMT K2": ["IMTK2"],
    "IMT 16116": ["IMT16116"],
    "TD7": [],
    "Nigerian": ["Nigeria"],
    "RF12": [],
    "S20": [],
    "CP286": [],
    "FCA 20/GHA": ["FCA10/GHA", "FCA 10/GHA", "FCA20/GHA"],
    "SKF58": [],
    "TM4/8.2": ["TM4-8.2", "TM 4/8.2", "wtTM4/8.2"],
    "T9/94 RC17": ["T9/94RC17"] }


lst_icm = {
    "WC4":  [],
    "XE7":  [],
    "JC3":  [],
    "JF11": [],
    "JH6":  [],
    "3C20": [],
    "XD8":  [],
    "KH7":  [],
    "JON":  [],
    "KC2":  [],
    "XB3":  [],
    "K1AM": ["K1-AM", "K1 AM"],
    "HB3-leuR1": [],
    "AUD":  [],
    "LA10": [],
    "NF10": [],
    "NIC":  [],
    "TF2":  [],
    "KA6":  [],
    "QF5":  [],
    "XG10": [],
    "DAN":  [],
    "DEV":  [],
    "JC9":  [],
    "JB8":  [],
    "KC5":  [],
    "K1HF": ["K1Hf"],
    "KB8":  [],
    "KT3":  [],
    "LC12": [],
    "Dd2-B2": ["DD2-B2"],
    "FcR3-A": ["FCR3-A"],
    "WE2":  [],
    "C2GC03": ["C2"],
    "MRA1240": ["MRA 1240"],
    "WF12": [],
    "106/1^76T": ["106/1'76T"],
    "AL2":  [],
    "34-1/E": [],
    "106/1^76I": ["106/1'76I"],
    "Smith": [],
    "JF6":  [],
    "FcR-3/A": ["FCR-3/A","FCR3/A", "FcR3/A"],
    "FcR-3/A2": ["FCR-3/A2", "FCR3/A2"],
    "GA3":  [],
    "106/1^76N": ["106/1'76N"],
    "CH12": [],
    "JB12": [],
    "AZ10011003": [],
    "K1H6/2": ["R1H6/2"],
    "F09N66": [],
    "F09N40": [],
    "F09N33": [],
    "TM3036": [],
    "F07-59": [],
    "F08B41": [],
    "AK144": [],
    "F09A55": [],
    "F07-56": [],
    "KC98": [],
    "F07-6":  [],
    "F09N68": [],
    "F07-13": [],
    "F08B53": [],
    "F07-10": [],
    "F09N35": [],
    "F09N6":  [],
    "F07-47": [],
    "F08B32": [],
    "F09N58": [],
    "AZ10011008": [],
    "F07-23": [],
    "F07-29": [],
    "SB":   [],
    "KC98 358": ["KC98358"],
    "3A7":  [],
    "F07-11": [],
    "F09A41": [],
    "F08B40": [],
    "F09N1":  [],
    "PF18":  [],
    "F08B63": [],
    "F09A61": [],
    "F09N22": [],
    "F07-7":  [],
    "F09N18": [],
    "F09A54": [],
    "F09A9":  [],
    "F07-3":  [],
    "F07-27": [],
    "F07-28": [],
    "Cam3.II^C580Y": [],
    "3A26": [],
    "F09N72": [],
    "F07-8":  [],
    "F07-35": [],
    "F07-40": [],
    "F09N44": [],
    "F09N64": [],
    "SM":   [],
    "F09A10": [],
    "Feng": [],
    "MB":   [],
    "3C18": [],
    "3A1":  [],
    "Kil-164": [],
    "F08B27": [],
    "F07-42": [],
    "AK127": [],
    "AK033 D-49": [],
    "W2AL80": [],
    "8968": [],
    "6816": [],
    "8973": [],
    "3A8":  [],
    "3A12": [],
    "3C1":  [],
    "AK158": [],
    "F07-58": [],
    "AK182": [],
    "XF12": [],
    "AK183": [],
    "AK022 D-28": [],
    "F07-31": [],
    "208432": [],
    "AK062": [],
    "F08B72": [],
    "8885": [],
    "9067": [],
    "AK033 D-0": [],
    "AZ10011017": [],
    "AK121": [],
    "AK152 D-56": [],
    "9070": [],
    "8977": [],
    "F08B9":  [],
    "AK222": [],
    "8948": [],
    "3C10": [],
    "F07-9":  [],
    "3A17": [],
    "F08B2":  [],
    "F08B61": [],
    "F09N29": [],
    "F09N78": [],
    "F07-34": [],
    "F08B7":  [],
    "F07-25": [],
    "AK150 D-42": [],
    "F08B60": [],
    "F07-50": [],
    "AK227": [],
    "AK167": [],
    "AK018": [],
    "8966": [],
    "TM90C6B": [],
    "AK150 D-0": [],
    "F09A21": []}

# -- Large dictionary C truncated for brevity in this snippet --
#   The full C is included in the notebook environment.
C={'NF54 3D7': '3D7', '3D7 NF54': '3D7', '3d7': '3D7', 'D7': '3D7', '3D7': '3D7', 'W-2': 'W2', 'W2': 'W2', 'KI': 'K1', 'Kl': 'K1', 'K1': 'K1', 'NF 54': 'NF54', 'NFS4': 'NF54', 'NF S4': 'NF54', 'NF-54': 'NF54', 'NF5': 'NF54', 'NF54-Mal8p1.16-GFP-Luc': 'NF54', 'NF54': 'NF54', 'Dd2': 'DD2', 'dD2': 'DD2', 'Dd2-luc': 'DD2', 'DD2': 'DD2', 'D-6': 'D6', 'D6': 'D6', 'D-10': 'D10', 'D10 CRT': 'D10', 'D10': 'D10', 'HB': 'HB3', 'HB-3': 'HB3', 'CDC/I/HB-3': 'HB3', 'HB3': 'HB3', '7GB': '7G8', '7G8': '7G8', 'FCB1': 'FcB1', 'FCB125': 'FcB1', 'Fcb1': 'FcB1', 'FcB1': 'FcB1', 'GB4': 'GB4', 'PC49': 'PC49', 'TM91C235': 'TM91-C235', 'TM91-C23': 'TM91-C235', 'TM91 c235': 'TM91-C235', 'TM91 C235': 'TM91-C235', 'C235': 'TM91-C235', 'TM1C235': 'TM91-C235', 'TM91-C235': 'TM91-C235', 'FCR-3': 'FcR3', 'FCR3': 'FcR3', 'FRC-3': 'FcR3', 'FRC3': 'FcR3', 'FcR3': 'FcR3', '3D7A': '3D7A', 'F-32': 'F32', 'FC32': 'F32', 'F32': 'F32', 'FcB1R': 'FcB1-R', 'FCB1-R': 'FcB1-R', 'FcBIR': 'FcB1-R', 'FcB1-R': 'FcB1-R', 'wtTM4': 'TM4', 'WTTM4': 'TM4', 'WT TM4': 'TM4', 'wt TM4': 'TM4', 'TM4': 'TM4', 'Indo': 'INDO', 'INDO': 'INDO', 'Ghana': 'GHA', 'ghana': 'GHA', '(Ghana)': 'GHA', 'GHA': 'GHA', 'FCM29': 'FcM29', 'FcM29': 'FcM29', 'K1CB1': 'K1CB1', 'K14': 'K14', 'FCB': 'FCB', 'GC03': 'GCO3', 'GCO3 C2': 'GCO3', 'GCO3': 'GCO3', 'FCK2': 'FCK2', 'RKL-9': 'RKL9', 'RKL09': 'RKL9', 'RKL9': 'RKL9', 'CsI-2': 'CsL-2', 'CSL-2': 'CsL-2', 'CsL-2': 'CsL-2', 'D-2': 'D2', 'D2': 'D2', 'FCA 20': 'FCA20', 'FCA': 'FCA20', 'FCA20': 'FCA20', 'T9-96': 'T9/96', 'T9.96': 'T9/96', 'T9 96': 'T9/96', 'T996': 'T9/96', 'T9/96': 'T9/96', 'T9-94': 'T9/94', 'T9.94': 'T9/94', 'T9 94': 'T9/94', 'T994': 'T9/94', 'T9/94': 'T9/94', 'VS/1': 'V1/S', 'V1S': 'V1/S', 'VS1': 'V1/S', 'V1/S': 'V1/S', 'Indochina I': 'Indochina I', 'TM6': 'TM6', 'Haiti 135': 'Haiti 135', 'W2MEF': 'W2mef', 'W2Mef': 'W2mef', 'W2-Mef': 'W2mef', 'W2-MEF': 'W2mef', 'W2mef': 'W2mef', 'NF54-R': 'NF54-R', 'Bre1': 'Bre1', 'MRC-02': 'MRC-2', 'MRC2': 'MRC-2', 'MRC-2': 'MRC-2', 'FSH 14': 'FSH14', 'FSH14': 'FSH14', 'FJB D9': 'FJB-D9', 'FJB-D9': 'FJB-D9', 'Thai': 'THAI', 'THAI/Thailand': 'THAI', 'THAI': 'THAI', 'FCBR': 'FcBR', 'FcBR': 'FcBR', 'PA': 'PA', 'SRIV35': 'SRIV35', '106/1': '106/1', 'PH 3': 'PH3', 'PH3': 'PH3', 'D-8': 'D8', 'D8': 'D8', 'IMT16332': 'IMT 16332', 'IMT 16332': 'IMT 16332', 'TM267TR': 'TM267TR', 'IMT A4': 'IMT A4', 'IMT Bres': 'IMT Bres', 'MRC20': 'MRC 20', 'MRC 20': 'MRC 20', 'NHP1337': 'NHP1337', 'IMT 9996': 'IMT 9996', 'IMT K4': 'IMT K4', 'IMT10336': 'IMT 10336', 'IMT 10336': 'IMT 10336', 'IMT 8425': 'IMT 8425', 'IMT 9881': 'IMT 9881', 'IMTL1': 'IMT L1', 'IMT L1': 'IMT L1', 'K-2': 'K2', 'K2': 'K2', 'FDL-B': 'FDL-B', 'K-39': 'K39', 'K39': 'K39', 'FLDB': 'FLD-B', 'FLD-B': 'FLD-B', 'FLD-NG': 'FLD-NG', 'MP14': 'MP-14', 'MP-14': 'MP-14', 'RKL2': 'RKL-2', 'RKL02': 'RKL-2', 'RKL-2': 'RKL-2', 'FAC8': 'FAC8', 'FDL-NG': 'FDL-NG', 'FDL-HD': 'FDL-HD', 'RCS': 'RCS', 'Cam3.II': 'Cam 3.II', 'Cam 3.II': 'Cam 3.II', 'SB1': 'SB1', 'IMT031': 'IMT 31', 'IMT 31': 'IMT 31', 'MAD20': 'Mad20', 'Mad20': 'Mad20', 'ItG': 'ItG', '6320': '6320', 'IMT Guy': 'IMT Guy', 'IMT10354': 'IMT 10354', 'IMT 10354': 'IMT 10354', 'PFB': 'PFB', 'CAM': 'CAM', 'Cam 3.1 R539': 'Cam 3.1^R539T', 'Cam 3.1^R539T': 'Cam 3.1^R539T', '6218': '6218', 'FcR 3F86': 'FcR3F86', 'FcR3F86': 'FcR3F86', 'F32-TEM': 'F32-TEM', 'F32A': 'F32A', 'CS2': 'CS2', 'L1': 'L1', 'IMTVo1': 'IMT Vol', 'IMT Vol': 'IMT Vol', 'PaloAlto': 'Palo Alto', 'Palo Alto': 'Palo Alto', 'NF54 HGL': 'NF54-HGL', 'NF54-HGL': 'NF54-HGL', 'FCC1': 'FCC1', 'GP1': 'GP1', 'VNS': 'VNS', 'CDC1': 'CDC1', 'IMT10500': 'IMT 10500', 'IMT 10500': 'IMT 10500', 'RSA11': 'RSA11', 'Voll': 'VOLL', 'VOLL': 'VOLL', 'Bres': 'BRES', 'BRES': 'BRES', '8425': '8425', 'J164': 'J164', '3BA6': '3BA6', 'C2B': 'TM90-C2B', 'Tm90-C2B': 'TM90-C2B', 'TM90 c2B': 'TM90-C2B', 'TM90 C2B': 'TM90-C2B', 'Tm90': 'TM90-C2B', 'Tm90C2b': 'TM90-C2B', 'TM90C2B': 'TM90-C2B', 'TM90-C2B': 'TM90-C2B', 'C2A': 'TM90-C2A', 'Tm90-C2A': 'TM90-C2A', 'TM90 c2A': 'TM90-C2A', 'TM90 C2A': 'TM90-C2A', 'Tm90C2a': 'TM90-C2A', 'TM90C2A': 'TM90-C2A', 'TM90-C2A': 'TM90-C2A', 'BHz26/86': 'BHz 26/86', 'BHz 26/86': 'BHz 26/86', 'IMT K14': 'IMT K14', 'TM93-C1088': 'TM93-C1088', 'SB1-A6': 'SB1-A6', 'Cam3.1': 'Cam 3.I', 'Cam 3.l': 'Cam 3.I', 'Cam 3.I': 'Cam 3.I', 'FVO': 'FVO', 'Camp': 'CAMP', 'CAMP': 'CAMP', 'E8B': 'E8B', 'FCC1/HN': 'FCC-1/HN', 'FCCI-HN': 'FCC-1/HN', 'FCC1-HN': 'FCC-1/HN', 'FCCI/HN': 'FCC-1/HN', 'FCC-1/HN': 'FCC-1/HN', 'Dd2048': 'Dd 2048', 'DD2048': 'Dd 2048', 'DD 2048': 'Dd 2048', 'Dd 2048': 'Dd 2048', 'Fab 9': 'Fab9', 'Fab9': 'Fab9', '2087': '2087', 'Cam3.II^rev': 'Cam 3.II^rev', 'Cam 3.II^rev': 'Cam 3.II^rev', 'Itg-2': 'Itg2', 'Itg 2': 'Itg2', 'Itg2': 'Itg2', 'MRA1239': 'MRA1239', 'FcB-2': 'FcB2', 'FCB2': 'FcB2', 'FCB-2': 'FcB2', 'FcB2': 'FcB2', 'KT-1': 'KT1', 'KT 1': 'KT1', 'KT1': 'KT1', 'FG 4': 'FG4', 'FG-4': 'FG4', 'FG4': 'FG4', 'FG3': 'FG3', 'FG1': 'FG1', 'FG2': 'FG2', 'FCM17': 'FCM 17', 'FCM 17': 'FCM 17', 'FCM 6': 'FCM6', 'FCM6': 'FCM6', 'SGE 2': 'SGE2', 'SGE2': 'SGE2', 'FRC 3TC': 'FcR3TC', 'FCR3TC': 'FcR3TC', 'FcR3TC': 'FcR3TC', 'D 3': 'D3', 'D-3': 'D3', 'D3': 'D3', 'Bres 1': 'Bres1', 'Bres1': 'Bres1', 'PoW': 'Pow', 'poW': 'Pow', 'Pow': 'Pow', 'FC 27': 'FC27', 'FC27': 'FC27', 'IMT6311': 'IMT 6311', 'IMT 6311': 'IMT 6311', 'K1MEF': 'K1Mef', 'K1mef': 'K1Mef', 'K1-mef': 'K1Mef', 'K1Mef': 'K1Mef', 'PH1263-C': 'PH1263-C', 'IMTK2': 'IMT K2', 'IMT K2': 'IMT K2', 'IMT16116': 'IMT 16116', 'IMT 16116': 'IMT 16116', 'TD7': 'TD7', 'Nigeria': 'Nigerian', 'Nigerian': 'Nigerian', 'RF12': 'RF12', 'S20': 'S20', 'CP286': 'CP286', 'FCA10/GHA': 'FCA 20/GHA', 'FCA 10/GHA': 'FCA 20/GHA', 'FCA20/GHA': 'FCA 20/GHA', 'FCA 20/GHA': 'FCA 20/GHA', 'SKF58': 'SKF58', 'TM4-8.2': 'TM4/8.2', 'TM 4/8.2': 'TM4/8.2', 'wtTM4/8.2': 'TM4/8.2', 'TM4/8.2': 'TM4/8.2', 'T9/94RC17': 'T9/94 RC17', 'T9/94 RC17': 'T9/94 RC17', 'WC4': 'WC4', 'XE7': 'XE7', 'JC3': 'JC3', 'JF11': 'JF11', 'JH6': 'JH6', '3C20': '3C20', 'XD8': 'XD8', 'KH7': 'KH7', 'JON': 'JON', 'KC2': 'KC2', 'XB3': 'XB3', 'K1-AM': 'K1AM', 'K1 AM': 'K1AM', 'K1AM': 'K1AM', 'HB3-leuR1': 'HB3-leuR1', 'AUD': 'AUD', 'LA10': 'LA10', 'NF10': 'NF10', 'NIC': 'NIC', 'TF2': 'TF2', 'KA6': 'KA6', 'QF5': 'QF5', 'XG10': 'XG10', 'DAN': 'DAN', 'DEV': 'DEV', 'JC9': 'JC9', 'JB8': 'JB8', 'KC5': 'KC5', 'K1Hf': 'K1HF', 'K1HF': 'K1HF', 'KB8': 'KB8', 'KT3': 'KT3', 'LC12': 'LC12', 'DD2-B2': 'Dd2-B2', 'Dd2-B2': 'Dd2-B2', 'FCR3-A': 'FcR3-A', 'FcR3-A': 'FcR3-A', 'WE2': 'WE2', 'C2': 'C2GC03', 'C2GC03': 'C2GC03', 'MRA 1240': 'MRA1240', 'MRA1240': 'MRA1240', 'WF12': 'WF12', "106/1'76T": '106/1^76T', '106/1^76T': '106/1^76T', 'AL2': 'AL2', '34-1/E': '34-1/E', "106/1'76I": '106/1^76I', '106/1^76I': '106/1^76I', 'Smith': 'Smith', 'JF6': 'JF6', 'FCR-3/A': 'FcR-3/A', 'FCR3/A': 'FcR-3/A', 'FcR3/A': 'FcR-3/A', 'FcR-3/A': 'FcR-3/A', 'FCR-3/A2': 'FcR-3/A2', 'FCR3/A2': 'FcR-3/A2', 'FcR-3/A2': 'FcR-3/A2', 'GA3': 'GA3', "106/1'76N": '106/1^76N', '106/1^76N': '106/1^76N', 'CH12': 'CH12', 'JB12': 'JB12', 'AZ10011003': 'AZ10011003', 'R1H6/2': 'K1H6/2', 'K1H6/2': 'K1H6/2', 'F09N66': 'F09N66', 'F09N40': 'F09N40', 'F09N33': 'F09N33', 'TM3036': 'TM3036', 'F07-59': 'F07-59', 'F08B41': 'F08B41', 'AK144': 'AK144', 'F09A55': 'F09A55', 'F07-56': 'F07-56', 'KC98': 'KC98', 'F07-6': 'F07-6', 'F09N68': 'F09N68', 'F07-13': 'F07-13', 'F08B53': 'F08B53', 'F07-10': 'F07-10', 'F09N35': 'F09N35', 'F09N6': 'F09N6', 'F07-47': 'F07-47', 'F08B32': 'F08B32', 'F09N58': 'F09N58', 'AZ10011008': 'AZ10011008', 'F07-23': 'F07-23', 'F07-29': 'F07-29', 'SB': 'SB', 'KC98358': 'KC98 358', 'KC98 358': 'KC98 358', '3A7': '3A7', 'F07-11': 'F07-11', 'F09A41': 'F09A41', 'F08B40': 'F08B40', 'F09N1': 'F09N1', 'PF18': 'PF18', 'F08B63': 'F08B63', 'F09A61': 'F09A61', 'F09N22': 'F09N22', 'F07-7': 'F07-7', 'F09N18': 'F09N18', 'F09A54': 'F09A54', 'F09A9': 'F09A9', 'F07-3': 'F07-3', 'F07-27': 'F07-27', 'F07-28': 'F07-28', 'Cam3.II^C580Y': 'Cam3.II^C580Y', '3A26': '3A26', 'F09N72': 'F09N72', 'F07-8': 'F07-8', 'F07-35': 'F07-35', 'F07-40': 'F07-40', 'F09N44': 'F09N44', 'F09N64': 'F09N64', 'SM': 'SM', 'F09A10': 'F09A10', 'Feng': 'Feng', 'MB': 'MB', '3C18': '3C18', '3A1': '3A1', 'Kil-164': 'Kil-164', 'F08B27': 'F08B27', 'F07-42': 'F07-42', 'AK127': 'AK127', 'AK033 D-49': 'AK033 D-49', 'W2AL80': 'W2AL80', '8968': '8968', '6816': '6816', '8973': '8973', '3A8': '3A8', '3A12': '3A12', '3C1': '3C1', 'AK158': 'AK158', 'F07-58': 'F07-58', 'AK182': 'AK182', 'XF12': 'XF12', 'AK183': 'AK183', 'AK022 D-28': 'AK022 D-28', 'F07-31': 'F07-31', '208432': '208432', 'AK062': 'AK062', 'F08B72': 'F08B72', '8885': '8885', '9067': '9067', 'AK033 D-0': 'AK033 D-0', 'AZ10011017': 'AZ10011017', 'AK121': 'AK121', 'AK152 D-56': 'AK152 D-56', '9070': '9070', '8977': '8977', 'F08B9': 'F08B9', 'AK222': 'AK222', '8948': '8948', '3C10': '3C10', 'F07-9': 'F07-9', '3A17': '3A17', 'F08B2': 'F08B2', 'F08B61': 'F08B61', 'F09N29': 'F09N29', 'F09N78': 'F09N78', 'F07-34': 'F07-34', 'F08B7': 'F08B7', 'F07-25': 'F07-25', 'AK150 D-42': 'AK150 D-42', 'F08B60': 'F08B60', 'F07-50': 'F07-50', 'AK227': 'AK227', 'AK167': 'AK167', 'AK018': 'AK018', '8966': '8966', 'TM90C6B': 'TM90C6B', 'AK150 D-0': 'AK150 D-0', 'F09A21': 'F09A21'}

# ------------------------
# Helper function
# ------------------------
def flatten_alias_dict(d):
    """
    Return a set containing the dictionary's keys and all strings in the value lists.
    """
    flattened = set(d.keys())
    for aliases in d.values():
        flattened.update(aliases)
    return flattened

flat_A = flatten_alias_dict(A)
flat_B = flatten_alias_dict(B)
flat_all = flat_A | flat_B

missing = sorted(s for s in flat_all if s not in C)

print("Summary")
print("-------")
print(f"Unique strings in A : {len(flat_A):>5}")
print(f"Unique strings in B : {len(flat_B):>5}")
print(f"Total unique strings: {len(flat_all):>5}")
print(f"Strings missing in C: {len(missing):>5}")
print()

print("List of missing strings")
print("-----------------------")
pprint(missing)


# x = '''K1
# FJB D9
# NF54
# D6
# W2
# Dd2
# FcB1R
# FSH 14
# FCB
# D2
# HB3
# Haiti 135
# W-2
# F09A41
# F08B40
# F09N1
# 3D7
# TM91C235
# D-8
# IMT Guy
# FcB1
# IMT 16332
# IMT K2
# SB1
# FCB1
# 7G8
# DD2
# FCR3
# Bre1
# D10
# TM90C2B
# GB4
# F32
# K1CB1
# CsI-2
# V1/S
# F-32
# wtTM4/8.2
# FcBIR
# FCA 20
# FCK2
# VS/1
# Indochina I
# 3D7A
# NHP1337
# Thai
# Dd2-luc
# 3C20
# PF18
# F09A21
# FCC1/HN
# Ghana
# IMT K14
# IMT 10500
# IMT K4
# IMT Bres
# TM1C235
# GCO3
# 3BA6
# Nigerian
# 106/1'76I
# FcM29
# FCM29
# VS1
# MRC-02
# C235
# (Ghana)
# NF 54
# GHA
# FCR-3
# K14
# KI
# IMT A4
# IMT 31
# IMT 8425
# IMT 10336
# IMT Vol
# FCR-3/A2
# THAI/Thailand
# FCB1-R
# T9/94
# F08B63
# F09A61
# F09N22
# FLD-B
# RKL9
# HB
# FCA
# T9-96
# E8B
# 3A26
# TM90-C2A
# TM90C2A
# V1S
# IMT 9881
# FLD-NG
# CDC1
# BHz26/86
# TM4/8.2
# Tm90C2b
# 3D7 K1
# FDL-HD
# KA6
# CDC/I/HB-3
# FCM17
# F07-7
# F09N18
# F09A54
# F09A9
# F07-3
# F07-27
# IMT 10354
# TM4
# Nigeria
# MRC2
# QF5
# DAN
# MP14
# TM6
# RKL-9
# FDL-B
# FCR3/A2
# IMT031
# Palo Alto
# CSL-2
# MRA1239
# IMT L1
# FCC1
# FcBR
# W2mef
# FCBR
# Indo
# VNS
# FAC8
# KT1
# dD2
# SKF58
# CP286
# PH1263-C
# NF54-R
# FCB125
# T9.94
# F07-11
# Mad20
# TM267TR
# 3D7 NF54
# INDO
# HB3-leuR1
# F09N6
# F07-47
# RSA11
# FCA10/GHA
# T9/94 RC17
# NF54-Mal8p1.16-GFP-Luc
# Cam 3.1 R539
# C2A
# FVO
# FCCI-HN
# PA
# THAI
# RKL-2
# Cam3.II
# W2AL80
# 3A17
# F08B32
# F09N58
# D-6
# D-10
# NF54 HGL
# K2
# Cam3.II^C580Y
# ItG
# RKL09
# 2087
# FG3
# AZ10011008
# SRIV35
# SGE2
# F32-TEM
# MRC-2
# PFB
# FRC-3
# Tm90-C2B
# 106/1
# GP1
# F07-23
# F07-29
# L1
# FG1
# PH3
# SB
# KC98358
# FC32
# FG2
# MRC 20
# FC27
# 6320
# BHz 26/86
# GC03
# F32A
# W2-Mef
# D10 CRT
# 3A7
# GA3
# CAMP
# AZ10011003
# AZ10011017
# Bres
# SB1-A6
# Cam3.II^rev
# 6816
# 8968
# 8973
# FDL-NG
# D10 7G8
# J164
# K39
# 208432
# 106/1'76N
# CH12
# TD7
# Kil-164
# 7GB
# KT3
# Dd2 CAM
# D-3
# Dd2-B2
# FCM6
# FCB-2
# PaloAlto
# RF12
# FG4
# JC9
# T996
# 34-1/E
# 3A8
# JB12
# XG10
# WE2
# TM90-C2B
# FcR3
# IMT10500
# Dd2048
# C2B
# 3A12
# 3C1
# 106/1'76T
# FCR3-A
# W2MEF
# Voll
# S20
# F07-42
# Tm90
# NFS4
# 8425
# 6218
# K1AM
# AK127
# AK033 D-49
# TM90C6B
# F08B27
# 3A1
# 3C18
# CS2
# Itg2
# FcB2
# TM93-C1088
# F09N72
# F07-8
# F07-35
# F07-40
# F09N44
# AK158
# Cam3.1
# F07-58
# Smith
# IMT 16116
# F09N64
# SM
# AK182
# PoW
# Fab9
# F09A10
# Feng
# R1H6/2
# XF12
# FCB2
# IMT 9996
# LC12
# MB
# JH6
# JON
# KH7
# JF11
# JC3
# AK183
# AK022 D-28
# K1Hf
# XE7
# F07-28
# WF12
# KC5
# FCR3TC
# F07-31
# PC49
# XD8
# D7
# AK062
# GCO3 C2
# F09N40
# F08B41
# IMT 6311
# C2
# F08B72
# AL2
# F07-59
# TM91-C235
# 8885
# F09N66
# IMTL1
# 9067
# AK033 D-0
# F09N33
# MRA1240
# DEV
# AK121
# AK152 D-56
# 9070
# 8977
# F09N68
# NIC
# NF10
# JF6
# LA10
# XB3
# 3C10
# F07-10
# AUD
# JB8
# 3d7
# KC2
# KB8
# WC4
# F08B9
# K1mef
# RKL2
# D3
# AK222
# 8948
# KC98
# TM91-C23
# TM3036
# poW
# F07-56
# F07-6
# F09N35
# F07-9
# F08B53
# F07-13
# F08B2
# F08B61
# F09N29
# F09N78
# F07-34
# F08B7
# F07-25
# AK150 D-42
# F08B60
# F07-50
# AK227
# AK167
# AK018
# 8966
# Bres1
# F09A55
# AK144
# AK150 D-0
# NF-54
# Fcb1'''

# list_of_strings = x.split('\n')
# # Remove leading and trailing whitespace from each string
# list_of_strings = [s.strip() for s in list_of_strings]
# # Print the duplicate strings
# duplicates = set([s for s in list_of_strings if list_of_strings.count(s) > 1])
# print("Duplicate strings:")
# for dup in duplicates:
#     print(dup)