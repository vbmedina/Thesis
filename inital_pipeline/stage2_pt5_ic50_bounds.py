import pandas as pd
import numpy as np

PATH_IN   = "./do_not_touch/postphase2_CorrectpIC50.csv"    # raw file you just uploaded
PATH_OUT  = "./postphase2_cleaned.csv"              # cleaned-and-labelled table

df = pd.read_csv(PATH_IN)

# ────────────────────────────────────────────────────────────────
# 2.  Hard-drop IC50 values that are chemically implausible
#     (0.01 nM  ≤  IC50  ≤  100 000 nM  ⇔  11 ≥ pIC50 ≥ 4)
# ────────────────────────────────────────────────────────────────
LOWER_NM = 1e-2          # 0.01 nM  → pIC50 = 11
UPPER_NM = 1e5           # 100 000 nM = 100 µM → pIC50 = 4

df = df[df["Standard_Value"].between(LOWER_NM, UPPER_NM)].copy()

# ────────────────────────────────────────────────────────────────
# 3.  Convert the surviving IC50s (in nM) to pIC50
# ────────────────────────────────────────────────────────────────
df["pIC50"] = -np.log10(df["Standard_Value"] * 1e-9)   #  nM → M → –log10

# ────────────────────────────────────────────────────────────────
# 4.  Assign potency labels required for classification
#     inactive  : pIC50 < 6
#     active    : 6 ≤ pIC50 < 7.5
#     high-potency : pIC50 ≥ 7.5
# ────────────────────────────────────────────────────────────────
df["potency_label"] = "inactive"
df.loc[df["pIC50"] >= 6.0,  "potency_label"] = "active"
df.loc[df["pIC50"] >= 7.5,  "potency_label"] = "high_potency"

# ────────────────────────────────────────────────────────────────
# 5.  Persist the cleaned table for downstream splits
# ────────────────────────────────────────────────────────────────
df.to_csv(PATH_OUT, index=False)
print(f"Cleaned dataset saved to {PATH_OUT}  (n = {len(df):,} rows)")