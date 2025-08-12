# make_violin_max_tanimoto.py
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# ---- paths / splits ----
BASE = Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/splits")
SPLITS = ["random", "scaffold", "butina", "umap"] 

# ---- RDKit helpers ----
def morgan_fp(smi, r=2, nBits=2048):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, r, nBits)

def max_sims_for_fold(train_csv: Path, test_csv: Path):
    tr = pd.read_csv(train_csv, usecols=["Smiles"])
    te = pd.read_csv(test_csv,  usecols=["Smiles"])
    tr_fps = [morgan_fp(s) for s in tr["Smiles"]]
    out = []
    for smi in te["Smiles"]:
        sims = DataStructs.BulkTanimotoSimilarity(morgan_fp(smi), tr_fps)
        out.append(max(sims))
    return out

def collect_split_max_sims(split: str):
    folder = BASE / split
    if not folder.exists():
        return []
    # find all train files; extract fold number; pair with the matching test file
    train_files = sorted(folder.glob(f"{split}_fold_*_train.csv"))
    results = []
    for train_csv in train_files:
        m = re.search(r"fold[_]?(\d+)_train\.csv", train_csv.name)
        if not m: 
            continue
        fnum = m.group(1)
        # try both test filename patterns
        candidates = [
            folder / f"{split}_fold{fnum}_test.csv",
            folder / f"{split}_fold_{fnum}_test.csv",
        ]
        test_csv = next((p for p in candidates if p.exists()), None)
        if test_csv is None:
            continue
        results.extend(max_sims_for_fold(train_csv, test_csv))
    return results

# ---- collect all max(test->train) sims ----
rows = []
for split in SPLITS:
    vals = collect_split_max_sims(split)
    rows += [{"split": split, "max_sim": v} for v in vals]

sim_df = pd.DataFrame(rows)
print(sim_df.groupby("split")["max_sim"].mean().round(3))  # quick sanity check

# ---- violin + box plot (all reds) ----
order = [s for s in SPLITS if s in sim_df["split"].unique()]
reds = sns.color_palette("Reds", n_colors=len(order))
palette = dict(zip(order, reds))

plt.figure(figsize=(9, 4.5))
ax = sns.violinplot(
    data=sim_df, x="split", y="max_sim",
    order=order, palette=palette, cut=0, inner=None
)
# overlay a box to mimic the paper's quartile box
sns.boxplot(
    data=sim_df, x="split", y="max_sim",
    order=order, width=0.22, showcaps=True, showfliers=False,
    boxprops={"facecolor": "#f6c7c7", "edgecolor": "#7f1d1d"},
    whiskerprops={"color": "#7f1d1d"}, medianprops={"color": "#7f1d1d"},
    ax=ax
)

ax.set_title("Train-Test Tanimoto Similarity")
ax.set_xlabel("Splitting Method")
ax.set_ylabel("Max Tanimoto Similarity to Training Set")
ax.set_ylim(0, 1)                    # Tanimoto is in [0,1]
plt.tight_layout()

out_png = BASE / "max_tanimoto_by_split.png"
plt.savefig(out_png, dpi=300)
print("Saved:", out_png)
plt.show()
