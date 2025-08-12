from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import umap, seaborn as sns, matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

BASE = Path("p2_models/data/splits")
METHODS = ["random","scaffold","butina","umap","umap_hdb"]
FOLDS = [1,2,3,4,5]
SMI_COL = "Smiles"
N_BITS, RADIUS = 2048, 2
OUT = Path("p2_models/data/figures/split_folds_embedding.png")

def find_file(dirp, method, fold, split):
    for pat in (f"{method}_fold_{fold}_{split}.csv", f"{method}_fold{fold}_{split}.csv"):
        hits = list((dirp/method).glob(pat))
        if hits: return hits[0]
    raise FileNotFoundError(f"{method} fold {fold} {split}")

def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m, canonical=True) if m is not None else None

# Collect all unique TEST smiles across methods/folds
all_smi = set()
present = []
for m in METHODS:
    if not (BASE/m).exists(): continue
    present.append(m)
    for f in FOLDS:
        try:
            p = find_file(BASE, m, f, "test")
        except FileNotFoundError:
            continue
        df = pd.read_csv(p)
        all_smi.update(df[SMI_COL].dropna().astype(str).tolist())
all_smi = [s for s in (canon(s) for s in all_smi) if s is not None]
if not all_smi: raise SystemExit("No molecules found.")

# ECFP4 bits and a single UMAP (Hamming ≈ Tanimoto on bits)
rows = []
for s in all_smi:
    bv = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=int)
    DataStructs.ConvertToNumpyArray(bv, arr)
    rows.append(arr > 0)
X = np.vstack(rows)
embed = umap.UMAP(n_components=2, metric="hamming", n_neighbors=30, min_dist=0.15, random_state=42)
XY = embed.fit_transform(X)
idx = {s:i for i,s in enumerate(all_smi)}

# For each method, label each point with its TEST fold id (or -1)
labels = {}
for m in present:
    lab = np.full(len(all_smi), -1, int)
    for f in FOLDS:
        try:
            p = find_file(BASE, m, f, "test")
        except FileNotFoundError:
            continue
        df = pd.read_csv(p)
        for s in df[SMI_COL].dropna().astype(str):
            cs = canon(s)
            if cs in idx: lab[idx[cs]] = f
    labels[m] = lab

# Plot panels
sns.set_style("whitegrid")
n = len(present)
fig, axes = plt.subplots(1, n, figsize=(5*n,5), sharex=True, sharey=True)
if n == 1: axes = [axes]
pad = 5
xlim = (XY[:,0].min()-pad, XY[:,0].max()+pad)
ylim = (XY[:,1].min()-pad, XY[:,1].max()+pad)
palette = sns.color_palette("Reds", len(FOLDS)+2)[2:]

for ax, m in zip(axes, present):
    lab = labels[m]
    ax.set_title(f"{m} split")
    for f in FOLDS:
        mask = lab == f
        if not np.any(mask): continue
        ax.scatter(XY[mask,0], XY[mask,1], s=6, alpha=0.45, color=palette[f-1], label=f"Fold {f}")
        pts = XY[mask]
        if pts.shape[0] >= 3:
            try:
                hull = ConvexHull(pts)
                ax.add_patch(plt.Polygon(pts[hull.vertices], fill=False, lw=1.2, ec=palette[f-1]))
            except Exception:
                pass
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)

handles, lbls = axes[-1].get_legend_handles_labels()
fig.legend(handles, lbls, loc="center left", bbox_to_anchor=(1.02, 0.5), title="Folds")
fig.suptitle("Molecular distribution for each fold under different splits", y=0.98, fontsize=14)
plt.tight_layout(rect=[0,0,0.9,0.95])
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300)
print(f"[ok] saved -> {OUT}")