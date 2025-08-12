from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import matplotlib.patches as mpatches

# Path
BASE = Path("p2_models/data/splits")

# Splitting methods
METHODS = ["random", "scaffold", "butina", "umap", "umap_hdb"]  # edit if needed

# Number of folds
FOLDS = [1, 2, 3, 4, 5]

# Parameters
SMI_COL = "Smiles"
N_BITS, RADIUS = 2048, 2  # ECFP4
OUT = Path("p2_models/data/figures/split_folds_embedding.png")
POINT_SIZE = 6.0
ALPHA = 0.45

# Finding files per splitting method and fold
def find_file(dirp: Path, method: str, fold: int, split: str) -> Path:
    """Handle both <m>_fold_1_split.csv and <m>_fold1_split.csv."""
    for pat in (f"{method}_fold_{fold}_{split}.csv", f"{method}_fold{fold}_{split}.csv"):
        hits = list((dirp / method).glob(pat))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"{method} fold {fold} {split}")

# Mols from Smiles
def canon(s: str) -> str | None:
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m, canonical=True) if m is not None else None

# 1) Collect ALL unique SMILES across train, validation, and test for every method and fold
all_smi = set()
present = []
for m in METHODS:
    mdir = BASE / m
    if not mdir.exists():
        continue
    present.append(m)
    for f in FOLDS:
        for split in ("test", "val", "train"):
            try:
                p = find_file(BASE, m, f, split)
            except FileNotFoundError:
                continue
            df = pd.read_csv(p)
            all_smi.update(df[SMI_COL].dropna().astype(str).tolist())

all_smi = [s for s in (canon(s) for s in all_smi) if s is not None]
if not all_smi:
    raise SystemExit("No molecules found. Check BASE/METHODS/SMI_COL.")

# ECFP4 bits and one global UMAP
rows = []
for s in all_smi:
    bv = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=int)
    DataStructs.ConvertToNumpyArray(bv, arr)
    rows.append(arr > 0)
X = np.vstack(rows)

# UMAP settings
embed = umap.UMAP(n_components=2, metric="hamming", n_neighbors=30, min_dist=0.15, random_state=42)
XY = embed.fit_transform(X)
idx = {s: i for i, s in enumerate(all_smi)}

# For each method, assign a fold to every molecule using priority: TEST > VAL > TRAIN
labels = {}
for m in present:
    lab = np.full(len(all_smi), -1, int)

    # set in priority order so 'test' wins over 'val', which wins over 'train'
    for split in ("test", "val", "train"):
        for f in FOLDS:
            try:
                p = find_file(BASE, m, f, split)
            except FileNotFoundError:
                continue
            df = pd.read_csv(p)
            for s in df[SMI_COL].dropna().astype(str):
                cs = canon(s)
                i = idx.get(cs)
                if i is None:
                    continue
                if lab[i] == -1:   # only set if not already labeled by a higher-priority split
                    lab[i] = f

    labels[m] = lab

# Colors
order = FOLDS[:]  # [1..5]
reds = sns.color_palette("Reds", n_colors=len(order))
color_map = dict(zip(order, reds))
legend_handles = [mpatches.Patch(color=color_map[f], label=f"Fold {f}") for f in order]

# Plot panels
sns.set_style("whitegrid")
n = len(present)
fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), sharex=True, sharey=True)
if n == 1:
    axes = [axes]

pad = 5
xlim = (XY[:, 0].min() - pad, XY[:, 0].max() + pad)
ylim = (XY[:, 1].min() - pad, XY[:, 1].max() + pad)

# Titles per panels
for ax, m in zip(axes, present):
    lab = labels[m]
    ax.set_title(f"{m} split")
    for f in order:
        mask = lab == f
        if not np.any(mask):
            continue
        col = color_map[f]
        ax.scatter(XY[mask, 0], XY[mask, 1], s=POINT_SIZE, alpha=ALPHA, color=col)
        pts = XY[mask]
        if pts.shape[0] >= 3:
            try:
                hull = ConvexHull(pts)
                ax.add_patch(plt.Polygon(pts[hull.vertices], fill=False, lw=1.2, ec=col))
            except Exception:
                pass
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

# Legend
leg = fig.legend(
    legend_handles, [f"Fold {f}" for f in order],
    loc="center left",
    bbox_to_anchor=(0.905, 0.5),
    title="Folds",
    frameon=False,
    borderaxespad=0.12,
    prop={"size": 13},
    title_fontsize=13,
    handlelength=1.4,
    handleheight=1.2,
    handletextpad=0.5,
    labelspacing=0.35
)
leg.get_title().set_fontweight("bold")
leg.get_title().set_fontsize(13)

# Title for whole figure
fig.suptitle("Molecular distribution for each fold under different splitting methods",
             y=0.94, fontsize=17, fontweight="bold")
plt.tight_layout(rect=[0, 0, 0.9, 0.95])

# Out
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight", pad_inches=0.15)
print(f"Saved to {OUT}")