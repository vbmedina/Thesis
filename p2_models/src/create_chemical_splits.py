#  Run after   conda activate molml
# Random, Random-Sim (fingerprint diversity), Scaffold, Butina, UMAP-HDBSCAN

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap, hdbscan
import argparse
import os
import time
import random
import pickle
import pathlib
import deepchem as dc
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Scaffolds
from rdkit.ML.Cluster import Butina

# Set random seed
SEED = 0
np.random.seed(SEED)
INPUTFILE = "pp_tidy.csv"

# Split object
def main():
    out_root = pathlib.Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/splits") 
    log = []

    df = pd.read_csv(INPUTFILE)
    if {"SMILES","pIC50"}.issubset(df.columns) is False:
        raise ValueError("CSV must contain SMILES and pIC50 columns.")
    active_threshold = 6.0

    # build DeepChem dataset (y only used for strat QC)
    y_arr = df["pIC50"].values.reshape(-1,1)
    ds = dc.data.DiskDataset.from_numpy(
            X = np.zeros((len(df),1)),
            y = y_arr,
            ids = df["SMILES"].values)

    # fingerprint cache for similarity stats
    fps = [morgan_fp(s) for s in df["SMILES"]]

    for split_name, fn in SPLITS.items():
        print(f"\n=== {split_name.upper()} split ===")
        split_dir = out_root / split_name; split_dir.mkdir(exist_ok=True)
        folds = fn(ds, k=5, seed=args.seed)

        for f, (train_ds, val_ds, test_ds) in enumerate(folds):
            foldnum = f+1
            # dataframe views
            tr, vl, te = (df.iloc[train_ds.indices],
                           df.iloc[val_ds.indices],
                           df.iloc[test_ds.indices])
            # save CSVs
            tr.to_csv(split_dir / f"fold{foldnum}_train.csv", index=False)
            vl.to_csv(split_dir / f"fold{foldnum}_val.csv"  , index=False)
            te.to_csv(split_dir / f"fold{foldnum}_test.csv" , index=False)
            # QC print
            qc_print(split_name, f+1, tr["pIC50"], vl["pIC50"], te["pIC50"],
                     active_threshold, log)
            # violin
            plot_df = pd.concat([tr.assign(split="train"),
                                 vl.assign(split="val"),
                                 te.assign(split="test")])
            save_violin(plot_df, "split", "pIC50",
                        f"{split_name} fold {foldnum}", split_dir / f"fold{foldnum}_violin.png")
            # tanimoto stats
            if split_name in {"random_sim","butina","umap"}:
                tr_fps = [fps[i] for i in train_ds.indices]
                mean_sim = []
                for idx in test_ds.indices:
                    sims = DataStructs.BulkTanimotoSimilarity(fps[idx], tr_fps)
                    mean_sim.append(max(sims))
                line = f"    ↳ mean max-Tanimoto(test-to-train) = {np.mean(mean_sim):.3f}"
                print(line); log.append(line)

    # write summary log
    with open(out_root / "split_stats.log","w") as f:
        f.write("\n".join(log))
    print("\nAll splits finished. QC in split_stats.log")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="pp_tidy.csv path")
    ap.add_argument("--out", required=True, help="output root dir for split files")
    ap.add_argument("--seed", type=int, default=0, help="random seed")
    args = ap.parse_args()
    main(args)


# ────────────────────────────────────────────────────────────────────────── helpers
def morgan_fp(smiles, r=2, nBits=1024):
    m = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(m, r, nBits)

def tanimoto_matrix(fps):
    n = len(fps)
    mat = np.empty((n, n), dtype=np.float32)
    for i, fp in enumerate(fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps[:i+1])
        mat[i, :i+1] = sims
        mat[:i+1, i] = sims
    return mat

def qc_print(split_name, fold, y_train, y_val, y_test, active_threshold, log):
    total_active = (y_test >= active_threshold).sum()
    line  = (f"[{split_name} | fold {fold}] "
             f"train {len(y_train):6,} | val {len(y_val):6,} | test {len(y_test):6,} "
             f"| actives≥{active_threshold} in test: {total_active:5,}")
    print(line); log.append(line)

def save_violin(df, label_col, value_col, title, out_png):
    plt.figure(figsize=(6,4))
    sns.violinplot(data=df, x=label_col, y=value_col, cut=0, inner="quartile")
    plt.title(title); plt.tight_layout(); plt.savefig(out_png, dpi=300); plt.close()

# ───────────────────────────────────────────────────────────────────── split blocks
def split_random(ds, k, seed):
    return dc.splits.RandomSplitter().k_fold_split(ds, k=k, seed=seed)

def split_random_sim(ds, k, seed):
    """FingerprintSplitter → 5 folds."""
    splitter = dc.splits.FingerprintSplitter()
    return splitter.k_fold_split(ds, k=k)          # deterministic; seed ignored

def split_scaffold(ds, k, seed):
    return dc.splits.ScaffoldSplitter().k_fold_split(ds, k=k, seed=seed)

def split_butina(ds, k, seed, cutoff=0.6):
    return dc.splits.ButinaSplitter(cutoff).k_fold_split(ds, k=k, seed=seed)

def split_umap_hdb(ds, k, seed):
    """Custom: UMAP(ECFP) → HDBSCAN clusters → greedy fold assignment."""
    smiles = ds.ids; fps = [morgan_fp(s) for s in smiles]
    arr    = np.asarray([np.frombuffer(fp.ToBitString().encode('utf-8'), 'S1').view('i1')
                         for fp in fps])           # binary array (n,1024)
    emb    = umap.UMAP(n_components=2, random_state=seed).fit_transform(arr)
    clus   = hdbscan.HDBSCAN(min_cluster_size=25).fit(emb).labels_
    # assign noise points (-1) to their own clusters
    uniq = np.unique(clus)
    clusters = {c: np.where(clus==c)[0].tolist() for c in uniq}
    if -1 in clusters:
        for i in clusters[-1]:
            clusters[i] = [i]                      # singleton clusters
        del clusters[-1]
    # sort clusters by size & drop into folds round-robin
    fold_inds = [[] for _ in range(k)]
    for c_idx in sorted(clusters, key=lambda c: len(clusters[c]), reverse=True):
        smallest = min(range(k), key=lambda f: len(fold_inds[f]))
        fold_inds[smallest].extend(clusters[c_idx])
    # build (train,val,test) tuples as DeepChem expects
    folds = []
    for f in range(k):
        test   = np.array(fold_inds[f])
        other  = np.setdiff1d(np.arange(len(ds)), test)
        val    = other[::10]                       # 10 % val stratified by order
        train  = np.setdiff1d(other, val)
        folds.append((train, val, test))
    return [ds.select(t)[0] for t in folds]        # convert index triplets to Datasets


# map names → splitter function
SPLITS = {
    "random"     : split_random,
    "random_sim" : split_random_sim,
    "scaffold"   : split_scaffold,
    "butina"     : split_butina,
    "umap"       : split_umap_hdb,
}