#  Run after   conda activate molml
# Random, Random-Sim (fingerprint diversity), Scaffold, Butina, UMAP-HDBSCAN

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap, hdbscan
import pathlib
import deepchem as dc
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Scaffolds
from rdkit.ML.Cluster import Butina
from sklearn.model_selection import train_test_split
RDLogger.DisableLog("rdApp.warning")

# Set random seed
SEED = 0
np.random.seed(SEED)
INPUTFILE = "/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/tidy/pp_tidy.csv"

# Generate Morgan fingerprints
def morgan_fp(smiles, r=2, nBits=1024):
    m = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(m, r, nBits)

# Tanimoto similarity matrix for fingerprints
def tanimoto_matrix(fps):
    n = len(fps)
    mat = np.empty((n, n), dtype=np.float32)
    for i, fp in enumerate(fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps[:i+1])
        mat[i, :i+1] = sims
        mat[:i+1, i] = sims
    return mat

# Quality control printout
def qc_print(split_name, fold, y_train, y_val, y_test, active_threshold, log):
    total_active = (y_test >= active_threshold).sum()
    line  = (f"[{split_name} | fold {fold}] "
             f"train {len(y_train):6,} | val {len(y_val):6,} | test {len(y_test):6,} "
             f"| actives≥{active_threshold} in test: {total_active:5,}")
    print(line); log.append(line)

# Generate violin plots
def save_violin(df, label_col, value_col, title, out_png):
    plt.figure(figsize=(6,4))
    sns.violinplot(data=df, x=label_col, y=value_col, cut=0, inner="quartile")
    plt.title(title); plt.tight_layout(); plt.savefig(out_png, dpi=300); plt.close()

# Split blocks
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
    # Smiles and Morgan fingerprints for Data Disk
    smiles = ds.ids
    fps    = [morgan_fp(s) for s in smiles]

    # Set array for UMAP
    arr = np.asarray([
        np.frombuffer(fp.ToBitString().encode("utf-8"), "S1").view("i1")
        for fp in fps]) 
    
    # UMAP embedding
    emb = umap.UMAP(
        n_components=2,
        init="random",
        random_state=seed
    ).fit_transform(arr)

    # HDBSCAN clustering 
    labels = hdbscan.HDBSCAN(min_cluster_size=25).fit(emb).labels_

    clusters = {}
    for idx, lab in enumerate(labels):
        clusters.setdefault(lab, []).append(idx)

    # Treat noise pointss as singles
    if -1 in clusters:
        for idx in clusters[-1]:
            clusters[idx] = [idx]
        del clusters[-1]

    #Round-robin cluster assignment
    fold_bins = [[] for _ in range(k)]
    for lab in sorted(clusters, key=lambda c: len(clusters[c]), reverse=True):
        target = min(range(k), key=lambda f: len(fold_bins[f]))
        fold_bins[target].extend(clusters[lab])

    # Build (train, validation, test) datasets
    fold_datasets = []
    n_all = len(ds)

    for f in range(k):
        test_idx  = np.array(fold_bins[f])
        other_idx = np.setdiff1d(np.arange(n_all), test_idx)

        # 10 % of 'other' for validation
        val_idx, train_idx = train_test_split(
            other_idx,
            test_size=0.90,
            random_state=seed + f,
            shuffle=True
        )

        train_ds = ds.select(train_idx)
        val_ds   = ds.select(val_idx)
        test_ds  = ds.select(test_idx)

        fold_datasets.append((train_ds, val_ds, test_ds))

    return fold_datasets

def ensure_val(train_ds, test_ds, val_frac=0.10, seed=0):

    # Case where custom splitter already returned (train, val, test)
    if isinstance(train_ds, tuple) and len(train_ds) == 3:
        return train_ds  

    n = len(train_ds)
    idx = np.arange(n)

    # stratify optionally on y if you like; here we ignore stratification
    train_idx, val_idx = train_test_split(
        idx,
        test_size=val_frac,
        random_state=seed,
        shuffle=True
    )

    # DeepChem .select() returns (new_dataset, selected_data_indices)
    new_train_ds = train_ds.select(train_idx)
    val_ds       = train_ds.select(val_idx)

    return new_train_ds, val_ds, test_ds

# Map names
SPLITS = {
    "random"     : split_random,
    "random_sim" : split_random_sim,
    "scaffold"   : split_scaffold,
    "butina"     : split_butina,
    "umap"       : split_umap_hdb,
}

# Split object
def main():
    out_root = pathlib.Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/splits") 
    log = []

    df = pd.read_csv(INPUTFILE)
    if {"Smiles","pIC50"}.issubset(df.columns) is False:
        raise ValueError("CSV must contain Smiles and pIC50 columns.")
    active_threshold = 6.0

    # Vuild DeepChem dataset
    y_arr = df["pIC50"].values.reshape(-1,1)
    ds = dc.data.DiskDataset.from_numpy(
            X = np.zeros((len(df),1)),
            y = y_arr,
            ids = df["Smiles"].values)

    # Fingerprint cache for similarity stats
    fps = [morgan_fp(s) for s in df["Smiles"]]

    id2row = dict(zip(df["Smiles"], df.index))
    def ds_to_df(dset):
        rows = [id2row[smi] for smi in dset.ids]
        return df.loc[rows]

    for split_name, fn in SPLITS.items():
        print(f"\n=== {split_name.upper()} split ===")
        split_dir = out_root / split_name; split_dir.mkdir(exist_ok=True)
        folds = fn(ds, k=5, seed=SEED)

        for f, fold_tuple in enumerate(folds):
            # Make sure we end up with (train_ds, val_ds, test_ds)
            if len(fold_tuple) == 2:
                train_ds, test_ds = fold_tuple
                train_ds, val_ds, test_ds = ensure_val(train_ds, test_ds, seed=SEED+f)
            else:
                train_ds, val_ds, test_ds = fold_tuple

            foldnum = f+1

            tr = ds_to_df(train_ds)
            vl = ds_to_df(val_ds)
            te = ds_to_df(test_ds)

            # Save CSVs
            tr.to_csv(split_dir / f"{split_name}_fold_{foldnum}_train.csv", index=False)
            vl.to_csv(split_dir / f"{split_name}_fold{foldnum}_val.csv"  , index=False)
            te.to_csv(split_dir / f"{split_name}_fold{foldnum}_test.csv" , index=False)

            # QC print
            qc_print(split_name, f+1, tr["pIC50"], vl["pIC50"], te["pIC50"],
                     active_threshold, log)
            
            # Violin plots
            plot_df = pd.concat([tr.assign(split="train"),
                                 vl.assign(split="val"),
                                 te.assign(split="test")])
            save_violin(plot_df, "split", "pIC50",
                        f"{split_name} fold {foldnum}", split_dir / f"{split_name}_fold{foldnum}_violin.png")
            # Tanimoto stats
            if split_name in {"random_sim","butina","umap"}:
                tr_fps = [fps[id2row[smi]] for smi in train_ds.ids]

                mean_sim = []
                for smi in test_ds.ids:
                    sims = DataStructs.BulkTanimotoSimilarity(
                        fps[id2row[smi]], tr_fps
                    )
                    mean_sim.append(max(sims))

                line = f"    mean max-Tanimoto(test→train) = {np.mean(mean_sim):.3f}"
                print(line); log.append(line)

    # Summary log
    with open(out_root / "split_stats.log","w") as f:
        f.write("\n".join(log))
    print("\nAll splits finished. QC in split_stats.log")

if __name__ == "__main__":
    main()