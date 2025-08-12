# conda activate molml
import pandas as pd
import numpy as np
from pathlib import Path
import deepchem as dc
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.cluster import AgglomerativeClustering
import umap
import hdbscan

RDLogger.DisableLog("rdApp.warning")

# Settings
SEED = 0
np.random.seed(SEED)
INPUTFILE = "/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/tidy/pp_tidy.csv"
OUTROOT  = Path("/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/training")
K = 5 
VAL_FRAC_TOTAL = 0.10
ACTIVE_THRESHOLD = 6.0

# Fingerprint helpers
def morgan_fp(smi, r=2, nBits=1024):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, r, nBits)

def fp_to_numpy(fp, nBits=1024):
    arr = np.zeros((nBits,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# Logging while splitting
def qc_print(split_name, fold, y_train, y_val, y_test, active_threshold, log):
    n_active_val  = (y_val  >= active_threshold).sum()
    n_active_test = (y_test >= active_threshold).sum()
    line = (f"[{split_name} | fold {fold}] "
            f"train {len(y_train):6,} | val {len(y_val):6,} | test {len(y_test):6,} "
            f"| actives≥{active_threshold} val: {n_active_val:5,} test: {n_active_test:5,}")
    print(line); log.append(line)

# Splits 
def split_random(ds, k, seed):
    return dc.splits.RandomSplitter().k_fold_split(ds, k=k, seed=seed)

def split_scaffold(ds, k, seed):
    return dc.splits.ScaffoldSplitter().k_fold_split(ds, k=k, seed=seed)

def split_butina(ds, k, seed, cutoff=0.6):
    return dc.splits.ButinaSplitter(cutoff).k_fold_split(ds, k=k, seed=seed)

# UMAP(k-fold) block
def split_umap_k(ds, k, seed):
    smiles = ds.ids
    fps = [morgan_fp(s) for s in smiles]
    arr = np.stack([fp_to_numpy(fp) for fp in fps]).astype(bool)
    emb = umap.UMAP(n_components=2, n_neighbors=25, min_dist=0.1,
                    metric="jaccard", random_state=seed).fit_transform(arr)
    labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(emb)

    folds = []
    n_all = len(ds)
    for f in range(k):
        test_idx  = np.where(labels == f)[0]
        train_idx = np.setdiff1d(np.arange(n_all), test_idx)
        folds.append((ds.select(train_idx), ds.select(test_idx)))
    return folds

# UMAP (HDBSCAN) block
def split_umap_hdb(ds, k, seed):
    smiles = ds.ids
    fps = [morgan_fp(s) for s in smiles]
    arr = np.stack([fp_to_numpy(fp) for fp in fps]).astype(bool)
    emb = umap.UMAP(n_components=2, n_neighbors=25, min_dist=0.1,
                    metric="jaccard", random_state=seed).fit_transform(arr)
    labels = hdbscan.HDBSCAN(min_cluster_size=25, min_samples=10).fit(emb).labels_

    clusters = {}
    for i, c in enumerate(labels):
        clusters.setdefault(c, []).append(i)
    if -1 in clusters:
        for i in clusters[-1]:
            clusters[i] = [i]
        del clusters[-1]

    fold_bins = [[] for _ in range(k)]
    for c in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
        tgt = min(range(k), key=lambda f: len(fold_bins[f]))
        fold_bins[tgt].extend(clusters[c])

    folds = []
    n_all = len(ds)
    for f in range(k):
        test_idx  = np.array(fold_bins[f], dtype=int)
        train_idx = np.setdiff1d(np.arange(n_all), test_idx)
        folds.append((ds.select(train_idx), ds.select(test_idx)))
    return folds

# Splits
SPLITS = {
    "random"   : split_random,
    "scaffold" : split_scaffold,
    "butina"   : split_butina,
    "umap"     : split_umap_k,
    "umap_hdb" : split_umap_hdb,
}

# 10% validation set (of the whole dataset) from train pool
def add_validation(ds_all, train_ds, test_ds, seed, val_frac_total=0.10):
    n_total = len(ds_all)
    desired_val = max(1, int(round(val_frac_total * n_total)))

    id2row = {smi: i for i, smi in enumerate(ds_all.ids)}

    train_idx_pool = np.array([id2row[s] for s in train_ds.ids], dtype=int)
    rng = np.random.RandomState(seed)
    val_idx = rng.choice(train_idx_pool, size=desired_val, replace=False)
    final_train_idx = np.setdiff1d(train_idx_pool, val_idx)

    return ds_all.select(final_train_idx), ds_all.select(val_idx), test_ds

# Main funct
def main():
    OUTROOT.mkdir(parents=True, exist_ok=True)
    log = []

    df = pd.read_csv(INPUTFILE)
    if not {"Smiles", "pIC50"}.issubset(df.columns):
        raise ValueError("CSV must contain Smiles and pIC50 columns.")

    # DeepChem dataset
    y_arr = df["pIC50"].values.reshape(-1, 1)
    ds = dc.data.DiskDataset.from_numpy(
        X=np.zeros((len(df), 1)),
        y=y_arr,
        ids=df["Smiles"].values
    )

    # Similarity logs MFs
    fps = [morgan_fp(s) for s in df["Smiles"]]
    id2row = dict(zip(df["Smiles"], df.index))
    def ds_to_df(dset):
        rows = [id2row[s] for s in dset.ids]
        return df.loc[rows]

    # Saving splits 
    for split_name, splitter in SPLITS.items():
        print(f"\n=== {split_name.upper()} (k={K}) ===")
        split_dir = OUTROOT / split_name
        split_dir.mkdir(exist_ok=True)

        base_folds = splitter(ds, k=K, seed=SEED)

        for f, (train_ds, test_ds) in enumerate(base_folds, start=1):
            train_ds, val_ds, test_ds = add_validation(ds, train_ds, test_ds, seed=SEED+f,
                                                       val_frac_total=VAL_FRAC_TOTAL)

            tr, vl, te = ds_to_df(train_ds), ds_to_df(val_ds), ds_to_df(test_ds)

            # Save CSVs
            (split_dir / f"{split_name}_fold_{f}_train.csv").write_text(tr.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_val.csv").write_text(vl.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_test.csv").write_text(te.to_csv(index=False))

            # QC sizes
            qc_print(split_name, f, tr["pIC50"], vl["pIC50"], te["pIC50"],
                     ACTIVE_THRESHOLD, log)

            # Mean max-Tanimoto(test to train)
            tr_fps = [fps[id2row[s]] for s in train_ds.ids]
            max_sims = []
            for s in test_ds.ids:
                sims = DataStructs.BulkTanimotoSimilarity(fps[id2row[s]], tr_fps)
                max_sims.append(max(sims))
            line = f"    mean max-Tanimoto(test→train) = {np.mean(max_sims):.3f}"
            print(line); log.append(line)

    (OUTROOT / "split_stats.log").write_text("\n".join(log))
    print("\nAll folds written. QC in split_stats.log")

if __name__ == "__main__":
    main()
