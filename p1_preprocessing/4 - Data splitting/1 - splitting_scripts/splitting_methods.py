"""
Description: This script uses random, scaffold, Butina, UMAP and UMAP with HDBSCAN as splitting methods to create training, 
validation, and test sets.
Dataset: ChEMBL 364 P. Falciparum IC50s (~40k pIC50 values for ~20k molecules)

Goal:
Build chemically realistic cross-validation folds for virtual screening style evaluation on the ChEMBL set. Prospective 
screens used on real-world data, encounter molecules unlike those found in training chemistry, therefore this experiment aims 
to study which splitting methods increase train-test dissimilarity.

K-fold scheme and set sizes
This script creates 5-fold cross-validation splits with train/val/test sets per fold. 
In each fold:
  • Test = 1/5 of molecules (≈ 20%) held out for evaluation as a baseline and final model testing.
  • Validation = 10% of molecules from the remaining pool (80%). This set will be used for model selection and early stopping.
  • Train = 70% of molecules for fitting the model.

This approach addresses the fundamental trade-off in machine learning with limited data: maximizing training data for model 
learning while maintaining rigorous evaluation standards. Each fold uses 70% of molecules for training while reserving 20% as 
a completely unseen test set. The 5-fold cross-validation framework reduces evaluation variance by ensuring every molecule 
serves as test data exactly once across folds. This method provides a comprehensive assessment of model performance across the 
entire dataset rather than relying on a single, potentially biased train-test partition.

Why five splitting methods?
Random — Computationally efficient baseline that assigns data points to splits purely at random, without regard to the 
    chemical similarity, diversity, or scaffold of molecules. As a results a carries higher probability of placing chemically 
    similar compounds across train/test sets when compared to other methods chemically aware splitting methods.
Scaffold (Bemis-Murcko) —  Extracts the Bemis-Murcko scaffold (core ring system) from each molecule and assigns all molecules 
    sharing identical scaffolds to the same split, preventing scaffold leakage between train and test sets. However, this 
    method treats structurally similar scaffolds as completely distinct and fails to account for molecules with different 
    scaffolds that share similar pharmacophores or binding modes, allowing chemical similarity to persist across splits. 
Butina (Tanimoto coefficient on Morgan fingerprints) — Creates clusters of molecules based on fingerprint similarity using a 
    distance threshold, keeping structurally related compounds within the same split to reduce near-neighbor leakage compared 
    to scaffold splitting. However, cluster quality depends heavily on the similarity threshold selection, and the method may 
    struggle with boundary cases where molecules fall near cluster edges or in sparse regions of chemical space.
UMAP + K-means clustering — Applies UMAP dimensionality reduction to Morgan fingerprints using Jaccard distance, then performs
    K-means clustering on the low-dimensional embedding to partition molecules into a predetermined number of groups. This 
    approach creates compact, spherical clusters in the embedding space and provides deterministic cluster assignments, but 
    requires pre-specifying the number of clusters and may struggle with irregular cluster shapes or varying cluster densities
    in chemical space.
UMAP + Ward clustering — Applies UMAP dimensionality reduction to Morgan fingerprints using Jaccard distance, then performs 
    hierarchical agglomerative clustering on the low-dimensional embedding to create chemically coherent groups. This two-step 
    process preserves both local and global chemical relationships in the embedding space, yielding larger distribution shifts 
    and more realistic test scenarios compared to simpler splitting methods. However, performance depends on UMAP 
    hyperparameter tuning (n_neighbors, min_dist) and the choice of cluster cutoff in the resulting dendrogram.

Dataset-specific constraints that guided this design
• Limited sample size: a single, large hold-out would waste data; 5-fold CV reuses data
  efficiently without test leakage.
• Strain imbalance: ~½ of measurements are concentrated in a few strains (e.g., 3D7/K1/W2/NF54),
  so we split **by molecule (SMILES)** and keep all measurements for the same molecule in the
  same fold to avoid “the same structure in train and test via different strains.”
• Class prevalence: real VS has tiny hit rates (many more negatives), but ChEMBL here is not
  extremely skewed (see the pIC50 histogram). Rather than force artificial imbalance, we rely on
  harder chemistry splits (UMAP-based) and, downstream, early-recognition metrics/top-k hit rate.

Outputs and QC
For each split method and fold we write train/val/test CSVs and log:
  • set sizes and active counts at pIC50 ≥ 6.0
  • mean max-Tanimoto(test→train) as a leakage/difficulty check (lower is harder)

References (code/data in paper):
1)“Comprehensive study showing UMAP > Butina > Scaffold > Random in realism and why ROC AUC can mislead virtual screenings: 
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00576-2; 
GitHub: https://github.com/Rong830/UMAP_split_for_VS archived in Zenodo: https://zenodo.org/records/14736486
2) https://github.com/rdkit
2) https://umap-learn.readthedocs.io/en/latest/faq.html
3) https://arxiv.org/pdf/2406.00873
"""

# conda activate molml
import pandas as pd
import numpy as np
from pathlib import Path
import deepchem as dc
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from sklearn.cluster import AgglomerativeClustering, KMeans
import umap

RDLogger.DisableLog("rdApp.warning")

# Settings
SEED = 0
np.random.seed(SEED)
INPUTFILE = "./p1_preprocessing/4 - Data splitting/1 - splitting_scripts/final_data_copy.csv"
OUTROOT  = Path("./p1_preprocessing/4 - Data splitting/2 - split_data")
K = 5
VAL_FRAC_TOTAL = 0.10
ACTIVE_THRESHOLD = 6.0

# Fingerprint helpers
def morgan_fp(smi, r=2, nBits=2048):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, r, nBits)

def fp_to_numpy(fp, nBits=2048):
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

# UMAP (+ k-means) block
def split_umap_kmeans(ds, k, seed):
    smiles = ds.ids
    fps = [morgan_fp(s) for s in smiles]
    arr = np.stack([fp_to_numpy(fp) for fp in fps]).astype(bool)

    emb = umap.UMAP(
        n_components=2, n_neighbors=25, min_dist=0.1,
        metric="jaccard", random_state=seed
    ).fit_transform(arr)

    labels = KMeans(n_clusters=k, random_state=seed, n_init=20).fit_predict(emb)

    folds = []
    n_all = len(ds)
    for f in range(k):
        test_idx  = np.where(labels == f)[0]
        train_idx = np.setdiff1d(np.arange(n_all), test_idx)
        folds.append((ds.select(train_idx), ds.select(test_idx)))
    return folds

# UMAP (+ Ward) block
def split_umap_ward(ds, k, seed):
    smiles = ds.ids
    fps = [morgan_fp(s) for s in smiles]
    arr = np.stack([fp_to_numpy(fp) for fp in fps]).astype(bool)

    emb = umap.UMAP(
        n_components=2, n_neighbors=25, min_dist=0.1,
        metric="jaccard", random_state=seed
    ).fit_transform(arr)

    labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(emb)

    folds = []
    n_all = len(ds)
    for f in range(k):
        test_idx  = np.where(labels == f)[0]
        train_idx = np.setdiff1d(np.arange(n_all), test_idx)
        folds.append((ds.select(train_idx), ds.select(test_idx)))
    return folds

# Split registry
SPLITS = {
    "random"      : split_random,
    "scaffold"    : split_scaffold,
    "butina"      : split_butina,
    "umap_kmeans" : split_umap_kmeans,
    "umap_ward"   : split_umap_ward,
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

    # DeepChem dataset (store SMILES in ids, pIC50 in y)
    y_arr = df["pIC50"].values.reshape(-1, 1)
    ds = dc.data.DiskDataset.from_numpy(
        X=np.zeros((len(df), 1)),
        y=y_arr,
        ids=df["Smiles"].values
    )

    # Precompute fingerprints for QC similarity
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
            train_ds, val_ds, test_ds = add_validation(
                ds, train_ds, test_ds, seed=SEED + f, val_frac_total=VAL_FRAC_TOTAL
            )

            tr, vl, te = ds_to_df(train_ds), ds_to_df(val_ds), ds_to_df(test_ds)

            # Save CSVs
            (split_dir / f"{split_name}_fold_{f}_train.csv").write_text(tr.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_val.csv").write_text(vl.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_test.csv").write_text(te.to_csv(index=False))

            # QC sizes
            qc_print(split_name, f, tr["pIC50"], vl["pIC50"], te["pIC50"],
                     ACTIVE_THRESHOLD, log)

            # Mean max-Tanimoto(test-train)
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