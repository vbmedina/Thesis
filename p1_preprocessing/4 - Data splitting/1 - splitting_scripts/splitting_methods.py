"""
Description: This script uses random, scaffold, Butina, UMAP (k-mean) and UMAP (ward) as splitting methods to create training, 
validation, and test sets. 

Goal:
Build chemically realistic cross-validation folds for virtual screening evaluation on ChEMBL P. falciparum IC50 data 
(~40k pIC50 values for ~20k molecules). Study which splitting methods increase train-test dissimilarity to better reflect 
real-world prospective screening scenarios where models encounter molecules unlike those in training data. (1)

K-fold scheme and set sizes
This script creates 5-fold cross-validation splits with train/val/test sets per fold. 
In each fold it aims to allocate:
  • Test = 1/5 of molecules (≈ 20%) held out for evaluation as a baseline and final model testing.
  • Validation = 10% of molecules from the remaining pool (80%), and is used for model selection and early stopping.
  • Train = 70% of molecules for fitting the model.

The fundamental trade-off between maximizing training data and maintaining rigorous evaluation standards. This approach 
prevents the same molecular structure from appearing in both training and test sets within each fold by keeping all strain 
measurements of a molecule together, while the 5-fold framework provides multiple independent train/test partitions to assess 
the consistency and robustness of each splitting method.

Why five splitting methods?
Random — Computationally efficient baseline that assigns data points to splits purely at random, without regard to the
    chemical similarity, diversity, or scaffold of molecules. As a result, it carries a higher probability of placing chemically
    similar compounds across train/test sets compared to other chemically aware splitting methods. (2)
Scaffold (Bemis-Murcko) — Extracts the Bemis-Murcko scaffold (core ring system) from each molecule and assigns all molecules
    sharing identical scaffolds to the same split, preventing scaffold leakage between train and test sets. However, this
    method treats structurally similar scaffolds as completely distinct and fails to account for molecules with different
    scaffolds that share similar pharmacophores or binding modes, allowing chemical similarity to persist across splits. (1)(3)
Butina (Tanimoto coefficient on Morgan fingerprints) — Creates clusters of molecules based on fingerprint similarity using a
    distance threshold, keeping structurally related compounds within the same split to reduce near-neighbor leakage compared
    to scaffold splitting. However, cluster quality depends heavily on the choice of similarity threshold, and the method may
    struggle with boundary cases where molecules fall near cluster edges or in sparse regions of chemical space. (1)(3)
UMAP + K-means clustering — Applies UMAP dimensionality reduction to Morgan fingerprints using Jaccard distance, then performs
    K-means clustering on the low-dimensional embedding to partition molecules into a predetermined number of groups. This
    approach creates compact, spherical clusters in the embedding space and provides deterministic cluster assignments, but
    requires specifying the number of clusters and may struggle with irregular cluster shapes or varying cluster densities
    in chemical space. (1)(2)(4)
UMAP + Ward clustering — Applies UMAP dimensionality reduction to Morgan fingerprints using Jaccard distance, then performs
    hierarchical agglomerative clustering on the low-dimensional embedding to create chemically coherent groups. This two-step
    process preserves both local and global chemical relationships in the embedding space, yielding larger distribution shifts.
     However, performance depends on UMAP hyperparameter tuning (n_neighbors, min_dist) and the method requires determining 
     the final number of clusters from the hierarchical structure. (2)(4)

All methods use established default parameters: Morgan fingerprints (radius=2, 2048 bits), Butina cutoff=0.6, 
UMAP (n_neighbors=25, min_dist=0.1, Jaccard metric), and pIC50 ≥ 6.0 activity threshold.

Dataset-specific constraints that guided this design:
- Limited sample size: The 5-fold cross-validation approach reuses all data without test leakage, maximizing both training 
    signal and evaluation coverage. 
- Multiple strain measurements: When the same molecule has been tested across multiple strains, all of their measurements are 
    kept in the same fold, allowing models to learn from the full biological activity profile of each compound and make robust 
    predictions despite inherent experimental variability. This reflects real-world scenarios where compounds show activity 
    variation across conditions.
- Class prevalence: While real virtual screening campaigns often have very low hit rates (many more inactives), this ChEMBL 
    dataset shows a relatively balanced pIC50 distribution (with slightly more inactive molecules). This study has not
    artificially force class imbalance in the splitting step. Instead, it approximates deployment difficulty via 
    chemistry-aware splits. (5)

Outputs and quality control for each split method and fold:
- Files: train/val/test CSVs per fold
- Logging: Set sizes, active counts (pIC50 ≥ 6.0), mean max-Tanimoto(test→train) as difficulty metric
- Expected ranking: Random < Scaffold < Butina < UMAP+K-means < UMAP+Ward

References (code/data in paper):
1)“Comprehensive study showing UMAP > Butina > Scaffold > Random in realism of train-test splits for virtual screening”
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-025-01039-8; 
2) https://github.com/rdkit
3) https://arxiv.org/pdf/2406.00873
4) https://umap-learn.readthedocs.io/en/latest/faq.html
5) https://pubs.acs.org/doi/10.1021/acsmedchemlett.4c00093
"""

# Conda activate molml
# Imports
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

import deepchem as dc
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from sklearn.cluster import AgglomerativeClustering, KMeans
import umap

RDLogger.DisableLog("rdApp.warning")

# Input/Output files
INPUTFILE = "./p1_preprocessing/4 - Data splitting/1 - splitting_scripts/final_data_copy.csv"
OUTROOT = Path("./p1_preprocessing/4 - Data splitting/2 - split_data")

# Settings
SEED = np.random.randint(1, 10000000000)
print(f"Random seed: {SEED}")
np.random.seed(SEED) # Seed used for study was 2858808528
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

# UMAP (+ k-means)
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
        test_idx = np.where(labels == f)[0]
        train_idx = np.setdiff1d(np.arange(n_all), test_idx)
        folds.append((ds.select(train_idx), ds.select(test_idx)))
    return folds

# UMAP (+ Ward)
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
    "random" : split_random,
    "scaffold" : split_scaffold,
    "butina" : split_butina,
    "umap_kmeans" : split_umap_kmeans,
    "umap_ward" : split_umap_ward,
}
# Helpers for non-unique SMILES
def build_id2rows(ids: List[str]) -> Dict[str, np.ndarray]:
    id2rows = defaultdict(list)
    for i, s in enumerate(ids):
        id2rows[s].append(i)
    return {s: np.asarray(ix, dtype=int) for s, ix in id2rows.items()}

def rows_for_subset(subset_ids: List[str],
                    id2rows: Dict[str, np.ndarray],
                    counters: Dict[str, int] = None) -> Tuple[np.ndarray, Dict[str, int]]:
    if counters is None:
        counters = defaultdict(int)
    out = np.empty(len(subset_ids), dtype=int)
    for j, s in enumerate(subset_ids):
        k = counters[s]
        if k >= len(id2rows[s]):
            raise IndexError( f"Subset expects more occurrences of SMILES {s!r} " 
                            f"({k+1}) than exist in the full dataset ({len(id2rows[s])}).")
        out[j] = id2rows[s][k]
        counters[s] = k + 1
    return out, counters

def choose_val_rows_by_smiles(train_ids: List[str], train_rows: np.ndarray, desired_val_rows: int,
                              seed: int) -> Tuple[np.ndarray, np.ndarray]:

    # Group assigned TRAIN rows by SMILES
    smi2rows = defaultdict(list)
    for s, r in zip(train_ids, train_rows):
        smi2rows[s].append(int(r))

    rng = np.random.RandomState(seed)
    smiles = list(smi2rows.keys())
    rng.shuffle(smiles)

    val_rows_list = []
    for s in smiles:
        val_rows_list.extend(smi2rows[s])
        if len(val_rows_list) >= desired_val_rows:
            break

    val_rows = set(val_rows_list)
    final_train_rows = [r for r in train_rows if int(r) not in val_rows]
    return np.asarray(final_train_rows, dtype=int), np.asarray(sorted(val_rows), dtype=int)

# Main loop
def main():
    OUTROOT.mkdir(parents=True, exist_ok=True)
    log = []

    df = pd.read_csv(INPUTFILE)
    if not {"Smiles", "pIC50"}.issubset(df.columns):
        raise ValueError("CSV must contain Smiles and pIC50 columns.")

    # DeepChem dataset:
    y_arr = df["pIC50"].values.reshape(-1, 1)
    ds = dc.data.DiskDataset.from_numpy(
        X=np.zeros((len(df), 1)),
        y=y_arr,
        ids=df["Smiles"].values
    )

    # Map each SMILES to its row indices
    id2rows_all = build_id2rows(list(ds.ids))

    # Dataset of UNIQUE SMILES for molecule-level splitting
    unique_smiles = np.array(sorted(id2rows_all.keys()))
    ds_u = dc.data.DiskDataset.from_numpy(
        X=np.zeros((len(unique_smiles), 1)),
        y=np.zeros((len(unique_smiles), 1)),
        ids=unique_smiles
    )

    # Precompute fingerprints for QC similarity (row-level)
    fps = [morgan_fp(s) for s in df["Smiles"]]

    # Per-split processing
    for split_name, splitter in SPLITS.items():
        print(f"\n=== {split_name.upper()} (k={K}) ===")
        split_dir = OUTROOT / split_name
        split_dir.mkdir(exist_ok=True, parents=True)

        # Split on UNIQUE SMILES rather than on rows
        base_folds_u = splitter(ds_u, k=K, seed=SEED)

        # Expand unique-SMILES folds back to ALL rows for those SMILES
        for f, (train_u, test_u) in enumerate(base_folds_u, start=1):

            # rows for all occurrences of each SMILES
            train_rows = np.concatenate([id2rows_all[s] for s in train_u.ids])
            test_rows  = np.concatenate([id2rows_all[s] for s in test_u.ids])

            # Validation: 10% of total rows, sampled by whole SMILES (keeps molecules intact)
            n_total = len(ds)
            desired_val = max(1, int(round(VAL_FRAC_TOTAL * n_total)))
            desired_val = min(desired_val, len(train_rows))

            # Pass unique train SMILES so validation respects molecule grouping 
            final_train_rows, val_rows = choose_val_rows_by_smiles(
                list(train_u.ids), train_rows, desired_val, seed=SEED + f
            )

            # DataFrames for each set
            tr = df.loc[final_train_rows]
            vl = df.loc[val_rows]
            te = df.loc[test_rows]

            # Save CSVs
            (split_dir / f"{split_name}_fold_{f}_train.csv").write_text(tr.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_val.csv").write_text(vl.to_csv(index=False))
            (split_dir / f"{split_name}_fold_{f}_test.csv").write_text(te.to_csv(index=False))

            # QC sizes
            qc_print(split_name, f, tr["pIC50"].values, vl["pIC50"].values, te["pIC50"].values,
                     ACTIVE_THRESHOLD, log)

            # Mean max-Tanimoto on test and train
            tr_fps = [fps[i] for i in final_train_rows]
            if len(tr_fps) == 0:
                mean_max = float("nan")
            else:
                max_sims = []
                for i in test_rows:
                    sims = DataStructs.BulkTanimotoSimilarity(fps[i], tr_fps)
                    max_sims.append(max(sims) if len(sims) else 0.0)
                mean_max = float(np.mean(max_sims))
            line = f" mean max-Tanimoto(test-train) = {mean_max:.3f}"
            print(line); log.append(line)

    log.append(f"Random seed: {SEED}")
    (OUTROOT / "split_stats.log").write_text("\n".join(log))
    print("\nAll folds written. QC in split_stats.log")

if __name__ == "__main__":
    main()