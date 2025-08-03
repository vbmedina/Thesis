# BEFORE RUNNING -------------------------------
# Conda activate molml
# Python -m src.create_splits

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import StratifiedKFold

# Config
root   = Path(__file__).resolve().parents[1] 
data   = root / "data"
tidy   = data / "tidy"
splits = data / "splits"
rng    = np.random.default_rng(seed=0)

# Murcko scaffold SMILES (if RDKit cant parse)
def murcko(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles or "")
    return MurckoScaffold.MurckoScaffoldSmiles(mol) if mol else ""

#  Write 80/10/10 CSVs for a single fold  
def write_fold(df_fold: pd.DataFrame, split_dir: Path, fold_id: int) -> None:
    test_idx  = df_fold.index
    other_idx = df.index.difference(test_idx)

    # stratified 10 % of “other” rows for validation
    val_idx = RNG.choice(other_idx, size=int(0.10 * len(df)), replace=False)
    train_idx = other_idx.difference(val_idx)

    fold_dir = split_dir / f"fold{fold_id}"
    fold_dir.mkdir(parents=True, exist_ok=True)

    df.loc[train_idx].to_csv(fold_dir / "train.csv", index=False)
    df.loc[val_idx].to_csv  (fold_dir / "val.csv",   index=False)
    df.loc[test_idx].to_csv (fold_dir / "test.csv",  index=False)

# main 
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate chemistry splits")
    parser.add_argument("--seed", type=int, default=0, help="RNG seed")
    args = parser.parse_args()
    RNG = np.random.default_rng(args.seed)

    # 1  load tidy table
    tidy_csv = tidy / "pp_tidy.csv"
    if not tidy_csv.exists():
        raise FileNotFoundError(f"{tidy_csv} not found.")

    df = pd.read_csv(tidy_csv)
    print(f"Loaded {len(df):,} rows")

    # 2  helper columns
    df["potency_bin"] = pd.cut(
        df.pIC50, bins=[-1, 5, 6, 7, 99], labels=["≤5", "5-6", "6-7", "≥7"]
    )
    mol_freq = df.groupby("compound_id").size()
    df["mol_freq"] = df.compound_id.map(mol_freq)

    df["scaffold"] = df.SMILES.map(murcko)

    # 3  Random-stratified split 
    strat_key = (
        df.strain + "_" +
        df.potency_bin.astype(str) + "_" +
        pd.qcut(df.mol_freq, [0, .8, .95, 1],
                labels=["common", "frequent", "hyper"]).astype(str)
    )

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.seed)
    split_dir = splits / "random"
    for fold, (_, test_idx) in enumerate(skf.split(df, strat_key)):
        write_fold(df.loc[test_idx], split_dir, fold)
    print("✓ random-stratified 5-fold saved to", split_dir)

    # 4  Scaffold split 

    scaf_groups = df["scaffold"].unique()
    RNG.shuffle(scaf_groups)
    chunks = np.array_split(scaf_groups, 5)

    split_dir = splits / "scaffold"
    for fold, chunk in enumerate(chunks):
        idx = df[df.scaffold.isin(chunk)].index
        write_fold(df.loc[idx], split_dir, fold)
    print("✓ scaffold 5-fold saved to", split_dir)

# # Run the following seperately if you have UMAP and HDBSCAN installed
# # Otherwise, it will fall back to Butina clustering
# # 5  UMAP-cluster split  (strictest chemistry novelty)
# try:
#     import umap, hdbscan
#     from rdkit.Chem import AllChem, DataStructs

#     print("→ generating UMAP-HDBSCAN clusters ...")

#     # 5-A  Morgan (ECFP4) fingerprints → dense NumPy
#     def ecfp_bits(smiles, n_bits=2048):
#         mol = Chem.MolFromSmiles(smiles or "")
#         fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits) if mol else None
#         arr = np.zeros((1,), dtype=np.uint8)
#         DataStructs.ConvertToNumpyArray(fp, arr)
#         return arr

#     fps = np.vstack(df.SMILES.map(ecfp_bits).tolist())            # (N, 2048)

#     # 5-B  2-D UMAP embedding with Jaccard metric
#     emb = umap.UMAP(metric="jaccard", n_neighbors=50,
#                     min_dist=0.1, random_state=args.seed).fit_transform(fps)

#     # 5-C  Density clustering
#     clusterer = hdbscan.HDBSCAN(min_cluster_size=25,
#                                 metric="euclidean",
#                                 cluster_selection_epsilon=0.5,
#                                 prediction_data=False).fit(emb)
#     df["umap_cluster"] = clusterer.labels_          # -1 = noise

#     # 5-D  Group-K-fold on cluster IDs
#     split_dir = SPLITS / "umap"
#     clusters = df.umap_cluster.unique()
#     RNG.shuffle(clusters)
#     chunks = np.array_split(clusters, 5)

#     for fold, chunk in enumerate(chunks):
#         idx = df[df.umap_cluster.isin(chunk)].index
#         write_fold(df.loc[idx], split_dir, fold)
#     print("✓ umap 5-fold saved to", split_dir)

# except ImportError:

#     # (Fallback) Butina Tanimoto clustering

#     from rdkit import DataStructs
#     from rdkit.Chem import AllChem
#     from rdkit.ML.Cluster import Butina

#     print("UMAP or hdbscan missing → falling back to Butina clustering ...")

#     n_bits = 2048
#     fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi or ""),
#                                                  radius=2, nBits=n_bits)
#            for smi in df.SMILES]

#     # pairwise similarity list required by Butina
#     dists = []
#     nfps = len(fps)
#     for i in range(1, nfps):
#         sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
#         dists.extend([1 - x for x in sims])

#     clusters = Butina.ClusterData(dists, nfps, distThresh=0.3, isDistData=True)
#     # Butina returns a list of tuples of indices
#     cluster_id = np.full(nfps, -1, dtype=int)
#     for cid, tup in enumerate(clusters):
#         cluster_id[list(tup)] = cid
#     df["butina_cluster"] = cluster_id

#     split_dir = splits / "butina"
#     clusters = np.unique(cluster_id)
#     RNG.shuffle(clusters)
#     chunks = np.array_split(clusters, 5)

#     for fold, chunk in enumerate(chunks):
#         idx = df[df.butina_cluster.isin(chunk)].index
#         write_fold(df.loc[idx], split_dir, fold)
#     print("✓ butina 5-fold saved to", split_dir)

# ---------------------------------------------------------------------------

    print("All splits complete.")
