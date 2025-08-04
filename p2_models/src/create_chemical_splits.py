# BEFORE RUNNING -------------------------------
# conda activate molml
# python -m src.create_splits
from pathlib import Path
import argparse, numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles
from sklearn.model_selection import (StratifiedKFold,
                                     StratifiedShuffleSplit,
                                     GroupKFold)

# Paths
ROOT   = Path(__file__).resolve().parents[1]
DATA   = ROOT / "data"
TIDY   = DATA / "tidy"
SPLITS = DATA / "splits"

# Helper function to get Murcko scaffold
def murcko(smiles: str) -> str:
    if not smiles:
        return ""
    try:
        return MurckoScaffoldSmilesFromSmiles(smiles, includeChirality=False)
    except Exception:
        return ""

# Show distribution of potency bins across train, val, test
def show_distribution(name, trn, val, tst):
    def c(d):
        return (d["potency_bin"].value_counts()
                  .reindex(["≤5","5-6","6-7","≥7"], fill_value=0))
    ct, cv, cs = map(c, (trn, val, tst))
    print(f"\n[{name}]")
    print("  set      rows   ≤5   5-6   6-7   ≥7")
    print("  -------------------------------------")
    print(f"  train {len(trn):7,} {ct['≤5']:5} {ct['5-6']:5} {ct['6-7']:5} {ct['≥7']:5}")
    print(f"  val   {len(val):7,} {cv['≤5']:5} {cv['5-6']:5} {cv['6-7']:5} {cv['≥7']:5}")
    print(f"  test  {len(tst):7,} {cs['≤5']:5} {cs['5-6']:5} {cs['6-7']:5} {cs['≥7']:5}")

# Write a fold to CSV
def write_fold(df_test, split_dir, fold_id, df_full, rng):
    tst_idx  = df_test.index
    oth_idx  = df_full.index.difference(tst_idx)
    if oth_idx.empty:
        raise ValueError("Fold has no train/val rows!")

    # 10 % validation - stratified by potency_bin
    try:
        sss = StratifiedShuffleSplit(
        n_splits=1,
        test_size=0.10,
        random_state=int(rng.integers(0, 1_000_000_000))
    )
        (_, v_arr) = next(
            sss.split(
                df_full.loc[oth_idx],
                df_full.loc[oth_idx, "potency_bin"]
            )
        )
        val_idx = pd.Index(oth_idx[v_arr])
    except ValueError:
        val_idx = pd.Index(
            rng.choice(oth_idx, int(0.10 * len(df_full)), replace=False)
        )

    # Remaining rows are training
    trn_idx = oth_idx.difference(val_idx)

    # Create directory for this fold
    fold_dir = split_dir / f"fold{fold_id}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    df_full.loc[trn_idx].to_csv(fold_dir/"train.csv", index=False)
    df_full.loc[val_idx ].to_csv(fold_dir/"val.csv",   index=False)
    df_full.loc[tst_idx ].to_csv(fold_dir/"test.csv",  index=False)

    # Show distribution
    return df_full.loc[trn_idx], df_full.loc[val_idx], df_full.loc[tst_idx]

# Main function
if __name__ == "__main__":
    ap  = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0, help="random seed")
    args = ap.parse_args()
    rng  = np.random.default_rng(args.seed)

    # Load tidy table
    df = pd.read_csv(TIDY / "pp_tidy.csv", low_memory=False)
    print(f"Loaded {len(df):,} rows")

    # Helper columns
    df["potency_bin"] = pd.cut(df.pIC50, [-1,5,6,7,99],
                               labels=["≤5","5-6","6-7","≥7"])
    freq = df.groupby("Molecule_ChEMBL_ID").size()
    df["mol_freq"] = df["Molecule_ChEMBL_ID"].map(freq)
    df["scaffold"] = df["Smiles"].map(murcko)

# # ----------------------------------------------------------------
# # Comment out the sections you run  5 ─ UMAP-HDBSCAN or Butina
#     # 3 ─ random-stratified
#     strat = (df["Strains"] + "_" +
#              df.potency_bin.astype(str) + "_" +
#              pd.qcut(df.mol_freq,[0,.8,.95,1],
#                      labels=["common","frequent","hyper"]).astype(str))
#     skf = StratifiedKFold(5, shuffle=True, random_state=args.seed)
#     rnd_dir = SPLITS/"random"
#     for f, (_, tst) in enumerate(skf.split(df, strat)):
#         trn, val, tst = write_fold(df.loc[tst], rnd_dir, f, df, rng)
#         show_distribution(f"random | fold{f}", trn, val, tst)
#     print("\n✓ random 5-fold CSVs →", rnd_dir)

#     # 4 ─ scaffold with GroupKFold
#     gkf = GroupKFold(5)
#     scaf_dir = SPLITS/"scaffold"
#     for f, (_, tst) in enumerate(gkf.split(df, groups=df["scaffold"])):
#         trn, val, tst = write_fold(df.loc[tst], scaf_dir, f, df, rng)
#         show_distribution(f"scaffold | fold{f}", trn, val, tst)
#     print("\n✓ scaffold 5-fold CSVs →", scaf_dir)
# # -----------------------------------------------------------------

# # Run the following seperately if you have UMAP and HDBSCAN installed
# # Otherwise, it will fall back to Butina clustering
    # 5 ─ UMAP-HDBSCAN or Butina
    try:
        import umap, hdbscan
        from rdkit.Chem import AllChem, DataStructs

# UMAP + HDBSCAN clustering
        print("\n→ generating UMAP-HDBSCAN clusters …")
        def ecfp_bits(smi, bits=2048):
            mol = Chem.MolFromSmiles(smi or "")
            fp  = AllChem.GetMorganFingerprintAsBitVect(mol, 2, bits) if mol else None
            arr = np.zeros((1,), np.uint8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            return arr

    # Generate ECFP fingerprints and UMAP embedding
        fps = np.vstack(df["Smiles"].map(ecfp_bits))
        emb = umap.UMAP(metric="jaccard", n_neighbors=50,
                        min_dist=0.1, random_state=args.seed).fit_transform(fps)
        df["umap_cluster"] = hdbscan.HDBSCAN(min_cluster_size=25).fit(emb).labels_

    # Write UMAP clusters to CSV
        umap_dir = SPLITS / "umap"
        for f, clusters in enumerate(np.array_split(df.umap_cluster.unique(), 5)):
            idx = df[df.umap_cluster.isin(clusters)].index
            trn, val, tst = write_fold(df.loc[idx], umap_dir, f, df, rng)
            show_distribution(f"umap | fold{f}", trn, val, tst)
        print("\n✓ umap 5-fold CSVs →", umap_dir)

    # Fallback to Butina clustering if UMAP or HDBSCAN is not available
    except ImportError:
        from rdkit import DataStructs
        from rdkit.Chem import AllChem
        from rdkit.ML.Cluster import Butina

        print("\nUMAP/HDBSCAN missing – Butina clustering …")
        fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s or ""),
                                                     2, 2048) for s in df["Smiles"]]
    # Butina clustering
        dists=[]
        for i in range(1,len(fps)):
            sims=DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend([1-x for x in sims])
        clusters = Butina.ClusterData(dists, len(fps), 0.3, True)
        cid = np.full(len(fps), -1, int)
        for n,tup in enumerate(clusters): cid[list(tup)] = n
        df["butina_cluster"]=cid
   
    # Write Butina clusters to CSV
        but_dir = SPLITS / "butina"
        for f, cls in enumerate(np.array_split(df.butina_cluster.unique(), 5)):
            idx = df[df.butina_cluster.isin(cls)].index
            trn, val, tst = write_fold(df.loc[idx], but_dir, f, df, rng)
            show_distribution(f"butina | fold{f}", trn, val, tst)
        print("\n✓ butina 5-fold CSVs →", but_dir)
# ---------------------------------------------------------------------------

    print("\nAll splits complete.")
