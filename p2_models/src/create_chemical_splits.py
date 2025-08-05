# BEFORE RUNNING:
# conda activate molml
# python -m src.create_chemical_splits
from pathlib import Path
import argparse, numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles
from sklearn.model_selection import StratifiedShuffleSplit

# Paths
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TIDY = DATA / "tidy"
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
                  .reindex(["<=5.0","5.0-6.0","6.0-7.5",">=7.5"], fill_value=0))
    ct, cv, cs = map(c, (trn, val, tst))
    print(f"\n[{name}]")
    print("  set      rows   <=5.0   5.0-6.0   6.0-7.5   >=7.5")
    print("  -------------------------------------")
    print(f"  train {len(trn):7,} {ct['<=5.0']:5} {ct['5.0-6.0']:5} {ct['6.0-7.5']:5} {ct['>=7.5']:5}")
    print(f"  val   {len(val):7,} {cv['<=5.0']:5} {cv['5.0-6.0']:5} {cv['6.0-7.5']:5} {cv['>=7.5']:5}")
    print(f"  test  {len(tst):7,} {cs['<=5.0']:5} {cs['5.0-6.0']:5} {cs['6.0-7.5']:5} {cs['>=7.5']:5}")

# Write a fold to CSV
def write_fold(df_test, split_dir, fold_id, df_full, rng):
    tst_idx  = df_test.index
    oth_idx  = df_full.index.difference(tst_idx)
    if oth_idx.empty:
        raise ValueError("Fold has no train/val rows!")

    # 10 % validation - stratified by potency_bin
    try:
        sss = StratifiedShuffleSplit(
            n_splits=1, test_size=0.10,
            random_state=int(rng.integers(0, 1_000_000_000))
        )
        (_, v_arr) = next(
            sss.split(df_full.loc[oth_idx],
                      df_full.loc[oth_idx, "potency_bin"])
        )
        val_idx = pd.Index(oth_idx[v_arr])
    except ValueError:
        val_idx = pd.Index(
            rng.choice(oth_idx, int(0.10 * len(df_full)), replace=False)
        )

    trn_idx = oth_idx.difference(val_idx)

    fold_dir = split_dir / f"fold{fold_id}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    df_full.loc[trn_idx].to_csv(fold_dir/"train.csv", index=False)
    df_full.loc[val_idx].to_csv(fold_dir/"val.csv",   index=False)
    df_full.loc[tst_idx].to_csv(fold_dir/"test.csv",  index=False)

    return df_full.loc[trn_idx], df_full.loc[val_idx], df_full.loc[tst_idx]

# Main function to create splits
if __name__ == "__main__":
    SEED = 0
    rng  = np.random.default_rng(SEED)

    df = pd.read_csv(TIDY / "pp_tidy.csv", dtype={13: "string"})
    print(f"Loaded {len(df):,} rows")

    df["potency_bin"] = pd.cut(df.pIC50, [-1,5,6,7.5,99],
                               labels=["<=5.0","5.0-6.0","6.0-7.5",">=7.5"])
    df["scaffold"] = df["Smiles"].map(murcko)
# # ----------------------------------------------------------------
# # 1. random-stratified baseline – even test rows, no grouping
#     rows = np.arange(len(df))
#     folds = np.array_split(rows, 5)          # five equal-size slices

#     rnd_dir = SPLITS / "random"
#     for f, idx in enumerate(folds):
#         trn, val, tst = write_fold(df.loc[idx], rnd_dir, f, df, rng)
#         show_distribution(f"random | fold{f}", trn, val, tst)
#     print("\n random 5-fold CSVs", rnd_dir)

# # ---------------------------------------------------------------------
# # 2. Scaffold split
#     scaf_sizes = df.groupby("scaffold").size().sort_values(ascending=False)
#     scaf_bins  = {i: [] for i in range(5)}
#     scaf_rows  = [0]*5
#     for scaf, n in scaf_sizes.items():
#         fid = scaf_rows.index(min(scaf_rows))
#         scaf_bins[fid].append(scaf)
#         scaf_rows[fid] += n

#     scaf_dir = SPLITS / "scaffold"
#     for f in range(5):
#         idx = df[df["scaffold"].isin(scaf_bins[f])].index
#         trn, val, tst = write_fold(df.loc[idx], scaf_dir, f, df, rng)
#         show_distribution(f"scaffold | fold{f}", trn, val, tst)
#     print("\n scaffold 5-fold CSVs", scaf_dir)


# ------------------------------------------------------------------
# 3. Chemistry-cluster split (UMAP-HDBSCAN or Butina)
    try:
        import umap,hdbscan
        from rdkit.Chem import AllChem,DataStructs

        print("\n→ UMAP-HDBSCAN clustering …")
        def ecfp_bits(smi,bits=2048):
            mol = Chem.MolFromSmiles(smi or "")
            fp  = AllChem.GetMorganFingerprintAsBitVect(mol,2,bits) if mol else None
            arr = np.zeros((1,),np.uint8); DataStructs.ConvertToNumpyArray(fp,arr)
            return arr

        fps  = np.vstack(df["Smiles"].map(ecfp_bits))
        emb  = umap.UMAP(metric="jaccard",n_neighbors=50,
                         min_dist=0.1,random_state=SEED).fit_transform(fps)
        df["cur_cluster"] = hdbscan.HDBSCAN(min_cluster_size=25)\
                               .fit(emb).labels_.astype(str)

    except ImportError:
        from rdkit import DataStructs,AllChem
        from rdkit.ML.Cluster import Butina

        print("\nUMAP/HDBSCAN missing – Butina clustering …")
        fps=[AllChem.GetMorganFingerprintAsBitVect(
             Chem.MolFromSmiles(s or ""),2,2048) for s in df["Smiles"]]
        dists=[]; 
        for i in range(1,len(fps)):
            dists.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])])
        clusters = Butina.ClusterData(dists,len(fps),0.3,True)
        df["cur_cluster"] = pd.Series(clusters,dtype="int").astype(str)

# Distrbute noise
    noise_idx  = df[df.cur_cluster == "-1"].index
    df.loc[noise_idx,"cur_cluster"] = "-1"
    noise_slices = np.array_split(noise_idx,5)

    bins, rows = {i:set() for i in range(5)}, [0]*5
    for f in range(5):
        nid = f"noise_{f}"
        df.loc[noise_slices[f],"cur_cluster"] = nid
        bins[f].add(nid); rows[f] += len(noise_slices[f])

# Greedy balancing
    clus_sizes = (df.cur_cluster.value_counts()
                  .drop([f"noise_{i}" for i in range(5)])).to_dict()

    for cid,n in sorted(clus_sizes.items(), key=lambda kv:-int(kv[0])):
        fid = rows.index(min(rows))
        bins[fid].add(cid); rows[fid] += n

# X^2 balancing
    global_cnt = df["potency_bin"].value_counts()
    def chi2(cnt):
        return sum((cnt.get(c,0)-global_cnt[c]*cnt.sum()/len(df))**2 /
                   (global_cnt[c]*cnt.sum()/len(df)) for c in global_cnt.index)

    for _ in range(1000):
        a,b = rng.integers(0,5,size=2)
        if not bins[a] or not bins[b]: continue
        ca,cb = rng.choice(list(bins[a])), rng.choice(list(bins[b]))
        if ca==cb: continue
        bins[a].remove(ca); bins[b].remove(cb)
        bins[a].add(cb);    bins[b].add(ca)

        before=chi2(df.loc[df.cur_cluster.isin(bins[a]),"potency_bin"].value_counts())+\
               chi2(df.loc[df.cur_cluster.isin(bins[b]),"potency_bin"].value_counts())
        after =chi2(df.loc[df.cur_cluster.isin(bins[a]),"potency_bin"].value_counts())+\
               chi2(df.loc[df.cur_cluster.isin(bins[b]),"potency_bin"].value_counts())
        if after>=before:
            bins[a].remove(cb); bins[b].remove(ca)
            bins[a].add(ca);   bins[b].add(cb)

# Row Count Swap
target = len(df) // 5
for _ in range(2000):
    big = rows.index(max(rows))
    small = rows.index(min(rows))
    if rows[big] - rows[small] <= 50:
        break

    # choose the smallest cluster in the big fold
    move = min((cid for cid in bins[big]),
               key=lambda cid: df.cur_cluster.eq(cid).sum())

    bins[big].remove(move); bins[small].add(move)
    n = df.cur_cluster.eq(move).sum()
    rows[big]  -= n
    rows[small]+= n
    
# Writing folds
    out_dir = SPLITS/"umap"
    for f in range(5):
        idx = df[df.cur_cluster.isin(bins[f])].index
        trn,val,tst = write_fold(df.loc[idx], out_dir, f, df, rng)
        show_distribution(f"umap | fold{f}", trn,val,tst)

    print("\n UMAP 5-fold CSVs", out_dir)

# ---------------------------------------------------------------------
    print("\n All splits complete.")