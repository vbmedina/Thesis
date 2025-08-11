import os, sys, math, glob, argparse
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
import matplotlib.pyplot as plt

RDLogger.DisableLog("rdApp.*")

def find_file("/Users/victoriamedina/Thesis_Project/thesis/p2_models/data/splits", method, fold, split):
    # robust to your naming (train has "_fold_#_", val/test have "_fold#_")
    pat = os.path.join(dir_path, f"{method}_fold*{fold}*_{split}.csv")
    matches = sorted(glob.glob(pat))
    if not matches:
        raise FileNotFoundError(f"Missing file for {method} fold {fold} {split}: {pat}")
    return matches[0]

def smiles_col_in(df, user_col):
    if user_col and user_col != "auto":
        if user_col not in df.columns:
            raise ValueError(f"SMILES column '{user_col}' not in {list(df.columns)}")
        return user_col
    for c in ["smiles","SMILES","canonical_smiles","Canonical_SMILES"]:
        if c in df.columns: return c
    raise ValueError("Could not find a SMILES column; pass --smiles_col")

def load_smiles(csv_path, smiles_col):
    df = pd.read_csv(csv_path)
    col = smiles_col_in(df, smiles_col)
    return df[col].astype(str).tolist()

def to_mols(smiles):
    mols, bad = [], 0
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            mols.append(None); bad += 1
        else:
            mols.append(m)
    if bad: print(f"[warn] {bad} invalid SMILES in file", file=sys.stderr)
    return mols

def to_fps(mols, n_bits=2048, radius=2):
    fps = []
    for m in mols:
        fps.append(None if m is None else AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits))
    return fps

def max_to_train(query_fps, train_fps, exclude_self=False):
    keep_train = [t for t in train_fps if t is not None]
    out = []
    for i, fp in enumerate(query_fps):
        if fp is None:
            out.append(float("nan")); continue
        if not exclude_self:
            sims = DataStructs.BulkTanimotoSimilarity(fp, keep_train)
            out.append(max(sims) if sims else float("nan"))
        else:
            # self vs full list; zero self if it is 1.0
            sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
            if i < len(sims) and not math.isnan(sims[i]) and sims[i] == 1.0:
                sims[i] = 0.0
            out.append(max([s for s in sims if not math.isnan(s)]))
    return out

def violin(stats, out_png):
    labels = list(stats.keys())
    data = [[v for v in stats[k] if pd.notna(v)] for k in labels]
    fig = plt.figure(figsize=(6,4.5), dpi=150)
    plt.violinplot(data, showmeans=True)
    plt.xticks(range(1, len(labels)+1), labels)
    plt.ylabel("Max Tanimoto to training set")
    plt.ylim(0,1.0)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[plot] {out_png}")

def run(base_dir, smiles_col, methods):
    all_rows = []
    for method in methods:
        mdir = os.path.join(base_dir, method)
        for fold in range(1, 6):
            train_csv = find_file(mdir, method, fold, "train")
            val_csv   = find_file(mdir, method, fold, "val")
            test_csv  = find_file(mdir, method, fold, "test")

            tr = load_smiles(train_csv, smiles_col)
            va = load_smiles(val_csv,   smiles_col)
            te = load_smiles(test_csv,  smiles_col)

            tr_fps = to_fps(to_mols(tr))
            va_fps = to_fps(to_mols(va))
            te_fps = to_fps(to_mols(te))

            tr_max = max_to_train(tr_fps, tr_fps, exclude_self=True)
            va_max = max_to_train(va_fps, tr_fps, exclude_self=False)
            te_max = max_to_train(te_fps, tr_fps, exclude_self=False)

            # save CSV per fold + accumulate summary
            out_dir_csv = os.path.join("reports","tanimoto",method)
            os.makedirs(out_dir_csv, exist_ok=True)
            df = pd.DataFrame({
                "split": ["train"]*len(tr_max) + ["val"]*len(va_max) + ["test"]*len(te_max),
                "max_sim": tr_max + va_max + te_max
            })
            out_csv = os.path.join(out_dir_csv, f"{method}_fold{fold}_max_tanimoto.csv")
            df.to_csv(out_csv, index=False); print(f"[csv]  {out_csv}")

            for s, arr in (("train", tr_max), ("val", va_max), ("test", te_max)):
                all_rows += [{"method":method,"fold":fold,"split":s,"max_sim":v} for v in arr]

            # plot violin per fold
            out_png = os.path.join("figures","tanimoto",method,f"{method}_fold{fold}_violin.png")
            violin({"train":tr_max, "val":va_max, "test":te_max}, out_png)

    # one big summary file (optional)
    out_all = os.path.join("reports","tanimoto","_summary_all_folds.csv")
    pd.DataFrame(all_rows).to_csv(out_all, index=False); print(f"[csv]  {out_all}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--base_dir", default="p2_models/data/splits",
                    help="folder that contains butina/random/random_sim/scaffold/umap")
    ap.add_argument("--smiles_col", default="auto", help="SMILES column name (or 'auto')")
    ap.add_argument("--methods", nargs="*", default=["butina","random","random_sim","scaffold","umap"])
    args = ap.parse_args()
    run(args.base_dir, args.smiles_col, args.methods)