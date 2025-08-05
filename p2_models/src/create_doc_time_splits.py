
#  Run:
#     conda activate molml
#     python -m src.create_doc_time_splits

from pathlib import Path
import argparse, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedShuffleSplit

# Paths
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TIDY_CSV = DATA / "tidy" / "pp_tidy.csv"
SPLITS = DATA / "splits"

# Constants
DOC_COL = "Document ChEMBL ID"
YEAR_COL = "Year"


# Print distribution of potency bins across train, val, test
def show(name, *frames):
    def cnt(f):
        return (f["potency_bin"].value_counts().reindex(["<=5.0","5.0-6.0","6.0-7.5",">=7.5"], fill_value=0))
    tags = ["train","val","test"]
    print(f"\n[{name}]")
    head = "  set      rows   <=5.0   5.0-6.0   6.0-7.5   >=7.5"
    print(head, "\n  " + "-"*(len(head)-2))
    for tag, df_, ct in zip(tags, frames, map(cnt, frames)):
        print(f"  {tag:<6} {len(df_):7,} {ct['<=5.0']:5} {ct['5.0-6.0']:5} "
              f"{ct['6.0-7.5']:5} {ct['>=7.5']:5}")

# Splirt writing helper
def write_fold(test_idx, out_dir, fold_id, df_full, rng):
    other = df_full.index.difference(test_idx)

    try:
        sss = StratifiedShuffleSplit(1, test_size= 0.10, random_state=int(rng.integers(0, 1e9)))
        (_, v_arr) = next(sss.split(df_full.loc[other], df_full.loc[other,"potency_bin"]))
        val_idx = pd.Index(other[v_arr])
    except ValueError:
        val_idx = pd.Index(rng.choice(other, int(0.10*len(df_full)), False))
    train_idx = other.difference(val_idx)

    fold_dir = out_dir / f"fold{fold_id}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    df_full.loc[train_idx].to_csv(fold_dir/"train.csv", index=False)
    df_full.loc[val_idx  ].to_csv(fold_dir/"val.csv",   index=False)
    df_full.loc[test_idx ].to_csv(fold_dir/"test.csv",  index=False)
    return (df_full.loc[train_idx], df_full.loc[val_idx], df_full.loc[test_idx])

# Main function
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng  = np.random.default_rng(args.seed)

    # Load tidy table
    df = pd.read_csv(TIDY_CSV, low_memory=False)
    df["potency_bin"] = pd.cut(df.pIC50, [-1,5,6,7,99], labels=["<=5.0","5.0-6.0","6.0-7.5",">=7.5"])
   
    if DOC_COL not in df.columns:
        raise KeyError(f"Document column “{DOC_COL}” not found in {TIDY_CSV.name}")

    # -----------------------------------------------------------------
    # 1. Document-group 5-fold
    doc_dir = SPLITS / "doc"
    gkf = GroupKFold(5)
    for fold, (_, test_idx) in enumerate(gkf.split(df, groups=df[DOC_COL])):
        trn, val, tst = write_fold(test_idx, doc_dir, fold, df, rng)
        show(f"doc | fold{fold}", trn, val, tst)
        # QC ≥ 50 actives
        assert (tst["pIC50"] >= 6).sum() >= 50, "doc fold test has < 50 actives"
    print("\n✓ document 5-fold CSVs →", doc_dir)

    # -----------------------------------------------------------------
    # 2. Chronological split  (only if we can get a year)
    # Ensure a year column exists or can be merged
    if YEAR_COL not in df.columns:
        helper = DATA / "doc_year_lookup.csv"
        if helper.exists():
            year_map = pd.read_csv(helper).set_index(DOC_COL)["year"]
            df[YEAR_COL] = df[DOC_COL].map(year_map)
            print(f"[info] merged year from {helper.name}")
            build_time = True
        else:
            print("[warning] no year column and no helper file : skipping time split")
    else:
        build_time = True

    # Build the split if possible
    if build_time:
        cut = 2018
        train_df = df[df[YEAR_COL] <=  cut]
        val_df   = df[df[YEAR_COL] == cut + 1]
        test_df  = df[df[YEAR_COL] >= cut + 2]

        time_dir = SPLITS / "time"
        time_dir.mkdir(parents=True, exist_ok=True)
        train_df.to_csv(time_dir/"train.csv", index=False)
        val_df  .to_csv(time_dir/"val.csv",   index=False)
        test_df .to_csv(time_dir/"test.csv",  index=False)
        show("time split", train_df, val_df, test_df)

        assert (test_df["pIC50"] >= 6).sum() >= 50, "time-split test has < 50 actives"
        print("\n time-split CSVs", time_dir)
    else:
        print(" time split skipped (no year information)")

    print("\n Doc/time split script finished.\n")


#!/usr/bin/env python
# """
# fetch_years.py  –  build doc_year_lookup.csv  (one-off)
# -------------------------------------------------------
# Run from project root:

#     conda activate molml
#     python fetch_years.py
# """
# import time, pathlib, warnings, traceback
# import pandas as pd, tqdm
# from chembl_webresource_client.new_client import new_client
# from requests.exceptions import RetryError, ConnectionError
# from urllib3.exceptions  import HTTPError

# # ─── project paths ─────────────────────────────────────────────────
# PROJECT = pathlib.Path(__file__).resolve().parents[2]            # …/thesis
# TIDY    = PROJECT / "p2_models/data/tidy/pp_tidy.csv"
# OUT_CSV = PROJECT / "p2_models/data/doc_year_lookup.csv"

# DOC_COL = "Document ChEMBL ID"          # header in pp_tidy.csv

# # ─── collect unique document IDs ──────────────────────────────────
# doc_ids = pd.read_csv(TIDY, usecols=[DOC_COL])[DOC_COL].unique()
# docs_api = new_client.document
# records  = []

# # ─── robust fetch loop ────────────────────────────────────────────
# for j, doc in enumerate(tqdm.tqdm(doc_ids, desc="fetching pub years")):
#     year = None
#     try:
#         rec  = docs_api.get(doc)            # single REST call
#         year = rec.get("year", None)        # None if missing
#     except (HTTPError, RetryError, ConnectionError):
#         # too many 404s or transient network issue
#         if j % 50 == 0:                     # log occasionally
#             warnings.warn(f"HTTP error on {doc}; continuing")
#     except Exception as e:
#         warnings.warn(f"unexpected error on {doc}: {e}")
#         traceback.print_exc(limit=1)

#     records.append({DOC_COL: doc, "year": year})

#     if (j + 1) % 100 == 0:                  # polite pause
#         time.sleep(0.5)

# # ─── save lookup table ────────────────────────────────────────────
# OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
# pd.DataFrame(records).to_csv(OUT_CSV, index=False)

# filled = sum(r["year"] is not None for r in records)
# print(f"\nSaved  {OUT_CSV}\n→ {filled:,} / {len(records):,} documents have a year.")