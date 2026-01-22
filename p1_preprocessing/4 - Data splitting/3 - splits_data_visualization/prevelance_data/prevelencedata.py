#!/usr/bin/env python3
"""
Compute prevalence of actives (pIC50 ≥ 6) for each data split and subset.

Input tree (example):
  p1_preprocessing/4 - Data splitting/2 - split_data/
    ├─ butina/
    │   ├─ butina_fold_1_train.csv
    │   ├─ butina_fold_1_val.csv
    │   └─ butina_fold_1_test.csv
    ├─ random/
    ├─ scaffold/
    ├─ umap_kmeans/
    └─ umap_ward/

Output (single CSV):
  p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/prevelance_data/prevelancedata.csv
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

# Configuring paths
SPLIT_ROOT = Path("p1_preprocessing/4 - Data splitting/2 - split_data")
OUT_CSV = Path("p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/prevelance_data/prevelancedata.csv")

# Activity definition
PIC50_THRESHOLD = 6.0
PIC50_CANDIDATES = ["pIC50", "pic50", "p_ic50", "PIC50"]

# Preferred reporting order & pretty labels
SPLIT_ORDER = ["random", "scaffold", "butina", "umap_kmeans", "umap_ward"]
SPLIT_DISPLAY = {
    "random": "Random (stratified)",
    "scaffold": "Scaffold",
    "butina": "Butina",
    "umap_kmeans": "UMAP (k-means)",
    "umap_ward": "UMAP (ward)",
}
SUBSET_ORDER = ["train", "val", "test"]

# Parse names like "<split>_fold_<n>_<subset>.csv"
FILE_RE = re.compile(r"(?P<split>[^/\\]+)_fold_(?P<fold>\d+)_(?P<subset>train|val|test)\.csv$", re.IGNORECASE)

# Map various tokens to canonical split keys above
EXACT_SPLIT_MAP = {
    "random": "random",
    "random_stratified": "random",
    "scaffold": "scaffold",
    "butina": "butina",
    "umap_kmeans": "umap_kmeans",
    "umap-kmeans": "umap_kmeans",
    "umap_k-means": "umap_kmeans",
    "umap_ward": "umap_ward",
}

def canon_split(token: str) -> str | None:
    s = str(token).strip().lower().replace(" ", "_")
    s = s.replace("-", "_")
    return EXACT_SPLIT_MAP.get(s, s if s in EXACT_SPLIT_MAP.values() else None)

def find_pic50_col(df: pd.DataFrame) -> str | None:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in PIC50_CANDIDATES:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None

def prevalence_from_pic50(series: pd.Series, thr: float = 6.0) -> dict:
    s = pd.to_numeric(series, errors="coerce")
    mask = s.notna()
    total = int(mask.sum())
    if total == 0:
        return {"n_total": 0, "n_active": 0, "n_inactive": 0, "prevalence_active": np.nan}
    actives = int((s[mask] >= thr).sum())
    inactives = total - actives
    return {
        "n_total": total,
        "n_active": actives,
        "n_inactive": inactives,
        "prevalence_active": actives / total,
    }

def main():
    SPLIT_ROOT.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for csv_path in SPLIT_ROOT.rglob("*.csv"):
        m = FILE_RE.search(csv_path.name)
        if not m:
            continue

        raw_split = m.group("split")
        split_key = canon_split(raw_split)
        subset = m.group("subset").lower()
        fold = int(m.group("fold"))

        if split_key is None:
            # skip unknown naming
            continue

        try:
            df = pd.read_csv(csv_path)
            col = find_pic50_col(df)
            if col is None:
                raise ValueError("No pIC50 column found.")
            stats = prevalence_from_pic50(df[col], thr=PIC50_THRESHOLD)
            status = "ok"
        except Exception as e:
            stats = {"n_total": np.nan, "n_active": np.nan, "n_inactive": np.nan, "prevalence_active": np.nan}
            status = f"ERROR: {e}"

        rows.append({
            "split_key": split_key,
            "subset": subset,
            "fold": fold,
            "file": str(csv_path),
            "status": status,
            **stats
        })

    if not rows:
        raise SystemExit(f"No matching CSVs found under {SPLIT_ROOT}")

    by_file = pd.DataFrame(rows)

    # Keep successful rows
    ok = by_file[by_file["status"] == "ok"].copy()
    if ok.empty:
        raise SystemExit("No valid rows to summarize (all failed).")

    # Pooled counts per split × subset, plus mean±SD across folds
    pooled = (ok.groupby(["split_key", "subset"], as_index=False)
                .agg(n_total=("n_total", "sum"),
                     n_active=("n_active", "sum"),
                     n_inactive=("n_inactive", "sum"),
                     prevalence_mean=("prevalence_active", "mean"),
                     prevalence_sd=("prevalence_active", "std"),
                     folds=("fold", "nunique")))

    pooled["prevalence_active_pooled"] = pooled["n_active"] / pooled["n_total"]

    # Add pretty labels + sort
    pooled["split"] = pooled["split_key"].map(SPLIT_DISPLAY).fillna(pooled["split_key"])
    pooled["subset"] = pd.Categorical(pooled["subset"], categories=SUBSET_ORDER, ordered=True)
    pooled["split_key"] = pd.Categorical(pooled["split_key"], categories=SPLIT_ORDER, ordered=True)
    pooled = pooled.sort_values(["split_key", "subset"]).reset_index(drop=True)

    # Final column order
    cols = [
        "split_key", "split", "subset", "folds",
        "n_total", "n_active", "n_inactive",
        "prevalence_active_pooled", "prevalence_mean", "prevalence_sd",
    ]
    pooled = pooled[cols]

    pooled.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV.resolve()}")

if __name__ == "__main__":
    main()
