import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Paths
inpath  = Path("./Do Not Touch/postphase4_validSmiles.csv")
outpath = inpath.with_name("postphase5.csv")

# ── LOAD DATA ───────────────────────────────────────────────────
df = pd.read_csv(inpath, low_memory=False)

# ── GENERATE SCAFFOLDS ─────────────────────────────────────────
def smiles_to_scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    return ""

df["Scaffold"] = df["Smiles"].apply(smiles_to_scaffold)

# ── FILTER TO NON-EMPTY SCAFFOLDS ───────────────────────────────
mask = df["Scaffold"].str.strip().astype(bool)
filtered = df[mask].reset_index(drop=True)

# ── SAVE CLEANED DATA ───────────────────────────────────────────
filtered.to_csv(outpath, index=False)

# ── REPORT ──────────────────────────────────────────────────────
print(f"Total rows before: {len(df)}")
print(f"Rows with a valid Murcko scaffold: {len(filtered)}")
print(f"Rows dropped (acyclic/no scaffold): {len(df) - len(filtered)}")
print(f"Filtered dataset written to: {outpath}")