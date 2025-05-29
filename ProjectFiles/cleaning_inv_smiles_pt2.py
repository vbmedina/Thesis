import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# ── CONFIGURE PATHS ────────────────────────────────────────────
IN_CSV  = Path("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv")
OUT_CSV = IN_CSV.with_name("chembl_with_scaffolds_only.csv")

# ── LOAD DATA ───────────────────────────────────────────────────
df = pd.read_csv(IN_CSV, low_memory=False)

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
filtered.to_csv(OUT_CSV, index=False)

# ── REPORT ──────────────────────────────────────────────────────
print(f"Total rows before: {len(df)}")
print(f"Rows with a valid Murcko scaffold: {len(filtered)}")
print(f"Rows dropped (acyclic/no scaffold): {len(df) - len(filtered)}")
print(f"Filtered dataset written to: {OUT_CSV}")