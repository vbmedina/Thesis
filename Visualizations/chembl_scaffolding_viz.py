from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def generate_scaffold(mol):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
  
scaffolds = {}
for _, row in data.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    try:
        scaffold = generate_scaffold(mol)
    except:
        scaffold = None
    if scaffold not in scaffolds:
        scaffolds[scaffold] = []
    scaffolds[scaffold].append(row)

scaffold_sets = list(scaffolds.values())