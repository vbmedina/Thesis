import pandas as pd
import numpy as np

path = "/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/chembl.csv"
df = pd.read_csv(path)

# delete rows with invalid SMILES
df = df[df['Smiles'].apply(lambda x: isinstance(x, str) and len(x) > 0)]

df.to_csv(path, index=False)