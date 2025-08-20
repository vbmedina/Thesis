''' Description: This script is the first step to the pipeline (strains). Starting with the CSV "prephase1.csv", which
was made using scripts from part 0 and manually, has a new column called "stand_strain". This has all the raw strains
and isolates before binning them into major strains and isolates based on parental similarity. The use of this script uses
pre-generated files list_of_strains and list_of_isolates manually made to map out associated strain or isolate. '''

# Imports
from pathlib import Path
import pandas as pd
import numpy as np

# Paths
MAPPING_CSV_STRAINS = Path(".Thesis_Project/thesis/p0_all_csvs//list_of_strains.csv")
MAPPING_CSV_ISOLATES = Path(".Thesis_Project/thesis/p0_all_csvs/list_of_isolates.csv")
IN_CSV = Path(".Thesis_Project/thesis/p0_all_csvs//prephase1.csv")

# Read in the CSV files
df = pd.read_csv(IN_CSV)
df_mapping_strains = pd.read_csv(MAPPING_CSV_STRAINS, encoding="latin1")
df_mapping_isolates = pd.read_csv(MAPPING_CSV_ISOLATES)

# Separate the alternative strain names by ;
df_mapping_strains["Alternative Names"] = df_mapping_strains["Alternative Names"].str.split(";")
df_mapping_strains = df_mapping_strains.explode("Alternative Names")
df_mapping_strains["Alternative Names"] = df_mapping_strains["Alternative Names"].str.strip()

# Separate the alternative isolate names by ;
df_mapping_isolates["Alternative Names"] = df_mapping_isolates["Alternative Names"].str.split(";")
df_mapping_isolates = df_mapping_isolates.explode("Alternative Names")
df_mapping_isolates["Alternative Names"] = df_mapping_isolates["Alternative Names"].str.strip()

# Filter out rows where "stand_strain" is NaN
print(f"{df['stand_strain'].isna().sum()} removed.")
df = df[~df["stand_strain"].isna()]

# Create the mapping dictionary from df_mapping_strains
strain_map = df_mapping_strains.set_index("Alternative Names")["Strains"].str.strip().to_dict()

# Add the isolates to the mapping dictionary
for isolate, altIsolate in zip(df_mapping_isolates["Suspected Isolates/Clones"].str.strip(), df_mapping_isolates["Alternative Names"].str.strip()):
    strain_map[isolate] = isolate
    if str(altIsolate) != "nan" and altIsolate != "":
        strain_map[altIsolate] = isolate

for strain in df_mapping_strains["Strains"].str.strip():
    strain_map[strain] = strain

# Map "stand_strain" in df to "Strain" using the mapping
df["stand_strain_bin"] = df["stand_strain"].str.strip().map(strain_map)

# print(df.loc[df["stand_strain_bin"].isna(), "stand_strain"].unique())
print(f"{df['stand_strain_bin'].isna().sum()} removed.")
df = df[~df["stand_strain_bin"].isna()]

df.to_csv(IN_CSV.with_name("postphase1_strainsMerged.csv"), index=False)