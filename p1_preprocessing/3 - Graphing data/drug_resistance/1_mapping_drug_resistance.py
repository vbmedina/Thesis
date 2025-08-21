''' Description: This script is used to scrape drug resistance flags from the 'Assay Description' for 
visualizations and analysis.

Requirements:
1) final_data.csv - from Step 2 if pipeline.
2) "mapped_resistance_flags.csv" - A pre-made mapping of flags to their groups. If using a different dataset,
or updated dataset this file will need to be updated.
'''

import pandas as pd

''' Section 1. Create personalized strain scraping CSV file with drug resistance flags.
1. Load the CSV file containing all chemical data.
2. We will scrape the 'Assay Description' column for drug resistance flags.
3. Create a new (temporary) column 'Flag resistance' to store the extracted flags.'''

# Load the csv
df = pd.read_csv("./p1_preprocessing/3 - Graphing data/drug_resistance/final_data_copy.csv")
def flag_drug_resistance(description):
    description = description.lower()
    description = description.replace("-", " ")
    description = description.replace(",", " ")
    words = description.split()

    indexAgainst = words.index("against") if "against" in words else -1
    indexActivity = words.index("activity") if "activity" in words else -1
    indexStart = max(indexAgainst, indexActivity)


    indexSensitive = words.index("sensitive") if "sensitive" in words else -1
    indexResistant = words.index("resistant") if "resistant" in words else -1
    indexSusceptible = words.index("susceptible") if "susceptible" in words else -1
    indexEnd = max(indexSensitive, indexResistant, indexSusceptible)     

    if indexEnd == -1:
        return ""

    return " ".join(words[indexStart+1:indexEnd+1])

df['Flag resistance'] = df['Assay Description'].apply(flag_drug_resistance)

unique_flags = df['Flag resistance'].unique()
print("Unique flags:", len(unique_flags))


''' Section 2. Used pre-made mapping in "mapped_resistance_flags.csv" to map the flags to their groups.'''
# Load the mapping CSV
flags = pd.read_csv("./p1_preprocessing/3 - Graphing data/drug_resistance/mapped_resistance_flags.csv")

# Map
mapping= dict(zip(flags["Unique Flags"], flags["Flag Groups"]))

# Map the flags to their groups
df['Flag Groups'] = df['Flag resistance'].map(mapping)

df = df.drop(columns=['Flag resistance'])

# Save the updated DataFrame with flag groups to a new CSV file
df.to_csv("./p1_preprocessing/3 - Graphing data/drug_resistance/final_data_with_flags.csv", index=False)

print("File saved as 'csv_with_flags.csv'")