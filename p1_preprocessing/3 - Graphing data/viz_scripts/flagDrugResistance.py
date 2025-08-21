import pandas as pd

# df = pd.read_csv("./pp.csv")

# def flag_drug_resistance(description):
#     description = description.lower()
#     description = description.replace("-", " ")
#     description = description.replace(",", " ")
#     words = description.split()

#     indexAgainst = words.index("against") if "against" in words else -1
#     indexActivity = words.index("activity") if "activity" in words else -1
#     indexStart = max(indexAgainst, indexActivity)


#     indexSensitive = words.index("sensitive") if "sensitive" in words else -1
#     indexResistant = words.index("resistant") if "resistant" in words else -1
#     indexSusceptible = words.index("susceptible") if "susceptible" in words else -1
#     indexEnd = max(indexSensitive, indexResistant, indexSusceptible)     

#     if indexEnd == -1:
#         return ""

#     return " ".join(words[indexStart+1:indexEnd+1])

# df['Flag resistance'] = df['Assay Description'].apply(flag_drug_resistance)

# unique_flags = df['Flag resistance'].unique()
# print("Unique flags:", len(unique_flags))

# # Save the updated DataFrame to a new CSV file
# output_path = "./pp_with_resistance_scraping.csv"
# df.to_csv(output_path, index=False)

# unique_output_path = "./unique_flags.csv"
# unique_flags_df = pd.DataFrame(unique_flags, columns=['Unique Flags'])
# unique_flags_df.to_csv(unique_output_path, index=False)

# Load the csv
csv_path = pd.read_csv("./pp_with_resistance_scraping.csv")
flags = pd.read_csv("./unique_flags.csv")

# Map
mapping= dict(zip(flags["Unique Flags"], flags["Flag Groups"]))

# Map the flags to their groups
csv_path['Flag Groups'] = csv_path['Flag resistance'].map(mapping)

# Save the updated DataFrame with flag groups to a new CSV file
csv_path.to_csv("./pp_with_flag_groups.csv", index=False)

print("File saved as 'pp_with_flag_groups_2.csv'")