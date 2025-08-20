# When downloading csv files from ChEMBL, they may be semicolon-delimited.
# This script converts them to comma-delimited CSV files.

# Import
import csv

# Paths
path= "/Users/victoriamedina/Thesis_Project/other.csv"
new = "/Users/victoriamedina/Thesis_Project/other_f.csv"

# Ensure the source file exists
with open(path, newline="", encoding="utf-8") as f_in, \
    open(new, "w", newline="", encoding="utf-8") as f_out:
    reader = csv.reader(f_in, delimiter=";", quotechar='"')
    writer = csv.writer(f_out)                               
    writer.writerows(reader)                                

print(f"Created comma-delimited CSV file: {new}")