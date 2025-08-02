import csv

# This script converts a semicolon-delimited CSV file to a comma-delimited CSV file.
path= "/Users/victoriamedina/Thesis_Project/other.csv"
dst = "/Users/victoriamedina/Thesis_Project/other_f.csv"

# Ensure the source file exists
with open(path, newline="", encoding="utf-8") as f_in, \
     open(dst, "w", newline="", encoding="utf-8") as f_out:
    # Read the semicolon-delimited file and write it as a comma-delimited file
    reader = csv.reader(f_in, delimiter=";", quotechar='"')
    writer = csv.writer(f_out)                               
    writer.writerows(reader)                                

print(f"✓ Created {dst}")