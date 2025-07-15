import csv
path= "/Users/victoriamedina/Thesis_Project/other.csv"
dst = "/Users/victoriamedina/Thesis_Project/other_f.csv"       # ← new comma-delimited file

with open(path, newline="", encoding="utf-8") as f_in, \
     open(dst, "w", newline="", encoding="utf-8") as f_out:

    reader = csv.reader(f_in, delimiter=";", quotechar='"')  # ← read with semicolons
    writer = csv.writer(f_out)                               # ← write with commas
    writer.writerows(reader)                                 # ← copy rows across

print(f"✓ Created {dst}")