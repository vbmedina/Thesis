''' 
Description: This script is not needed for pipeline. Only use when freshly downloading CSV files from ChEMBL. Data directly
downloaded from ChEMBL as a CSV is automatically downloaded with semicolon-delimiters. This script convert the CSV to a 
comma-delimited file.
''' 

# Import
import csv

# Paths
path= "./Thesis_Project/thesis/p0_all_csvs/insert_downloaded_csv_here.csv"
new = "./Thesis_Project/thesis/p0_all_csvs/chembl_data_original.csv"

# Ensure the source file exists
with open(path, newline="", encoding="utf-8") as f_in, \
    open(new, "w", newline="", encoding="utf-8") as f_out:
    reader = csv.reader(f_in, delimiter=";", quotechar='"')
    writer = csv.writer(f_out)                               
    writer.writerows(reader)                                

print(f"Created comma-delimited CSV file: {new}")