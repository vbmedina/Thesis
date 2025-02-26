import requests
import json
import os
import pandas as pd
import time
import re

# Directory to save downloaded data
SAVE_DIR = "chembl_data"
os.makedirs(SAVE_DIR, exist_ok=True)

# Base URL for ChEMBL API
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

# List of known Plasmodium falciparum strains to look for
KNOWN_STRAINS = ["W2", "NF54", "D6", "HB3", "K1", "D10", "3D7", "Dd2", "7G8", "FCR3"]

def fetch_data(endpoint, params=None):
    """Fetch data from ChEMBL API with retries and backoff."""
    url = f"{BASE_URL}/{endpoint}.json"
    retries = 5
    delay = 5

    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.Timeout:
            print(f"Timeout. Retrying in {delay} seconds...")
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}")
            if response.status_code == 500:
                print(f"Server Error (500). Retrying in {delay} seconds...")
            return None  

        time.sleep(delay)
        delay *= 2

    print(f"API request failed after {retries} retries: {url}")
    return None

def fetch_multiple_molecules(chembl_ids, batch_size=50):
    """Fetch multiple molecule properties in batches to prevent API overload."""
    molecules = {}
    
    for i in range(0, len(chembl_ids), batch_size):
        batch_ids = chembl_ids[i:i+batch_size]
        ids_query = ",".join(batch_ids)
        data = fetch_data("molecule", params={"molecule_chembl_id__in": ids_query})

        if data and "molecules" in data:
            for mol in data["molecules"]:
                molecule_structures = mol.get("molecule_structures") or {}
                molecule_properties = mol.get("molecule_properties") or {}

                molecules[mol["molecule_chembl_id"]] = {
                    "Canonical_SMILES": molecule_structures.get("canonical_smiles"),
                    "pChEMBL_Value": molecule_properties.get("pchembl_value"),
                    "Molecule_Max_Phase": mol.get("max_phase"),
                    "RO5_Violations": molecule_properties.get("num_ro5_violations"),
                    "AlogP": molecule_properties.get("alogp"),
                    "Molecular_Weight": molecule_properties.get("full_mwt"),
                }
        time.sleep(1)

    return molecules

def fetch_multiple_assays(assay_ids, batch_size=50):
    """Fetch multiple assay descriptions in batches to prevent API overload."""
    assays = {}
    
    for i in range(0, len(assay_ids), batch_size):
        batch_ids = assay_ids[i:i+batch_size]
        ids_query = ",".join(batch_ids)
        data = fetch_data("assay", params={"assay_chembl_id__in": ids_query})

        if data and "assays" in data:
            for assay in data["assays"]:
                description = assay.get("description", "")
                assays[assay["assay_chembl_id"]] = {
                    "Assay_Description": description,
                    "Strain_Tested": extract_strain(description),
                }
        time.sleep(1)

    return assays

def extract_strain(description):
    """Extract strain name from assay description."""
    if not description:
        return "Unknown"

    for strain in KNOWN_STRAINS:
        if re.search(rf"\b{strain}\b", description, re.IGNORECASE):
            return strain
    return "Unknown"

def save_data_chunk(data, filename, mode="a"):
    """Save a chunk of data to a CSV file to prevent memory overload."""
    df = pd.DataFrame(data)

    # Sort by strain tested
    df = df.sort_values(by="Strain_Tested")

    filepath = os.path.join(SAVE_DIR, f"{filename}.csv")

    # Append mode if the file already exists
    if mode == "a" and os.path.exists(filepath):
        df.to_csv(filepath, mode="a", header=False, index=False)
    else:
        df.to_csv(filepath, index=False)

    print(f"Saved {len(data)} records to {filename}.csv")

def download_ic50_for_plasmodium_falciparum():
    """Download all IC50 bioactivity data for Plasmodium falciparum with optimized pagination."""
    
    params = {
        "target_chembl_id": "CHEMBL364",
        "standard_type": "IC50",
        "limit": 500,
        "offset": 0
    }

    all_data = []
    batch_size = 1000  # Save every 1000 records to avoid crashes
    filename = "ic50_plasmodium_falciparum_sorted"

    while True:
        print(f"Fetching records {params['offset']} to {params['offset'] + params['limit']}...")

        data = fetch_data("activity", params=params)
        
        if not data or "activities" not in data:
            print("No data found or API failed.")
            break

        activities = data["activities"]

        # Collect ChEMBL IDs for batch fetching
        chembl_ids = list({act.get("molecule_chembl_id") for act in activities if act.get("molecule_chembl_id")})
        assay_ids = list({act.get("assay_chembl_id") for act in activities if act.get("assay_chembl_id")})

        molecule_data = fetch_multiple_molecules(chembl_ids, batch_size=50)
        assay_data = fetch_multiple_assays(assay_ids, batch_size=50)

        batch_data = []
        for activity in activities:
            chembl_id = activity.get("molecule_chembl_id")
            assay_id = activity.get("assay_chembl_id")

            formatted_entry = {
                "Compound_ID": chembl_id,
                "IC50_Value": activity.get("standard_value"),
                "Units": activity.get("standard_units"),
                "BAO_Format": activity.get("bao_format"),
                "Assay_ID": assay_id,
                "Document_ID": activity.get("document_chembl_id"),
                "Relation": activity.get("standard_relation"),
                "Pubmed_ID": activity.get("document_id"),
                "Target_ID": activity.get("target_chembl_id"),
                "Organism": "Plasmodium falciparum",
                **molecule_data.get(chembl_id, {}),
                **assay_data.get(assay_id, {}),
            }

            batch_data.append(formatted_entry)

        all_data.extend(batch_data)

        # Save every batch_size records to avoid memory overload
        if len(all_data) >= batch_size:
            save_data_chunk(all_data, filename)
            all_data = []

        params["offset"] += params["limit"]
        print(f"{params['offset']} records fetched so far...")

        if len(activities) < params["limit"]:
            print("All available records fetched.")
            break

        time.sleep(2)

    # Save any remaining data
    if all_data:
        save_data_chunk(all_data, filename)

if __name__ == "__main__":
    download_ic50_for_plasmodium_falciparum()