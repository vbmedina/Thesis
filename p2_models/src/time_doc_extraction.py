# Get publication years from documents
import time, pathlib, warnings, traceback
import pandas as pd, tqdm
from chembl_webresource_client.new_client import new_client
from requests.exceptions import RetryError, ConnectionError
from urllib3.exceptions  import HTTPError

# Paths
PROJECT = pathlib.Path(__file__).resolve().parents[2]
TIDY = PROJECT / "p2_models/data/tidy/pp_tidy.csv"
OUT_CSV = PROJECT / "p2_models/data/doc_year_lookup.csv"

# Document ID column in CSV
DOC_COL = "Document ChEMBL ID" 

# Get unique document IDs 
doc_ids = pd.read_csv(TIDY, usecols=[DOC_COL])[DOC_COL].unique()
docs_api = new_client.document
records = []

# Fetch loop
for j, doc in enumerate(tqdm.tqdm(doc_ids, desc="Getting publication years")):
    year = None
    try:
        rec  = docs_api.get(doc) # Single REST call
        year = rec.get("year", None) # None if missing
    except (HTTPError, RetryError, ConnectionError):
        # Too many 404s or transient network issue
        if j % 50 == 0:
            warnings.warn(f"HTTP error on {doc}; continuing")
    except Exception as e:
        warnings.warn(f"Unexpected error on {doc}: {e}")
        traceback.print_exc(limit=1)

    records.append({DOC_COL: doc, "year": year})

    # Throttle requests to avoid rate limits
    if (j + 1) % 100 == 0: 
        time.sleep(0.5)

# Save results
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(records).to_csv(OUT_CSV, index=False)

# Count filled years
filled = sum(r["year"] is not None for r in records)
print(f"\n Saved  {OUT_CSV}\n - {filled:,} / {len(records):,} documents have a year.")