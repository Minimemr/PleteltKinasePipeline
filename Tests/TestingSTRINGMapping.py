
import requests
import pandas as pd
from io import StringIO
from pathlib import Path

# 1. STRING API endpoint
url = "https://version-12-5.string-db.org/api/tsv/get_string_ids"

# 2. Test UniProt accessions
accessions = [
    "P31749",  # AKT1
    "P28482",  # MAPK1
    "P04637",  # TP53
    "P00533",  # EGFR
    "P42336",  # PIK3CA
    "NOT_A_REAL_UNIPROT_ID"
]

# 3. Prepare API parameters
parameters = {
    "identifiers": "\r".join(accessions),
    "species": 9606,
    "echo_query": 1,
    "caller_identity": "platelet_kinase_mapping_test"
}

# 4. Submit request
response = requests.post(
    url,
    data=parameters,
    timeout=30
)

# 5. Check response
response.raise_for_status()

# 6. Convert returned TSV into a DataFrame
results = pd.read_csv(
    StringIO(response.text),
    sep="\t"
)

# 7. Display results
print("\nSTRING MAPPING RESULTS\n")
print(results.to_string(index=False))

# 8. Save results
# Find the folder this Python file is in (the Tests folder).
tests_folder = Path(__file__).resolve().parent
# Build the full path to the CSV file inside that folder.
output_path = tests_folder / "test_STRING_results.csv"
results.to_csv(
    output_path,
    index=False
)

print(f"\nResults saved successfully to: {output_path}")
