"""Compare platelet protein accessions with entries in KinaseDB."""
                            # 1. the Import section: 
from pathlib import Path    #imports in the built in python path tool, helping with locating files and folders (without manual long windows paths)
import re                   #brings in text pattern tool, the script utilises it to check if the IDs in either files match up to uniProt accession format

import pandas as pd         #loads the pandas library, enable the code to work with spreadsheet data like and excel sheet or csv files

                                                        # 2. projects and file paths - defines where the input files like the two datasets and the output file for potential candidate kinases are located
PROJECT_DIR = Path(__file__).resolve().parents[1]       # Project directory is the main folder, path file is the current script, parents[1] means it goes up one level to the main project folder, so it can access the data and output folders
DATA_DIR = PROJECT_DIR / "Data"                         # Data directory is a subfolder of the project directory, its just identifying paths so it knows where to pull things from, the / is used to define a subfolder within the main project folder, in this case the data folder
OUTPUT_DIR = PROJECT_DIR / "Outputs"                    # Output directory is a subfolder of the project directory, identifying things so it knows where to deposit the results 
PLATELET_FILE = DATA_DIR / "All Proteins Confidently Identified In The Platelet Proteome.xls"       # Platelet file path, locating the excel file with the platelet proteome data, the script will read this file to get the protein accessions
KINASE_FILE = DATA_DIR / "KinaseDB.csv"                                                             # Kinase file path, locating the csv file with the kinase database, the script will read this file to get the kinase accessions
OUTPUT_FILE = OUTPUT_DIR / "platelet_kinase_matches.csv"                                            # Output file path, where the results of the comparison will be saved as a csv file
                                                        # 3. 
PLATELET_SHEET = "sorted according to accession"        # Sheet name in the platelet file
PLATELET_HEADER_ROW = 6
PLATELET_ID_COLUMN = "Protein accession"
KINASE_ID_COLUMN = "uniprot_id"

KINASE_ANNOTATION_COLUMNS = [
    "gene names (primary)",
    "HGNC_ID",
    "ensembl_gene_id",
    "ensembl_protein_id",
    "protein names",
    "kinase name",
    "group",
    "species",
    "description",
    "data_sources",
    "substrates",
]


def normalise_accessions(values: pd.Series) -> pd.Series:
    """Trim spaces and standardise letter case before comparing identifiers."""
    return values.fillna("").astype(str).str.strip().str.upper()


def combine_unique_values(values: pd.Series) -> str:
    """Keep distinct non-empty values when several source rows share an ID."""
    unique_values = []
    for value in values.dropna().astype(str):
        value = value.strip()
        if value and value not in unique_values:
            unique_values.append(value)
    return "; ".join(unique_values)


def require_columns(data: pd.DataFrame, required: list[str], source_name: str) -> None:
    missing_columns = [column for column in required if column not in data.columns]
    if missing_columns:
        raise ValueError(
            f"{source_name} is missing required column(s): {', '.join(missing_columns)}. "
            f"Columns found: {', '.join(map(str, data.columns))}"
        )


def main() -> None:
    for file_path in (PLATELET_FILE, KINASE_FILE):
        if not file_path.is_file():
            raise FileNotFoundError(f"Could not find the input file: {file_path}")

    # The workbook has six introductory rows before its table headings.
    platelet_data = pd.read_excel(
        PLATELET_FILE,
        sheet_name=PLATELET_SHEET,
        header=PLATELET_HEADER_ROW,
        dtype=str,
    )
    kinase_data = pd.read_csv(KINASE_FILE, dtype=str, low_memory=False)
    require_columns(platelet_data, [PLATELET_ID_COLUMN], "Platelet workbook")
    require_columns(kinase_data, [KINASE_ID_COLUMN], "KinaseDB CSV")

    platelet_data["_accession"] = normalise_accessions(
        platelet_data[PLATELET_ID_COLUMN]
    )
    kinase_data["_accession"] = normalise_accessions(kinase_data[KINASE_ID_COLUMN])

    platelet_with_ids = platelet_data[platelet_data["_accession"] != ""].copy()
    kinase_with_ids = kinase_data[kinase_data["_accession"] != ""].copy()

    platelet_id_counts = platelet_with_ids["_accession"].value_counts()
    kinase_id_counts = kinase_with_ids["_accession"].value_counts()
    platelet_duplicate_rows = int((platelet_id_counts - 1).clip(lower=0).sum())
    platelet_duplicate_ids = int((platelet_id_counts > 1).sum())
    kinase_duplicate_rows = int((kinase_id_counts - 1).clip(lower=0).sum())

    accession_pattern = re.compile(r"[A-Z0-9]{6,10}(?:-[0-9]+)?")
    unusual_kinase_ids = sorted(
        accession
        for accession in kinase_with_ids["_accession"].unique()
        if not accession_pattern.fullmatch(accession)
    )

    platelet_columns = {
        "Other proteins with same peptides": (
            "other proteins with same peptides (protein group)"
        ),
        "Platelet confidence": "Confidence",
        "Estimated copy numbers per platelet": "Estimated copy numbers per platelet",
        "Platelet description": "Description",
    }
    platelet_aggregation = {
        "Platelet rows for accession": ("_accession", "size"),
    }
    for output_column, source_column in platelet_columns.items():
        if source_column in platelet_with_ids.columns:
            platelet_aggregation[output_column] = (
                source_column,
                combine_unique_values,
            )

    # Grouping first ensures repeated platelet rows do not inflate protein counts.
    platelet_summary = platelet_with_ids.groupby(
        "_accession", as_index=False, sort=False
    ).agg(**platelet_aggregation)

    available_annotations = [
        column
        for column in KINASE_ANNOTATION_COLUMNS
        if column in kinase_with_ids.columns
    ]
    kinase_aggregation = {
        column: (column, combine_unique_values) for column in available_annotations
    }
    kinase_summary = kinase_with_ids.groupby(
        "_accession", as_index=False, sort=False
    ).agg(**kinase_aggregation)

    matches = platelet_summary.merge(
        kinase_summary,
        on="_accession",
        how="inner",
        validate="one_to_one",
    )
    matches = matches.rename(columns={"_accession": "UniProt accession"})
    matches = matches.sort_values("UniProt accession").reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    matches.to_csv(OUTPUT_FILE, index=False)

    unique_platelet_ids = int(platelet_with_ids["_accession"].nunique())
    unique_kinase_ids = int(kinase_with_ids["_accession"].nunique())
    print("Platelet proteome vs KinaseDB")
    print(f"Platelet table rows after its headings: {len(platelet_data):,}")
    print(f"Platelet rows without an accession (not compared): {len(platelet_data) - len(platelet_with_ids):,}")
    print(f"Platelet rows with an accession: {len(platelet_with_ids):,}")
    print(f"Unique platelet accessions checked: {unique_platelet_ids:,}")
    print(
        f"Repeated platelet rows: {platelet_duplicate_rows:,} extra rows "
        f"across {platelet_duplicate_ids:,} accessions"
    )
    print(f"KinaseDB rows: {len(kinase_data):,}")
    print(f"KinaseDB rows without an accession: {len(kinase_data) - len(kinase_with_ids):,}")
    print(f"Unique nonblank KinaseDB IDs: {unique_kinase_ids:,}")
    print(f"Repeated KinaseDB rows: {kinase_duplicate_rows:,} extra rows")
    if unusual_kinase_ids:
        print(
            "KinaseDB ID values outside the observed UniProt-like format "
            f"(still considered in exact matching): {', '.join(unusual_kinase_ids)}"
        )
    print(f"Unique platelet accessions matched: {len(matches):,}")
    print(f"Unique platelet accessions without a match: {unique_platelet_ids - len(matches):,}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(
        "A match means the accession occurs in both files; it does not show "
        "that the kinase is active or that an interaction is occurring."
    )


if __name__ == "__main__":
    main()