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
                                                        # 3. set of instructions telling the script where to look in berkhart file, creating the settings for later comparisons
PLATELET_SHEET = "sorted according to accession"        # Sheet name in the platelet file, in the excel file look at the sheet sorted according to accesion numbers
PLATELET_HEADER_ROW = 6                                 # at the top of the file there is introductory text and formating, therefore pandas must be told the correct row at the begining of the data table
PLATELET_ID_COLUMN = "Protein accession"                # Column name in the platelet file that contains the protein accessions, this is the column that will be used to compare with the kinase database
KINASE_ID_COLUMN = "uniprot_id"                     # Column name in the kinase file that contains the kinase accessions, this is the column that will be used to compare with the platelet proteome

KINASE_ANNOTATION_COLUMNS = [       # 4. list of columns to keep from the kinase database, these are the columns that will be kept in the output file if a match is found
    "gene names (primary)",         # Column containing the primary gene names
    "HGNC_ID",                      # Column containing the HGNC ID
    "ensembl_gene_id",          # Column containing the Ensembl gene ID
    "ensembl_protein_id",       # Column containing the Ensembl protein ID
    "protein names",            # Column containing the protein names
    "kinase name",              # Column containing the kinase name
    "group",                    # Column containing the kinase group
    "species",                  # Column containing the species
    "description",              # Column containing the description
    "data_sources",             # Column containing the data sources
    "substrates",               # Column containing the substrates
]

                                                                    # 5. functions to normalise and combine values, and check for required columns
def normalise_accessions(values: pd.Series) -> pd.Series:                   # def starts the function and gives it a name normalise accessions for it to be recalled, and the values are a pandas series, which is a one-dimensional array-like structure that can hold any data type, in this case the protein accessions
    """Trim spaces and standardise letter case before comparing identifiers."""         # description of function, also values is the functions input e.g. the accession column from the platelet and kinase files, and the output is a pandas series with the normalised accessions
    return values.fillna("").astype(str).str.strip().str.upper()                 # return sends the results back to the code which calls the function, the code after return is for cleaning the data
                                                                                # fillna("") replaces any missing values (NaN) with an empty string, ensuring that the subsequent operations can be performed without errors. astype(str) converts all values to strings, allowing for consistent string operations. str.strip() removes any leading or trailing whitespace from the strings, which is important for accurate comparisons. str.upper() converts all characters in the strings to uppercase, standardizing the letter case for comparison purposes.
                                                                                # so this cleans up by removing missing values, converting to string, removing whitespace, and converting to uppercase for consistent comparison of protein accessions

                                                            # 6. function to combine unique values from multiple rows sharing the same ID
def combine_unique_values(values: pd.Series) -> str:            # def to define a function called combine_unique_values, which takes a pandas Series as input and returns a string. 
    """Keep distinct non-empty values when several source rows share an ID."""      # description of function
    unique_values = []                                          # Initialize an empty list to store unique values
    for value in values.dropna().astype(str):            # Iterate over the values in the Series, dropping any NaN values and converting them to strings, during repeat the current item is labelled value
        value = value.strip()                       # Remove leading and trailing whitespace from the value
        if value and value not in unique_values:    # Check if the value is non-empty and not already in the list of unique values ???
            unique_values.append(value)             # If the value is valid, append it to the list of unique values
    return "; ".join(unique_values)                 # Return a single string containing the unique values, separated by "; ". This creates a consolidated representation of the distinct values associated with a shared ID.

                                                                                                # 7. function to check for required columns in the input data, if there are missing columns it raises an error and stops the script
def require_columns(data: pd.DataFrame, required: list[str], source_name: str) -> None:         # def to define a function called require_columns, which takes three parameters: data (a pandas DataFrame), required (a list of strings representing required column names), and source_name (a string representing the name of the data source). The function does not return any value (None).
    missing_columns = [column for column in required if column not in data.columns]             # Create a list of missing columns by checking if each required column is present in the DataFrame's columns. If a required column is not found, it is added to the missing_columns list.
    if missing_columns:                                                                         # If there are any missing columns
        raise ValueError(                                                                       # Raise a value error, stops the program 
            f"{source_name} is missing required column(s): {', '.join(missing_columns)}. "      # Provide a descriptive error message, tell which input file is missing which columns
            f"Columns found: {', '.join(map(str, data.columns))}"                               # Lists the columns present in the data, so the user can see what is available and compare it to what is required
        )

                                                                                # 8. main function to orchestrate the comparison of platelet protein accessions with KinaseDB entries
def main() -> None:                                                             # def to define the main function, which serves as the entry point for the script. It does not return any value (None).
    for file_path in (PLATELET_FILE, KINASE_FILE):                              # for loop to chceck each input path for the two files, PLATELET_FILE and KINASE_FILE, to ensure they exist before proceeding with the comparison. The loop iterates over each file path in the tuple (PLATELET_FILE, KINASE_FILE). ???
        if not file_path.is_file():                                             # means if it is not a file 
            raise FileNotFoundError(f"Could not find the input file: {file_path}")  # if the file does not exist, raise a FileNotFoundError with a descriptive message indicating which input file could not be found. This prevents the script from proceeding with missing input files and provides clear feedback to the user.

                                        # 8. part 1: read the platelet and kinase data from their respective files, normalise the accessions, and perform the comparison
    platelet_data = pd.read_excel(      # asks pandas to read the workbook
        PLATELET_FILE,                  # says which file to read through the path defined earlier
        sheet_name=PLATELET_SHEET,      # says which sheet to read through the sheet name defined earlier
        header=PLATELET_HEADER_ROW,     # says which row to use as the header through the header row defined earlier
        dtype=str,                      # ask pandas to read the cell values as strings, important as some variables need to stay text and some values are not purely numbers
    )                                                                           # 8. part 2: read the kinase data from the KinaseDB CSV file, and perform the comparison
    kinase_data = pd.read_csv(KINASE_FILE, dtype=str, low_memory=False)         # identifies the file with kinase_data, asks pandas to read it, then details its location, to read as a text and finally asks pandas not to use its low memory type guessing behavior????
    require_columns(platelet_data, [PLATELET_ID_COLUMN], "Platelet workbook")   # details the required column is the platelet ID column, and if it is missing it will raise an error with the message "Platelet workbook"
    require_columns(kinase_data, [KINASE_ID_COLUMN], "KinaseDB CSV")            # details the required column is the kinase ID column, and if it is missing it will raise an error with the message "KinaseDB CSV"
                                                                    # 8. part 3: normalise the accessions in both datasets, filter out rows without accessions, and prepare for comparison
    platelet_data["_accession"] = normalise_accessions(             # creates a new column in the platelet_data DataFrame called "_accession" and populates it with the normalised accessions obtained by calling the normalise_accessions function on the values in the PLATELET_ID_COLUMN of the
        platelet_data[PLATELET_ID_COLUMN]                           # detailing its for the platelet_data DataFrame, and the column to be normalised is the PLATELET_ID_COLUMN
    )
    kinase_data["_accession"] = normalise_accessions(kinase_data[KINASE_ID_COLUMN])     # creates a new column in the kinase_data DataFrame called "_accession" and populates it with the normalised accessions obtained by calling the normalise_accessions function on the values in the KINASE_ID_COLUMN of the kinase_data DataFrame. This ensures that the kinase accessions are cleaned and standardised for accurate comparison with the platelet proteome.
                                                                                        # 8. part 4: filter out rows without accessions, count duplicates, and prepare for aggregation
    platelet_with_ids = platelet_data[platelet_data["_accession"] != ""].copy()         # creates a new DataFrame called platelet_with_ids that contains only the rows from platelet_data where the "_accession" column is not empty. This filters out any rows without valid accessions, ensuring that only relevant data is considered for comparison with the kinase database. The .copy() method is used to create a copy of the filtered DataFrame, preventing potential issues with chained indexing and allowing for safe modifications to the new DataFrame without affecting the original platelet_data.
    kinase_with_ids = kinase_data[kinase_data["_accession"] != ""].copy()               # creates another new dataframe called kinase_with_ids that contains only the rows from kinase_data where the "_accession" column is not empty. The != means not equal to so filters out rows with nothing

    platelet_id_counts = platelet_with_ids["_accession"].value_counts()                 # counts how many times each unique accession appears in the platelet data, helps identify duplicates and understand the distribution of accessions in the dataset
    kinase_id_counts = kinase_with_ids["_accession"].value_counts()                     # same as the line above but for the kinase data
    platelet_duplicate_rows = int((platelet_id_counts - 1).clip(lower=0).sum())         # calculates how many rows are duplicates, by subtracting 1 from the count of each unique accession (to account for the first occurrence), clipping any negative values to 0 (to avoid counting non-duplicates), and summing the results. This gives the total number of extra rows that are duplicates across all accessions in the platelet data.
    platelet_duplicate_ids = int((platelet_id_counts > 1).sum())                        # .clip(lower=0) ensures that any negative values resulting from the subtraction are set to 0, preventing them from affecting the sum. The final result is converted to an integer for clarity and consistency in reporting. .sum() adds up all these occurences together
    kinase_duplicate_rows = int((kinase_id_counts - 1).clip(lower=0).sum())             # duplicate IDs just tells you how many accessions have duplicated not how many extra rows there are
                                                                                        # line above is just repeating it but for kinase data
    accession_pattern = re.compile(r"[A-Z0-9]{6,10}(?:-[0-9]+)?")           #prepearing a text pattern ro check if IDs are in the correct uniprot format, re.compile() compiles a regular expression pattern, used for matching and searching text. The pattern r"[A-Z0-9]{6,10}(?:-[0-9]+)?" is a regular expression that defines the expected format of UniProt accession numbers. It matches strings that consist of 6 to 10 uppercase letters or digits, optionally followed by a hyphen and one or more digits. This pattern is used to validate the format of kinase accessions in the dataset.
                                                                            # r in the brackets means raw string so certain characters are not treated as escape sequences. [A-Z0-9] matches any uppercase letter or digit, {6,10} specifies that the preceding character class should occur between 6 and 10 times, (?:-[0-9]+)? is a non-capturing group that matches an optional hyphen followed by one or more digits, and the ? at the end indicates that this group is optional.
    unusual_kinase_ids = sorted(                                            # creates a sorted list of kinase accessions that do not match the expected UniProt format, helping to identify any anomalies or unexpected values in the kinase dataset. The sorted() function is used to arrange the accessions in ascending order for easier review and analysis.
        accession
        for accession in kinase_with_ids["_accession"].unique()             # .unique() returns distinct values avoiding checking the same value, for loop for iteration of accessions
        if not accession_pattern.fullmatch(accession)                       # if not for if it doesnt comply with the pattern its added to unusual kinase IDs
    )
                                                                            # 8. part 5: aggregate the platelet and kinase data, merge the summaries, and save the results to a CSV file
    platelet_columns = {                                                    # creates a dictionary with {} then defines it as platelet_columns, a dictionary is a collection of key-value pairs, where each key is a string representing the output column name, and the corresponding value is a string representing the source column name in the platelet data. This dictionary is used to specify which columns from the platelet data should be included in the aggregated summary and how they should be named in the final output.
        "Other proteins with same peptides": (                              # first key value pair, key is other proteins with same peptides
            "other proteins with same peptides (protein group)"             # value is the column name in the platelet data that contains information about other proteins that share the same peptides
        ),
        "Platelet confidence": "Confidence",                                                # second key value pair, key is platelet confidence, value is name in platelet data about the confidence of protein identification
        "Estimated copy numbers per platelet": "Estimated copy numbers per platelet",       # third key value pair, key is estimated copy numbers per platelet, value is name in platelet data about the estimated copy numbers
        "Platelet description": "Description",                                              # fourth key value pair, key is platelet description, value is name in platelet data about the description of the protein
    }
    platelet_aggregation = {                                            # creates a dictionary called platelet_aggregation
        "Platelet rows for accession": ("_accession", "size"),          # first key-value pair, key is "Platelet rows for accession", value is a tuple containing the column name "_accession" and the aggregation function "size". This means that when aggregating the platelet data, it will count the number of rows for each unique accession.
    }
    for output_column, source_column in platelet_columns.items():       # for loop to iterate over the items in the platelet_columns dictionary, where output_column is the key (the desired output column name) and source_column is the value (the corresponding source column name in the platelet data), ,items returns a view object that displays a list of a dictionary's key-value tuple pairs, allowing the loop to access both the output column name and the source column name for each entry in the platelet_columns dictionary.
        if source_column in platelet_with_ids.columns:                  # if step which checks if the source_column exists in the columns of the platelet_with_ids DataFrame. This ensures that only valid columns present in the platelet data are included ???
            platelet_aggregation[output_column] = (                     # creates a new entry in the platelet_aggregation dictionary, where the key is the output_column (the desired output column name) and the value is a tuple containing the
                source_column,                                          # the source column name in the platelet data, and the aggregation function
                combine_unique_values,                                  # function from earlier which combs columns and if repeated row is viewed it is not included in the list its building, also cleans up text
            )

                                                                # 8. part 6: groups platelet rows by accession once cleaned and calculates summary for each group
    platelet_summary = platelet_with_ids.groupby(               # creates the platelet_summary, from the platelet with ids column, the .groupby puts rows together with the same accession
        "_accession", as_index=False, sort=False                # as_index=False keep the accessions as a regular column in the results, sort=False avoids asking pandas to sort the groups at this step
    ).agg(**platelet_aggregation)                               # .agg applies aggregation instructions from the dicitionary, **platelet_aggregation passes the dictionarys entries into .agg as instructions (** is python for unpacking a dictionary into keywords)

    available_annotations = [                                   # 8. Part 7: preparing KinaseDB annotations
        column                                                  # list of the columns listed early on, title of each column from the CSV file
        for column in KINASE_ANNOTATION_COLUMNS                 # builds a list of the requested annotation columns, check items one by one to build the list with if and for
        if column in kinase_with_ids.columns                    
    ]
    kinase_aggregation = {                                                                  # creates a aggregation rule for every annotation column
        column: (column, combine_unique_values) for column in available_annotations         # utilises the combine function from earlier in the selected columns
    }
    kinase_summary = kinase_with_ids.groupby(               # groups all kinaseDB records now with cleaned accessions then combines the selected annoatation into one row per accession
        "_accession", as_index=False, sort=False            # 
    ).agg(**kinase_aggregation)
                                                    # 8. Part 8: matching up the datasets
    matches = platelet_summary.merge(               # combines two summary tables, matching accession values, .merge is pandas table combining operation
        kinase_summary,                             
        on="_accession",                            # tells to compare that column in both tables
        how="inner",                                # keeps only rows with and accession present in both tables
        validate="one_to_one",                      # asks pandas to cehck that each accession appears at most once in each summary, error raised if not met
    )
    matches = matches.rename(columns={"_accession": "UniProt accession"})           # renames the internal helper column to a better name for the file
    matches = matches.sort_values("UniProt accession").reset_index(drop=True)       # orders the matched results by accession with .sort_values, .reset_index(drop=True), drop=True means dont keep old index as a column
                                                                # 9. Saving the results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)               # Saving the results, making sure folder exists, .mkdir is used to create folders in the correct location given by the path, Parents=True allows for python to make missing folder if needed, exist_ok=True error avoided if there is a folder
    matches.to_csv(OUTPUT_FILE, index=False)                    # output file details which folder to save it in, index=False stops pandas row numbers from being written into the file as a column
                                                                                # 10. Printing
    unique_platelet_ids = int(platelet_with_ids["_accession"].nunique())        # .nunique() counts distinct values
    unique_kinase_ids = int(kinase_with_ids["_accession"].nunique())            # these two lines print how many unique non blank accessions were available
    print("Platelet proteome vs KinaseDB")                                      # print heading
    print(f"Platelet table rows after its headings: {len(platelet_data):,}")                #these next couple lines report, table rows, blank IDs, duplicate rows and matches found
    print(f"Platelet rows without an accession (not compared): {len(platelet_data) - len(platelet_with_ids):,}")            # Subtracts rows with IDs from rows read to give number of rows excluded for having blank accessions
    print(f"Platelet rows with an accession: {len(platelet_with_ids):,}")       # f strings used to insert values for each variables
    print(f"Unique platelet accessions checked: {unique_platelet_ids:,}")       # :, adds thousands seperators e.g. 4,000 instead of 4000
    print(
        f"Repeated platelet rows: {platelet_duplicate_rows:,} extra rows "
        f"across {platelet_duplicate_ids:,} accessions"
    )
    print(f"KinaseDB rows: {len(kinase_data):,}")
    print(f"KinaseDB rows without an accession: {len(kinase_data) - len(kinase_with_ids):,}")
    print(f"Unique nonblank KinaseDB IDs: {unique_kinase_ids:,}")
    print(f"Repeated KinaseDB rows: {kinase_duplicate_rows:,} extra rows")
    if unusual_kinase_ids:                                                              # prints a list of unusual formatted IDs only if it is non-empty
        print(
            "KinaseDB ID values outside the observed UniProt-like format "
            f"(still considered in exact matching): {', '.join(unusual_kinase_ids)}"
        )
    print(f"Unique platelet accessions matched: {len(matches):,}")
    print(f"Unique platelet accessions without a match: {unique_platelet_ids - len(matches):,}")        # subtracts the num of matched accessions from total num of unique accessions to find unique platelet accessions without a match
    print(f"Results saved to: {OUTPUT_FILE}")
    print(
        "A match means the accession occurs in both files; it does not show "
        "that the kinase is active or that an interaction is occurring."
    )

                                                    # 10. Starting the script
if __name__ == "__main__":                          # for if functions within the file needed to be imported, it wont run the whole file just the requested function
                                                    # instead of just pulling the function it will run all the code within this file first then the one its imported to
                                                    # if something has been imported the name variable within the imported script is equal to its title not __main__  so it wont execute the programme under the if