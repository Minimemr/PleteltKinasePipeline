# Path builds filesystem paths relative to this coordinator script.
from pathlib import Path
# pandas reads the mapper's CSV so this coordinator can extract its STRING IDs.
import pandas as pd

# Import the identification stage's output path and entry point directly.
from step_1_platelet_kinase_pipeline_identify_and_annotate import (
    OUTPUT_FILE as IDENTIFICATION_OUTPUT_FILE,
    main as identify_main,
)
# Import the mapper's reusable entry point, configuration, and output path.
from step_2_platelet_kinase_pipeline_string_mapping import (
    DEFAULT_OUTPUT_PATH as DEFAULT_STRING_MAPPING_OUTPUT_PATH,
    StringMappingConfig,
    run_mapping,
)
# Import the reusable STRING network retrieval function.
from step_3_platelet_kinase_pipeline_string_relationship_retriever_with_variable_flavor import (
    retrieve_string_network,
)

# Find the root folder by moving two levels up from this script's resolved path.
PROJECT_DIR = Path(__file__).resolve().parents[1]
# Store the stage scripts in their existing project subfolder.
PYTHON_SCRIPT_DIR = PROJECT_DIR / "PythonScriptKinase"
# Identify the current platelet-protein identification stage.
IDENTIFICATION_SCRIPT = (
    PYTHON_SCRIPT_DIR / "step_1_platelet_kinase_pipeline_identify_and_annotate.py"
)
# Identify the current UniProt-to-STRING mapping stage.
STRING_MAPPING_SCRIPT = (
    PYTHON_SCRIPT_DIR / "step_2_platelet_kinase_pipeline_string_mapping.py"
)
# Identify the reusable interaction retriever module.
STRING_RETRIEVAL_SCRIPT = (
    PYTHON_SCRIPT_DIR
    / (
        "step_3_platelet_kinase_pipeline_string_relationship_retriever_"
        "with_variable_flavor.py"
    )
)

# Match the accession header written by the identification stage.
STRING_INPUT_COLUMN = "UniProt accession"
# Match the identifier header written by the mapping stage.
MAPPING_STRING_ID_COLUMN = "STRING_ID"
# Match the status header used to distinguish successful mappings.
MAPPING_STATUS_COLUMN = "Mapping_Status"
# Use the same medium-confidence score threshold for both network flavors.
REQUIRED_SCORE = 400


def run_pipeline() -> None:
    """Run identification, mapping, then retrieve confidence and typed networks."""
    # Check every required stage before doing any work, avoiding a partial run
    # when one of the later scripts is absent.
    for script_path in (
        IDENTIFICATION_SCRIPT,
        STRING_MAPPING_SCRIPT,
        STRING_RETRIEVAL_SCRIPT,
    ):
        # Stop immediately and name the missing file if a pipeline stage is absent.
        if not script_path.is_file():
            raise FileNotFoundError(f"Could not find pipeline stage: {script_path}")

    # Announce and run the first imported stage; its main guard prevents import-time execution.
    print("Stage 1 of 3: identify platelet proteins listed in KinaseDB.")
    identify_main()

    # Reuse the output path imported from the identification module.
    kinase_matches_file = Path(IDENTIFICATION_OUTPUT_FILE)
    # Verify that identification produced the file required by the mapping stage.
    if not kinase_matches_file.is_file():
        raise FileNotFoundError(
            "The identification stage finished without creating its expected "
            f"output file: {kinase_matches_file}"
        )

    # Announce the second imported stage before mapping the identified proteins.
    print("Stage 2 of 3: map the kinase accessions to STRING IDs.")
    # Instantiate the mapper's default human STRING settings.
    mapping_config = StringMappingConfig()
    # Use the mapper's declared output path instead of duplicating its location.
    string_output_file = Path(DEFAULT_STRING_MAPPING_OUTPUT_PATH)
    # Map the UniProt accessions from the just-created identification output.
    run_mapping(
        source=kinase_matches_file,  # Read the file produced by identification.
        input_column=STRING_INPUT_COLUMN,  # Select the actual UniProt column header.
        output_path=string_output_file,  # Save the mapping CSV at its configured path.
        config=mapping_config,  # Apply the mapper's existing default configuration.
    )

    # Read mapping output as text and retain empty CSV cells as empty strings.
    mapping_data = pd.read_csv(string_output_file, dtype=str, keep_default_na=False)
    # Define the two headers needed to select successfully mapped identifiers.
    required_mapping_columns = {MAPPING_STRING_ID_COLUMN, MAPPING_STATUS_COLUMN}
    # Find any expected headers missing from the mapping file.
    missing_mapping_columns = required_mapping_columns.difference(mapping_data.columns)
    # Refuse to guess column names if the mapper's output format has changed.
    if missing_mapping_columns:
        raise ValueError(
            "The STRING mapping output is missing required column(s): "
            + ", ".join(sorted(missing_mapping_columns))
        )
    # Select only IDs whose mapping status says they were successfully mapped.
    string_ids = mapping_data.loc[
        mapping_data[MAPPING_STATUS_COLUMN] == "Mapped",  # Exclude failed mappings.
        MAPPING_STRING_ID_COLUMN,  # Pass only STRING IDs to the network API.
    ].tolist()  # Convert the selected pandas column to the list expected by the tool.

    # Announce the third stage, which retrieves and saves both network flavors.
    print("Stage 3 of 3: retrieve STRING interaction networks.")
    # Store both network files in the project's shared Outputs directory.
    outputs_directory = PROJECT_DIR / "Outputs"
    # Create that directory if necessary; leave it unchanged if it already exists.
    outputs_directory.mkdir(parents=True, exist_ok=True)
    # Pair each requested API flavor with its distinct output filename.
    for flavour, filename in (
        ("confidence", "Kinase_Network_Confidence.csv"),
        ("typed", "Kinase_Network_Typed.csv"),
    ):
        # Report which flavor is being requested before making its API call.
        print(f"Retrieving STRING {flavour} network.")
        # Request the same input IDs and score threshold for a fair comparison.
        network = retrieve_string_network(
            string_ids,  # Reuse the exact mapped-ID list for both flavors.
            flavour=flavour,  # Let this call choose the network flavor.
            required_score=REQUIRED_SCORE,  # Keep the confidence cutoff fixed.
        )
        # Combine the output directory and this flavor's distinct CSV filename.
        output_path = outputs_directory / filename
        # Save every returned interaction column, without an extra pandas index.
        network.to_csv(output_path, index=False)
        # Tell the user the path where this network was saved.
        print(f"STRING {flavour} network: {output_path}")

    # Report successful completion and the location of the mapping table.
    print("Pipeline complete.")
    print(f"STRING mapping results: {string_output_file}")


# Run the workflow only when this file is executed as a script, not imported.
if __name__ == "__main__":
    # Start the complete identify-map-retrieve pipeline.
    run_pipeline()
