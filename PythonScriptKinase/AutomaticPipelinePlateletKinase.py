# Path builds paths relative to this project, rather than relying on the
# terminal's current folder.
from pathlib import Path
# runpy loads and runs Python files whose names cannot be used in a normal
# import statement, such as filenames that start with "1." or "2.".
import runpy

# Find the project folder from this coordinator script's location.
PROJECT_DIR = Path(__file__).resolve().parents[1]
# Keep the two processing stages as separate scripts; this file only coordinates them.
PYTHON_SCRIPT_DIR = PROJECT_DIR / "PythonScriptKinase"
IDENTIFICATION_SCRIPT = (
    PYTHON_SCRIPT_DIR / "1. PlateletKinasePipelineIdentifyAndAnnotate.py"
)
STRING_MAPPING_SCRIPT = (
    PYTHON_SCRIPT_DIR / "2. PlateletKinasePipelineStringMapping.py"
)

# The identification stage writes this exact heading for its UniProt accession.
STRING_INPUT_COLUMN = "UniProt accession"


def run_pipeline() -> None:
    """Run the kinase identification stage, then map its accessions with STRING."""
    # Check both stage files before starting. This avoids completing step one and
    # only then discovering that the STRING mapping script is missing.
    for script_path in (IDENTIFICATION_SCRIPT, STRING_MAPPING_SCRIPT):
        if not script_path.is_file():
            raise FileNotFoundError(f"Could not find pipeline stage: {script_path}")

    # Load the identification script under a name other than "__main__".
    # This defines its functions and settings without automatically running it.
    print("Stage 1 of 2: identify platelet proteins listed in KinaseDB.")
    identification = runpy.run_path(
        str(IDENTIFICATION_SCRIPT),
        run_name="platelet_identification_stage",
    )

    # Get and call the existing stage's main() function. This reuses the working
    # identification code instead of copying it into this coordinator.
    identify_main = identification.get("main")
    if not callable(identify_main):
        raise RuntimeError(
            f"The identification script has no callable main(): "
            f"{IDENTIFICATION_SCRIPT}"
        )
    identify_main()

    # Reuse the identification script's own output path so the two scripts do
    # not need separate, potentially inconsistent copies of that setting.
    kinase_matches_file = Path(identification["OUTPUT_FILE"])
    if not kinase_matches_file.is_file():
        raise FileNotFoundError(
            "The identification stage finished without creating its expected "
            f"output file: {kinase_matches_file}"
        )

    # Load the mapping script in the same way. Its main() function is not called,
    # so the mapper's command-line interface does not run during this import.
    print("Stage 2 of 2: map the kinase accessions to STRING IDs.")
    string_mapping = runpy.run_path(
        str(STRING_MAPPING_SCRIPT),
        run_name="string_mapping_stage",
    )

    # Retrieve the reusable mapping function and its settings class from the
    # loaded script. This avoids duplicating the STRING API code here.
    map_accessions = string_mapping.get("run_mapping")
    mapping_config_type = string_mapping.get("StringMappingConfig")
    if not callable(map_accessions) or not callable(mapping_config_type):
        raise RuntimeError(
            "The STRING mapping script must define run_mapping() and "
            f"StringMappingConfig: {STRING_MAPPING_SCRIPT}"
        )

    # Use the mapper's configured output path (Outputs/STRING_Mapping.csv) and
    # its default human STRING settings. Running the pipeline again regenerates
    # the two output CSVs; the input datasets are not changed.
    mapping_config = mapping_config_type()
    string_output_file = Path(string_mapping["DEFAULT_OUTPUT_PATH"])
    map_accessions(
        source=kinase_matches_file,
        input_column=STRING_INPUT_COLUMN,
        output_path=string_output_file,
        config=mapping_config,
    )

    # Tell the user where to find the final mapping results.
    print("Pipeline complete.")
    print(f"STRING mapping results: {string_output_file}")


if __name__ == "__main__":
    # Start the two-stage workflow only when this coordinator is run directly.
    run_pipeline()