"""Map protein accessions to STRING identifiers and save a CSV report."""

# Tell Python to evaluate type hints later; this helps when a type is named
# before Python has finished loading the whole module.
from __future__ import annotations

# argparse is Python's built-in helper for reading options typed in a terminal.
import argparse
# dataclass makes a small class for holding related settings.
from dataclasses import dataclass
# datetime and timezone are used to record when the mappings were retrieved.
from datetime import datetime, timezone
# StringIO lets pandas read response text as though it came from a file.
from io import StringIO
# Path provides a safe way to build and check file paths on different systems.
from pathlib import Path
# re provides regular-expression tools for checking accession shapes.
import re
# time lets the tool pause between API requests and retries.
import time
# These typing names make the kinds of accepted inputs clearer to readers.
from typing import Any, Optional, Sequence, Union

# pandas reads CSV files and stores data in table-like DataFrames.
import pandas as pd
# requests sends HTTP requests to the STRING web API.
import requests


# __file__ is this script's path; parents[1] moves up to the project folder.
PROJECT_DIR = Path(__file__).resolve().parents[1]
# The default results file is kept in the project's Outputs folder.
DEFAULT_OUTPUT_PATH = PROJECT_DIR / "Outputs" / "STRING_Mapping.csv"

# These are the column names that STRING documents for the TSV mapping response.
# The program checks for all of them before trusting the response.
REQUIRED_RESPONSE_COLUMNS = {
    "queryItem",
    "queryIndex",
    "stringId",
    "ncbiTaxonId",
    "taxonName",
    "preferredName",
    "annotation",
}
# These are the columns written to the final CSV, in the order they will appear.
OUTPUT_COLUMNS = [
    "UniProt_Accession",
    "STRING_ID",
    "Preferred_Name",
    "Organism_ID",
    "Organism",
    "Annotation",
    "Mapping_Status",
    "STRING_Version",
    "Retrieved_At_UTC",
]

# This checks the general shape of a UniProt accession, including isoform suffixes.
# It checks the format, not whether UniProt knows the accession.
ACCESSION_PATTERN = re.compile(r"[A-Z0-9]{6,10}(?:-[1-9][0-9]*)?")


@dataclass(frozen=True)
class StringMappingConfig:
    """Settings that can be changed without editing the mapping workflow."""

    # Use a version-specific STRING URL so later database updates do not silently
    # change results from this run.
    api_url: str = "https://version-12-5.string-db.org/api/tsv/get_string_ids"
    # Record the same version in the result CSV for reproducibility.
    string_version: str = "12.5"
    # 9606 is the NCBI taxonomic identifier for Homo sapiens.
    organism_id: int = 9606
    # STRING asks API users to identify the tool making the request.
    caller_identity: str = "reading_platelet_kinase_pipeline"
    # Send at most this many distinct accessions in one API request.
    # This is configurable; it is not a maximum specified by STRING's docs.
    batch_size: int = 100
    # Stop waiting for a single web request after this many seconds.
    timeout_seconds: float = 30.0
    # Number of additional attempts after the first request fails temporarily.
    max_retries: int = 3
    # Wait at least this long between separate requests to STRING.
    request_interval_seconds: float = 1.0
    # Increase the wait between retries to avoid repeatedly hitting a busy server.
    retry_backoff_seconds: float = 2.0

    def __post_init__(self) -> None:
        # Check settings immediately when this configuration object is created.
        # Invalid settings should produce a clear error instead of failing later.
        if self.batch_size < 1:
            raise ValueError("batch_size must be at least 1.")
        if self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be greater than 0.")
        if self.max_retries < 0:
            raise ValueError("max_retries cannot be negative.")
        if self.request_interval_seconds < 0:
            raise ValueError("request_interval_seconds cannot be negative.")
        if self.retry_backoff_seconds < 0:
            raise ValueError("retry_backoff_seconds cannot be negative.")
        if not self.caller_identity.strip():
            raise ValueError("caller_identity must not be blank.")


class StringApiError(RuntimeError):
    """An API request failed or STRING returned an unexpected response."""


def load_input(
    source: Union[str, Path, Sequence[Any]],
    input_column: Optional[str] = None,
) -> pd.DataFrame:
    """Load accessions from a CSV column or a Python list, preserving their order."""
    # A string or Path means the input is a filename; other accepted inputs are
    # sequences such as Python lists or tuples.
    if isinstance(source, (str, Path)):
        # Convert the supplied path to a Path object so file checks are available.
        input_path = Path(source)
        # Stop with the missing path in the message rather than continue silently.
        if not input_path.is_file():
            raise FileNotFoundError(f"Could not find the input CSV: {input_path}")

        # Read values as text to preserve accessions exactly. Keep blank cells as
        # blanks and keep fully blank CSV rows so they can be marked Invalid_Input.
        data = pd.read_csv(
            input_path,
            dtype=str,
            keep_default_na=False,
            skip_blank_lines=False,
        )
        # If the user did not choose a column, infer it only when the CSV has
        # exactly one column. Multiple columns would make guessing unsafe.
        if input_column is None:
            if len(data.columns) != 1:
                raise ValueError(
                    "This CSV has more than one column. Pass input_column with "
                    "the column containing UniProt accessions."
                )
            input_column = str(data.columns[0])
        # A misspelled column name should produce a clear error listing the headers.
        if input_column not in data.columns:
            raise ValueError(
                f"Input column {input_column!r} was not found. "
                f"Columns found: {', '.join(map(str, data.columns))}"
            )
        # Pull out only the selected identifier column, preserving its row order.
        values = data[input_column].tolist()
    else:
        # Bytes are sequences too, but are not a list of protein identifiers.
        # Reject them and any other unsupported input with an understandable error.
        if isinstance(source, (bytes, bytearray)) or not isinstance(source, Sequence):
            raise TypeError("source must be a CSV path or a sequence of accessions.")
        # Make a list copy so later processing has a consistent, ordered collection.
        values = list(source)

    # Store each original value and its position. The original value will be used
    # in the output, while the position helps preserve the original row order.
    return pd.DataFrame(
        {
            "UniProt_Accession": values,
            "_input_order": range(len(values)),
        }
    )


def validate_accessions(accessions: Sequence[Any]) -> pd.DataFrame:
    """Return cleaned values and validity flags without changing the original inputs."""
    # Keep cleaned values and validity results in separate lists, one item per input.
    normalised_values = []
    valid_values = []

    # Check every original accession in turn without changing the input list.
    for accession in accessions:
        # pandas represents missing values in several possible ways; replace them
        # with an empty string for checking.
        if pd.isna(accession):
            normalised = ""
        else:
            # Convert to text, remove surrounding spaces, and use uppercase so
            # minor formatting differences do not prevent an otherwise exact match.
            normalised = str(accession).strip().upper()

        # Keep the cleaned value for matching against STRING's echoed query.
        normalised_values.append(normalised)
        # fullmatch requires the entire value to fit the pattern, not just part of it.
        valid_values.append(bool(ACCESSION_PATTERN.fullmatch(normalised)))

    # Return two aligned columns: cleaned accession and True/False validity flag.
    return pd.DataFrame(
        {
            "_normalised_accession": normalised_values,
            "_is_valid": valid_values,
        }
    )


def remove_duplicate_accessions(accessions: Sequence[str]) -> list[str]:
    """Keep the first occurrence of each accession, preserving input order."""
    # dict keys are unique and keep insertion order in supported Python versions.
    # dict.fromkeys therefore removes duplicates without alphabetically sorting IDs.
    return list(dict.fromkeys(accessions))


def split_into_batches(
    accessions: Sequence[str], batch_size: int
) -> list[list[str]]:
    """Split accessions into request-sized groups."""
    # A zero or negative batch size cannot make progress through the input.
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1.")
    # range advances by batch_size; slicing takes the current group.
    # The final slice may be shorter when the list does not divide evenly.
    return [
        list(accessions[start : start + batch_size])
        for start in range(0, len(accessions), batch_size)
    ]


def _retry_delay(
    attempt: int,
    config: StringMappingConfig,
    response: Optional[requests.Response] = None,
) -> float:
    """Choose a polite wait before retrying a failed request."""
    # If STRING tells us how long to wait, respect that instruction when its
    # Retry-After value is a number of seconds.
    if response is not None:
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                # Do not wait less than the configured minimum request interval.
                return max(config.request_interval_seconds, float(retry_after))
            except ValueError:
                # If the header is not a number, use the normal backoff instead.
                pass
    # Exponential backoff grows the wait after each failed attempt:
    # 2 seconds, then 4, then 8 with the default settings.
    backoff = config.retry_backoff_seconds * (2**attempt)
    # Also keep the wait no shorter than STRING's configured request interval.
    return max(config.request_interval_seconds, backoff)


def query_string(
    accessions: Sequence[str],
    config: StringMappingConfig,
    session: requests.Session,
    wait_before_request: bool = False,
) -> Optional[str]:
    """POST one batch to STRING; return TSV text, or None for STRING's no-match 404."""
    # Pause before every batch after the first. Retry pauses are handled below.
    if wait_before_request:
        time.sleep(config.request_interval_seconds)

    # These values are sent in the body of the POST request, not placed in the URL.
    parameters = {
        # STRING accepts multiple identifiers separated by carriage returns.
        "identifiers": "\r".join(accessions),
        # Restrict lookups to the requested organism.
        "species": config.organism_id,
        # Ask STRING to echo each submitted ID in the response as queryItem.
        "echo_query": 1,
        # Included to follow the requested setting. STRING's docs say this method
        # returns a single best match, but do not document the limit parameter.
        "limit": 1,
        # Identify this script to STRING.
        "caller_identity": config.caller_identity,
    }
    # These temporary server/rate-limit responses are worth retrying.
    retryable_statuses = {429, 500, 502, 503, 504}

    # max_retries counts extra attempts, so +1 includes the first request.
    for attempt in range(config.max_retries + 1):
        try:
            # Send the batch. timeout prevents a request hanging indefinitely.
            response = session.post(
                config.api_url,
                data=parameters,
                timeout=config.timeout_seconds,
            )
        # Connection errors and timeouts may be temporary, so retry them.
        except (requests.Timeout, requests.ConnectionError) as error:
            # After the final allowed attempt, report a clear API error.
            if attempt >= config.max_retries:
                raise StringApiError(
                    f"STRING request failed after {attempt + 1} attempt(s): {error}"
                ) from error
            # Pause before retrying and then continue to the next loop attempt.
            time.sleep(_retry_delay(attempt, config))
            continue
        # Other request setup errors are not assumed to be temporary.
        except requests.RequestException as error:
            raise StringApiError(f"Could not send STRING request: {error}") from error

        # STRING documents 404 when none of the submitted identifiers resolve.
        # Returning None lets run_mapping label that batch Unmapped.
        if response.status_code == 404:
            return None

        # Retry server errors and rate limiting, but do not retry permanent
        # client errors such as a malformed request.
        if response.status_code in retryable_statuses:
            # Stop retrying when the configured attempts have all been used.
            if attempt >= config.max_retries:
                raise StringApiError(
                    f"STRING returned HTTP {response.status_code} after "
                    f"{attempt + 1} attempt(s): {response.text[:500]}"
                )
            # Respect Retry-After if present; otherwise use exponential backoff.
            time.sleep(_retry_delay(attempt, config, response))
            continue

        # For other non-success HTTP statuses, raise_for_status creates an
        # HTTPError, which is converted to the mapper's clearer StringApiError.
        try:
            response.raise_for_status()
        except requests.HTTPError as error:
            raise StringApiError(
                f"STRING returned HTTP {response.status_code}: "
                f"{response.text[:500]}"
            ) from error

        # A successful response body is TSV text for process_api_response().
        return response.text

    # Defensive fallback: normally the loop returns a response or raises above.
    raise StringApiError("STRING request ended without a response.")


def process_api_response(
    response_text: str, requested_accessions: Sequence[str]
) -> dict[str, dict[str, str]]:
    """Parse STRING's TSV and associate each result with its submitted accession."""
    # An empty successful response is unexpected; do not silently call it Unmapped.
    if not response_text.strip():
        raise StringApiError("STRING returned an empty response with HTTP 200.")

    try:
        # StringIO wraps text in a file-like object; sep="\t" means tab-separated.
        response_data = pd.read_csv(StringIO(response_text), sep="\t", dtype=str)
    # Malformed or headerless TSV should be reported as an API response problem.
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise StringApiError(f"Could not parse STRING's TSV response: {error}") from error

    # Only use a response that has all fields this tool needs.
    missing_columns = REQUIRED_RESPONSE_COLUMNS.difference(response_data.columns)
    if missing_columns:
        raise StringApiError(
            "STRING response is missing expected column(s): "
            + ", ".join(sorted(missing_columns))
        )
    # A header without any result rows is not the documented no-match response;
    # STRING documents 404 for a query where no identifiers resolve.
    if response_data.empty:
        raise StringApiError("STRING returned a TSV header but no mapping rows.")

    # Build a dictionary keyed by our cleaned accession for quick matching later.
    mappings = {}
    # Keep track of positions already used; each input should have no more than
    # one best-match row according to STRING's mapping endpoint documentation.
    seen_query_indexes = set()
    # Iterate through each returned protein mapping row.
    for _, row in response_data.iterrows():
        try:
            # queryIndex is STRING's zero-based position of this query in the batch.
            query_index = int(row["queryIndex"])
        # If STRING did not return a usable numeric position, we cannot safely
        # connect that row to the correct input.
        except (TypeError, ValueError) as error:
            raise StringApiError(
                f"STRING returned an invalid queryIndex: {row['queryIndex']!r}"
            ) from error

        # Check that the returned position points to an item in this batch.
        if not 0 <= query_index < len(requested_accessions):
            raise StringApiError(
                f"STRING returned queryIndex {query_index}, outside this batch."
            )
        # Repeated positions would mean ambiguous or unexpected API results.
        if query_index in seen_query_indexes:
            raise StringApiError(
                f"STRING returned more than one result for queryIndex {query_index}."
            )
        # Mark this position as handled before checking the echoed identifier.
        seen_query_indexes.add(query_index)

        # Use the request position to get the corresponding cleaned accession.
        accession = requested_accessions[query_index]
        # queryItem is the identifier STRING echoed back because echo_query=1.
        echoed_accession = str(row["queryItem"]).strip().upper()
        # Verify the echo agrees with the position; this prevents attaching the
        # correct-looking mapping to the wrong original protein by mistake.
        if echoed_accession != accession:
            raise StringApiError(
                f"STRING echoed {row['queryItem']!r} for submitted accession "
                f"{accession!r}."
            )

        # Rename STRING's field names to the clearer column names used in our CSV.
        mappings[accession] = {
            "STRING_ID": row["stringId"],
            "Preferred_Name": row["preferredName"],
            "Organism_ID": row["ncbiTaxonId"],
            "Organism": row["taxonName"],
            "Annotation": row["annotation"],
        }

    # Return only successfully mapped identifiers. Missing valid IDs are labelled
    # Unmapped later, after we know their request batch succeeded.
    return mappings


def match_results(
    input_data: pd.DataFrame,
    validation: pd.DataFrame,
    mappings: dict[str, dict[str, str]],
    api_error_accessions: set[str],
    string_version: str,
    retrieved_at_utc: str,
) -> pd.DataFrame:
    """Restore original rows and order, assigning a status to every input."""
    # Append results in input order; do not reorder or discard duplicate inputs.
    results = []
    for position, row in input_data.iterrows():
        # Retain exactly what the user supplied for the output accession column.
        original_accession = row["UniProt_Accession"]
        # The cleaned value is only for matching, not for replacing the original.
        normalised_accession = validation.at[position, "_normalised_accession"]
        # Read whether this input passed the simple accession-format check.
        is_valid = validation.at[position, "_is_valid"]

        # Start with blank output fields; fill in the available mapping or status.
        result = {column: pd.NA for column in OUTPUT_COLUMNS}
        # Preserve the original value, including original case or surrounding spaces.
        result["UniProt_Accession"] = original_accession
        # Record which STRING database version and retrieval time apply to this row.
        result["STRING_Version"] = string_version
        result["Retrieved_At_UTC"] = retrieved_at_utc

        # Invalid or missing values are never sent to STRING.
        if not is_valid:
            result["Mapping_Status"] = "Invalid_Input"
        # A valid ID in a failed batch is not evidence that it is absent from STRING.
        elif normalised_accession in api_error_accessions:
            result["Mapping_Status"] = "API_Error"
        # A result returned by STRING is a successful mapping.
        elif normalised_accession in mappings:
            result.update(mappings[normalised_accession])
            result["Mapping_Status"] = "Mapped"
        # No returned row after a successful query means no mapping was found.
        else:
            result["Mapping_Status"] = "Unmapped"

        # Save this completed row before processing the next input.
        results.append(result)

    # Use the declared column order, including when there are no input rows.
    return pd.DataFrame(results, columns=OUTPUT_COLUMNS)


def export_results(results: pd.DataFrame, output_path: Union[str, Path]) -> Path:
    """Write the mapping table to a CSV without changing any input file."""
    # Accept either a path string or Path and make sure parent folders exist.
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # index=False avoids adding pandas row numbers; na_rep keeps missing values blank.
    results.to_csv(destination, index=False, na_rep="")
    # Return the final path so the caller can report where it saved the file.
    return destination


def run_mapping(
    source: Union[str, Path, Sequence[Any]],
    input_column: Optional[str] = None,
    output_path: Union[str, Path] = DEFAULT_OUTPUT_PATH,
    config: StringMappingConfig = StringMappingConfig(),
) -> pd.DataFrame:
    """Run input loading, validation, STRING lookup, matching, and CSV export."""
    # Load the CSV column or Python sequence into a consistent table format.
    input_data = load_input(source, input_column)

    # Protect the source CSV from accidental replacement by the output file.
    if isinstance(source, (str, Path)):
        input_path = Path(source).resolve()
        if Path(output_path).resolve() == input_path:
            raise ValueError("The output CSV must not overwrite the input CSV.")

    # Validate cleaned identifiers while keeping the original column unchanged.
    validation = validate_accessions(input_data["UniProt_Accession"].tolist())
    # Exclude invalid values from API requests; these will still get output rows.
    valid_accessions = [
        accession
        for accession, is_valid in zip(
            validation["_normalised_accession"], validation["_is_valid"]
        )
        if is_valid
    ]
    # Count duplicates so the progress message can report how many API calls
    # were avoided by querying each valid identifier only once.
    valid_accession_counts = pd.Series(valid_accessions, dtype="string").value_counts()
    duplicate_accessions = valid_accession_counts[valid_accession_counts > 1]
    duplicate_rows = int((duplicate_accessions - 1).sum())
    # Deduplicate in first-seen order, then divide the list into configured batches.
    unique_accessions = remove_duplicate_accessions(valid_accessions)
    batches = split_into_batches(unique_accessions, config.batch_size)

    # Give a summary before making network requests so the user can see what
    # will be sent and how duplicate rows are being handled.
    print(
        f"Input rows: {len(input_data):,}; valid accession rows: "
        f"{len(valid_accessions):,}; unique accessions to query: "
        f"{len(unique_accessions):,}; repeated rows: {duplicate_rows:,} "
        f"across {len(duplicate_accessions):,} accession(s)."
    )

    # These dictionaries collect successful mappings and failed batches separately.
    mappings = {}
    api_error_accessions = set()
    # A Session reuses the HTTP connection for multiple batches where possible.
    with requests.Session() as session:
        # Send each batch in turn; enumerate also gives a human-friendly batch number.
        for batch_number, batch in enumerate(batches, start=1):
            print(
                f"Querying STRING batch {batch_number} of {len(batches)} "
                f"({len(batch)} unique accession(s))..."
            )
            try:
                # Wait between batches; query_string handles waits between retries.
                response_text = query_string(
                    batch,
                    config,
                    session,
                    wait_before_request=batch_number > 1,
                )
                # None means documented no-match HTTP 404; leave these unmapped.
                if response_text is not None:
                    # Parse a successful TSV and add its matches to the results map.
                    mappings.update(process_api_response(response_text, batch))
            # Mark every ID in this batch API_Error if requesting or parsing failed.
            except StringApiError as error:
                api_error_accessions.update(batch)
                # Report the reason instead of silently turning an API failure into
                # an Unmapped biological result.
                print(f"API error for this batch: {error}")

    # Use one timezone-aware timestamp for every row in this mapping run.
    retrieved_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    # Restore original rows, duplicates, values, and order; assign final statuses.
    results = match_results(
        input_data,
        validation,
        mappings,
        api_error_accessions,
        config.string_version,
        retrieved_at_utc,
    )
    # Write the completed table only after all batches have been processed.
    saved_path = export_results(results, output_path)

    # Print the file location and counts so the user can check the run in the terminal.
    print(f"Results saved to: {saved_path}")
    print("Mapping status counts:")
    print(results["Mapping_Status"].value_counts().to_string())
    print(
        f"STRING version recorded: {config.string_version}; "
        f"retrieved at (UTC): {retrieved_at_utc}"
    )
    # Return the DataFrame too, so another Python script can use the results directly.
    return results


def main() -> None:
    """Run the CSV workflow from the command line."""
    # argparse creates a beginner-friendly help page and checks the supplied options.
    parser = argparse.ArgumentParser(
        description="Map UniProt accessions in a CSV file to STRING protein IDs."
    )
    # The CSV path is required because the command-line interface handles CSV input.
    parser.add_argument("input_csv", type=Path, help="CSV file containing accessions.")
    # This option is needed for CSVs with more than one column.
    parser.add_argument(
        "--column",
        dest="input_column",
        help="CSV column containing accessions (required for multi-column CSV files).",
    )
    # Let the user choose a different output location; otherwise use the default.
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT_PATH}).",
    )
    # Parse the options typed after the script name in the terminal.
    arguments = parser.parse_args()
    # Hand the parsed values to the same workflow used by Python callers.
    run_mapping(
        arguments.input_csv,
        input_column=arguments.input_column,
        output_path=arguments.output,
    )


if __name__ == "__main__":
    # Run the command-line entry point only when this file itself is executed.
    # Importing its functions from another script will not start a mapping run.
    main()