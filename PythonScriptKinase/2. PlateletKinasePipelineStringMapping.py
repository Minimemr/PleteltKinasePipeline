"""Map protein accessions to STRING identifiers and save a CSV report."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
import re
import time
from typing import Any, Optional, Sequence, Union

import pandas as pd
import requests


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_PATH = PROJECT_DIR / "Outputs" / "STRING_Mapping.csv"

REQUIRED_RESPONSE_COLUMNS = {
    "queryItem",
    "queryIndex",
    "stringId",
    "ncbiTaxonId",
    "taxonName",
    "preferredName",
    "annotation",
}
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
ACCESSION_PATTERN = re.compile(r"[A-Z0-9]{6,10}(?:-[1-9][0-9]*)?")


@dataclass(frozen=True)
class StringMappingConfig:
    """Settings that can be changed without editing the mapping workflow."""

    api_url: str = "https://version-12-5.string-db.org/api/tsv/get_string_ids"
    string_version: str = "12.5"
    organism_id: int = 9606
    caller_identity: str = "reading_platelet_kinase_pipeline"
    batch_size: int = 100
    timeout_seconds: float = 30.0
    max_retries: int = 3
    request_interval_seconds: float = 1.0
    retry_backoff_seconds: float = 2.0

    def __post_init__(self) -> None:
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
    if isinstance(source, (str, Path)):
        input_path = Path(source)
        if not input_path.is_file():
            raise FileNotFoundError(f"Could not find the input CSV: {input_path}")

        data = pd.read_csv(
            input_path,
            dtype=str,
            keep_default_na=False,
            skip_blank_lines=False,
        )
        if input_column is None:
            if len(data.columns) != 1:
                raise ValueError(
                    "This CSV has more than one column. Pass input_column with "
                    "the column containing UniProt accessions."
                )
            input_column = str(data.columns[0])
        if input_column not in data.columns:
            raise ValueError(
                f"Input column {input_column!r} was not found. "
                f"Columns found: {', '.join(map(str, data.columns))}"
            )
        values = data[input_column].tolist()
    else:
        if isinstance(source, (bytes, bytearray)) or not isinstance(source, Sequence):
            raise TypeError("source must be a CSV path or a sequence of accessions.")
        values = list(source)

    return pd.DataFrame(
        {
            "UniProt_Accession": values,
            "_input_order": range(len(values)),
        }
    )


def validate_accessions(accessions: Sequence[Any]) -> pd.DataFrame:
    """Return cleaned values and validity flags without changing the original inputs."""
    normalised_values = []
    valid_values = []

    for accession in accessions:
        if pd.isna(accession):
            normalised = ""
        else:
            normalised = str(accession).strip().upper()

        normalised_values.append(normalised)
        valid_values.append(bool(ACCESSION_PATTERN.fullmatch(normalised)))

    return pd.DataFrame(
        {
            "_normalised_accession": normalised_values,
            "_is_valid": valid_values,
        }
    )


def remove_duplicate_accessions(accessions: Sequence[str]) -> list[str]:
    """Keep the first occurrence of each accession, preserving input order."""
    return list(dict.fromkeys(accessions))


def split_into_batches(
    accessions: Sequence[str], batch_size: int
) -> list[list[str]]:
    """Split accessions into request-sized groups."""
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1.")
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
    if response is not None:
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                return max(config.request_interval_seconds, float(retry_after))
            except ValueError:
                pass
    backoff = config.retry_backoff_seconds * (2**attempt)
    return max(config.request_interval_seconds, backoff)


def query_string(
    accessions: Sequence[str],
    config: StringMappingConfig,
    session: requests.Session,
    wait_before_request: bool = False,
) -> Optional[str]:
    """POST one batch to STRING; return TSV text, or None for STRING's no-match 404."""
    if wait_before_request:
        time.sleep(config.request_interval_seconds)

    parameters = {
        "identifiers": "\r".join(accessions),
        "species": config.organism_id,
        "echo_query": 1,
        "limit": 1,
        "caller_identity": config.caller_identity,
    }
    retryable_statuses = {429, 500, 502, 503, 504}

    for attempt in range(config.max_retries + 1):
        try:
            response = session.post(
                config.api_url,
                data=parameters,
                timeout=config.timeout_seconds,
            )
        except (requests.Timeout, requests.ConnectionError) as error:
            if attempt >= config.max_retries:
                raise StringApiError(
                    f"STRING request failed after {attempt + 1} attempt(s): {error}"
                ) from error
            time.sleep(_retry_delay(attempt, config))
            continue
        except requests.RequestException as error:
            raise StringApiError(f"Could not send STRING request: {error}") from error

        if response.status_code == 404:
            # STRING documents 404 when none of the submitted identifiers resolve.
            return None

        if response.status_code in retryable_statuses:
            if attempt >= config.max_retries:
                raise StringApiError(
                    f"STRING returned HTTP {response.status_code} after "
                    f"{attempt + 1} attempt(s): {response.text[:500]}"
                )
            time.sleep(_retry_delay(attempt, config, response))
            continue

        try:
            response.raise_for_status()
        except requests.HTTPError as error:
            raise StringApiError(
                f"STRING returned HTTP {response.status_code}: "
                f"{response.text[:500]}"
            ) from error

        return response.text

    raise StringApiError("STRING request ended without a response.")


def process_api_response(
    response_text: str, requested_accessions: Sequence[str]
) -> dict[str, dict[str, str]]:
    """Parse STRING's TSV and associate each result with its submitted accession."""
    if not response_text.strip():
        raise StringApiError("STRING returned an empty response with HTTP 200.")

    try:
        response_data = pd.read_csv(StringIO(response_text), sep="\t", dtype=str)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise StringApiError(f"Could not parse STRING's TSV response: {error}") from error

    missing_columns = REQUIRED_RESPONSE_COLUMNS.difference(response_data.columns)
    if missing_columns:
        raise StringApiError(
            "STRING response is missing expected column(s): "
            + ", ".join(sorted(missing_columns))
        )
    if response_data.empty:
        raise StringApiError("STRING returned a TSV header but no mapping rows.")

    mappings = {}
    seen_query_indexes = set()
    for _, row in response_data.iterrows():
        try:
            query_index = int(row["queryIndex"])
        except (TypeError, ValueError) as error:
            raise StringApiError(
                f"STRING returned an invalid queryIndex: {row['queryIndex']!r}"
            ) from error

        if not 0 <= query_index < len(requested_accessions):
            raise StringApiError(
                f"STRING returned queryIndex {query_index}, outside this batch."
            )
        if query_index in seen_query_indexes:
            raise StringApiError(
                f"STRING returned more than one result for queryIndex {query_index}."
            )
        seen_query_indexes.add(query_index)

        accession = requested_accessions[query_index]
        echoed_accession = str(row["queryItem"]).strip().upper()
        if echoed_accession != accession:
            raise StringApiError(
                f"STRING echoed {row['queryItem']!r} for submitted accession "
                f"{accession!r}."
            )

        mappings[accession] = {
            "STRING_ID": row["stringId"],
            "Preferred_Name": row["preferredName"],
            "Organism_ID": row["ncbiTaxonId"],
            "Organism": row["taxonName"],
            "Annotation": row["annotation"],
        }

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
    results = []
    for position, row in input_data.iterrows():
        original_accession = row["UniProt_Accession"]
        normalised_accession = validation.at[position, "_normalised_accession"]
        is_valid = validation.at[position, "_is_valid"]

        result = {column: pd.NA for column in OUTPUT_COLUMNS}
        result["UniProt_Accession"] = original_accession
        result["STRING_Version"] = string_version
        result["Retrieved_At_UTC"] = retrieved_at_utc

        if not is_valid:
            result["Mapping_Status"] = "Invalid_Input"
        elif normalised_accession in api_error_accessions:
            result["Mapping_Status"] = "API_Error"
        elif normalised_accession in mappings:
            result.update(mappings[normalised_accession])
            result["Mapping_Status"] = "Mapped"
        else:
            result["Mapping_Status"] = "Unmapped"

        results.append(result)

    return pd.DataFrame(results, columns=OUTPUT_COLUMNS)


def export_results(results: pd.DataFrame, output_path: Union[str, Path]) -> Path:
    """Write the mapping table to a CSV without changing any input file."""
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(destination, index=False, na_rep="")
    return destination


def run_mapping(
    source: Union[str, Path, Sequence[Any]],
    input_column: Optional[str] = None,
    output_path: Union[str, Path] = DEFAULT_OUTPUT_PATH,
    config: StringMappingConfig = StringMappingConfig(),
) -> pd.DataFrame:
    """Run input loading, validation, STRING lookup, matching, and CSV export."""
    input_data = load_input(source, input_column)

    if isinstance(source, (str, Path)):
        input_path = Path(source).resolve()
        if Path(output_path).resolve() == input_path:
            raise ValueError("The output CSV must not overwrite the input CSV.")

    validation = validate_accessions(input_data["UniProt_Accession"].tolist())
    valid_accessions = [
        accession
        for accession, is_valid in zip(
            validation["_normalised_accession"], validation["_is_valid"]
        )
        if is_valid
    ]
    valid_accession_counts = pd.Series(valid_accessions, dtype="string").value_counts()
    duplicate_accessions = valid_accession_counts[valid_accession_counts > 1]
    duplicate_rows = int((duplicate_accessions - 1).sum())
    unique_accessions = remove_duplicate_accessions(valid_accessions)
    batches = split_into_batches(unique_accessions, config.batch_size)

    print(
        f"Input rows: {len(input_data):,}; valid accession rows: "
        f"{len(valid_accessions):,}; unique accessions to query: "
        f"{len(unique_accessions):,}; repeated rows: {duplicate_rows:,} "
        f"across {len(duplicate_accessions):,} accession(s)."
    )

    mappings = {}
    api_error_accessions = set()
    with requests.Session() as session:
        for batch_number, batch in enumerate(batches, start=1):
            print(
                f"Querying STRING batch {batch_number} of {len(batches)} "
                f"({len(batch)} unique accession(s))..."
            )
            try:
                response_text = query_string(
                    batch,
                    config,
                    session,
                    wait_before_request=batch_number > 1,
                )
                if response_text is not None:
                    mappings.update(process_api_response(response_text, batch))
            except StringApiError as error:
                api_error_accessions.update(batch)
                print(f"API error for this batch: {error}")

    retrieved_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    results = match_results(
        input_data,
        validation,
        mappings,
        api_error_accessions,
        config.string_version,
        retrieved_at_utc,
    )
    saved_path = export_results(results, output_path)

    print(f"Results saved to: {saved_path}")
    print("Mapping status counts:")
    print(results["Mapping_Status"].value_counts().to_string())
    print(
        f"STRING version recorded: {config.string_version}; "
        f"retrieved at (UTC): {retrieved_at_utc}"
    )
    return results


def main() -> None:
    """Run the CSV workflow from the command line."""
    parser = argparse.ArgumentParser(
        description="Map UniProt accessions in a CSV file to STRING protein IDs."
    )
    parser.add_argument("input_csv", type=Path, help="CSV file containing accessions.")
    parser.add_argument(
        "--column",
        dest="input_column",
        help="CSV column containing accessions (required for multi-column CSV files).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT_PATH}).",
    )
    arguments = parser.parse_args()
    run_mapping(
        arguments.input_csv,
        input_column=arguments.input_column,
        output_path=arguments.output,
    )


if __name__ == "__main__":
    main()