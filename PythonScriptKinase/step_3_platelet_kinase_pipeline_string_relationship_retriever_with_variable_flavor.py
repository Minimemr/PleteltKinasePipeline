"""Retrieve STRING interaction networks for a set of protein identifiers."""

# Postpone evaluation of annotations so type references remain straightforward.
from __future__ import annotations

# StringIO presents the API's response text as a file-like object for pandas.
from io import StringIO
# Any describes identifier values before they are normalized to strings.
# Sequence documents that the function accepts an ordered collection of IDs.
from typing import Any, Sequence

# pandas parses TSV responses and detects missing identifier values.
import pandas as pd
# requests sends the API POST and exposes HTTP, connection, and timeout errors.
import requests


# Pin requests to STRING 12.5 so future database releases do not change results.
STRING_NETWORK_URL = "https://version-12-5.string-db.org/api/tsv/network"
# Use the NCBI taxonomy identifier for Homo sapiens.
SPECIES = 9606
# Ask STRING for functional associations for both requested network flavors.
NETWORK_TYPE = "functional"
# Identify this dissertation pipeline to the STRING service.
CALLER_IDENTITY = "reading_platelet_kinase_pipeline"
# Stop waiting if the API request takes longer than this number of seconds.
REQUEST_TIMEOUT_SECONDS = 30
# Accept only the two network flavors requested by the pipeline.
VALID_FLAVOURS = {"confidence", "typed"}
# Require the two IDs and preferred names shared by both flavor response schemas.
REQUIRED_NETWORK_COLUMNS = {
    "stringId_A",
    "stringId_B",
    "preferredName_A",
    "preferredName_B",
}
# STRING uses a different combined-score header for the typed flavor.
SCORE_COLUMN_BY_FLAVOUR = {
    "confidence": "score",
    "typed": "functional_combined_score",
}


class StringNetworkError(RuntimeError):
    """Represent a failed request or an invalid STRING network response."""


def retrieve_string_network(
    string_ids: Sequence[Any],  # Mapped STRING identifiers supplied by the caller.
    flavour: str,  # Selects "confidence" or "typed" for this individual request.
    required_score: int = 400,  # Default to STRING's medium confidence threshold.
) -> pd.DataFrame:  # Return all parsed interaction rows and response columns.
    """Return a STRING functional-association network in the requested flavor.

    STRING associations alone do not establish direct binding or
    kinase-substrate relationships.
    """
    # Reject misspellings and unsupported flavors before contacting the API.
    if not isinstance(flavour, str) or flavour not in VALID_FLAVOURS:
        raise ValueError(
            f"flavour must be one of {', '.join(sorted(VALID_FLAVOURS))}; "
            f"received {flavour!r}."
        )
    # Reject a single string (which is itself a sequence) and non-sequence inputs.
    if isinstance(string_ids, (str, bytes)) or not isinstance(string_ids, Sequence):
        raise TypeError("string_ids must be a sequence of STRING identifiers.")
    # Require an integer score while excluding bool, which is an int subclass.
    if isinstance(required_score, bool) or not isinstance(required_score, int):
        raise TypeError("required_score must be an integer from 0 to 1000.")
    # STRING confidence scores use the inclusive range from 0 to 1000.
    if not 0 <= required_score <= 1000:
        raise ValueError("required_score must be between 0 and 1000.")

    # Keep cleaned identifiers in input order and use a set for duplicate checks.
    identifiers = []
    seen = set()
    # Examine each original value so missing and duplicate values can be removed.
    for value in string_ids:
        # Skip None and pandas-style missing values before converting to text.
        if value is None or pd.isna(value):
            continue
        # Strip whitespace so accidental spaces do not become part of a STRING ID.
        identifier = str(value).strip()
        # Keep each non-empty identifier only on its first appearance.
        if identifier and identifier not in seen:
            identifiers.append(identifier)
            seen.add(identifier)

    # A network query requires at least two distinct, non-missing identifiers.
    if len(identifiers) < 2:
        raise ValueError(
            "At least two distinct, non-missing STRING identifiers are required."
        )

    # Build the form fields sent in the POST body to the STRING network endpoint.
    parameters = {
        # Separate multiple IDs with carriage returns as accepted by STRING POST.
        "identifiers": "\r".join(identifiers),
        # Restrict the search to human proteins.
        "species": SPECIES,
        # Apply the caller's threshold to both flavors.
        "required_score": required_score,
        # Keep functional associations distinct from the flavor selection below.
        "network_type": NETWORK_TYPE,
        # Pass the caller-selected flavor rather than hardcoding one flavor.
        "network_flavor": flavour,
        # Do not expand the input set with additional interaction partners.
        "add_nodes": 0,
        # Identify the client making this API request.
        "caller_identity": CALLER_IDENTITY,
    }

    # Send the request and convert network/HTTP problems into a clear tool error.
    try:
        # Use POST so the multi-protein identifier list is sent in the request body.
        response = requests.post(
            STRING_NETWORK_URL,  # Use the version-pinned TSV network endpoint.
            data=parameters,  # Submit all STRING options as form data.
            timeout=REQUEST_TIMEOUT_SECONDS,  # Prevent an indefinite wait.
        )
        # Raise an HTTPError for non-success status codes before parsing the body.
        response.raise_for_status()
    # Catch requests-specific failures, including HTTP errors, timeouts, and connections.
    except requests.RequestException as error:
        # Extract the status when this exception was caused by an HTTP response.
        status = getattr(getattr(error, "response", None), "status_code", None)
        # Include the status and a short server response when available.
        if status is not None:
            detail = error.response.text[:500]
            raise StringNetworkError(
                f"STRING network request failed with HTTP {status}: {detail}"
            ) from error
        # Report connection and timeout failures without hiding their original cause.
        raise StringNetworkError(
            f"Could not request the STRING network: {error}"
        ) from error

    # Treat a successful HTTP response with no body as an invalid API response.
    if not response.text.strip():
        raise StringNetworkError("STRING returned an empty network response.")

    # Parse the tab-separated response while retaining all columns STRING returned.
    try:
        network = pd.read_csv(StringIO(response.text), sep="\t")
    # Convert empty or malformed TSV parser errors into a STRING-specific message.
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise StringNetworkError(
            f"Could not parse STRING's TSV network response: {error}"
        ) from error

    # Require the score header specific to the requested response flavor.
    required_columns = REQUIRED_NETWORK_COLUMNS | {
        SCORE_COLUMN_BY_FLAVOUR[flavour]
    }
    # Find fields that would indicate an unexpected response schema.
    missing_columns = required_columns.difference(network.columns)
    # Report the exact expected network fields that were absent.
    if missing_columns:
        raise StringNetworkError(
            "STRING returned TSV data missing expected network column(s): "
            + ", ".join(sorted(missing_columns))
        )
    # A valid header without rows means STRING returned no interactions for this query.
    if network.empty:
        raise StringNetworkError("STRING returned valid TSV data but no interactions.")

    # Return the complete DataFrame, including combined and any evidence scores.
    return network
