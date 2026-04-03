#!/usr/bin/env python3
"""download_gdc_data.py — Download TCGA splice-junction quantification files
from the GDC Data Portal API.

This script is called by two Snakemake rules:
  * download_gdc_manifest  – queries the API and writes a manifest TSV
  * download_gdc_files     – downloads every file listed in the manifest

The Snakemake rule that invokes this script sets ``snakemake.rule`` to either
``download_gdc_manifest`` or ``download_gdc_files`` so the script knows which
action to perform.  When run outside Snakemake (e.g., for testing), supply
``--mode manifest|files`` on the command line.

**Authentication**: TCGA data on GDC is controlled-access and requires a valid
GDC token. Obtain one from https://portal.gdc.cancer.gov/ (login → Token →
Download). Set the path to the downloaded token file in config/config.yaml
under ``gdc.token_file``.

GDC API reference: https://docs.gdc.cancer.gov/API/Users_Guide/
"""

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GDC query helpers
# ---------------------------------------------------------------------------
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

_FIELDS = [
    "file_id",
    "file_name",
    "cases.samples.sample_type",
    "cases.project.project_id",
    "md5sum",
    "file_size",
]


def _load_gdc_token(token_file: str | Path | None) -> str | None:
    """Load GDC authentication token from file.

    Args:
        token_file: Path to the GDC token file, or None.

    Returns:
        The token string (stripped of whitespace), or None if no file provided.

    Raises:
        FileNotFoundError: If the token file path is provided but does not exist.
    """
    if token_file is None:
        return None
    token_path = Path(token_file).expanduser()
    if not token_path.exists():
        raise FileNotFoundError(
            f"GDC token file not found: {token_path}\n"
            "Obtain a token from https://portal.gdc.cancer.gov/ (login → Token → Download)\n"
            "and set the path in config/config.yaml under gdc.token_file"
        )
    token = token_path.read_text().strip()
    if not token:
        raise ValueError(f"GDC token file is empty: {token_path}")
    log.info("Loaded GDC authentication token from %s", token_path)
    return token


def _build_filters(project_id: str) -> dict:
    """Return GDC API filter JSON for splice-junction quantification files."""
    return {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": ["Splice Junction Quantification"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [project_id],
                },
            },
        ],
    }


def fetch_manifest(
    project_id: str,
    files_endpoint: str = GDC_FILES_ENDPOINT,
    max_files: int | None = None,
    page_size: int = 100,
) -> list[dict]:
    """Query the GDC API and return a list of file metadata dicts.

    Args:
        project_id:      TCGA project id, e.g. ``TCGA-BRCA``.
        files_endpoint:  GDC files API endpoint URL.
        max_files:       Cap on the number of files returned (``None`` for all).
        page_size:       Number of files to request per API page.

    Returns:
        List of file metadata dicts with keys: file_id, file_name,
        sample_type, project_id.
    """
    filters = _build_filters(project_id)
    records: list[dict] = []
    from_offset = 0

    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(_FIELDS),
            "format": "JSON",
            "size": page_size,
            "from": from_offset,
        }
        log.info(
            "Querying GDC files endpoint: project=%s offset=%d", project_id, from_offset
        )
        try:
            response = requests.get(files_endpoint, params=params, timeout=60)
            response.raise_for_status()
        except requests.RequestException as exc:
            log.error("GDC API request failed: %s", exc)
            raise

        payload = response.json()
        hits = payload.get("data", {}).get("hits", [])
        if not hits:
            break

        for hit in hits:
            file_id = hit.get("file_id", "")
            file_name = hit.get("file_name", "")
            cases = hit.get("cases", [{}])
            project = cases[0].get("project", {}).get("project_id", project_id)
            samples = cases[0].get("samples", [{}])
            sample_type = samples[0].get("sample_type", "Unknown") if samples else "Unknown"

            records.append(
                {
                    "file_id": file_id,
                    "file_name": file_name,
                    "sample_type": sample_type,
                    "project_id": project,
                }
            )

        from_offset += len(hits)
        pagination = payload.get("data", {}).get("pagination", {})
        total = pagination.get("total", 0)
        log.info("Fetched %d / %d records for %s", from_offset, total, project_id)

        if from_offset >= total:
            break
        if max_files and from_offset >= max_files:
            records = records[:max_files]
            break

        time.sleep(0.2)  # polite rate-limiting

    log.info("Total records for %s: %d", project_id, len(records))
    return records


def write_manifest(records: list[dict], manifest_path: str | Path) -> None:
    """Write manifest records to a TSV file.

    Args:
        records:       List of dicts returned by :func:`fetch_manifest`.
        manifest_path: Destination path for the TSV.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w") as fh:
        fh.write("file_id\tfile_name\tsample_type\tproject_id\n")
        for rec in records:
            fh.write(
                f"{rec['file_id']}\t{rec['file_name']}\t{rec['sample_type']}\t{rec['project_id']}\n"
            )
    log.info("Manifest written to %s (%d files)", manifest_path, len(records))


def download_files(
    manifest_path: str | Path,
    output_dir: str | Path,
    data_endpoint: str = GDC_DATA_ENDPOINT,
    token_file: str | Path | None = None,
    chunk_size: int = 1 << 20,
) -> None:
    """Download every file listed in the manifest using the GDC data endpoint.

    Files that already exist on disk (by file_id sub-directory) are skipped.

    **Authentication**: TCGA Splice Junction Quantification files are controlled-
    access data and require a valid GDC authentication token. Obtain one from
    https://portal.gdc.cancer.gov/ (login with eRA Commons → Token → Download).

    Args:
        manifest_path:  Path to the manifest TSV.
        output_dir:     Directory in which to save downloaded files.
        data_endpoint:  GDC data API endpoint URL.
        token_file:     Path to GDC token file (required for controlled-access data).
        chunk_size:     HTTP streaming chunk size in bytes.

    Raises:
        FileNotFoundError: If token_file is set but does not exist.
    """
    import csv

    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load authentication token
    token = _load_gdc_token(token_file)
    headers = {}
    if token:
        headers["X-Auth-Token"] = token
        log.info("Using GDC authentication token for controlled-access downloads")
    else:
        log.warning(
            "No GDC token provided. Downloads will fail for controlled-access data.\n"
            "To fix: obtain a token from https://portal.gdc.cancer.gov/ and set\n"
            "gdc.token_file in config/config.yaml"
        )

    with manifest_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        records = list(reader)

    log.info("Downloading %d files to %s", len(records), output_dir)
    failed_count = 0
    for i, rec in enumerate(records, 1):
        file_id = rec["file_id"]
        file_name = rec.get("file_name", file_id + ".tsv")
        dest = output_dir / (file_id + ".tsv")

        if dest.exists():
            log.debug("[%d/%d] Skipping %s (already exists)", i, len(records), file_id)
            continue

        url = f"{data_endpoint}/{file_id}"
        log.info("[%d/%d] Downloading %s → %s", i, len(records), file_id, dest)
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                r.raise_for_status()
                with dest.open("wb") as out:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        out.write(chunk)
        except requests.RequestException as exc:
            log.error("Failed to download %s: %s", file_id, exc)
            failed_count += 1
            if "403" in str(exc) and not token:
                log.error(
                    "403 Forbidden — this is controlled-access data. "
                    "Provide a GDC token via gdc.token_file in config/config.yaml"
                )
            if dest.exists():
                dest.unlink()
            continue

        time.sleep(0.1)

    if failed_count > 0:
        log.warning(
            "Download completed with %d failures out of %d files.",
            failed_count, len(records)
        )
    else:
        log.info("Download complete — all %d files successful.", len(records))


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------
def _snakemake_main() -> None:
    """Entry point when called from a Snakemake rule via ``script:``."""
    rule_name = snakemake.rule  # type: ignore[name-defined]  # noqa: F821

    if rule_name == "download_gdc_manifest":
        cancer_type = snakemake.wildcards.cancer_type  # type: ignore[name-defined]  # noqa: F821
        manifest_out = snakemake.output.manifest  # type: ignore[name-defined]  # noqa: F821
        files_endpoint = snakemake.params.gdc_files_endpoint  # type: ignore[name-defined]  # noqa: F821
        max_files = snakemake.params.max_files  # type: ignore[name-defined]  # noqa: F821

        records = fetch_manifest(
            project_id=cancer_type,
            files_endpoint=files_endpoint,
            max_files=max_files,
        )
        write_manifest(records, manifest_out)

    elif rule_name == "download_gdc_files":
        manifest_in = snakemake.input.manifest  # type: ignore[name-defined]  # noqa: F821
        data_dir = snakemake.output.data_dir  # type: ignore[name-defined]  # noqa: F821
        done_file = snakemake.output.done  # type: ignore[name-defined]  # noqa: F821
        data_endpoint = snakemake.params.gdc_data_endpoint  # type: ignore[name-defined]  # noqa: F821
        token_file = snakemake.params.get("gdc_token_file", None)  # type: ignore[name-defined]  # noqa: F821

        download_files(
            manifest_path=manifest_in,
            output_dir=data_dir,
            data_endpoint=data_endpoint,
            token_file=token_file,
        )
        Path(done_file).touch()

    else:
        raise ValueError(f"Unexpected rule name: {rule_name!r}")


def _cli_main() -> None:
    """Command-line entry point for standalone use."""
    parser = argparse.ArgumentParser(
        description="Download TCGA splice-junction data from the GDC API."
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    # manifest sub-command
    m = sub.add_parser("manifest", help="Fetch file manifest from GDC API")
    m.add_argument("--project-id", required=True, help="TCGA project ID (e.g. TCGA-BRCA)")
    m.add_argument("--output", required=True, help="Output manifest TSV path")
    m.add_argument(
        "--files-endpoint",
        default=GDC_FILES_ENDPOINT,
        help="GDC files API endpoint",
    )
    m.add_argument("--max-files", type=int, default=None, help="Cap on file count")

    # download sub-command
    d = sub.add_parser("files", help="Download files listed in a manifest")
    d.add_argument("--manifest", required=True, help="Manifest TSV path")
    d.add_argument("--output-dir", required=True, help="Directory for downloaded files")
    d.add_argument(
        "--data-endpoint",
        default=GDC_DATA_ENDPOINT,
        help="GDC data API endpoint",
    )
    d.add_argument(
        "--token-file",
        default=None,
        help="Path to GDC authentication token file (required for controlled-access data)",
    )

    args = parser.parse_args()

    if args.mode == "manifest":
        records = fetch_manifest(
            project_id=args.project_id,
            files_endpoint=args.files_endpoint,
            max_files=args.max_files,
        )
        write_manifest(records, args.output)
    elif args.mode == "files":
        download_files(
            manifest_path=args.manifest,
            output_dir=args.output_dir,
            data_endpoint=args.data_endpoint,
            token_file=args.token_file,
        )


if __name__ == "__main__":
    # Detect whether we are inside a Snakemake execution
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
