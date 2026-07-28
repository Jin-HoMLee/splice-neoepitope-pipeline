#!/usr/bin/env python3
"""Pull the human-source T-cell-assayed epitope set from the IEDB IQ-API.

Stage 1 of the Issue #1290 reverse-discovery pipeline: take *every* T-cell-assayed
epitope with a human source antigen, regardless of how its splice origin was
annotated, so later stages can remap the peptide sequence and test junction-crossing
independently. This is the "ceiling" complement to the annotation-based mine in
Issue #734, which can only find epitopes whose authors already knew were splice-derived.

The API is PostgREST, public, unauthenticated:
    https://query-api.iedb.org/tcell_search

Two API behaviours worth knowing before changing this:

  - `offset` without `order` is REFUSED (HTTP 400), because paging without a sort key
    is not stable. We always send `order=tcell_id.asc`.
  - the server caps a page at 10,000 rows.

Verified against the live API 2026-07-28: 66,358 human-source assays, 66,353 carrying
a `parent_source_antigen_name` anchor, split 33,674 positive / 32,684 negative.
"""

import argparse
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

API_ROOT = "https://query-api.iedb.org"
ENDPOINT = "tcell_search"

# The tumour/self universe. IEDB is mostly pathogen data, so without this filter the
# pool is 576,931 assays of which the overwhelming majority are irrelevant here.
HUMAN_SOURCE_FILTER = ("parent_source_antigen_source_org_name", "ilike.*Homo sapiens*")

PAGE_SIZE = 10_000
ORDER_KEY = "tcell_id"

# Only the columns later stages need: the peptide, the anchor that makes
# peptide -> transcript mapping tractable, the assay outcome, and provenance.
SELECT_COLUMNS = [
    "tcell_id",
    "linear_sequence",
    "linear_sequence_length",
    "structure_type",
    "qualitative_measure",
    "parent_source_antigen_name",
    "parent_source_antigen_iri",
    "parent_source_antigen_source_org_name",
    "mhc_allele_name",
    "mhc_class",
    "mhc_allele_resolution",
    "assay_names",
    "host_organism_name",
    "disease_names",
    "pubmed_id",
    "reference_id",
]

# IEDB's full `qualitative_measure` vocabulary for this pool. Confirmed exhaustive on
# 2026-07-28: the five values sum to exactly the unfiltered total (66,358), so there
# are no nulls and no sixth category hiding in the pool.
OUTCOME_BY_MEASURE = {
    "Positive": "positive",
    "Positive-High": "positive",
    "Positive-Intermediate": "positive",
    "Positive-Low": "positive",
    "Negative": "negative",
}

LIST_COLUMNS = {"assay_names", "disease_names"}
LIST_JOIN = "|"


class IntegrityError(RuntimeError):
    """A fetched result set does not reconcile with what the server reported."""


def normalise_outcome(qualitative_measure):
    """Map an IEDB `qualitative_measure` onto `positive` / `negative`.

    Raises on anything unrecognised. A permissive default here would silently move
    rows between the two classes, and the negative class is the scarce one this whole
    Issue exists to grow.
    """
    try:
        return OUTCOME_BY_MEASURE[qualitative_measure]
    except KeyError:
        raise IntegrityError(
            f"unrecognised qualitative_measure {qualitative_measure!r}; "
            f"known values are {sorted(OUTCOME_BY_MEASURE)}. IEDB may have added a "
            f"category - decide its class explicitly rather than defaulting."
        ) from None


def plan_pages(total, page_size):
    """Offsets covering `total` rows, with no empty trailing page."""
    return list(range(0, total, page_size))


def assert_no_truncation(fetched, expected_total):
    """Fetched row count must equal the server's own reported total.

    Short means a truncated walk (silent under-report); long means the pages
    overlapped (silent double-count). Both corrupt the yield figure.
    """
    if fetched != expected_total:
        raise IntegrityError(
            f"fetched {fetched} rows but server reported {expected_total}; "
            f"page walk was {'truncated' if fetched < expected_total else 'overlapping'}"
        )


def assert_unique(rows, key):
    """No duplicate `key` across the assembled result set."""
    seen = set()
    duplicates = set()
    for row in rows:
        value = row.get(key)
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    if duplicates:
        raise IntegrityError(
            f"duplicate {key} values in result set: {sorted(duplicates)[:10]}"
        )


def _build_url(offset=None, select=None, count_only=False):
    params = [HUMAN_SOURCE_FILTER]
    if select:
        params.append(("select", ",".join(select)))
    if count_only:
        params.append(("limit", "1"))
    else:
        params.append(("order", f"{ORDER_KEY}.asc"))
        params.append(("limit", str(PAGE_SIZE)))
        params.append(("offset", str(offset)))
    query = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    return f"{API_ROOT}/{ENDPOINT}?{query}"


def fetch_total():
    """Server-reported row count, read off the Content-Range header."""
    url = _build_url(select=["tcell_id"], count_only=True)
    request = urllib.request.Request(url, headers={"Prefer": "count=exact"})
    with urllib.request.urlopen(request, timeout=120) as response:
        content_range = response.headers.get("Content-Range", "")
    if "/" not in content_range:
        raise IntegrityError(f"no Content-Range total in response: {content_range!r}")
    return int(content_range.rsplit("/", 1)[1])


def fetch_page(offset):
    url = _build_url(offset=offset, select=SELECT_COLUMNS)
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=300) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_all(verbose=True):
    """Walk every page, then prove the result set reconciles before returning it."""
    total = fetch_total()
    rows = []
    for offset in plan_pages(total, PAGE_SIZE):
        page = fetch_page(offset)
        rows.extend(page)
        if verbose:
            print(f"  offset {offset:>6}: {len(page):>6} rows (running {len(rows)})", file=sys.stderr)

    assert_no_truncation(len(rows), total)
    assert_unique(rows, key=ORDER_KEY)
    return rows, total


def _flatten(value, column):
    if value is None:
        return ""
    if column in LIST_COLUMNS and isinstance(value, list):
        return LIST_JOIN.join(str(v) for v in value)
    if isinstance(value, list):
        return LIST_JOIN.join(str(v) for v in value)
    return str(value)


def write_tsv(rows, path):
    """Deterministic TSV: sorted by the order key, fixed column order."""
    header = SELECT_COLUMNS + ["outcome"]
    lines = ["\t".join(header)]
    for row in sorted(rows, key=lambda r: r[ORDER_KEY]):
        fields = [_flatten(row.get(c), c).replace("\t", " ").replace("\n", " ") for c in SELECT_COLUMNS]
        fields.append(normalise_outcome(row.get("qualitative_measure")))
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarise(rows):
    counts = {"positive": 0, "negative": 0}
    anchored = 0
    for row in rows:
        counts[normalise_outcome(row.get("qualitative_measure"))] += 1
        if row.get("parent_source_antigen_name"):
            anchored += 1
    return {
        "total": len(rows),
        "positive": counts["positive"],
        "negative": counts["negative"],
        "with_source_antigen_anchor": anchored,
    }


def write_sidecar(path, tsv_path, server_total, summary):
    digest = hashlib.sha256(tsv_path.read_bytes()).hexdigest()
    sidecar = {
        "issue": 1290,
        "stage": "1-fetch",
        "fetched_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "api_root": API_ROOT,
        "endpoint": ENDPOINT,
        "filter": {HUMAN_SOURCE_FILTER[0]: HUMAN_SOURCE_FILTER[1]},
        "order": f"{ORDER_KEY}.asc",
        "page_size": PAGE_SIZE,
        "columns": SELECT_COLUMNS,
        "server_reported_total": server_total,
        "counts": summary,
        "output_tsv": tsv_path.name,
        "output_sha256": digest,
    }
    path.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")
    return sidecar


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs",
        help="directory for the TSV and its provenance sidecar",
    )
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.outdir / "iedb_human_source_tcell_assays.tsv"
    sidecar_path = args.outdir / "iedb_human_source_tcell_assays.provenance.json"

    print(f"Fetching {ENDPOINT} (human source antigens) from {API_ROOT} ...", file=sys.stderr)
    rows, server_total = fetch_all()

    write_tsv(rows, tsv_path)
    summary = summarise(rows)
    write_sidecar(sidecar_path, tsv_path, server_total, summary)

    print(
        f"\n{summary['total']} assays "
        f"({summary['positive']} positive / {summary['negative']} negative; "
        f"{summary['with_source_antigen_anchor']} source-antigen-anchored)\n"
        f"  -> {tsv_path}\n"
        f"  -> {sidecar_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
