#!/usr/bin/env python3
"""Derive the venue_type column (Issue #1001).

venue_type records the *publication venue class* of each row's source, so a
scoring run (and any audit) can query and optionally down-weight non-peer-
reviewed evidence without parsing free-text `source` / PROVENANCE prose.
Standard evidence-grading practice (GRADE, systematic-review methodology)
down-weights non-peer-reviewed sources; a functional ground-truth set used to
benchmark callers should be able to flag preprint-sourced rows at a glance.

The assignment is **source-keyed**: the `source` column is the stable curation
key (same convention as derive_assay_context.py), matched on lowercased
substrings so the minor per-source string variants (e.g. "IRIS" vs
"IRIS (Pan/Xing, PNAS)") map to one venue. The per-source venue was resolved
first-hand in the 2026-07-04 venue audit recorded in PROVENANCE.md; all 11
current studies (14 distinct `source` strings) are peer-reviewed journal
articles (0 preprints folded to date - the three preprints encountered were
each deferred for sequence-unavailability, not because they are preprints).

This script OWNS the venue_type column: re-running overwrites every value from
source. When you fold a new source, add it to VENUE_BY_SOURCE_SUBSTR
(labeling_constants.py) with its audited venue (never edit registry.tsv directly,
or the next run reverts it).

Two guards catch a mis-marked venue (both in validate_registry.py):
  1. A source matching NO key derives to the out-of-vocab sentinel `unclassified`,
     which the validator rejects - so a genuinely new source cannot slip in
     venue-unmarked.
  2. Because a study-substring key CANNOT self-distinguish a preprint from the
     journal version of an *already-mapped* study (a future "SNAF ... bioRxiv"
     source would first-match "snaf" -> journal), the validator separately fails
     any row whose `source` names a preprint server (PREPRINT_MARKERS) but whose
     venue_type is not `preprint`. The offline-optional Zotero cross-check is a
     third, independent backstop when creds are present.

Run: research/.venv/bin/python derive_venue_type.py
"""
import pandas as pd
from pathlib import Path

from labeling_constants import VENUE_BY_SOURCE_SUBSTR, VENUE_UNCLASSIFIED

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"


def venue_type(r):
    """Source-keyed venue_type. VENUE_BY_SOURCE_SUBSTR (labeling_constants.py) is
    the single source of truth, matched on a lowercased substring of `source`."""
    src = str(r["source"]).lower()
    for substr, spec in VENUE_BY_SOURCE_SUBSTR.items():
        if substr in src:
            return spec["venue"]
    # unmapped source: force classification rather than silently defaulting.
    return VENUE_UNCLASSIFIED


def main():
    df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")
    df["venue_type"] = df.apply(venue_type, axis=1)
    df.to_csv(REG, sep="\t", index=False)

    print(df["venue_type"].value_counts().to_string())
    unclassified = int((df["venue_type"] == VENUE_UNCLASSIFIED).sum())
    if unclassified:
        print(f"WARNING: {unclassified} row(s) with unclassified venue_type "
              f"(source not in VENUE_BY_SOURCE_SUBSTR) - validate_registry.py will fail.")


if __name__ == "__main__":
    main()
