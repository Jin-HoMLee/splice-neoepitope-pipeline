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
first-hand in the 2026-07-04 venue audit recorded in PROVENANCE.md; all 14
current sources are peer-reviewed journal articles (0 preprints folded to date -
the three preprints encountered were each deferred for sequence-unavailability,
not because they are preprints).

This script OWNS the venue_type column: re-running overwrites every value from
source. A source NOT in the map below derives to the out-of-vocab sentinel
`unclassified`, which validate_registry.py rejects - so a newly-folded source
cannot slip in venue-unmarked. When you fold a new source, add it here with its
audited venue (never edit registry.tsv directly, or the next run reverts it).

Run: research/.venv/bin/python derive_venue_type.py
"""
import pandas as pd
from pathlib import Path

from labeling_constants import VENUE_UNCLASSIFIED

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"

# source-substring -> venue_type. All current sources are peer-reviewed journals
# (2026-07-04 audit, PROVENANCE.md). Keyed on a lowercased substring of `source`.
VENUE_BY_SOURCE_SUBSTR = {
    "bigot": "journal",       # Bigot 2021, JCI Insight (SF3B1 uveal melanoma)
    "kim 2025": "journal",    # Kim 2025, Cell (SF-mutant leukemia)
    "manoharan": "journal",   # Manoharan 2026, OA journal (IR-CRC)
    "merlotti": "journal",    # Merlotti 2023, Sci Immunol (NSCLC exon-TE)
    "snaf": "journal",        # SNAF / Li 2024, Sci Transl Med (the folded rows; the SNAF preprint was deferred)
    "long-read": "journal",   # Long-read UM 2023, Cancer Immunol Res
    "iris": "journal",        # IRIS (Pan/Xing), PNAS
    "fisher": "journal",      # Fisher 2026 (CoREST)
    "kwok": "journal",        # Kwok 2024, Nature
    "xiong": "journal",       # Xiong 2025 (GBM)
    "postn": "journal",       # POSTN-203 study
}


def venue_type(r):
    src = str(r["source"]).lower()
    for substr, venue in VENUE_BY_SOURCE_SUBSTR.items():
        if substr in src:
            return venue
    # unmapped source: force classification rather than silently defaulting.
    return VENUE_UNCLASSIFIED


df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")
df["venue_type"] = df.apply(venue_type, axis=1)
df.to_csv(REG, sep="\t", index=False)

print(df["venue_type"].value_counts().to_string())
unclassified = int((df["venue_type"] == VENUE_UNCLASSIFIED).sum())
if unclassified:
    print(f"WARNING: {unclassified} row(s) with unclassified venue_type "
          f"(source not in VENUE_BY_SOURCE_SUBSTR) - validate_registry.py will fail.")
