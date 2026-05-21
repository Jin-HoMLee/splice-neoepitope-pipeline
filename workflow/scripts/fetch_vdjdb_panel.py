"""Fetch a per-patient VDJdb TCR panel (Issue #204).

Reads a pinned VDJdb release, filters to HomoSapiens + MHCI + paired α/β +
vdjdb.score >= threshold, exact-matches by 4-digit HLA allele to a patient's
alleles.tsv, picks the top-N TCRs per allele (ranked by score, deterministic
tiebreak by donor ID), reconstructs full α/β chains via stitchr, and writes:

    results/{patient_id}/tcr_panel/vdjdb/panel.tsv    — TCRs with full sequences
    results/{patient_id}/tcr_panel/vdjdb/panel_qc.tsv — per-allele coverage status

Lazy-imports pandas + stitchr at first use (per PR #428 pattern) so pytest
collection stays fast.
"""

from __future__ import annotations

import logging
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Allele normalization
# ---------------------------------------------------------------------------

def normalize_allele_to_4digit(allele: str) -> Optional[str]:
    """Truncate an HLA allele string to its 4-digit form (`HLA-X*GG:PP`).

    Anything shorter than 4-digit returns None — we require exact 4-digit
    matching, and shorter forms (e.g. `HLA-A*02`) are excluded.

    Examples:
        HLA-A*02:01      -> HLA-A*02:01
        HLA-A*02:01:110  -> HLA-A*02:01
        HLA-A*02         -> None
        ""               -> None
    """
    if not allele or ":" not in allele:
        return None
    parts = allele.split(":")
    if len(parts) < 2:
        return None
    return f"{parts[0]}:{parts[1]}"


# ---------------------------------------------------------------------------
# VDJdb load + filter
# ---------------------------------------------------------------------------

def load_and_filter_vdjdb(vdjdb_full_tsv, min_score: int):
    """Load vdjdb_full.txt and apply HomoSapiens + MHCI + score filters.

    Adds a `mhc.a_4digit` column with 4-digit-normalized allele.
    Drops rows where the allele cannot be normalized to 4-digit.
    Returns a pandas DataFrame.
    """
    import pandas as pd  # lazy import

    df = pd.read_csv(vdjdb_full_tsv, sep="\t", dtype=str)
    df["vdjdb.score"] = pd.to_numeric(df["vdjdb.score"], errors="coerce")
    df = df[
        (df["species"] == "HomoSapiens")
        & (df["mhc.class"] == "MHCI")
        & (df["vdjdb.score"] >= min_score)
    ].copy()
    df["mhc.a_4digit"] = df["mhc.a"].apply(normalize_allele_to_4digit)
    df = df[df["mhc.a_4digit"].notna()].copy()
    log.info(
        "Loaded VDJdb (%s): %d rows after filters (HomoSapiens + MHCI + score>=%d)",
        vdjdb_full_tsv, len(df), min_score,
    )
    return df
