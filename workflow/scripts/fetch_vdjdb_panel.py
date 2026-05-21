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


# ---------------------------------------------------------------------------
# Per-allele top-N selection
# ---------------------------------------------------------------------------

def select_top_n_for_allele(df, allele: str, n: int):
    """Filter `df` to exact 4-digit matches for `allele`, sort by
    (vdjdb.score DESC, meta.subject.id ASC) for deterministic tiebreak,
    return top `n` rows.
    """
    matched = df[df["mhc.a_4digit"] == allele].copy()
    matched = matched.sort_values(
        ["vdjdb.score", "meta.subject.id"],
        ascending=[False, True],
        kind="mergesort",  # stable
    )
    return matched.head(n).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Panel status classification
# ---------------------------------------------------------------------------

def classify_panel_status(n_in_panel: int, target_size: int) -> str:
    """Return 'ok' if panel reached target_size, 'low_coverage' if 1..target_size-1,
    'empty' if 0.
    """
    if n_in_panel <= 0:
        return "empty"
    if n_in_panel < target_size:
        return "low_coverage"
    return "ok"


# ---------------------------------------------------------------------------
# stitchr wrapper
# ---------------------------------------------------------------------------

def stitch_chain(v_gene: str, j_gene: str, cdr3: str, chain: str) -> Optional[str]:
    """Reconstruct a full TCR chain (V + CDR3 + J + framework) via stitchr.

    Returns the AA sequence string, or None if stitchr fails for any reason
    (caller is expected to log and skip).

    `chain` is 'A' for alpha (TRA) or 'B' for beta (TRB).
    """
    try:
        from Stitchr import stitchrfunctions, stitchr  # noqa: F401  presence-check
    except ImportError:
        log.error("stitchr not installed in this environment")
        return None

    import subprocess

    try:
        cmd = [
            "stitchr",
            "-v", v_gene,
            "-j", j_gene,
            "-cdr3", cdr3,
            "-species", "HUMAN",
            "-m", "AA",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if proc.returncode != 0:
            log.error("stitchr failed (chain=%s, V=%s, J=%s, CDR3=%s): %s",
                      chain, v_gene, j_gene, cdr3, proc.stderr.strip())
            return None
        seq_lines = [ln.strip() for ln in proc.stdout.splitlines() if ln and not ln.startswith(">")]
        return "".join(seq_lines) if seq_lines else None
    except Exception as exc:
        log.error("stitchr crashed (chain=%s, V=%s, J=%s, CDR3=%s): %s",
                  chain, v_gene, j_gene, cdr3, exc)
        return None
