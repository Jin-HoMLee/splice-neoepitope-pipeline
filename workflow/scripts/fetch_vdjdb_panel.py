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


# ---------------------------------------------------------------------------
# Build panel — orchestrate + write outputs
# ---------------------------------------------------------------------------

PANEL_COLUMNS = [
    "allele", "va_gene", "ja_gene", "cdr3a",
    "vb_gene", "jb_gene", "cdr3b",
    "alpha_seq", "beta_seq",
    "vdjdb_score", "vdjdb_donor_id",
]

QC_COLUMNS = ["allele", "n_exact_matches", "n_in_panel", "panel_status"]


def build_panel(
    vdjdb_full_tsv,
    alleles: list,
    output_panel,
    output_qc,
    min_score: int,
    panel_size: int,
) -> None:
    """Build a per-patient VDJdb panel and write panel.tsv + panel_qc.tsv.

    For each allele in `alleles`:
      - exact-match filter (4-digit normalized)
      - sort by score DESC + donor_id ASC
      - iterate top-down: call stitchr per row, accumulate successes,
        skip+log on stitch failures, stop at panel_size or exhaustion
      - record n_exact_matches + n_in_panel + panel_status in QC
    """
    import pandas as pd  # lazy import
    from pathlib import Path

    output_panel = Path(output_panel)
    output_qc = Path(output_qc)
    output_panel.parent.mkdir(parents=True, exist_ok=True)
    output_qc.parent.mkdir(parents=True, exist_ok=True)

    df = load_and_filter_vdjdb(vdjdb_full_tsv, min_score=min_score)

    panel_rows = []
    qc_rows = []
    for allele in alleles:
        candidates = select_top_n_for_allele(df, allele=allele, n=len(df))
        n_exact = len(candidates)
        n_in_panel = 0
        for _, row in candidates.iterrows():
            if n_in_panel >= panel_size:
                break
            alpha = stitch_chain(
                v_gene=row["v.alpha"], j_gene=row["j.alpha"],
                cdr3=row["cdr3.alpha"], chain="A",
            )
            if alpha is None:
                log.error("Skipping VDJdb row for allele %s (donor %s) — alpha stitch failed",
                          allele, row["meta.subject.id"])
                continue
            beta = stitch_chain(
                v_gene=row["v.beta"], j_gene=row["j.beta"],
                cdr3=row["cdr3.beta"], chain="B",
            )
            if beta is None:
                log.error("Skipping VDJdb row for allele %s (donor %s) — beta stitch failed",
                          allele, row["meta.subject.id"])
                continue
            panel_rows.append({
                "allele": allele,
                "va_gene": row["v.alpha"], "ja_gene": row["j.alpha"], "cdr3a": row["cdr3.alpha"],
                "vb_gene": row["v.beta"],  "jb_gene": row["j.beta"],  "cdr3b": row["cdr3.beta"],
                "alpha_seq": alpha, "beta_seq": beta,
                "vdjdb_score": int(row["vdjdb.score"]),
                "vdjdb_donor_id": row["meta.subject.id"],
            })
            n_in_panel += 1

        status = classify_panel_status(n_in_panel, target_size=panel_size)
        qc_rows.append({
            "allele": allele,
            "n_exact_matches": n_exact,
            "n_in_panel": n_in_panel,
            "panel_status": status,
        })
        if status == "empty":
            log.warning("Allele %s: zero VDJdb entries (panel_status=empty)", allele)
        elif status == "low_coverage":
            log.warning("Allele %s: only %d entries available (panel_status=low_coverage)",
                        allele, n_in_panel)

    pd.DataFrame(panel_rows, columns=PANEL_COLUMNS).to_csv(output_panel, sep="\t", index=False)
    pd.DataFrame(qc_rows, columns=QC_COLUMNS).to_csv(output_qc, sep="\t", index=False)

    n_ok = sum(1 for r in qc_rows if r["panel_status"] == "ok")
    n_low = sum(1 for r in qc_rows if r["panel_status"] == "low_coverage")
    n_empty = sum(1 for r in qc_rows if r["panel_status"] == "empty")
    log.info(
        "VDJdb panel built — %d alleles: %d ok, %d low_coverage, %d empty. Wrote %s and %s.",
        len(alleles), n_ok, n_low, n_empty, output_panel, output_qc,
    )


# ---------------------------------------------------------------------------
# Allele loading from alleles.tsv (HLA typing output)
# ---------------------------------------------------------------------------

def load_alleles_tsv(alleles_tsv) -> list:
    """Load unique 4-digit alleles from alleles.tsv produced by aggregate_hla_alleles.

    alleles.tsv schema (existing): rows per locus (A/B/C) with allele1, allele2 columns.
    Returns a deduplicated list of allele strings, all 4-digit.
    """
    import pandas as pd  # lazy import
    df = pd.read_csv(alleles_tsv, sep="\t", dtype=str)
    alleles: set = set()
    for col in ("allele1", "allele2"):
        if col not in df.columns:
            continue
        for a in df[col].dropna().unique():
            normalized = normalize_allele_to_4digit(a)
            if normalized:
                alleles.add(normalized)
    return sorted(alleles)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    """Entry point when called via Snakemake `script:` directive."""
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s — %(message)s")
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    alleles = load_alleles_tsv(snakemake.input.alleles_tsv)  # type: ignore[name-defined]  # noqa: F821
    build_panel(
        vdjdb_full_tsv=snakemake.input.vdjdb_tsv,  # type: ignore[name-defined]  # noqa: F821
        alleles=alleles,
        output_panel=snakemake.output.panel,  # type: ignore[name-defined]  # noqa: F821
        output_qc=snakemake.output.qc,  # type: ignore[name-defined]  # noqa: F821
        min_score=snakemake.params.min_score,  # type: ignore[name-defined]  # noqa: F821
        panel_size=snakemake.params.panel_size,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    """Entry point when called from the command line (for development/debug)."""
    import argparse

    parser = argparse.ArgumentParser(description="Build a VDJdb TCR panel for one patient.")
    parser.add_argument("--vdjdb-tsv", required=True, help="Path to vdjdb_full.txt")
    parser.add_argument("--alleles-tsv", required=True, help="Path to alleles.tsv")
    parser.add_argument("--output-panel", required=True, help="Output panel.tsv path")
    parser.add_argument("--output-qc", required=True, help="Output panel_qc.tsv path")
    parser.add_argument("--min-score", type=int, default=2)
    parser.add_argument("--panel-size", type=int, default=10)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s — %(message)s")
    alleles = load_alleles_tsv(args.alleles_tsv)
    build_panel(
        vdjdb_full_tsv=args.vdjdb_tsv, alleles=alleles,
        output_panel=args.output_panel, output_qc=args.output_qc,
        min_score=args.min_score, panel_size=args.panel_size,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
    except NameError:
        _cli_main()
    else:
        _snakemake_main()
