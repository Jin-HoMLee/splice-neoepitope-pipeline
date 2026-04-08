#!/usr/bin/env python3
"""run_mhcflurry.py — Wrapper to run MHCflurry 2.x and parse the output.

MHCflurry is an open-source MHC-I binding predictor that achieves
state-of-the-art performance.  Unlike NetMHCPan, it does not require
academic registration and can be installed via pip.

Reference:
  O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
  of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
  Cell Systems, 11(1), 42-48.e7.

The script:
  1. Reads the 16-mer peptide FASTA as input.
  2. Generates all 9-mer sub-peptides (sliding window).
  3. Runs MHCflurry affinity prediction for the specified HLA allele.
  4. Classifies each 9-mer as a strong binder (IC50 < 50 nM), weak binder
     (IC50 < 500 nM), or non-binder.

Output TSV columns:
  peptide_16mer  position  peptide_9mer  allele  ic50_nM  percentile_rank
  binder_class  source_header

Usage (standalone):
  python run_mhcflurry.py \\
      --peptides-fasta results/peptides/TCGA-BRCA/peptides.fa \\
      --output results/predictions/TCGA-BRCA/predictions.tsv \\
      --allele HLA-A*02:01 \\
      --peptide-length 9 \\
      --ic50-strong 50 \\
      --ic50-weak 500

Usage (Snakemake):
  Called automatically by the ``run_mhcflurry`` rule.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Peptide extraction
# ---------------------------------------------------------------------------

def _read_peptides_fasta(fasta_path: str | Path) -> Iterator[tuple[str, str]]:
    """Read 16-mer peptides from a FASTA file.

    Yields:
        Tuples of (header, sequence).
    """
    for record in SeqIO.parse(fasta_path, "fasta"):
        yield record.description, str(record.seq)


def _parse_frame_offset(header: str) -> int | None:
    """Extract the 0-based reading frame offset from a peptide FASTA header.

    Headers from translate_peptides.py end with ``|frame1``, ``|frame2``, or
    ``|frame3``, corresponding to 0-based offsets 0, 1, 2.

    Returns None if the frame token is absent (no junction filter applied).
    """
    for part in header.split("|"):
        if part.startswith("frame") and part[5:].isdigit():
            return int(part[5:]) - 1  # frame1 → 0, frame2 → 1, frame3 → 2
    return None


def _generate_9mers(
    peptide_16mer: str,
    window_size: int = 9,
    frame_offset: int | None = None,
    upstream_nt: int = 26,
) -> list[tuple[int, str]]:
    """Generate 9-mer sub-peptides from a 16-mer that span the splice junction.

    Only 9-mers containing at least one complete amino acid from each side of
    the junction are returned.  A 9-mer at 0-indexed position i with frame
    offset f starts at nucleotide f + i*3 in the 50 nt contig.  The junction
    falls between nucleotides upstream_nt-1 and upstream_nt (default 25/26).

    Spanning condition: 2 <= f + i*3 <= 23  (with defaults upstream_nt=26,
    window_size=9).

    Args:
        peptide_16mer: The 16-mer amino acid sequence.
        window_size:   Length of each sub-peptide (default 9).
        frame_offset:  0-based reading frame offset used to translate the contig.
                       If None, all 9-mers are returned (no junction filter).
        upstream_nt:   Number of upstream nucleotides in the contig (default 26).

    Returns:
        List of (position, 9-mer) tuples where position is 1-indexed.
    """
    min_start = upstream_nt - (window_size - 1) * 3      # = 2 for defaults
    max_start = upstream_nt - 3                           # = 23 for defaults

    nmers = []
    for i in range(len(peptide_16mer) - window_size + 1):
        nmer = peptide_16mer[i : i + window_size]
        # Skip peptides containing stop codons or invalid characters
        if "*" in nmer or "X" in nmer:
            continue
        # Keep only 9-mers that span the junction: must include at least one
        # complete codon from the upstream exon and one from the downstream exon.
        # Without this, purely exonic 9-mers match normal proteins (false positives).
        if frame_offset is not None:
            start_nt = frame_offset + i * 3
            if not (min_start <= start_nt <= max_start):
                continue
        nmers.append((i + 1, nmer))  # 1-indexed position
    return nmers


# ---------------------------------------------------------------------------
# MHCflurry prediction
# ---------------------------------------------------------------------------

def _run_mhcflurry_predictions(
    peptides: list[str],
    allele: str,
) -> pd.DataFrame:
    """Run MHCflurry predictions on a list of peptides.

    Args:
        peptides: List of peptide sequences.
        allele:   HLA allele in MHCflurry format (e.g., 'HLA-A*02:01').

    Returns:
        DataFrame with columns: peptide, allele, mhcflurry_affinity,
        mhcflurry_affinity_percentile.
    """
    # Import MHCflurry here to avoid import errors if not installed
    try:
        from mhcflurry import Class1PresentationPredictor, Class1AffinityPredictor
    except ImportError:
        log.error(
            "MHCflurry is not installed. Install it with: "
            "pip install mhcflurry && mhcflurry-downloads fetch"
        )
        raise

    # Normalise allele format (MHCflurry uses HLA-A*02:01 format)
    normalised_allele = allele.replace("HLA-", "").replace("*", "").replace(":", "")
    # Convert back to MHCflurry format: HLA-A*02:01
    if len(normalised_allele) >= 4:
        mhcflurry_allele = f"HLA-{normalised_allele[0]}*{normalised_allele[1:3]}:{normalised_allele[3:5]}"
    else:
        mhcflurry_allele = allele

    log.info("Using MHCflurry with allele: %s (normalised: %s)", allele, mhcflurry_allele)

    # Use affinity predictor for IC50 values
    try:
        predictor = Class1AffinityPredictor.load()
    except Exception as exc:
        log.error(
            "Failed to load MHCflurry models. Run: mhcflurry-downloads fetch\n%s", exc
        )
        raise

    # Check if allele is supported
    supported_alleles = predictor.supported_alleles
    if mhcflurry_allele not in supported_alleles:
        # Try alternative format
        alt_allele = allele.replace("-", "").replace("*", "").replace(":", "")
        log.warning(
            "Allele %s not directly supported, trying alternatives. "
            "Supported alleles: %d total",
            mhcflurry_allele,
            len(supported_alleles),
        )
        # Find closest match
        matching = [a for a in supported_alleles if alt_allele[:5] in a.replace("*", "").replace(":", "")]
        if matching:
            mhcflurry_allele = matching[0]
            log.info("Using matched allele: %s", mhcflurry_allele)
        else:
            log.warning("No matching allele found, using original: %s", mhcflurry_allele)

    # Run predictions in batches for efficiency
    log.info("Running MHCflurry predictions for %d peptides...", len(peptides))

    # predict_to_dataframe() returns a DataFrame with columns:
    # peptide, allele, prediction, prediction_low, prediction_high, prediction_percentile
    # (mhcflurry 2.2.x renamed mhcflurry_affinity → prediction)
    results = predictor.predict_to_dataframe(
        peptides=peptides,
        alleles=[mhcflurry_allele] * len(peptides),
    )

    return results


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_fasta: str | Path,
    output_tsv: str | Path,
    allele: str = "HLA-A*02:01",
    peptide_length: int = 9,
    ic50_strong: float = 50.0,
    ic50_weak: float = 500.0,
) -> None:
    """Run MHCflurry and write parsed results to a TSV file.

    Args:
        peptides_fasta: FASTA of 16-mer peptides.
        output_tsv:     Destination TSV file.
        allele:         HLA allele (MHCflurry format, e.g., HLA-A*02:01).
        peptide_length: Sliding-window peptide length (9 for 9-mers).
        ic50_strong:    Strong-binder IC50 threshold (nM).
        ic50_weak:      Weak-binder IC50 threshold (nM).
    """
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Read 16-mers and generate 9-mers
    records: list[dict] = []
    peptide_to_source: dict[str, list[dict]] = {}

    for header, seq_16mer in _read_peptides_fasta(peptides_fasta):
        frame_offset = _parse_frame_offset(header)
        for pos, nmer in _generate_9mers(seq_16mer, peptide_length, frame_offset):
            record = {
                "source_header": header,
                "peptide_16mer": seq_16mer,
                "position": pos,
                "peptide_9mer": nmer,
            }
            records.append(record)
            peptide_to_source.setdefault(nmer, []).append(record)

    if not records:
        log.warning("No valid peptides found in %s", peptides_fasta)
        empty_df = pd.DataFrame(
            columns=[
                "source_header", "peptide_16mer", "position", "peptide_9mer",
                "allele", "ic50_nM", "percentile_rank", "binder_class",
            ]
        )
        empty_df.to_csv(output_tsv, sep="\t", index=False)
        return

    # Get unique peptides for prediction
    unique_peptides = list({r["peptide_9mer"] for r in records})
    log.info(
        "Extracted %d 9-mers (%d unique) from %d 16-mers",
        len(records), len(unique_peptides),
        sum(1 for _ in _read_peptides_fasta(peptides_fasta)),
    )

    # Step 2: Run MHCflurry predictions
    predictions_df = _run_mhcflurry_predictions(unique_peptides, allele)

    # Step 3: Map predictions back to original records
    peptide_to_prediction = {}
    for _, row in predictions_df.iterrows():
        pep = row["peptide"]
        peptide_to_prediction[pep] = {
            "ic50_nM": row["prediction"],
            "percentile_rank": row.get("prediction_percentile", float("nan")),
        }

    # Step 4: Build output records
    output_records: list[dict] = []
    for rec in records:
        pep = rec["peptide_9mer"]
        pred = peptide_to_prediction.get(pep, {"ic50_nM": float("inf"), "percentile_rank": float("nan")})
        ic50 = pred["ic50_nM"]

        if ic50 <= ic50_strong:
            binder_class = "strong"
        elif ic50 <= ic50_weak:
            binder_class = "weak"
        else:
            binder_class = "non"

        output_records.append({
            "source_header": rec["source_header"],
            "peptide_16mer": rec["peptide_16mer"],
            "position": rec["position"],
            "peptide_9mer": pep,
            "allele": allele,
            "ic50_nM": ic50,
            "percentile_rank": pred["percentile_rank"],
            "binder_class": binder_class,
        })

    # Step 5: Write output
    df = pd.DataFrame(
        output_records,
        columns=[
            "source_header", "peptide_16mer", "position", "peptide_9mer",
            "allele", "ic50_nM", "percentile_rank", "binder_class",
        ],
    )
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Predictions: %d total, %d strong, %d weak binders → %s",
        len(df),
        (df["binder_class"] == "strong").sum(),
        (df["binder_class"] == "weak").sum(),
        output_tsv,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    run_prediction(
        peptides_fasta=snakemake.input.peptides_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        allele=snakemake.params.hla_allele,  # type: ignore[name-defined]  # noqa: F821
        peptide_length=snakemake.params.peptide_length,  # type: ignore[name-defined]  # noqa: F821
        ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
        ic50_weak=float(snakemake.params.ic50_weak),  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run MHCflurry epitope prediction on 16-mer peptides."
    )
    parser.add_argument("--peptides-fasta", required=True, help="Input peptides FASTA")
    parser.add_argument("--output", required=True, help="Output predictions TSV")
    parser.add_argument("--allele", default="HLA-A*02:01", help="HLA allele")
    parser.add_argument("--peptide-length", type=int, default=9)
    parser.add_argument("--ic50-strong", type=float, default=50.0)
    parser.add_argument("--ic50-weak", type=float, default=500.0)
    args = parser.parse_args()

    run_prediction(
        peptides_fasta=args.peptides_fasta,
        output_tsv=args.output,
        allele=args.allele,
        peptide_length=args.peptide_length,
        ic50_strong=args.ic50_strong,
        ic50_weak=args.ic50_weak,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
