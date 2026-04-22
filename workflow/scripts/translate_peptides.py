#!/usr/bin/env python3
"""translate_peptides.py — Extract junction-spanning peptides from splice-junction contigs.

For each contig, the junction breakpoint sits at nucleotide ``upstream_nt``.
For each configured peptide length L and each of the three reading frames
(offsets 0, 1, 2), the script identifies all L*3 nt windows whose first codon
is fully upstream of the junction and whose last codon is fully downstream.

Spanning condition for a window of length L starting at nucleotide ``start``:
    upstream_nt - (L-1)*3  <=  start  <=  upstream_nt - 3

Each valid window is translated directly to an L-mer.  Windows containing a
stop codon or ambiguous amino acid (X) are discarded.

Output: TSV with columns  contig_key | start_nt | peptide

Usage (standalone):
  python translate_peptides.py \\
      --contigs-fasta results/contigs/patient_001/contigs.fa \\
      --output results/peptides/patient_001/peptides.tsv \\
      --upstream-nt 27 --peptide-lengths 8 9 10

Usage (Snakemake):
  Called automatically by the ``translate_peptides`` rule.
"""

import argparse
import csv
import logging
from pathlib import Path

from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_CODON_SIZE = 3


# ---------------------------------------------------------------------------
# FASTA helper
# ---------------------------------------------------------------------------

def _parse_fasta(fasta_path: str | Path) -> list[tuple[str, str]]:
    """Parse a FASTA file and return a list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    current_header: str | None = None
    seq_parts: list[str] = []
    with Path(fasta_path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(seq_parts)))
                current_header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_header is not None:
        records.append((current_header, "".join(seq_parts)))
    return records


# ---------------------------------------------------------------------------
# Junction-spanning peptide extraction
# ---------------------------------------------------------------------------

def extract_spanning_peptides(
    contig_seq: str,
    upstream_nt: int = 27,
    peptide_lengths: list[int] | None = None,
    reading_frames: list[int] | None = None,
) -> list[tuple[int, str]]:
    """Extract junction-spanning peptides from a junction contig.

    For each peptide length and reading frame, finds all windows satisfying
    the spanning condition: first codon fully upstream, last codon fully
    downstream.

    Args:
        contig_seq:      Nucleotide sequence (upper-case).
        upstream_nt:     Position of the junction breakpoint. Must equal
                         3 * (max(peptide_lengths) - 1).
        peptide_lengths: Amino acid lengths to extract. Defaults to [8, 9, 10].
        reading_frames:  0-based frame offsets to consider. Defaults to [0, 1, 2].

    Returns:
        List of (start_nt, peptide) tuples.
    """
    if peptide_lengths is None:
        peptide_lengths = [8, 9, 10]
    if reading_frames is None:
        reading_frames = [0, 1, 2]

    expected_upstream = _CODON_SIZE * (max(peptide_lengths) - 1)
    if upstream_nt != expected_upstream:
        log.warning(
            "upstream_nt=%d does not match 3*(max_length-1)=%d; "
            "downstream coverage may be incomplete for the longest peptide length",
            upstream_nt, expected_upstream,
        )

    results: list[tuple[int, str]] = []
    max_start = upstream_nt - _CODON_SIZE
    for length in peptide_lengths:
        window_nt = length * _CODON_SIZE
        min_start = upstream_nt - (length - 1) * _CODON_SIZE
        for frame in reading_frames:
            start = frame
            while start + window_nt <= len(contig_seq):
                if min_start <= start <= max_start:
                    window = contig_seq[start : start + window_nt]
                    peptide = str(Seq(window).translate())
                    if "*" not in peptide and "X" not in peptide:
                        results.append((start, peptide))
                start += _CODON_SIZE

    return results


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def translate_all(
    contigs_fasta: str | Path,
    output_tsv: str | Path,
    upstream_nt: int = 27,
    peptide_lengths: list[int] | None = None,
    reading_frames: list[int] | None = None,
) -> None:
    """Extract junction-spanning peptides from all contigs and write to TSV.

    Args:
        contigs_fasta:   Input FASTA of junction contigs.
        output_tsv:      Output TSV (columns: contig_key, start_nt, peptide).
        upstream_nt:     Junction breakpoint position. Must equal
                         3 * (max(peptide_lengths) - 1).
        peptide_lengths: Amino acid lengths to extract. Defaults to [8, 9, 10].
        reading_frames:  0-based frame offsets. Defaults to [0, 1, 2].
    """
    if peptide_lengths is None:
        peptide_lengths = [8, 9, 10]
    if reading_frames is None:
        reading_frames = [0, 1, 2]

    contigs_fasta = Path(contigs_fasta)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    records = _parse_fasta(contigs_fasta)
    log.info(
        "Extracting junction-spanning %s-mers from %d contigs (upstream_nt=%d)",
        "/".join(str(plen) for plen in peptide_lengths), len(records), upstream_nt,
    )

    n_peptides = 0
    with output_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["contig_key", "start_nt", "peptide"])
        for header, seq in records:
            if not seq:
                continue
            for start_nt, peptide in extract_spanning_peptides(
                seq, upstream_nt, peptide_lengths, reading_frames
            ):
                writer.writerow([header, start_nt, peptide])
                n_peptides += 1

    log.info(
        "Wrote %d junction-spanning peptides from %d contigs to %s",
        n_peptides, len(records), output_tsv,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    translate_all(
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.peptides_tsv,  # type: ignore[name-defined]  # noqa: F821
        upstream_nt=snakemake.params.upstream_nt,  # type: ignore[name-defined]  # noqa: F821
        peptide_lengths=snakemake.params.peptide_lengths,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract junction-spanning peptides from splice-junction contigs."
    )
    parser.add_argument("--contigs-fasta", required=True, help="Input contigs FASTA")
    parser.add_argument("--output", required=True, help="Output peptides TSV")
    parser.add_argument(
        "--upstream-nt", type=int, default=27,
        help="Junction breakpoint position (default: 27)",
    )
    parser.add_argument(
        "--peptide-lengths", type=int, nargs="+", default=[8, 9, 10],
        help="Peptide lengths to extract (default: 8 9 10)",
    )
    args = parser.parse_args()

    translate_all(
        contigs_fasta=args.contigs_fasta,
        output_tsv=args.output,
        upstream_nt=args.upstream_nt,
        peptide_lengths=args.peptide_lengths,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
