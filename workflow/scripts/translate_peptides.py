#!/usr/bin/env python3
"""translate_peptides.py — Extract junction-spanning 9-mer peptides from 50 nt
splice-junction contigs.

For each contig, the junction breakpoint sits at nucleotide ``upstream_nt``
(default 26).  For each of the three reading frames (offsets 0, 1, 2), the
script identifies all 27 nt windows (= 9 codons) whose first codon is fully
upstream of the junction and whose last codon is fully downstream.

Spanning condition for a window starting at nucleotide ``start``:
    upstream_nt - 24  <=  start  <=  upstream_nt - 3

With defaults (upstream_nt = 26):  2 <= start <= 23.

Each valid window is translated directly to a 9-mer.  Windows containing a
stop codon or ambiguous amino acid (X) are discarded.

Output: TSV with columns  contig_key | start_nt | peptide

Usage (standalone):
  python translate_peptides.py \\
      --contigs-fasta results/contigs/patient_001/contigs.fa \\
      --output results/peptides/patient_001/peptides.tsv \\
      --upstream-nt 26

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

_PEPTIDE_LENGTH = 9
_CODON_SIZE = 3
_WINDOW_NT = _PEPTIDE_LENGTH * _CODON_SIZE  # 27 nt


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
# Junction-spanning 9-mer extraction
# ---------------------------------------------------------------------------

def extract_spanning_9mers(
    contig_seq: str,
    upstream_nt: int = 26,
    reading_frames: list[int] | None = None,
) -> list[tuple[int, str]]:
    """Extract junction-spanning 9-mers from a 50 nt contig.

    For each reading frame, finds all 27 nt windows satisfying the spanning
    condition: first codon fully upstream, last codon fully downstream.

    Args:
        contig_seq:     Nucleotide sequence (upper-case, typically 50 nt).
        upstream_nt:    Position of the junction breakpoint (number of upstream
                        nucleotides). Must match config[assembly][upstream_nt].
        reading_frames: 0-based frame offsets to consider. Defaults to [0, 1, 2].

    Returns:
        List of (start_nt, peptide) tuples. start_nt is the 0-based nucleotide
        start of the 27 nt window in the contig.
    """
    if reading_frames is None:
        reading_frames = [0, 1, 2]

    min_start = upstream_nt - (_PEPTIDE_LENGTH - 1) * _CODON_SIZE  # = 2 for defaults
    max_start = upstream_nt - _CODON_SIZE                           # = 23 for defaults

    results: list[tuple[int, str]] = []
    for frame in reading_frames:
        start = frame
        while start + _WINDOW_NT <= len(contig_seq):
            if min_start <= start <= max_start:
                window = contig_seq[start : start + _WINDOW_NT]
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
    upstream_nt: int = 26,
    reading_frames: list[int] | None = None,
) -> None:
    """Extract junction-spanning 9-mers from all contigs and write to TSV.

    Args:
        contigs_fasta:  Input FASTA of 50 nt contigs.
        output_tsv:     Output TSV (columns: contig_key, start_nt, peptide).
        upstream_nt:    Junction breakpoint position. Must match
                        config[assembly][upstream_nt].
        reading_frames: 0-based frame offsets. Defaults to [0, 1, 2].
    """
    if reading_frames is None:
        reading_frames = [0, 1, 2]

    contigs_fasta = Path(contigs_fasta)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    records = _parse_fasta(contigs_fasta)
    log.info(
        "Extracting junction-spanning 9-mers from %d contigs (upstream_nt=%d)",
        len(records), upstream_nt,
    )

    n_peptides = 0
    with output_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["contig_key", "start_nt", "peptide"])
        for header, seq in records:
            if not seq:
                continue
            for start_nt, peptide in extract_spanning_9mers(seq, upstream_nt, reading_frames):
                writer.writerow([header, start_nt, peptide])
                n_peptides += 1

    log.info(
        "Wrote %d junction-spanning 9-mers from %d contigs to %s",
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
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract junction-spanning 9-mers from 50 nt splice-junction contigs."
    )
    parser.add_argument("--contigs-fasta", required=True, help="Input contigs FASTA")
    parser.add_argument("--output", required=True, help="Output peptides TSV")
    parser.add_argument(
        "--upstream-nt", type=int, default=26,
        help="Junction breakpoint position (default: 26)",
    )
    args = parser.parse_args()

    translate_all(
        contigs_fasta=args.contigs_fasta,
        output_tsv=args.output,
        upstream_nt=args.upstream_nt,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
