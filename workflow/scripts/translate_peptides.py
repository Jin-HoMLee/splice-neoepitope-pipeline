#!/usr/bin/env python3
"""translate_peptides.py — In-silico translation of 50 nt contigs into 16-mer
peptides using modern Biopython (≥1.78).

**Important**: This implementation deliberately avoids ``Bio.Alphabet``, which
was removed in Biopython 1.78.  All translation is performed via
``Bio.Seq.Seq.translate()`` without any alphabet argument.

For each 50 nt contig, three reading frames are used:
  * Frame 1: translation begins at nucleotide position 1 (0-indexed: 0)
  * Frame 2: translation begins at nucleotide position 2 (0-indexed: 1)
  * Frame 3: translation begins at nucleotide position 3 (0-indexed: 2)

Each frame yields at most ``floor((50 - frame_offset) / 3)`` codons, which is
16 codons (amino acids) for frame offsets 0 and 1, and 16 for frame offset 2
as well (48 nt / 3 = 16).

Peptide post-processing:
  * Truncate at the first stop codon (``*``) if present.
  * Record the position of the first methionine (M) if present (but do not
    truncate — the original paper notes start/stop truncation for M only for
    display; we include the full peptide up to the stop).
  * Peptides shorter than 8 aa after truncation are discarded.

Output: FASTA file of 16-mer (or truncated) peptides.

Usage (standalone):
  python translate_peptides.py \\
      --contigs-fasta results/contigs/TCGA-BRCA/contigs.fa \\
      --output results/peptides/TCGA-BRCA/peptides.fa \\
      --reading-frames 0 1 2 --peptide-length 16

Usage (Snakemake):
  Called automatically by the ``translate_peptides`` rule.
"""

import argparse
import logging
import sys
from pathlib import Path

from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Minimum peptide length to keep after truncation
_MIN_PEPTIDE_LENGTH = 8


# ---------------------------------------------------------------------------
# FASTA helpers
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
# Translation logic
# ---------------------------------------------------------------------------

def translate_contig(
    contig_seq: str,
    reading_frames: list[int],
    peptide_length: int,
) -> list[tuple[int, str]]:
    """Translate a nucleotide contig into peptides in the specified reading frames.

    Uses the modern ``Bio.Seq`` API (no ``Bio.Alphabet``).

    Args:
        contig_seq:    Nucleotide sequence string (upper-case).
        reading_frames: List of 0-based frame offsets (e.g. [0, 1, 2]).
        peptide_length: Target amino acid length (default 16).

    Returns:
        List of (frame_offset, peptide_string) tuples.  Only non-empty
        peptides that pass the minimum length filter are returned.
    """
    results: list[tuple[int, str]] = []
    for frame in reading_frames:
        # Slice the contig to the correct reading frame
        sub = contig_seq[frame:]
        # Trim to a multiple of 3
        trimmed = sub[: len(sub) - (len(sub) % 3)]
        if len(trimmed) < 3:
            continue

        # Translate using Bio.Seq (no alphabet argument — modern API)
        aa_seq = str(Seq(trimmed).translate())

        # Truncate at first stop codon
        stop_pos = aa_seq.find("*")
        if stop_pos == 0:
            continue  # stop codon immediately — no peptide
        if stop_pos > 0:
            aa_seq = aa_seq[:stop_pos]

        # Trim to target peptide length
        peptide = aa_seq[:peptide_length]

        if len(peptide) < _MIN_PEPTIDE_LENGTH:
            continue

        results.append((frame, peptide))
    return results


def translate_all(
    contigs_fasta: str | Path,
    output_fasta: str | Path,
    reading_frames: list[int] | None = None,
    peptide_length: int = 16,
) -> None:
    """Translate all contigs in a FASTA file and write peptides to output FASTA.

    Args:
        contigs_fasta:  Input FASTA of 50 nt contigs.
        output_fasta:   Output FASTA of translated peptides.
        reading_frames: List of 0-based frame offsets.  Defaults to [0, 1, 2].
        peptide_length: Target peptide length in amino acids.
    """
    if reading_frames is None:
        reading_frames = [0, 1, 2]

    contigs_fasta = Path(contigs_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    records = _parse_fasta(contigs_fasta)
    log.info("Translating %d contigs in frames %s", len(records), reading_frames)

    n_peptides = 0
    with output_fasta.open("w") as out:
        for header, seq in records:
            if not seq:
                continue
            peptides = translate_contig(seq, reading_frames, peptide_length)
            for frame, peptide in peptides:
                pep_header = f">{header}|frame{frame + 1}"
                out.write(f"{pep_header}\n{peptide}\n")
                n_peptides += 1

    log.info(
        "Wrote %d peptide sequences from %d contigs to %s",
        n_peptides, len(records), output_fasta,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    translate_all(
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_fasta=snakemake.output.peptides_fasta,  # type: ignore[name-defined]  # noqa: F821
        reading_frames=snakemake.params.reading_frames,  # type: ignore[name-defined]  # noqa: F821
        peptide_length=snakemake.params.peptide_length,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Translate 50 nt contigs to 16-mer peptides (3 reading frames)."
    )
    parser.add_argument("--contigs-fasta", required=True, help="Input contigs FASTA")
    parser.add_argument("--output", required=True, help="Output peptides FASTA")
    parser.add_argument(
        "--reading-frames", type=int, nargs="+", default=[0, 1, 2],
        help="0-based reading frame offsets (default: 0 1 2)",
    )
    parser.add_argument(
        "--peptide-length", type=int, default=16,
        help="Target peptide length in amino acids (default: 16)",
    )
    args = parser.parse_args()

    translate_all(
        contigs_fasta=args.contigs_fasta,
        output_fasta=args.output,
        reading_frames=args.reading_frames,
        peptide_length=args.peptide_length,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
