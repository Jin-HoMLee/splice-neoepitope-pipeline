#!/usr/bin/env python3
"""proteome_filter.py — Filter junction-spanning peptides against the human reference proteome.

Parses the UniProt Swiss-Prot human FASTA once and builds a k-mer index:
  {k_mer: first_accession}  for each configured peptide length.

Each query peptide is then checked against the index in O(1).  Peptides
present in the canonical proteome are self-peptides the immune system is
tolerized to and cannot be neoantigens.

Output:
  peptides_novel.tsv    — same schema as peptides.tsv, exact matches removed
  peptides_excluded.tsv — excluded peptides with all matching UniProt accessions (semicolon-separated)

# Future extension — approximate matching (≤N mismatches)
# To also exclude near-self peptides (cross-reactivity / tolerisation):
# generate all single-substitution variants of each query and check each
# against the exact k-mer index.  For a 9-mer: 9 × 19 = 171 variants per
# query; O(query_count × L × 20) total — still fast.
# The key insight: expand neighbours on the query side (small: ~424K peptides
# × 171 variants = 72M lookups) rather than the proteome side (huge: ~8M
# proteome k-mers × 171 variants each ≈ 1.4B index entries).

Usage (Snakemake):
  Called automatically by the proteome_filter_peptides rule.

Usage (standalone):
  python proteome_filter.py \\
      --peptides-tsv results/.../peptides/peptides.tsv \\
      --proteome-fasta resources/human_proteome.fasta \\
      --novel-tsv results/.../peptides/peptides_novel.tsv \\
      --excluded-tsv results/.../peptides/peptides_excluded.tsv \\
      --peptide-lengths 8 9 10
"""

import argparse
import csv
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def _parse_accession(header: str) -> str:
    """Extract accession from a UniProt FASTA header.

    Swiss-Prot format: >sp|P12345|PROT_HUMAN ...
    Falls back to the first whitespace-delimited token for non-standard headers.
    """
    if header.startswith("sp|") or header.startswith("tr|"):
        return header.split("|")[1]
    return header.split()[0]


def _build_kmer_index(
    proteome_fasta: Path,
    peptide_lengths: list[int],
) -> dict[str, list[str]]:
    """Parse Swiss-Prot FASTA once; return {k_mer: [accession, ...]}.

    Slides a window of each configured peptide length over every protein
    sequence.  All accessions containing a given k-mer are recorded in order
    of first appearance.
    """
    index: dict[str, list[str]] = {}
    n_proteins = 0
    current_accession: str | None = None
    seq_parts: list[str] = []

    def _index_sequence(accession: str, seq: str) -> None:
        for length in peptide_lengths:
            for i in range(len(seq) - length + 1):
                kmer = seq[i : i + length]
                if kmer not in index:
                    index[kmer] = [accession]
                elif accession not in index[kmer]:
                    index[kmer].append(accession)

    with proteome_fasta.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_accession is not None:
                    _index_sequence(current_accession, "".join(seq_parts))
                    n_proteins += 1
                current_accession = _parse_accession(line[1:])
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_accession is not None:
            _index_sequence(current_accession, "".join(seq_parts))
            n_proteins += 1

    log.info(
        "Built k-mer index: %d unique k-mers from %d proteins (lengths %s)",
        len(index), n_proteins, "/".join(str(l) for l in peptide_lengths),
    )
    return index


def proteome_filter(
    peptides_tsv: Path,
    proteome_fasta: Path,
    novel_tsv: Path,
    excluded_tsv: Path,
    peptide_lengths: list[int] | None = None,
    stats_output_path: str | Path | None = None,
) -> None:
    """Filter translated peptides against the human proteome via k-mer lookup.

    When ``stats_output_path`` is supplied, also writes a 2-column funnel
    slice (``category, count``) for the proteome-filter step (Issue #215):
    categories ``novel`` and ``excluded``.
    """
    if peptide_lengths is None:
        peptide_lengths = [8, 9, 10]

    peptides_tsv = Path(peptides_tsv)
    proteome_fasta = Path(proteome_fasta)
    novel_tsv = Path(novel_tsv)
    excluded_tsv = Path(excluded_tsv)

    for p in (novel_tsv, excluded_tsv):
        p.parent.mkdir(parents=True, exist_ok=True)

    with peptides_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        peptides = list(reader)
        fieldnames = reader.fieldnames

    log.info("Loaded %d peptides from %s", len(peptides), peptides_tsv)

    kmer_index = _build_kmer_index(proteome_fasta, peptide_lengths)

    excluded_fieldnames = list(fieldnames) + ["matched_accessions"]

    n_novel = n_excluded = 0
    with novel_tsv.open("w", newline="") as novel_fh, \
         excluded_tsv.open("w", newline="") as excl_fh:

        novel_writer = csv.DictWriter(novel_fh, fieldnames=fieldnames, delimiter="\t")
        excl_writer = csv.DictWriter(excl_fh, fieldnames=excluded_fieldnames, delimiter="\t")
        novel_writer.writeheader()
        excl_writer.writeheader()

        for row in peptides:
            pep = row["peptide"]
            if pep in kmer_index:
                accessions = ";".join(kmer_index[pep])
                excl_writer.writerow({**row, "matched_accessions": accessions})
                n_excluded += 1
            else:
                novel_writer.writerow(row)
                n_novel += 1

    log.info(
        "Result: %d novel peptides passed, %d excluded (%.1f%% reduction)",
        n_novel, n_excluded,
        100 * n_excluded / len(peptides) if peptides else 0,
    )

    if stats_output_path is not None:
        stats_output_path = Path(stats_output_path)
        stats_output_path.parent.mkdir(parents=True, exist_ok=True)
        with stats_output_path.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["category", "count"])
            writer.writerow(["novel", n_novel])
            writer.writerow(["excluded", n_excluded])
        log.info("Proteome-filter stats written to %s", stats_output_path)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry points
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    sm = snakemake  # type: ignore[name-defined]  # noqa: F821
    log_file = sm.log[0]
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    proteome_filter(
        peptides_tsv=sm.input.peptides_tsv,
        proteome_fasta=sm.input.proteome_fasta,
        novel_tsv=sm.output.novel_tsv,
        excluded_tsv=sm.output.excluded_tsv,
        peptide_lengths=sm.params.peptide_lengths,
        stats_output_path=getattr(sm.output, "stats", None),
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter translated peptides against the human proteome via k-mer set lookup."
    )
    parser.add_argument("--peptides-tsv", required=True)
    parser.add_argument("--proteome-fasta", required=True)
    parser.add_argument("--novel-tsv", required=True)
    parser.add_argument("--excluded-tsv", required=True)
    parser.add_argument(
        "--peptide-lengths", type=int, nargs="+", default=[8, 9, 10],
    )
    parser.add_argument(
        "--stats-output", default=None,
        help="Optional proteome-filter stats TSV (Issue #215)",
    )
    args = parser.parse_args()

    proteome_filter(
        peptides_tsv=args.peptides_tsv,
        proteome_fasta=args.proteome_fasta,
        novel_tsv=args.novel_tsv,
        excluded_tsv=args.excluded_tsv,
        peptide_lengths=args.peptide_lengths,
        stats_output_path=args.stats_output,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
