#!/usr/bin/env python3
"""assemble_contigs.py — Assemble 50 nt nucleotide contigs around novel splice
junctions.

For each novel junction, a contig is constructed by joining:
  * 26 nt immediately **upstream**  of the junction (last 26 nt of the upstream exon)
  * 24 nt immediately **downstream** of the junction (first 24 nt of the downstream exon)

This yields a 50 nt contig that spans the novel splice site.

The script uses ``bedtools getfasta`` to retrieve sequences from the reference
genome.  Contigs with soft-clipped bases (lower-case nucleotides) are excluded.

Output: FASTA file where each sequence header encodes the junction coordinates.

Usage (standalone):
  python assemble_contigs.py \\
      --novel-junctions results/junctions/patient_001/novel_junctions.tsv \\
      --genome-fasta resources/GRCh38.primary_assembly.genome.fa \\
      --output results/contigs/patient_001/contigs.fa \\
      --upstream-nt 26 --downstream-nt 24

Usage (Snakemake):
  Called automatically by the ``assemble_contigs`` rule.
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _build_upstream_bed(
    junctions: pd.DataFrame,
    upstream_nt: int,
) -> pd.DataFrame:
    """Return a BED DataFrame for the upstream arm of each junction."""
    bed = pd.DataFrame(
        {
            "chrom": junctions["chrom"],
            "start": (junctions["start"] - upstream_nt).clip(lower=0),
            "end": junctions["start"],
            "name": junctions["junction_id"],
            "score": 0,
            "strand": junctions["strand"],
        }
    )
    return bed


def _build_downstream_bed(
    junctions: pd.DataFrame,
    downstream_nt: int,
) -> pd.DataFrame:
    """Return a BED DataFrame for the downstream arm of each junction."""
    bed = pd.DataFrame(
        {
            "chrom": junctions["chrom"],
            "start": junctions["end"],
            "end": junctions["end"] + downstream_nt,
            "name": junctions["junction_id"],
            "score": 0,
            "strand": junctions["strand"],
        }
    )
    return bed


def _write_bed(df: pd.DataFrame, path: str | Path) -> None:
    """Write a BED DataFrame to a file (no header)."""
    df.to_csv(path, sep="\t", header=False, index=False)


def _run_bedtools_getfasta(
    bed_path: str | Path,
    genome_fasta: str | Path,
    output_fasta: str | Path,
    strand: bool = True,
) -> None:
    """Run ``bedtools getfasta`` to retrieve sequences for BED intervals.

    Args:
        bed_path:     BED file of intervals to retrieve.
        genome_fasta: Indexed reference genome FASTA.
        output_fasta: Destination FASTA file.
        strand:       If True, reverse-complement minus-strand intervals.

    Raises:
        subprocess.CalledProcessError: if bedtools exits non-zero.
    """
    cmd = [
        "bedtools", "getfasta",
        "-fi", str(genome_fasta),
        "-bed", str(bed_path),
        "-fo", str(output_fasta),
        "-name",
    ]
    if strand:
        cmd.append("-s")
    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("bedtools getfasta stderr:\n%s", result.stderr)
        result.check_returncode()


def _parse_fasta(fasta_path: str | Path) -> dict[str, str]:
    """Parse a FASTA file and return a dict of {header: sequence}."""
    sequences: dict[str, str] = {}
    current_header: str | None = None
    seq_parts: list[str] = []
    with Path(fasta_path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(seq_parts)
                current_header = line[1:].split("::")[0].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_header is not None:
        sequences[current_header] = "".join(seq_parts)
    return sequences


def _has_soft_clip(sequence: str) -> bool:
    """Return True if the sequence contains any lower-case (soft-clipped) bases."""
    return any(c.islower() for c in sequence)


def assemble_contigs(
    novel_junctions_tsv: str | Path,
    genome_fasta: str | Path,
    output_fasta: str | Path,
    upstream_nt: int = 26,
    downstream_nt: int = 24,
) -> None:
    """Assemble 50 nt contigs for all novel junctions and write to FASTA.

    Args:
        novel_junctions_tsv: Path to the novel junctions TSV.
        genome_fasta:        Path to the GRCh38 reference FASTA (must be indexed).
        output_fasta:        Destination FASTA file.
        upstream_nt:         Nucleotides upstream of the junction to include.
        downstream_nt:       Nucleotides downstream of the junction to include.
    """
    novel_junctions_tsv = Path(novel_junctions_tsv)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(novel_junctions_tsv, sep="\t")
    if df.empty:
        log.warning("No novel junctions found in %s", novel_junctions_tsv)
        output_fasta.touch()
        return

    # Keep only tumor_exclusive junctions for neoepitope prediction.
    # normal_shared junctions (also present in normal) are excluded here;
    # they remain in the TSV for reporting but are not prediction candidates.
    if "junction_origin" in df.columns:
        n_total = len(df)
        df = df[df["junction_origin"] == "tumor_exclusive"].copy()
        log.info(
            "Filtered to tumor_exclusive junctions: %d / %d",
            len(df), n_total,
        )

    if df.empty:
        log.warning("No tumor_exclusive junctions remain after origin filter")
        output_fasta.touch()
        return

    # Deduplicate junctions (same coordinates may appear in multiple samples)
    junc_df = df.drop_duplicates(subset=["junction_id"]).copy()
    log.info("Assembling contigs for %d unique junctions", len(junc_df))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        upstream_bed = tmpdir / "upstream.bed"
        downstream_bed = tmpdir / "downstream.bed"
        upstream_fa = tmpdir / "upstream.fa"
        downstream_fa = tmpdir / "downstream.fa"

        _write_bed(_build_upstream_bed(junc_df, upstream_nt), upstream_bed)
        _write_bed(_build_downstream_bed(junc_df, downstream_nt), downstream_bed)

        _run_bedtools_getfasta(upstream_bed, genome_fasta, upstream_fa)
        _run_bedtools_getfasta(downstream_bed, genome_fasta, downstream_fa)

        upstream_seqs = _parse_fasta(upstream_fa)
        downstream_seqs = _parse_fasta(downstream_fa)

    n_written = 0
    n_skipped_softclip = 0
    n_skipped_length = 0

    with output_fasta.open("w") as out:
        for _, row in junc_df.iterrows():
            junc_id = row["junction_id"]
            up_seq = upstream_seqs.get(junc_id, "")
            dn_seq = downstream_seqs.get(junc_id, "")

            # Validate lengths
            if len(up_seq) < upstream_nt or len(dn_seq) < downstream_nt:
                log.debug(
                    "Skipping %s: insufficient sequence length (up=%d dn=%d)",
                    junc_id, len(up_seq), len(dn_seq),
                )
                n_skipped_length += 1
                continue

            # Take exactly the required number of nucleotides
            up_seq = up_seq[-upstream_nt:]
            dn_seq = dn_seq[:downstream_nt]
            contig = up_seq + dn_seq

            # Exclude contigs with soft-clipped regions
            if _has_soft_clip(contig):
                n_skipped_softclip += 1
                continue

            contig = contig.upper()
            header = (
                f">{junc_id}|{row['chrom']}:{row['start']}-{row['end']}"
                f":{row['strand']}|{row['sample_type']}"
            )
            out.write(f"{header}\n{contig}\n")
            n_written += 1

    log.info(
        "Contigs: %d written, %d skipped (soft-clip), %d skipped (length)",
        n_written, n_skipped_softclip, n_skipped_length,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    assemble_contigs(
        novel_junctions_tsv=snakemake.input.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        genome_fasta=snakemake.input.genome_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_fasta=snakemake.output.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        upstream_nt=snakemake.params.upstream_nt,  # type: ignore[name-defined]  # noqa: F821
        downstream_nt=snakemake.params.downstream_nt,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Assemble 50 nt contigs around novel splice junctions."
    )
    parser.add_argument(
        "--novel-junctions", required=True,
        help="Novel junctions TSV file",
    )
    parser.add_argument("--genome-fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--output", required=True, help="Output contigs FASTA")
    parser.add_argument("--upstream-nt", type=int, default=26)
    parser.add_argument("--downstream-nt", type=int, default=24)
    args = parser.parse_args()

    assemble_contigs(
        novel_junctions_tsv=args.novel_junctions,
        genome_fasta=args.genome_fasta,
        output_fasta=args.output,
        upstream_nt=args.upstream_nt,
        downstream_nt=args.downstream_nt,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
