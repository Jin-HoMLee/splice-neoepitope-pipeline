#!/usr/bin/env python3
"""build_reference_junctions.py — Build a reference splice-junction list from a
GENCODE GTF annotation file for GRCh38.

For each annotated transcript, consecutive exon boundaries define a splice
junction.  The junction is recorded as the half-open interval
``(exon_end, next_exon_start)`` on the chromosome (0-based, BED-style).

Output BED columns (tab-separated, no header):
  chrom  start  end  name  score  strand

where ``name`` is ``<chrom>:<start>-<end>:<strand>``.

Usage (standalone):
  python build_reference_junctions.py \\
      --gtf resources/gencode.v47.annotation.gtf.gz \\
      --output resources/reference_junctions.bed

Usage (Snakemake):
  Called automatically by the ``build_reference_junctions`` rule.
"""

import argparse
import gzip
import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GTF parsing helpers
# ---------------------------------------------------------------------------

def _open_gtf(path: str | Path):
    """Return a file handle for a plain or gzip-compressed GTF."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _parse_attribute(attributes: str, key: str) -> str | None:
    """Extract the value of a GTF attribute field by key."""
    match = re.search(rf'{key}\s+"([^"]+)"', attributes)
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def extract_junctions(gtf_path: str | Path) -> set[tuple[str, int, int, str]]:
    """Parse a GENCODE GTF and return a set of splice junctions.

    Each junction is represented as ``(chrom, start, end, strand)`` where
    ``start`` and ``end`` follow 0-based half-open BED convention:
    ``start`` = last nucleotide position of the upstream exon (0-based),
    ``end``   = first nucleotide position of the downstream exon (0-based).

    Args:
        gtf_path: Path to a GENCODE GTF file (plain or gzip-compressed).

    Returns:
        Set of ``(chrom, start, end, strand)`` tuples.
    """
    # transcript_id → list of (start0, end0) exon intervals (0-based, half-open)
    transcript_exons: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)

    log.info("Parsing GTF: %s", gtf_path)
    n_exons = 0
    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature != "exon":
                continue

            chrom = fields[0]
            # GTF is 1-based inclusive; convert to 0-based half-open
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])       # already exclusive after -1 + 1
            strand = fields[6]
            attributes = fields[8]
            transcript_id = _parse_attribute(attributes, "transcript_id")
            if not transcript_id:
                continue

            transcript_exons[transcript_id].append((start0, end0, chrom, strand))
            n_exons += 1

    log.info("Parsed %d exon records across %d transcripts", n_exons, len(transcript_exons))

    junctions: set[tuple[str, int, int, str]] = set()
    for transcript_id, exons in transcript_exons.items():
        # Validate consistency
        chroms = {e[2] for e in exons}
        strands = {e[3] for e in exons}
        if len(chroms) > 1 or len(strands) > 1:
            continue  # malformed transcript

        chrom = exons[0][2]
        strand = exons[0][3]
        # Sort exons by start position
        sorted_exons = sorted(exons, key=lambda e: e[0])

        for i in range(len(sorted_exons) - 1):
            upstream_end = sorted_exons[i][1]    # end of exon i (exclusive)
            downstream_start = sorted_exons[i + 1][0]  # start of exon i+1
            if downstream_start <= upstream_end:
                continue  # overlapping exons — skip
            junctions.add((chrom, upstream_end, downstream_start, strand))

    log.info("Extracted %d unique splice junctions", len(junctions))
    return junctions


def write_bed(
    junctions: set[tuple[str, int, int, str]],
    output_path: str | Path,
) -> None:
    """Write junctions to a sorted BED file.

    Args:
        junctions:   Set of ``(chrom, start, end, strand)`` tuples.
        output_path: Destination BED file path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_junctions = sorted(junctions, key=lambda j: (j[0], j[1], j[2]))
    with output_path.open("w") as fh:
        for chrom, start, end, strand in sorted_junctions:
            name = f"{chrom}:{start}-{end}:{strand}"
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    log.info("Wrote %d junctions to %s", len(junctions), output_path)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    gtf_path = snakemake.input.gtf  # type: ignore[name-defined]  # noqa: F821
    output_bed = snakemake.output.bed  # type: ignore[name-defined]  # noqa: F821

    junctions = extract_junctions(gtf_path)
    write_bed(junctions, output_bed)


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Build reference splice-junction BED from GENCODE GTF."
    )
    parser.add_argument("--gtf", required=True, help="GENCODE GTF file (plain or .gz)")
    parser.add_argument("--output", required=True, help="Output BED file path")
    args = parser.parse_args()

    junctions = extract_junctions(args.gtf)
    write_bed(junctions, args.output)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
