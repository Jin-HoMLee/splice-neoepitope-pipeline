#!/usr/bin/env python3
"""filter_junctions.py — Filter raw splice-junction quantification files for
novel (non-reference) junctions.

Two sequential filtering steps are applied:

1. **Read-count filter**: within each individual junction quantification file,
   keep only junctions whose mapped read count exceeds the *mean* mapped read
   count for that file.  This removes low-read background junctions.

2. **Novelty filter**: remove any junction that appears in the reference
   junction BED file derived from GENCODE annotations.

Input junction quantification files are tab-separated with columns:
  junction_id  mapped_reads  [additional columns ignored]

where ``junction_id`` has the form ``chr:start:end:strand`` (STAR output).

Output TSV columns (tab-separated, with header):
  junction_id  chrom  start  end  strand  mapped_reads  sample_id  sample_type

Usage (standalone):
  python filter_junctions.py \\
      --junction-files results/raw_data/TCGA-BRCA/files/*.tsv \\
      --manifest results/raw_data/TCGA-BRCA/manifest.tsv \\
      --reference-junctions resources/reference_junctions.bed \\
      --output results/junctions/TCGA-BRCA/novel_junctions.tsv

Usage (Snakemake):
  Called automatically by the ``filter_junctions`` rule.
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_reference_junctions(bed_path: str | Path) -> frozenset[tuple[str, int, int, str]]:
    """Load reference junctions from a BED file.

    Returns a frozenset of ``(chrom, start, end, strand)`` tuples.
    """
    ref: set[tuple[str, int, int, str]] = set()
    with Path(bed_path).open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]
            ref.add((chrom, start, end, strand))
    log.info("Loaded %d reference junctions from %s", len(ref), bed_path)
    return frozenset(ref)


def _load_manifest(manifest_path: str | Path) -> dict[str, str]:
    """Return a mapping of file_id → sample_type from the manifest TSV."""
    mapping: dict[str, str] = {}
    with Path(manifest_path).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["file_id"]] = row.get("sample_type", "Unknown")
    return mapping


def _parse_junction_id(junction_id: str) -> tuple[str, int, int, str] | None:
    """Parse a junction identifier of the form ``chr:start:end:strand``.

    STAR junction quantification files use 1-based coordinates; we convert to
    0-based half-open (BED) for consistency with the reference.

    Returns ``(chrom, start, end, strand)`` or ``None`` on parse error.
    """
    parts = junction_id.split(":")
    if len(parts) < 4:
        # Try alternative format: chr:start-end:strand
        parts = junction_id.replace("-", ":").split(":")
    if len(parts) < 4:
        return None
    chrom = parts[0]
    try:
        start = int(parts[1]) - 1  # 1-based → 0-based
        end = int(parts[2])         # end is already exclusive in STAR output
        strand = parts[3] if parts[3] in ("+", "-") else "."
    except ValueError:
        return None
    return chrom, start, end, strand


def filter_single_file(
    file_path: str | Path,
    reference_junctions: frozenset[tuple[str, int, int, str]],
    sample_id: str,
    sample_type: str,
    strategy: str = "mean",
) -> list[dict]:
    """Filter one junction quantification file and return novel junctions.

    Args:
        file_path:            Path to the raw TSV junction quantification file.
        reference_junctions:  Frozenset of known junction tuples to exclude.
        sample_id:            File/sample identifier (used in output rows).
        sample_type:          e.g., "Primary Tumor" or "Solid Tissue Normal".
        strategy:             Read-count filter strategy; currently only "mean".

    Returns:
        List of dicts with keys: junction_id, chrom, start, end, strand,
        mapped_reads, sample_id, sample_type.
    """
    file_path = Path(file_path)
    rows = []
    with file_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            junc_id = parts[0]
            try:
                reads = float(parts[1])
            except ValueError:
                continue
            rows.append((junc_id, reads))

    if not rows:
        log.warning("No data rows found in %s", file_path)
        return []

    # Step 1: read-count filter
    if strategy == "mean":
        mean_reads = sum(r for _, r in rows) / len(rows)
        rows = [(jid, r) for jid, r in rows if r > mean_reads]
    else:
        raise ValueError(f"Unknown filtering strategy: {strategy!r}")

    log.debug(
        "File %s: %d junctions pass read-count filter (mean=%.2f)",
        file_path.name,
        len(rows),
        mean_reads if strategy == "mean" else 0.0,
    )

    # Step 2: novelty filter
    novel = []
    for junc_id, reads in rows:
        parsed = _parse_junction_id(junc_id)
        if parsed is None:
            log.warning("Could not parse junction ID: %r", junc_id)
            continue
        chrom, start, end, strand = parsed
        if (chrom, start, end, strand) in reference_junctions:
            continue
        novel.append(
            {
                "junction_id": junc_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "mapped_reads": int(reads),
                "sample_id": sample_id,
                "sample_type": sample_type,
            }
        )

    return novel


def filter_all_junctions(
    junction_files: list[str | Path],
    manifest_path: str | Path,
    reference_bed: str | Path,
    output_path: str | Path,
    strategy: str = "mean",
) -> None:
    """Filter all junction quantification files for a single cancer type.

    Args:
        junction_files:  List of raw junction quantification TSV paths.
        manifest_path:   Manifest TSV mapping file_id → sample metadata.
        reference_bed:   Path to reference junction BED file.
        output_path:     Destination TSV for novel junctions.
        strategy:        Read-count filtering strategy.
    """
    ref_junctions = _load_reference_junctions(reference_bed)
    manifest = _load_manifest(manifest_path)

    all_novel: list[dict] = []
    for fp in junction_files:
        fp = Path(fp)
        # The file is named <file_id>.tsv
        sample_id = fp.stem
        sample_type = manifest.get(sample_id, "Unknown")
        novel = filter_single_file(
            fp, ref_junctions, sample_id=sample_id, sample_type=sample_type,
            strategy=strategy,
        )
        all_novel.extend(novel)
        log.info(
            "File %s (%s): %d novel junctions", sample_id, sample_type, len(novel)
        )

    df = pd.DataFrame(
        all_novel,
        columns=["junction_id", "chrom", "start", "end", "strand",
                 "mapped_reads", "sample_id", "sample_type"],
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log.info(
        "Wrote %d novel junction records to %s", len(df), output_path
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    filter_all_junctions(
        junction_files=snakemake.input.junction_files,  # type: ignore[name-defined]  # noqa: F821
        manifest_path=snakemake.input.manifest,  # type: ignore[name-defined]  # noqa: F821
        reference_bed=snakemake.input.reference_junctions,  # type: ignore[name-defined]  # noqa: F821
        output_path=snakemake.output.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        strategy=snakemake.params.strategy,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter raw TCGA junction files for novel splice junctions."
    )
    parser.add_argument(
        "--junction-files", nargs="+", required=True,
        help="Raw junction quantification TSV files",
    )
    parser.add_argument("--manifest", required=True, help="Manifest TSV")
    parser.add_argument(
        "--reference-junctions", required=True,
        help="Reference junction BED file",
    )
    parser.add_argument("--output", required=True, help="Output novel junctions TSV")
    parser.add_argument(
        "--strategy", default="mean",
        help="Read-count filter strategy (default: mean)",
    )
    args = parser.parse_args()

    filter_all_junctions(
        junction_files=args.junction_files,
        manifest_path=args.manifest,
        reference_bed=args.reference_junctions,
        output_path=args.output,
        strategy=args.strategy,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
