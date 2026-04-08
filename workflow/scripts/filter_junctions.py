#!/usr/bin/env python3
"""filter_junctions.py — Classify splice junctions by origin.

Junction origin hierarchy:
  all junctions
    └─ annotated         (in GENCODE)              → discarded
    └─ unannotated       (not in GENCODE)
         ├─ patient_specific  (also in normal)     → kept with label, excluded downstream
         └─ tumor_specific    (absent in normal)   → neoepitope prediction candidates

When no normal sample is present, all unannotated tumor junctions are labeled
`tumor_specific` with a warning.

Output TSV columns:
  junction_id  chrom  start  end  strand  mapped_reads  sample_id  sample_type  junction_origin

Usage (standalone):
  python filter_junctions.py \\
      --junction-files results/raw_data/local/files/*.tsv \\
      --manifest results/raw_data/local/manifest.tsv \\
      --reference-junctions resources/reference_junctions.bed \\
      --output results/junctions/local/novel_junctions.tsv

Usage (Snakemake):
  Called automatically by the ``filter_junctions`` rule.
"""

import argparse
import csv
import logging
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

    STAR/regtools junction files use 1-based coordinates; we convert to
    0-based half-open (BED) for consistency with the reference.

    Returns ``(chrom, start, end, strand)`` or ``None`` on parse error.
    """
    parts = junction_id.split(":")
    if len(parts) < 4:
        parts = junction_id.replace("-", ":").split(":")
    if len(parts) < 4:
        return None
    chrom = parts[0]
    try:
        start = int(parts[1]) - 1  # 1-based → 0-based
        end = int(parts[2])
        strand = parts[3] if parts[3] in ("+", "-") else "."
    except ValueError:
        return None
    return chrom, start, end, strand


def _read_junction_file(file_path: str | Path) -> list[tuple[str, float]]:
    """Read a raw junction quantification TSV.

    Returns a list of ``(junction_id, read_count)`` tuples.
    """
    rows: list[tuple[str, float]] = []
    with Path(file_path).open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                rows.append((parts[0], float(parts[1])))
            except ValueError:
                continue
    return rows


# ---------------------------------------------------------------------------
# Core classification logic
# ---------------------------------------------------------------------------

def _build_normal_junction_set(
    normal_files: list[tuple[Path, str]],
    reference_junctions: frozenset[tuple[str, int, int, str]],
    min_reads: int,
) -> frozenset[tuple[str, int, int, str]]:
    """Build the set of unannotated junctions reliably seen in normal samples.

    A junction is included if it has >= min_reads in at least one normal sample
    and is not in the reference annotation.

    Args:
        normal_files:         List of (file_path, sample_id) for normal samples.
        reference_junctions:  Frozenset of annotated junction tuples to exclude.
        min_reads:            Minimum read count to trust a normal junction.

    Returns:
        Frozenset of ``(chrom, start, end, strand)`` tuples.
    """
    normal_set: set[tuple[str, int, int, str]] = set()
    for path, sample_id in normal_files:
        rows = _read_junction_file(path)
        n_added = 0
        for junc_id, reads in rows:
            if reads < min_reads:
                continue
            parsed = _parse_junction_id(junc_id)
            if parsed is None or parsed in reference_junctions:
                continue
            normal_set.add(parsed)
            n_added += 1
        log.info(
            "Normal sample %s: %d unannotated junctions (min_reads=%d)",
            sample_id, n_added, min_reads,
        )
    log.info("Normal junction set: %d unique unannotated junctions", len(normal_set))
    return frozenset(normal_set)


def classify_junctions(
    junction_files: list[str | Path],
    manifest_path: str | Path,
    reference_bed: str | Path,
    output_path: str | Path,
    min_normal_reads: int = 2,
) -> None:
    """Classify all junction quantification files by origin.

    Tumor junctions that survive the reference filter are labeled:
      - ``tumor_specific``  — absent in the matched normal
      - ``patient_specific`` — also present in the matched normal

    Normal samples are used only to build the exclusion set and do not
    appear in the output TSV.

    Args:
        junction_files:   List of raw junction quantification TSV paths.
        manifest_path:    Manifest TSV mapping file_id → sample_type.
        reference_bed:    Path to reference junction BED file.
        output_path:      Destination TSV for classified junctions.
        min_normal_reads: Minimum reads required to trust a normal junction.
    """
    ref_junctions = _load_reference_junctions(reference_bed)
    manifest = _load_manifest(manifest_path)

    # Split files into tumor and normal
    tumor_files: list[tuple[Path, str, str]] = []
    normal_files: list[tuple[Path, str]] = []

    for fp in junction_files:
        fp = Path(fp)
        sample_id = fp.stem
        sample_type = manifest.get(sample_id, "Unknown")
        if "normal" in sample_type.lower():
            normal_files.append((fp, sample_id))
        else:
            tumor_files.append((fp, sample_id, sample_type))

    if not normal_files:
        log.warning(
            "No normal samples found — all unannotated tumor junctions will be "
            "labeled 'tumor_specific'. Add a 'Solid Tissue Normal' sample to the "
            "manifest for patient_specific filtering."
        )

    # Build normal junction set
    normal_junction_set = _build_normal_junction_set(
        normal_files, ref_junctions, min_reads=min_normal_reads
    )

    # Classify tumor junctions
    all_classified: list[dict] = []

    for path, sample_id, sample_type in tumor_files:
        rows = _read_junction_file(path)
        if not rows:
            log.warning("No data rows found in %s", path)
            continue

        # Read-count filter: keep junctions above mean read count in this file
        mean_reads = sum(r for _, r in rows) / len(rows)
        rows = [(jid, r) for jid, r in rows if r > mean_reads]

        n_annotated = n_patient_specific = n_tumor_specific = 0

        for junc_id, reads in rows:
            parsed = _parse_junction_id(junc_id)
            if parsed is None:
                log.warning("Could not parse junction ID: %r", junc_id)
                continue

            chrom, start, end, strand = parsed

            # Discard annotated junctions
            if parsed in ref_junctions:
                n_annotated += 1
                continue

            # Classify unannotated junctions
            if parsed in normal_junction_set:
                origin = "patient_specific"
                n_patient_specific += 1
            else:
                origin = "tumor_specific"
                n_tumor_specific += 1

            all_classified.append(
                {
                    "junction_id": junc_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "mapped_reads": int(reads),
                    "sample_id": sample_id,
                    "sample_type": sample_type,
                    "junction_origin": origin,
                }
            )

        log.info(
            "Tumor sample %s: %d annotated (discarded), %d patient_specific, %d tumor_specific",
            sample_id, n_annotated, n_patient_specific, n_tumor_specific,
        )

    df = pd.DataFrame(
        all_classified,
        columns=[
            "junction_id", "chrom", "start", "end", "strand",
            "mapped_reads", "sample_id", "sample_type", "junction_origin",
        ],
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log.info(
        "Wrote %d classified junction records to %s "
        "(%d tumor_specific, %d patient_specific)",
        len(df),
        output_path,
        (df["junction_origin"] == "tumor_specific").sum() if not df.empty else 0,
        (df["junction_origin"] == "patient_specific").sum() if not df.empty else 0,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    classify_junctions(
        junction_files=snakemake.input.junction_files,  # type: ignore[name-defined]  # noqa: F821
        manifest_path=snakemake.input.manifest,  # type: ignore[name-defined]  # noqa: F821
        reference_bed=snakemake.input.reference_junctions,  # type: ignore[name-defined]  # noqa: F821
        output_path=snakemake.output.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        min_normal_reads=snakemake.params.min_normal_reads,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify splice junctions by origin (tumor_specific / patient_specific)."
    )
    parser.add_argument(
        "--junction-files", nargs="+", required=True,
        help="Raw junction quantification TSV files (tumor and normal)",
    )
    parser.add_argument("--manifest", required=True, help="Manifest TSV (file_id → sample_type)")
    parser.add_argument(
        "--reference-junctions", required=True,
        help="Reference junction BED file (GENCODE-derived)",
    )
    parser.add_argument("--output", required=True, help="Output classified junctions TSV")
    parser.add_argument(
        "--min-normal-reads", type=int, default=2,
        help="Minimum reads to trust a junction in the normal sample (default: 2)",
    )
    args = parser.parse_args()

    classify_junctions(
        junction_files=args.junction_files,
        manifest_path=args.manifest,
        reference_bed=args.reference_junctions,
        output_path=args.output,
        min_normal_reads=args.min_normal_reads,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
