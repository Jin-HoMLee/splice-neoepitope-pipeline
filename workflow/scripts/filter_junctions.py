#!/usr/bin/env python3
"""filter_junctions.py — Classify splice junctions by origin.

Junction origin hierarchy:
  all junctions
    └─ annotated         (in GENCODE)              → discarded
    └─ unannotated       (not in GENCODE)
         ├─ normal_shared          (also in matched normal)        → kept with label, excluded downstream
         ├─ gtex_pantissue_shared  (in GTEx pan-tissue blacklist)  → kept with label, excluded downstream
         └─ tumor_exclusive        (absent in normal AND GTEx)     → neoepitope prediction candidates

When no normal sample is present, all unannotated tumor junctions are labeled
`tumor_exclusive` (subject to the GTEx filter) with a warning.

The GTEx pan-tissue blacklist (Issue #211/#212) is an optional population-normal
filter: a junction seen in GTEx normal tissue is not tumor-specific even if it is
absent from the matched normal. Membership is checked AFTER the matched-normal
step, so the gtex_pantissue_shared count is the filter's *marginal* contribution
beyond the matched normal (no double-counting). When ``gtex_bed`` is None the
filter is a no-op and no junction is labeled gtex_pantissue_shared.

Output TSV columns:
  junction_id  chrom  start  end  strand  mapped_reads  sample_id  sample_type
  junction_origin  normal_source  reading_frame

normal_source (Issue #940) records which normal sample_type(s) removed a
`normal_shared` junction (e.g. ``Solid Tissue Normal``, ``Blood Derived Normal``,
or both, comma-joined), so the exclusion's provenance is explicit rather than
collapsed into one bucket. It is empty for `gtex_pantissue_shared` and
`tumor_exclusive` rows. All available normals still contribute to subtraction
(a blood-present junction is a self-antigen, a valid conservative filter); the
column only makes *which* normal did the removing visible.

reading_frame (0, 1, or 2) is the canonical CDS reading frame at the splice donor,
derived from the GENCODE GTF when gencode_gtf is supplied. Empty string means the frame
could not be determined (novel donor, no protein-coding CDS match, or ambiguous). This
annotation is informational only; all three frames are still translated downstream.

Usage (standalone):
  python filter_junctions.py \\
      --junction-files results/{patient_id}/alignment/{sample_id}/raw_junctions.tsv ... \\
      --manifest results/{patient_id}/alignment/manifest.tsv \\
      --reference-junctions references/reference_junctions.bed \\
      --output results/{patient_id}/junctions/novel_junctions.tsv \\
      [--gencode-gtf references/gencode.v47.annotation.gtf.gz]

Usage (Snakemake):
  Called automatically by the ``filter_junctions`` rule.
"""

import argparse
import csv
import gzip
import logging
import re
import statistics
from collections import defaultdict
from pathlib import Path

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


def _read_junction_file(file_path: str | Path) -> list[tuple[str, float, int | None]]:
    """Read a raw junction quantification TSV.

    Returns a list of ``(junction_id, read_count, annotated)`` tuples.

    The optional 3rd column is STAR's col-6 annotated flag (0=novel, 1=annotated),
    surfaced by ``star_sj_to_junctions.py`` (Issue #375). The HISAT2 path
    (``bed12_to_junctions.py``) emits only 2 columns, so ``annotated`` is ``None``
    for legacy 2-column rows — back-compat is preserved. A 3rd column that is not
    a valid integer is treated as absent (``None``).
    """
    rows: list[tuple[str, float, int | None]] = []
    with Path(file_path).open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                reads = float(parts[1])
            except ValueError:
                continue
            annotated: int | None = None
            if len(parts) >= 3:
                try:
                    annotated = int(parts[2])
                except ValueError:
                    annotated = None
            rows.append((parts[0], reads, annotated))
    return rows


# ---------------------------------------------------------------------------
# CDS reading frame lookup
# ---------------------------------------------------------------------------

def _open_gtf(path: str | Path):
    """Return a file handle for a plain or gzip-compressed GTF."""
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _parse_gtf_attribute(attributes: str, key: str) -> str | None:
    """Extract the value of a GTF attribute field by key."""
    match = re.search(rf'{key}\s+"([^"]+)"', attributes)
    return match.group(1) if match else None


def _build_cds_donor_lookup(
    gtf_path: str | Path,
) -> dict[tuple[str, int, str], set[int]]:
    """Parse protein-coding CDS records from a GENCODE GTF.

    For each CDS exon end that acts as a splice donor, computes the canonical
    reading frame offset within a junction contig (where upstream_nt is always
    divisible by 3):

        phase_at_donor = (exon_length - gtf_frame) % 3
        frame_offset   = (-phase_at_donor) % 3   →  0, 1, or 2

    Donor site coordinates (0-based):
        + strand: junction.start == cds_exon_end0
        − strand: junction.end   == cds_exon_start0

    Returns:
        dict mapping (chrom, donor_coord, strand) → set of frame offsets across
        all protein-coding transcripts at that donor.  Only entries with at least
        one frame offset are included.
    """
    donor_frames: dict[tuple[str, int, str], set[int]] = defaultdict(set)

    log.info("Parsing protein-coding CDS records from %s", gtf_path)
    n_records = 0
    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = fields[8]
            if _parse_gtf_attribute(attrs, "transcript_type") != "protein_coding":
                continue
            chrom = fields[0]
            start0 = int(fields[3]) - 1   # GTF 1-based → 0-based
            end0 = int(fields[4])          # GTF end is 1-based inclusive = 0-based exclusive (no adjustment needed)
            strand = fields[6]
            try:
                frame = int(fields[7])
            except ValueError:
                continue
            exon_len = end0 - start0
            phase = (exon_len - frame) % 3
            frame_offset = (-phase) % 3
            donor_coord = end0 if strand == "+" else start0
            donor_frames[(chrom, donor_coord, strand)].add(frame_offset)
            n_records += 1

    log.info("Parsed %d protein-coding CDS records, built %d unique donor sites", n_records, len(donor_frames))
    return dict(donor_frames)


# ---------------------------------------------------------------------------
# Core classification logic
# ---------------------------------------------------------------------------

def _build_normal_junction_sources(
    normal_files: list[tuple[Path, str, str]],
    reference_junctions: frozenset[tuple[str, int, int, str]],
    min_reads: int,
) -> dict[tuple[str, int, int, str], set[str]]:
    """Map each unannotated junction reliably seen in a normal sample to the set
    of normal ``sample_type``\\ s that contain it.

    A junction is included if it has >= min_reads in at least one normal sample
    and is not in the reference annotation. The value set records *which* normal
    types (e.g. ``Solid Tissue Normal``, ``Blood Derived Normal``) it was seen in,
    so the downstream ``normal_shared`` label carries its provenance rather than
    collapsing every normal into one undifferentiated bucket (Issue #940).

    Both normal types are used for subtraction on purpose: a junction present in
    the patient's blood is a self-antigen, not tumor-specific, so subtracting it is
    a valid conservative filter. The labeling stays honest about which normal did
    the removing; it does not drop blood.

    Args:
        normal_files:         List of (file_path, sample_id, sample_type) for normals.
        reference_junctions:  Frozenset of annotated junction tuples to exclude.
        min_reads:            Minimum read count to trust a normal junction.

    Returns:
        dict mapping ``(chrom, start, end, strand)`` to a set of ``sample_type``
        strings. Membership (``junction in result``) preserves the old set semantics.
    """
    sources: dict[tuple[str, int, int, str], set[str]] = defaultdict(set)
    for path, sample_id, sample_type in normal_files:
        rows = _read_junction_file(path)
        n_added = 0
        for junc_id, reads, _annotated in rows:  # annotated flag unused here
            if reads < min_reads:
                continue
            parsed = _parse_junction_id(junc_id)
            if parsed is None or parsed in reference_junctions:
                continue
            sources[parsed].add(sample_type)
            n_added += 1
        log.info(
            "Normal sample %s (%s): %d unannotated junctions (min_reads=%d)",
            sample_id, sample_type, n_added, min_reads,
        )
    log.info("Normal junction set: %d unique unannotated junctions across %d normal type(s)",
             len(sources), len({st for s in sources.values() for st in s}))
    return dict(sources)


def classify_junctions(
    junction_files: list[str | Path],
    manifest_path: str | Path,
    reference_bed: str | Path,
    output_path: str | Path,
    stats_output_path: str | Path | None = None,
    min_normal_reads: int = 2,
    gencode_gtf: str | Path | None = None,
    gtex_bed: str | Path | None = None,
) -> None:
    """Classify all junction quantification files by origin.

    Tumor junctions that survive the reference filter are labeled:
      - ``normal_shared``          — also present in the matched normal
      - ``gtex_pantissue_shared``  — absent in matched normal but present in the
                                     GTEx pan-tissue population blacklist (Issue #212)
      - ``tumor_exclusive``        — absent in both the matched normal AND GTEx

    Normal samples are used only to build the exclusion set and do not
    appear in the output TSV.

    When ``gencode_gtf`` is supplied, a ``reading_frame`` column (0/1/2 or
    empty string) is added to the output TSV recording the canonical CDS
    reading frame at the splice donor of each junction.  See
    ``_build_cds_donor_lookup`` for details.

    When ``stats_output_path`` is supplied, also writes a long-format
    per-tumor-sample stats TSV (Issue #214) with columns:
    ``sample_id, sample_type, category, count``. Categories form a funnel
    that reconciles arithmetically — ``junctions_raw`` (pre-filter raw
    regtools output) equals the sum of ``mean_reads_filtered`` (noise removed
    by the per-sample mean-reads filter) + ``annotated_discarded`` (in
    GENCODE) + ``normal_shared`` + ``gtex_pantissue_shared`` +
    ``tumor_exclusive`` (the three unannotated classes). The
    ``gtex_pantissue_shared`` row is always emitted (count 0 when no GTEx
    blacklist is supplied) so the funnel schema is stable. Normal samples are
    omitted.

    Issue #215 adds four descriptive (non-funnel) rows per tumor sample —
    ``min_reads``, ``mean_reads``, ``median_reads``, ``max_reads`` — that
    summarise the raw read-count distribution before the per-file mean
    threshold is applied. Use these to sanity-check whether the silent
    mean threshold was appropriate for the sample's depth.

    Args:
        junction_files:    List of raw junction quantification TSV paths.
        manifest_path:     Manifest TSV mapping file_id → sample_type.
        reference_bed:     Path to reference junction BED file.
        output_path:       Destination TSV for classified junctions.
        stats_output_path: Optional destination TSV for per-tumor-sample funnel
                           counts. When None, no stats file is written.
        min_normal_reads:  Minimum reads required to trust a normal junction.
        gencode_gtf:       Optional path to GENCODE GTF for reading frame
                           annotation.  When None, reading_frame is always "".
        gtex_bed:          Optional path to the GTEx pan-tissue population-normal
                           junction blacklist (BED6, Issue #211/#212). When supplied,
                           unannotated tumor junctions absent from the matched normal
                           but present in this set are labeled gtex_pantissue_shared
                           and excluded downstream. When None, the filter is a no-op.
    """
    import pandas as pd

    ref_junctions = _load_reference_junctions(reference_bed)
    gtex_junctions = _load_reference_junctions(gtex_bed) if gtex_bed else frozenset()
    if gtex_bed:
        log.info("GTEx pan-tissue blacklist active: %d junctions from %s",
                 len(gtex_junctions), gtex_bed)
    donor_frames = _build_cds_donor_lookup(gencode_gtf) if gencode_gtf else {}
    manifest = _load_manifest(manifest_path)

    # Split files into tumor and normal
    tumor_files: list[tuple[Path, str, str]] = []
    normal_files: list[tuple[Path, str, str]] = []

    for fp in junction_files:
        fp = Path(fp)
        sample_id = fp.parent.name
        sample_type = manifest.get(sample_id, "Unknown")
        if "normal" in sample_type.lower():
            normal_files.append((fp, sample_id, sample_type))
        else:
            tumor_files.append((fp, sample_id, sample_type))

    if not normal_files:
        log.warning(
            "No normal samples found — all unannotated tumor junctions will be "
            "labeled 'tumor_exclusive'. Add a 'Solid Tissue Normal' sample to the "
            "manifest for normal_shared filtering."
        )

    # Build the normal-junction provenance map: junction → set of normal sample_types
    # that contain it. Membership drives the normal_shared classification (unchanged);
    # the value set drives the normal_source provenance column (Issue #940).
    normal_junction_sources = _build_normal_junction_sources(
        normal_files, ref_junctions, min_reads=min_normal_reads
    )

    # Classify tumor junctions
    all_classified: list[dict] = []
    stats_rows: list[dict] = []

    for path, sample_id, sample_type in tumor_files:
        rows = _read_junction_file(path)
        if not rows:
            log.warning("No data rows found in %s", path)
            continue

        n_raw = len(rows)

        # Capture the raw read-count distribution before the per-file mean
        # filter is applied — supports the "was the threshold appropriate?"
        # diagnostic in the filtering audit trail (Issue #215).
        raw_reads = [r for _, r, _ in rows]
        dist_min = float(min(raw_reads))
        dist_max = float(max(raw_reads))
        dist_median = float(statistics.median(raw_reads))

        # Keep only junctions with read count above the per-file mean,
        # reducing noise from low-evidence junctions.
        mean_reads = sum(raw_reads) / len(raw_reads)
        rows = [(jid, r, ann) for jid, r, ann in rows if r > mean_reads]
        n_mean_filtered = n_raw - len(rows)

        n_annotated = n_normal_shared = n_gtex_shared = n_tumor_exclusive = 0
        # Per-normal-source breakdown of normal_shared (Issue #940). A junction seen
        # in >1 normal type counts once per type it matched, so these are descriptive
        # (they can sum to >= n_normal_shared) and are NOT part of the funnel partition.
        normal_source_counts: dict[str, int] = defaultdict(int)
        # STAR col-6 annotated-flag cross-check counters (Issue #375).
        n_star_annot_not_in_ref = 0
        n_ref_not_star_annot = 0

        for junc_id, reads, annotated in rows:
            parsed = _parse_junction_id(junc_id)
            if parsed is None:
                log.warning("Could not parse junction ID: %r", junc_id)
                continue

            chrom, start, end, strand = parsed

            in_ref = parsed in ref_junctions

            # STAR annotated-flag cross-check (Issue #375). Diagnostic ONLY —
            # GENCODE/ref_junctions membership stays authoritative for the
            # classification below; this never alters the decision. When STAR's
            # col-6 flag is present (3-column STAR input) and disagrees with our
            # GENCODE-BED membership, warn: a mismatch flags a chromosome-naming
            # bug (cf. Issue #148) or an index built from a different GTF.
            if annotated is not None:
                star_annotated = annotated == 1
                if star_annotated and not in_ref:
                    n_star_annot_not_in_ref += 1
                    log.warning(
                        "Annotated-flag mismatch for %s: STAR=annotated but absent "
                        "from reference BED (possible chr-naming mismatch or stale "
                        "reference)",
                        junc_id,
                    )
                elif not star_annotated and in_ref:
                    n_ref_not_star_annot += 1
                    log.warning(
                        "Annotated-flag mismatch for %s: in reference BED but "
                        "STAR=novel (STAR index likely built from a different GTF)",
                        junc_id,
                    )

            # Discard annotated junctions
            if in_ref:
                n_annotated += 1
                continue

            # Classify unannotated junctions. Priority order: matched normal
            # (patient-specific, highest confidence) → GTEx pan-tissue population
            # blacklist → tumor-exclusive. Checking GTEx only after the normal step
            # makes gtex_pantissue_shared the filter's *marginal* contribution beyond
            # the matched normal (a junction already removed as normal_shared is not
            # re-counted), keeping the funnel a clean partition.
            normal_source = ""
            if parsed in normal_junction_sources:
                origin = "normal_shared"
                n_normal_shared += 1
                # Provenance: which normal sample_type(s) removed this junction
                # (Issue #940). Sorted + comma-joined so a junction seen in both a
                # solid and a blood normal reads e.g. "Blood Derived Normal,Solid Tissue Normal".
                matched_types = normal_junction_sources[parsed]
                normal_source = ",".join(sorted(matched_types))
                for st in matched_types:
                    normal_source_counts[st] += 1
            elif parsed in gtex_junctions:
                origin = "gtex_pantissue_shared"
                n_gtex_shared += 1
            else:
                origin = "tumor_exclusive"
                n_tumor_exclusive += 1

            donor_coord = start if strand == "+" else end
            frames = donor_frames.get((chrom, donor_coord, strand), set())
            reading_frame = ",".join(str(f) for f in sorted(frames)) if frames else ""

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
                    "normal_source": normal_source,
                    "reading_frame": reading_frame,
                }
            )

        log.info(
            "Tumor sample %s: %d annotated (discarded), %d normal_shared, "
            "%d gtex_pantissue_shared, %d tumor_exclusive",
            sample_id, n_annotated, n_normal_shared, n_gtex_shared, n_tumor_exclusive,
        )

        # Summary of the STAR annotated-flag cross-check (Issue #375). Only
        # meaningful for 3-column STAR input; both counts are 0 for the HISAT2
        # 2-column path (annotated is None there, so no comparison is made).
        if n_star_annot_not_in_ref or n_ref_not_star_annot:
            log.warning(
                "Tumor sample %s: STAR annotated-flag cross-check found %d "
                "disagreements (STAR-annotated/not-in-BED=%d, in-BED/STAR-novel=%d)",
                sample_id, n_star_annot_not_in_ref + n_ref_not_star_annot,
                n_star_annot_not_in_ref, n_ref_not_star_annot,
            )

        for category, count in (
            ("junctions_raw", n_raw),
            ("mean_reads_filtered", n_mean_filtered),
            ("annotated_discarded", n_annotated),
            ("normal_shared", n_normal_shared),
            # Always emitted (0 when no GTEx blacklist) so the funnel schema is stable.
            ("gtex_pantissue_shared", n_gtex_shared),
            ("tumor_exclusive", n_tumor_exclusive),
            # Distribution summaries — descriptive, not part of the funnel sum.
            ("min_reads", dist_min),
            ("mean_reads", mean_reads),
            ("median_reads", dist_median),
            ("max_reads", dist_max),
        ):
            stats_rows.append({
                "sample_id": sample_id,
                "sample_type": sample_type,
                "category": category,
                "count": count,
            })

        # Per-normal-source breakdown of normal_shared (Issue #940), one descriptive
        # row per normal type that removed at least one junction. Category is
        # ``normal_shared:<normal sample_type>``. Descriptive, not a funnel partition
        # (a junction in >1 normal type is counted under each), so these are kept
        # separate from the reconciling funnel rows above.
        for normal_type, count in sorted(normal_source_counts.items()):
            stats_rows.append({
                "sample_id": sample_id,
                "sample_type": sample_type,
                "category": f"normal_shared:{normal_type}",
                "count": count,
            })

    df = pd.DataFrame(
        all_classified,
        columns=[
            "junction_id", "chrom", "start", "end", "strand",
            "mapped_reads", "sample_id", "sample_type", "junction_origin",
            "normal_source", "reading_frame",
        ],
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log.info(
        "Wrote %d classified junction records to %s "
        "(%d tumor_exclusive, %d normal_shared, %d gtex_pantissue_shared)",
        len(df),
        output_path,
        (df["junction_origin"] == "tumor_exclusive").sum() if not df.empty else 0,
        (df["junction_origin"] == "normal_shared").sum() if not df.empty else 0,
        (df["junction_origin"] == "gtex_pantissue_shared").sum() if not df.empty else 0,
    )

    if stats_output_path is not None:
        stats_df = pd.DataFrame(
            stats_rows,
            columns=["sample_id", "sample_type", "category", "count"],
        )
        stats_output_path = Path(stats_output_path)
        stats_output_path.parent.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(stats_output_path, sep="\t", index=False)
        log.info("Wrote junction filter stats to %s (%d rows)", stats_output_path, len(stats_df))


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    # Optional named input: the gtex_bed input function returns [] when the GTEx
    # filter is disabled, so .get() yields a falsy value. When enabled it returns a
    # single path, which Snakemake may surface as a bare string OR a 1-element list —
    # handle both robustly (indexing a bare string would slice the first character).
    gtex_in = snakemake.input.get("gtex_bed")  # type: ignore[name-defined]  # noqa: F821
    if not gtex_in:
        gtex_bed = None
    elif isinstance(gtex_in, str):
        gtex_bed = gtex_in
    else:
        gtex_bed = gtex_in[0]

    classify_junctions(
        junction_files=snakemake.input.junction_files,  # type: ignore[name-defined]  # noqa: F821
        manifest_path=snakemake.input.manifest,  # type: ignore[name-defined]  # noqa: F821
        reference_bed=snakemake.input.reference_junctions,  # type: ignore[name-defined]  # noqa: F821
        output_path=snakemake.output.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        stats_output_path=snakemake.output.stats,  # type: ignore[name-defined]  # noqa: F821
        min_normal_reads=snakemake.params.min_normal_reads,  # type: ignore[name-defined]  # noqa: F821
        gencode_gtf=snakemake.input.gencode_gtf,  # type: ignore[name-defined]  # noqa: F821
        gtex_bed=gtex_bed,
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify splice junctions by origin "
                    "(tumor_exclusive / normal_shared / gtex_pantissue_shared)."
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
        "--stats-output", default=None,
        help="Optional output TSV for per-tumor-sample funnel counts",
    )
    parser.add_argument(
        "--min-normal-reads", type=int, default=2,
        help="Minimum reads to trust a junction in the normal sample (default: 2)",
    )
    parser.add_argument(
        "--gencode-gtf", default=None,
        help="GENCODE GTF file for reading frame annotation (plain or .gz)",
    )
    parser.add_argument(
        "--gtex-bed", default=None,
        help="GTEx pan-tissue population-normal junction blacklist BED6 (Issue #211/#212). "
             "When supplied, unannotated junctions present in this set (and absent from the "
             "matched normal) are labeled gtex_pantissue_shared and excluded downstream.",
    )
    args = parser.parse_args()

    classify_junctions(
        junction_files=args.junction_files,
        manifest_path=args.manifest,
        reference_bed=args.reference_junctions,
        output_path=args.output,
        stats_output_path=args.stats_output,
        min_normal_reads=args.min_normal_reads,
        gencode_gtf=args.gencode_gtf,
        gtex_bed=args.gtex_bed,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
