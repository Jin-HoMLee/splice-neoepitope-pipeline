#!/usr/bin/env python3
"""star_sj_to_junctions.py — Convert STAR ``SJ.out.tab`` to raw_junctions.tsv.

STAR emits 9 tab-separated columns in ``SJ.out.tab``:

    1: chromosome
    2: intron first base (1-based)
    3: intron last base  (1-based, inclusive)
    4: strand    (0=undefined, 1=+, 2=-)
    5: intron motif (0=non-canonical, 1=GT/AG +, 2=CT/AC -, 3=GC/AG +,
                     4=CT/GC -, 5=AT/AC +, 6=GT/AT -)
    6: annotated (0=novel, 1=annotated)
    7: # uniquely-mapping reads crossing the junction
    8: # multi-mapping reads crossing the junction
    9: maximum spliced alignment overhang

This script emits the same 2-column TSV that ``bed12_to_junctions.py`` produces
on the HISAT2 path, so downstream rules are aligner-agnostic:

    <chrom>:<1-based donor>:<1-based inclusive end>:<strand>\\t<reads>

Numerically, STAR's 1-based inclusive end equals the 0-based half-open exclusive
end consumed by ``filter_junctions._parse_junction_id`` (``end = int(parts[2])``).

Strand resolution (Issue #374):
    Previously inline awk emitted strand `.` whenever col 4 = 0, silently
    contaminating the candidate set — ``assemble_contigs.py`` then took
    forward-orientation flanking sequence for what may have been minus-strand
    introns. The rescue uses col 5 (intron motif) when col 4 = 0:
      motif ∈ {1,3,5} → '+' (GT/AG, GC/AG, AT/AC)
      motif ∈ {2,4,6} → '-' (CT/AC, CT/GC, GT/AT)
      motif = 0 or unknown → drop (truly non-canonical; cannot infer)

See: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


_MOTIF_TO_STRAND = {1: "+", 2: "-", 3: "+", 4: "-", 5: "+", 6: "-"}


def _resolve_strand(strand_code: int, motif_code: int) -> str | None:
    """Return '+'/'-' for this junction, or None if it should be dropped.

    STAR's col 4 takes priority — if STAR inferred strand directly, use it.
    Only ``strand_code == 0`` (STAR could not infer) triggers the motif rescue.
    Motif 0 (truly non-canonical), unknown motif codes, and unknown strand codes
    all return None so the caller drops the record rather than emitting strand
    '.', which downstream `bedtools getfasta` would treat as forward orientation.
    """
    if strand_code == 1:
        return "+"
    if strand_code == 2:
        return "-"
    if strand_code == 0:
        return _MOTIF_TO_STRAND.get(motif_code)
    return None


def convert_sj_to_junctions(input_path: str | Path, output_path: str | Path) -> int:
    """Convert STAR ``SJ.out.tab`` to a 2-column junctions TSV.

    Args:
        input_path:  STAR ``SJ.out.tab`` (9-column TSV).
        output_path: Destination TSV (``<chrom>:<donor>:<end>:<strand>\\treads``).

    Returns:
        Number of junctions written.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_direct = 0
    n_rescued = 0
    n_dropped_motif = 0
    n_dropped_zero_reads = 0
    n_dropped_malformed = 0

    with input_path.open() as fh_in, output_path.open("w") as fh_out:
        for line in fh_in:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                n_dropped_malformed += 1
                continue

            try:
                strand_code = int(fields[3])
                motif_code = int(fields[4])
                reads = int(fields[6])
            except ValueError:
                n_dropped_malformed += 1
                continue

            if reads <= 0:
                n_dropped_zero_reads += 1
                continue

            strand = _resolve_strand(strand_code, motif_code)
            if strand is None:
                n_dropped_motif += 1
                continue

            if strand_code in (1, 2):
                n_direct += 1
            else:
                n_rescued += 1

            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            fh_out.write(f"{chrom}:{start}:{end}:{strand}\t{reads}\n")
            n_written += 1

    log.info(
        "Wrote %d junctions to %s "
        "(direct=%d, rescued_via_motif=%d, "
        "dropped_motif_unknown=%d, dropped_zero_reads=%d, dropped_malformed=%d)",
        n_written,
        output_path,
        n_direct,
        n_rescued,
        n_dropped_motif,
        n_dropped_zero_reads,
        n_dropped_malformed,
    )
    return n_written


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert STAR SJ.out.tab to a 2-column junctions TSV."
    )
    parser.add_argument("--input", required=True, help="STAR SJ.out.tab input")
    parser.add_argument("--output", required=True, help="junctions TSV output")
    args = parser.parse_args()
    convert_sj_to_junctions(args.input, args.output)


if __name__ == "__main__":
    _cli_main()
