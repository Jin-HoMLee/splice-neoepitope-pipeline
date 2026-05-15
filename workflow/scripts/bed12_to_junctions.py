#!/usr/bin/env python3
"""bed12_to_junctions.py — Convert a regtools BED12 file to junctions.tsv.

regtools' ``junctions extract`` output is BED12 where columns 2–3 are the
**anchor outer boundaries** (chromStart, chromEnd), NOT the intron donor/
acceptor coordinates. The actual intron donor/acceptor is derived from
``blockSizes`` and ``blockStarts``:

    donor    (0-based)            = chromStart + blockSizes[0]
    acceptor (0-based exclusive)  = chromStart + blockStarts[1]

This script emits a 2-column TSV:

    <chrom>:<donor_1based>:<acceptor_0based_exclusive>:<strand>\\t<reads>

The format matches what ``filter_junctions._parse_junction_id`` consumes:
``start = int(parts[1]) - 1`` (1-based donor → 0-based intron start) and
``end = int(parts[2])`` (0-based half-open intron end).

See: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370
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


def convert_bed12_to_junctions(input_path: str | Path, output_path: str | Path) -> int:
    """Convert a regtools BED12 file to a 2-column junctions TSV.

    Args:
        input_path:  BED12 file produced by ``regtools junctions extract``.
        output_path: Destination TSV (``<chrom>:<donor>:<end>:<strand>\\treads``).

    Returns:
        Number of junctions written.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with input_path.open() as fh_in, output_path.open("w") as fh_out:
        for line in fh_in:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            try:
                reads = int(fields[4])
            except ValueError:
                continue
            if reads <= 0:
                continue

            chrom = fields[0]
            chrom_start = int(fields[1])
            strand = fields[5]
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            if len(block_sizes) < 2 or len(block_starts) < 2:
                continue

            donor_0based = chrom_start + block_sizes[0]
            acceptor_0based_exclusive = chrom_start + block_starts[1]

            fh_out.write(
                f"{chrom}:{donor_0based + 1}:{acceptor_0based_exclusive}:{strand}"
                f"\t{reads}\n"
            )
            n_written += 1

    log.info("Wrote %d junctions to %s", n_written, output_path)
    return n_written


def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    convert_bed12_to_junctions(
        snakemake.input.bed,  # type: ignore[name-defined]  # noqa: F821
        snakemake.output.junctions,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert regtools BED12 to a 2-column junctions TSV."
    )
    parser.add_argument("--input", required=True, help="regtools BED12 input")
    parser.add_argument("--output", required=True, help="junctions TSV output")
    args = parser.parse_args()
    convert_bed12_to_junctions(args.input, args.output)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
