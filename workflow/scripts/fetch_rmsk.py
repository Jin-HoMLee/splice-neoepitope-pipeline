#!/usr/bin/env python3
"""fetch_rmsk.py - Fetch a chromosome slice of the UCSC RepeatMasker track as BED.

Used to categorize junctions by repeat overlap (Issue #919): the NH-uniqueness
prefilter is meant to remove multimapper-driven junction calls at repeat copies,
so "are the junctions it removes actually in repeats?" is the check that tells us
whether it does what it claims.

Source is the UCSC REST API (``/getData/track``) rather than the goldenPath
``rmsk.txt.gz`` dump, because the API takes a ``chrom=`` filter server-side. The
dump is genome-wide (all 24 chromosomes plus alts), so a chr22 slice off it would
mean downloading and discarding ~98% of the table.

Truncation is the hazard worth guarding
---------------------------------------
The API silently caps its response and reports the cap in a ``maxItemsLimit``
key, which is absent on a complete response. A truncated pull looks exactly like
a complete one - a valid BED, just short - and would understate repeat overlap
without any error. So a ``maxItemsLimit`` in the payload is a hard failure, not
a warning. (chr22 carries 79,521 repeats, well under the default cap; a future
cap change, or a larger chromosome, is what this guards.)

Output is BED6 + 2 (``repClass``, ``repFamily``), coordinate-sorted:

    chrom  start  end  repName  swScore  strand  repClass  repFamily

UCSC ``genoStart``/``genoEnd`` are already 0-based half-open, i.e. BED
convention, so no coordinate arithmetic is needed here (contrast the regtools
BED12 anchor-outer trap in ``bed12_to_junctions.py``).
"""

import argparse
import json
import logging
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

UCSC_API = "https://api.genome.ucsc.edu/getData/track"

# Columns of the emitted BED. BED6 plus the two rmsk fields the #919 analysis
# groups on.
BED_COLUMNS = (
    "chrom",
    "start",
    "end",
    "repName",
    "swScore",
    "strand",
    "repClass",
    "repFamily",
)


class RmskTruncatedError(RuntimeError):
    """The UCSC API capped its response, so the track slice is incomplete."""


def build_track_url(genome: str, chrom: str, track: str = "rmsk") -> str:
    """URL for a whole-chromosome slice of a UCSC track.

    UCSC's API wants ``;``-separated params, not ``&``.
    """
    return f"{UCSC_API}?genome={genome};track={track};chrom={chrom}"


def extract_items(payload: Dict[str, Any], chrom: str, track: str = "rmsk") -> List[Dict[str, Any]]:
    """Pull the track's feature list out of an API payload, refusing a truncated one.

    Raises:
        RmskTruncatedError: the payload carries ``maxItemsLimit`` (the API's
            truncation marker), or its item count disagrees with the
            ``itemsReturned`` the API self-reports.
        KeyError: the payload carries no feature list under ``track`` - i.e. the
            request succeeded but named a track UCSC does not serve.
    """
    if "maxItemsLimit" in payload:
        raise RmskTruncatedError(
            f"UCSC capped the {track} response for {chrom} at "
            f"{payload['maxItemsLimit']} items - the slice is incomplete and would "
            f"silently understate repeat overlap. Page the request or use the "
            f"goldenPath dump."
        )

    if track not in payload:
        raise KeyError(
            f"UCSC payload for {chrom} has no '{track}' feature list "
            f"(keys: {sorted(payload)})"
        )

    items = payload[track]

    # The API self-reports its count; a disagreement means we are reading a
    # payload shape we do not understand, which is not something to paper over.
    reported = payload.get("itemsReturned")
    if reported is not None and reported != len(items):
        raise RmskTruncatedError(
            f"UCSC reported itemsReturned={reported} for {chrom} but the payload "
            f"carries {len(items)} items"
        )

    return items


def items_to_bed_rows(items: List[Dict[str, Any]], chrom: str) -> List[List[str]]:
    """Convert rmsk API records to coordinate-sorted BED rows.

    Records on a chromosome other than ``chrom`` are dropped - the API filters
    server-side, so this only fires if UCSC ever widens the response.
    """
    rows = []
    for item in items:
        if item.get("genoName") != chrom:
            continue
        rows.append(
            [
                item["genoName"],
                str(item["genoStart"]),
                str(item["genoEnd"]),
                item.get("repName", "."),
                str(item.get("swScore", 0)),
                item.get("strand", "."),
                item.get("repClass", "."),
                item.get("repFamily", "."),
            ]
        )
    rows.sort(key=lambda r: (int(r[1]), int(r[2])))
    return rows


def fetch_payload(url: str, timeout: int = 120) -> Dict[str, Any]:
    """GET a UCSC API URL and parse the JSON body."""
    log.info("Fetching %s", url)
    with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310
        return json.loads(response.read().decode("utf-8"))


def write_bed(rows: List[List[str]], output_path: Path) -> None:
    """Write BED rows, creating the parent directory."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", default="hg38", help="UCSC genome (default: hg38)")
    parser.add_argument("--chrom", required=True, help="Chromosome to fetch, e.g. chr22")
    parser.add_argument("--output", required=True, type=Path, help="Destination BED path")
    parser.add_argument("--track", default="rmsk", help="UCSC track name (default: rmsk)")
    args = parser.parse_args(argv)

    payload = fetch_payload(build_track_url(args.genome, args.chrom, args.track))
    items = extract_items(payload, args.chrom, args.track)
    rows = items_to_bed_rows(items, args.chrom)

    if not rows:
        # An empty repeat track for a real chromosome means the fetch is wrong,
        # not that the chromosome has no repeats.
        log.error("No %s features returned for %s - refusing to write an empty BED", args.track, args.chrom)
        return 1

    write_bed(rows, args.output)
    log.info(
        "Wrote %d %s features for %s (%s assembly) to %s",
        len(rows),
        args.track,
        args.chrom,
        payload.get("genome", args.genome),
        args.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
