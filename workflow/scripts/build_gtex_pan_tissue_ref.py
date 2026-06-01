#!/usr/bin/env python3
"""build_gtex_pan_tissue_ref.py — Build a pan-tissue novel-junction blacklist BED
from the Snaptron gtexv2 endpoint (recount3 / GTEx v8 / hg38 / 19,214 samples).

One-shot reference builder (NOT part of the per-patient Snakemake DAG). See
docs/gtex_pan_tissue_build.md for provenance + the coordinate transform, pinned
against real data and validated 259/259 via Issue #225's helpers (Issue #211).

Output BED is BED6 (chrom start end name score strand), drop-in with
references/reference_junctions.bed and filter_junctions.py's 4-tuple reference
reader. Snaptron carries strand, so (chrom,start,end,strand) is used end-to-end.
"""

import argparse
import json
import logging
import shutil
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SNAPTRON_GTEXV2_URL = "https://snaptron.cs.jhu.edu/gtexv2/snaptron"

# Columns this build addresses by name (Snaptron gtexv2 TSV carries a header row).
_REQUIRED_COLS = ("chromosome", "start", "end", "strand", "samples_count")

# samples_count thresholds reported in the QC sensitivity sweep.
QC_SWEEP_THRESHOLDS = (1, 2, 5, 10, 20)

# hg38 (UCSC) primary chromosome sizes — the genome-wide region set (one query each).
# chr22 = 50,818,468 matches Issue #225's SNAPTRON_CHR22_REGION exactly.
HG38_CHROM_SIZES: Dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def build_col_index(header_line: str) -> Dict[str, int]:
    """Map Snaptron column names -> positional indices.

    Lookup is by name so a future Snaptron column reorder does not silently
    shift the parse. Raises ValueError if any required column is absent.
    """
    cols = header_line.rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(cols)}
    missing = [c for c in _REQUIRED_COLS if c not in idx]
    if missing:
        raise ValueError(
            f"Snaptron header missing columns {missing}; got {cols}"
        )
    return idx


def parse_snaptron_line(
    fields: List[str], col_idx: Dict[str, int]
) -> Optional[Tuple[str, int, int, str, int]]:
    """Parse one Snaptron data row into (chrom, bed_start, bed_end, strand, samples_count).

    Snaptron start/end are 1-based inclusive intron donor/acceptor; convert to
    BED 0-based half-open via bed_start = start - 1, bed_end = end (the transform
    pinned + validated 259/259 in Issue #225's snaptron_to_key_set). Returns None
    on a short or non-integer row.
    """
    try:
        chrom = fields[col_idx["chromosome"]]
        start = int(fields[col_idx["start"]])
        end = int(fields[col_idx["end"]])
        strand = fields[col_idx["strand"]]
        samples_count = int(fields[col_idx["samples_count"]])
    except (IndexError, ValueError):
        return None
    return chrom, start - 1, end, strand, samples_count


def accumulate_union(
    lines: Iterable[str],
    min_samples: int,
    restrict_chrom: Optional[str] = None,
) -> Tuple[Set[Tuple[str, int, int, str]], Counter]:
    """Consume Snaptron TSV lines (first line = header) into a junction-key set.

    Returns (keys, sweep):
      - keys: {(chrom, bed_start, bed_end, strand)} for samples_count >= min_samples
        (and chrom == restrict_chrom, if given).
      - sweep: Counter mapping each QC_SWEEP_THRESHOLDS value t -> number of
        (restrict-filtered) junctions with samples_count >= t. Drives the QC sidecar.

    Each Snaptron junction is one row, so within a single region there is no
    dedup; the set still guards against any Snaptron-side duplicate.
    """
    it = iter(lines)
    try:
        header = next(it)
    except StopIteration:
        return set(), Counter()
    col_idx = build_col_index(header)

    keys: Set[Tuple[str, int, int, str]] = set()
    sweep: Counter = Counter()
    n_rows = 0
    for line in it:
        if not line:
            continue
        parsed = parse_snaptron_line(line.split("\t"), col_idx)
        if parsed is None:
            continue
        chrom, bstart, bend, strand, sc = parsed
        if restrict_chrom is not None and chrom != restrict_chrom:
            continue
        n_rows += 1
        for t in QC_SWEEP_THRESHOLDS:
            if sc >= t:
                sweep[t] += 1
        if sc >= min_samples:
            keys.add((chrom, bstart, bend, strand))
    log.info(
        "Accumulated %d junctions (min_samples=%d) from %d parsed rows",
        len(keys), min_samples, n_rows,
    )
    return keys, sweep


def fetch_snaptron_region(
    region: str,
    endpoint: str = SNAPTRON_GTEXV2_URL,
    timeout_s: int = 180,
    retries: int = 2,
) -> Iterator[str]:
    """Stream a Snaptron region query response as decoded, newline-stripped lines.

    Ported from Issue #225's fetch_snaptron_chr22: connect with a bounded retry
    (transient URLError -> 5s backoff), then yield lines lazily so a whole
    chromosome's TSV is never fully held in memory. The connection is retried
    only at open time; a mid-stream failure surfaces to the caller.
    """
    url = f"{endpoint}?regions={region}"
    attempt = 0
    resp = None
    while True:
        attempt += 1
        try:
            log.info("Snaptron query (attempt %d/%d): %s", attempt, retries, url)
            resp = urllib.request.urlopen(url, timeout=timeout_s)
            break
        except urllib.error.URLError as e:
            if attempt >= retries:
                raise
            log.warning("Snaptron query failed (%s); retrying in 5s", e)
            time.sleep(5)
    with resp:
        for raw in resp:
            yield raw.decode("utf-8").rstrip("\n")


def write_bed6(keys: Set[Tuple[str, int, int, str]], path: str) -> None:
    """Write the pan-tissue blacklist as sorted BED6 (chrom,start,end,name,0,strand).

    name = 'chrom:start-end:strand' and score 0 — identical in shape to
    build_reference_junctions.py, so filter_junctions.py's 4-tuple reference
    reader consumes it unchanged.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(keys, key=lambda k: (k[0], k[1], k[2], k[3]))
    with out.open("w") as fh:
        for chrom, start, end, strand in ordered:
            name = f"{chrom}:{start}-{end}:{strand}"
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    log.info("Wrote %d junctions to %s", len(ordered), out)


def write_qc_sidecar(sweep: Counter, n_junctions: int, path: str) -> None:
    """Write a QC TSV: total union size + the samples_count sensitivity sweep."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"n_junctions\t{n_junctions}\n")
        fh.write("min_samples_count\tn_junctions\n")
        for t in QC_SWEEP_THRESHOLDS:
            fh.write(f"{t}\t{sweep.get(t, 0)}\n")
    log.info("Wrote QC sidecar to %s", out)


def _write_chrom_sweep(sweep: Counter, path: Path) -> None:
    """Persist one chromosome's QC sweep next to its part, so a resumed run can
    fold in already-built chromosomes without re-querying them."""
    path.write_text(json.dumps({str(k): v for k, v in sweep.items()}))


def _read_chrom_sweep(path: Path) -> Counter:
    return Counter({int(k): v for k, v in json.loads(path.read_text()).items()})


def _merge_parts(part_dir: Path, chroms: List[str], bed_path: str) -> int:
    """Concatenate per-chromosome BED parts into the final BED in chrom-name
    order, streaming line-by-line. Each part is internally sorted and holds a
    single (disjoint) chromosome, so this reproduces a global ``sorted(union)``
    without ever materialising the union. Returns the total junction count."""
    out = Path(bed_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with out.open("w") as dst:
        for chrom in sorted(set(chroms)):
            part = part_dir / f"{chrom}.bed"
            if not part.exists():
                continue
            with part.open() as src:
                for line in src:
                    dst.write(line)
                    total += 1
    log.info("Merged %d junctions into %s", total, out)
    return total


def build(
    regions: List[str],
    bed_path: str,
    qc_path: str,
    min_samples: int = 1,
    restrict_chrom: Optional[str] = None,
    endpoint: str = SNAPTRON_GTEXV2_URL,
    line_source=None,
    resume: bool = False,
) -> int:
    """Query each region, stream the kept junctions to a per-chromosome part on
    disk, then merge the parts into the final sorted BED6 + QC. Returns the
    junction count.

    Regions are disjoint per-chromosome queries, so the union never needs to be
    held in memory: peak memory is one chromosome's key set (genome-wide that is
    ~1 GB, vs ~tens of GB for the full union — Issue #211). ``line_source(region)
    -> iterable[str]`` is injectable for tests; default is the live Snaptron
    fetch.

    Parts live under ``<bed_path>.parts/`` and are removed on success. With
    ``resume=True`` a chromosome whose ``.done`` marker already exists is folded
    in without re-querying — so an interrupted multi-hour build picks up where it
    stopped. ``resume=False`` (default) clears any stale parts and rebuilds fresh.

    Assumes one region per chromosome (the default genome-wide region set);
    windowed multi-region-per-chromosome queries are not supported here.
    """
    if line_source is None:
        line_source = lambda region: fetch_snaptron_region(region, endpoint=endpoint)

    # Parts key on chromosome only, so two regions on one chromosome would
    # collide and silently truncate a window (the second write_bed6 overwrites
    # the first). Refuse rather than corrupt the BED/QC.
    region_chroms = [r.split(":")[0] for r in regions]
    dups = sorted({c for c in region_chroms if region_chroms.count(c) > 1})
    if dups:
        raise ValueError(
            f"Multiple regions target the same chromosome(s) {dups}; the "
            f"per-chromosome streaming build requires one region per chromosome. "
            f"Merge the windows into a single region per chromosome."
        )

    part_dir = Path(bed_path + ".parts")
    if not resume and part_dir.exists():
        shutil.rmtree(part_dir)
    part_dir.mkdir(parents=True, exist_ok=True)

    # Stamp the parameters that define a part's contents and refuse to --resume
    # parts built under different parameters — otherwise the merge would silently
    # mix regimes (e.g. min_samples=10 early chromosomes + a min_samples=1 resumed
    # tail). endpoint and restrict_chrom change part contents too.
    params = {"min_samples": min_samples, "endpoint": endpoint,
              "restrict_chrom": restrict_chrom}
    params_path = part_dir / "params.json"
    if resume and params_path.exists():
        prior = json.loads(params_path.read_text())
        if prior != params:
            raise ValueError(
                f"Cannot --resume: parts in {part_dir} were built with {prior} but "
                f"this run uses {params}. Delete {part_dir} to rebuild, or re-run "
                f"with matching parameters."
            )
    elif resume and any(part_dir.glob("*.done")):
        raise ValueError(
            f"Cannot --resume: parts in {part_dir} have no params stamp, so their "
            f"build parameters cannot be verified. Delete {part_dir} to rebuild."
        )
    params_path.write_text(json.dumps(params))

    chroms: List[str] = []
    sweeps: List[Counter] = []
    for region in regions:
        chrom = region.split(":")[0]
        chroms.append(chrom)
        done = part_dir / f"{chrom}.done"
        sweep_path = part_dir / f"{chrom}.sweep.json"
        if resume and done.exists():
            sweeps.append(_read_chrom_sweep(sweep_path))
            log.info("Resume: %s already built, skipping query", chrom)
            continue
        keys, region_sweep = accumulate_union(
            line_source(region), min_samples=min_samples, restrict_chrom=restrict_chrom
        )
        # Write part + sweep first, touch the .done marker last: a crash between
        # the two leaves an incomplete part with no marker, so resume re-queries.
        write_bed6(keys, str(part_dir / f"{chrom}.bed"))
        _write_chrom_sweep(region_sweep, sweep_path)
        done.touch()
        sweeps.append(region_sweep)
        log.info("Region %s: +%d junctions (chromosome part written)", region, len(keys))

    total = _merge_parts(part_dir, chroms, bed_path)
    total_sweep: Counter = Counter()
    for s in sweeps:
        total_sweep += s
    write_qc_sidecar(total_sweep, n_junctions=total, path=qc_path)
    shutil.rmtree(part_dir, ignore_errors=True)
    return total


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a Snaptron gtexv2 pan-tissue novel-junction blacklist BED."
    )
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-qc", required=True)
    parser.add_argument("--min-samples", type=int, default=1,
                        help="Keep junctions seen in >= this many samples (default 1).")
    parser.add_argument("--endpoint", default=SNAPTRON_GTEXV2_URL)
    parser.add_argument(
        "--region", action="append", default=None,
        help="Snaptron region (e.g. chr22:1-50818468). Repeatable. "
             "Default: all hg38 primary chromosomes (genome-wide).")
    parser.add_argument("--restrict-chrom", default=None,
                        help="Emit only this chromosome (chr22 fixture build).")
    parser.add_argument("--resume", action="store_true",
                        help="Skip chromosomes already built in a prior interrupted run "
                             "(reuses <output-bed>.parts/). Use for the multi-hour "
                             "genome-wide build so a dropped connection doesn't restart it.")
    args = parser.parse_args()

    if args.region:
        regions = args.region
    else:
        regions = [f"{c}:1-{size}" for c, size in HG38_CHROM_SIZES.items()]

    n = build(
        regions=regions, bed_path=args.output_bed, qc_path=args.output_qc,
        min_samples=args.min_samples, restrict_chrom=args.restrict_chrom,
        endpoint=args.endpoint, resume=args.resume,
    )
    log.info("Done: %d junctions in the pan-tissue union", n)


if __name__ == "__main__":
    _cli_main()
