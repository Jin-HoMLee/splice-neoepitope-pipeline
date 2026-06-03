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
import os
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Bulk gtexv2 junction dump (recount3 / GTEx v8 / hg38), the source the region API itself
# serves from. Fetched via remote tabix per chromosome. Served static with a real
# Content-Length (NOT chunked), so a short read raises in htslib/tabix instead of silently
# truncating — the region API's chunked+no-Content-Length transport was the root cause of the
# genome-wide silent truncation (Issue #211). BGZF, tabix-indexed (.bgz.tbi alongside).
SNAPTRON_GTEXV2_BULK_BGZF = "https://snaptron.cs.jhu.edu/data/gtexv2/junctions.bgz"

# The bulk file's 17-column header (junctions.header.tsv). Prepended to tabix output (which
# carries no header) so build_col_index maps columns BY NAME: the bulk file's column order
# differs from the region API, but the needed names are all present, so the parser is unchanged.
BULK_HEADER = (
    "snaptron_id\tchromosome\tstart\tend\tlength\tstrand\tannotated\t"
    "left_motif\tright_motif\tleft_annotated\tright_annotated\t"
    "samples\tsamples_count\tcoverage_sum\tcoverage_avg\tcoverage_median\tsource_dataset_id"
)

# Columns this build addresses by name (the prepended header carries a header row).
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


def fetch_chromosome_bgz(
    chrom: str,
    bgz_url: str = SNAPTRON_GTEXV2_BULK_BGZF,
    retries: int = 3,
    work_dir: Optional[str] = None,
) -> Iterator[str]:
    """Yield BULK_HEADER then the bulk file's data rows for one chromosome, via tabix.

    ``chrom`` is a contig name (``chr22``) or any tabix region string
    (``chr22:1-50818468``) — both are passed to ``tabix`` verbatim; the CLI uses bare
    contig names, the tests exercise the windowed form.

    Runs ``tabix <bgz_url> <chrom>`` — htslib byte-range-streams just that contig over
    HTTPS. Three fail-loud truncation defenses (Issue #211); none may silently undercount:
      * transport short read — the static file's real Content-Length makes htslib error and
        tabix exit non-zero -> retried, then RuntimeError;
      * contig absent from the index / name mismatch (UCSC ``chr22`` vs ENSEMBL ``22``) /
        empty response — tabix exits 0 with NO rows -> the zero-row guard raises (a primary
        contig is never legitimately empty), so no empty ``.done`` part is ever written;
      * tabix not installed — a PATH preflight raises an actionable error (not FileNotFoundError).
    The contig TSV and the remote ``.tbi`` index htslib downloads both land in a private temp
    dir (passed as ``cwd``) co-located with the build via ``work_dir``, and are removed
    together. The contig is fetched with whole-operation retry (a failed run is re-run from
    scratch — never a partial yield). BULK_HEADER is prepended so accumulate_union's by-name
    column mapping is unchanged.
    """
    if shutil.which("tabix") is None:
        raise RuntimeError(
            "tabix not found on PATH. Install htslib/tabix (Ubuntu: 'apt-get install tabix'; "
            "conda: 'conda install -c bioconda htslib') — the pan-tissue builder fetches the "
            "bulk junctions.bgz via remote tabix (Issue #211)."
        )
    yield BULK_HEADER
    n_rows = 0
    with tempfile.TemporaryDirectory(dir=work_dir, prefix=f"gtex_{chrom}_") as td:
        tsv = os.path.join(td, f"{chrom}.tsv")
        attempt = 0
        while True:
            attempt += 1
            log.info("tabix fetch (attempt %d/%d): %s %s", attempt, retries, bgz_url, chrom)
            with open(tsv, "wb") as out:
                # cwd=td so the remote .tbi index htslib downloads lands in the temp dir
                # (not CWD, where it would risk an accidental commit) and is cleaned up here.
                proc = subprocess.run(
                    ["tabix", bgz_url, chrom], stdout=out, stderr=subprocess.PIPE, cwd=td
                )
            if proc.returncode == 0:
                break
            err = (proc.stderr or b"").decode("utf-8", "replace").strip()
            if "No space left on device" in err or "ENOSPC" in err:
                raise RuntimeError(
                    f"tabix ran out of disk staging {chrom} under {td}: {err}. Point "
                    f"TMPDIR / the build's output partition at >~10 GB free (chr1 ~6 GB)."
                )
            if attempt >= retries:
                raise RuntimeError(
                    f"tabix failed for {chrom} after {retries} attempts "
                    f"(exit {proc.returncode}): {err}"
                )
            log.warning("tabix %s failed (exit %d), retrying in 5s: %s",
                        chrom, proc.returncode, err[:200])
            time.sleep(5)
        with open(tsv, "rt") as inp:
            for line in inp:
                n_rows += 1
                yield line.rstrip("\n")
    if n_rows == 0:
        raise RuntimeError(
            f"tabix returned 0 rows for {chrom} (exit 0) — contig absent from the index, a "
            f"name mismatch (UCSC chr22 vs ENSEMBL 22?), or an empty response. A primary "
            f"contig is never legitimately empty (Issue #211 silent-undercount guard)."
        )


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
        # Format: line 1 is a key-value header (n_junctions\t<total>); the remaining
        # lines are a 2-column table (min_samples_count\tn_junctions) for the sweep.
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
    bgz_url: str = SNAPTRON_GTEXV2_BULK_BGZF,
    line_source=None,
    resume: bool = False,
) -> int:
    """Fetch each region (a chromosome), stream the kept junctions to a
    per-chromosome part on disk, then merge the parts into the final sorted BED6
    + QC. Returns the junction count.

    Regions are disjoint per-chromosome fetches, so the union never needs to be
    held in memory: peak memory is one chromosome's key set (genome-wide that is
    ~1 GB, vs ~tens of GB for the full union — Issue #211). ``line_source(region)
    -> iterable[str]`` is injectable for tests; default is a remote-tabix fetch of
    that contig from the bulk junctions.bgz (BULK_HEADER prepended).

    Parts live under ``<bed_path>.parts/`` and are removed on success. With
    ``resume=True`` a chromosome whose ``.done`` marker already exists is folded
    in without re-querying — so an interrupted multi-hour build picks up where it
    stopped. ``resume=False`` (default) clears any stale parts and rebuilds fresh.

    Assumes one region per chromosome (the default genome-wide region set);
    windowed multi-region-per-chromosome queries are not supported here.
    """
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

    if line_source is None:
        # Stage each contig's TSV + the remote .tbi under part_dir (same partition as the
        # output) rather than $TMPDIR, so disk headroom follows the operator's chosen path.
        line_source = lambda region: fetch_chromosome_bgz(
            region, bgz_url=bgz_url, work_dir=str(part_dir)
        )

    # Stamp the parameters that define a part's contents and refuse to --resume
    # parts built under different parameters — otherwise the merge would silently
    # mix regimes (e.g. min_samples=10 early chromosomes + a min_samples=1 resumed
    # tail). bgz_url and restrict_chrom change part contents too.
    params = {"min_samples": min_samples, "bgz_url": bgz_url,
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
            # .done implies write_bed6 ran (write order: bed -> sweep -> done.touch),
            # but a post-write deletion of the .bed part leaves a stale marker that
            # _merge_parts silently skips — dropping this chromosome from the final BED
            # while its sweep is still counted (silent undercount, Issue #211). Refuse.
            bed_part = part_dir / f"{chrom}.bed"
            if not bed_part.exists():
                raise RuntimeError(
                    f"Cannot --resume: {done.name} marks {chrom} as built, but its BED "
                    f"part {bed_part} is missing — _merge_parts would silently omit "
                    f"{chrom} from the final BED while its sweep is still counted, "
                    f"undercounting the union (Issue #211 silent-undercount guard). "
                    f"Delete {done} (and {sweep_path}) to re-query {chrom}, or delete "
                    f"{part_dir} to rebuild."
                )
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
    parser.add_argument("--bgz-url", default=SNAPTRON_GTEXV2_BULK_BGZF,
                        help="Snaptron bulk junctions.bgz URL fetched via remote tabix.")
    parser.add_argument(
        "--region", action="append", default=None,
        help="Contig to fetch via tabix (e.g. chr22). Repeatable. "
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
        regions = list(HG38_CHROM_SIZES.keys())  # bare contig names for tabix

    n = build(
        regions=regions, bed_path=args.output_bed, qc_path=args.output_qc,
        min_samples=args.min_samples, restrict_chrom=args.restrict_chrom,
        bgz_url=args.bgz_url, resume=args.resume,
    )
    log.info("Done: %d junctions in the pan-tissue union", n)


if __name__ == "__main__":
    _cli_main()
