"""Evaluate the anchor-vs-intron repeat-embedding score on the chr22 fixture.

[Issue #1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116).

Answers three questions, in increasing order of how much they can embarrass us:

1. **Does the score separate GENCODE-annotated from unannotated junctions?**
   Annotated junctions are the only local ground truth for "real". If the score
   is measuring anything, real junctions should be *less* repeat-embedded
   (higher Hamming) than the unannotated pool.

2. **Does it agree with the NH gate?** Score the junctions the `[NH]==1` filter
   removes (the #919 differential). Agreement means the gate was groping for a
   real signal; disagreement is independent evidence for
   [#1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)'s
   conclusion that the gate removes the wrong things. **Both outcomes inform.**

3. **Is the score measuring repeats at all?** This is the falsifier, and it is
   the reason this script can fail rather than merely produce numbers. The #919
   run already annotated every junction with **RepeatMasker** overlap
   (`donor_in_repeat` / `acceptor_in_repeat`) - an *independent* source of truth,
   built from a completely different method than Hamming distance. If low-Hamming
   junctions are NOT enriched for RepeatMasker overlap, then this score is not a
   repeat detector and the whole premise of #1116 is wrong. Say so if so.

Coordinate check runs first: if the 1-based-inclusive reading of the junction id
is off by even one base, the introns will not begin `GT` and end `AG`. A score
computed on shifted coordinates would be confident garbage.
"""

import argparse
import gzip
import statistics
import sys
from pathlib import Path

from repeat_score import (
    DEFAULT_ANCHOR,
    JunctionTooCloseToContigEnd,
    repeat_embedding_score,
)

REPO = Path(__file__).resolve().parents[3]
CHR22_FA = REPO / "resources/test/chr22.fa"
CHR22_GTF = REPO / "resources/test/chr22.gtf.gz"
EXP919 = REPO / "research/experiments/issue_919_nh_uniqueness_filter/outputs"


def load_chr22(path: Path) -> str:
    """Load the single-contig chr22 FASTA into one string (case preserved)."""
    seq = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def parse_junction_id(jid: str):
    """`chr22:<intron_start>:<intron_end>:<strand>`, 1-based inclusive."""
    chrom, start, end, strand = jid.split(":")
    return chrom, int(start), int(end), strand


def read_junctions(path: Path) -> dict:
    """junction_id -> read count."""
    out = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out[parts[0]] = int(parts[1])
    return out


def annotated_introns(gtf: Path) -> set:
    """GENCODE intron set, derived from consecutive exons of each transcript.

    Built independently of the pipeline's own `annotated` flag on purpose: an
    annotation check that reuses the pipeline's annotation would be circular.
    """
    exons = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = None
            for field in f[8].split(";"):
                field = field.strip()
                if field.startswith("transcript_id"):
                    tid = field.split('"')[1]
                    break
            if tid is None:
                continue
            exons.setdefault((tid, f[0], f[6]), []).append((int(f[3]), int(f[4])))

    introns = set()
    for (_tid, chrom, strand), blocks in exons.items():
        blocks.sort()
        for (_s1, e1), (s2, _e2) in zip(blocks, blocks[1:]):
            # intron spans the bases between two exons, 1-based inclusive
            introns.add(f"{chrom}:{e1 + 1}:{s2 - 1}:{strand}")
    return introns


def motif_of(ref: str, start: int, end: int) -> str:
    """Dinucleotides at the intron ends, as they appear on the + strand."""
    return (ref[start - 1 : start + 1] + ".." + ref[end - 2 : end]).upper()


def check_coordinates(ref: str, junction_ids, limit=400) -> dict:
    """Falsifier: canonical introns must read GT..AG (+) or CT..AC (-).

    An off-by-one in the coordinate convention silently shifts every window and
    produces a score that is precise and meaningless. If this check fails, stop.
    """
    counts = {}
    for jid in list(junction_ids)[:limit]:
        _c, s, e, _strand = parse_junction_id(jid)
        if s - 1 < 0 or e > len(ref):
            continue
        counts[motif_of(ref, s, e)] = counts.get(motif_of(ref, s, e), 0) + 1
    return counts


def summarize(name, values):
    if not values:
        return f"  {name:<34} n=0"
    return (
        f"  {name:<34} n={len(values):<5} "
        f"median={statistics.median(values):>5.1f}  mean={statistics.mean(values):>5.2f}  "
        f"frac_min<=2={sum(1 for v in values if v <= 2) / len(values):.3f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchor", type=int, default=DEFAULT_ANCHOR)
    ap.add_argument("--sample", default="tumor", choices=["tumor", "normal"])
    args = ap.parse_args()

    ref = load_chr22(CHR22_FA)
    print(f"chr22 loaded: {len(ref):,} bp\n")

    off = read_junctions(EXP919 / f"raw_junctions.{args.sample}.filter_off.tsv")
    on = read_junctions(EXP919 / f"raw_junctions.{args.sample}.filter_on.tsv")
    removed_by_nh = set(off) - set(on)
    print(f"junctions (filter off): {len(off):,}")
    print(f"junctions (filter on):  {len(on):,}")
    print(f"removed by the NH gate: {len(removed_by_nh):,}\n")

    # --- Falsifier 1: are we even reading the coordinates right? ---
    motifs = check_coordinates(ref, off)
    top = sorted(motifs.items(), key=lambda kv: -kv[1])[:4]
    canonical = sum(v for k, v in motifs.items() if k in ("GT..AG", "CT..AC"))
    total = sum(motifs.values())
    print("COORDINATE CHECK (canonical motif must dominate, else the offsets are wrong)")
    for k, v in top:
        print(f"  {k}  {v}")
    frac = canonical / total if total else 0
    print(f"  canonical GT..AG / CT..AC: {canonical}/{total} = {frac:.3f}")
    if frac < 0.9:
        print("\n  ABORT: coordinate convention is wrong; every score below would be garbage.")
        return 1
    print("  OK - 1-based inclusive intron coords confirmed.\n")

    annotated = annotated_introns(CHR22_GTF)
    print(f"GENCODE chr22 introns: {len(annotated):,}\n")

    scored, skipped = {}, 0
    for jid in off:
        _c, s, e, _strand = parse_junction_id(jid)
        try:
            scored[jid] = repeat_embedding_score(ref, s, e, anchor=args.anchor)
        except (ValueError, JunctionTooCloseToContigEnd):
            skipped += 1
    print(f"scored {len(scored):,} junctions (skipped {skipped}: short intron / contig edge)\n")

    ann = [sc.min_hamming for j, sc in scored.items() if j in annotated]
    unann = [sc.min_hamming for j, sc in scored.items() if j not in annotated]
    nh_removed = [sc.min_hamming for j, sc in scored.items() if j in removed_by_nh]
    nh_kept = [sc.min_hamming for j, sc in scored.items() if j not in removed_by_nh]

    print("Q1  Does the score separate annotated from unannotated? (higher = less repeat-like)")
    print(summarize("GENCODE-annotated", ann))
    print(summarize("unannotated", unann))
    print()
    print("Q2  Does it agree with the NH gate? (agreement => gate removed repeat-embedded ones)")
    print(summarize("removed by NH gate", nh_removed))
    print(summarize("kept by NH gate", nh_kept))
    print()

    # --- Falsifier 2: is this a repeat detector at all? ---
    # Independent ground truth: RepeatMasker overlap, from the #919 run.
    cat = EXP919 / f"junction_repeat_categorization.{args.sample}.tsv"
    rm_flag = {}
    with cat.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        di, ai = header.index("donor_in_repeat"), header.index("acceptor_in_repeat")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rm_flag[f[0]] = (f[di] == "True") or (f[ai] == "True")

    in_rm = [sc.min_hamming for j, sc in scored.items() if rm_flag.get(j) is True]
    out_rm = [sc.min_hamming for j, sc in scored.items() if rm_flag.get(j) is False]
    print("Q3  FALSIFIER: is the score measuring repeats at all?")
    print("    Independent truth = RepeatMasker overlap (#919). If these two look the")
    print("    same, the Hamming score is NOT a repeat detector and #1116's premise fails.")
    print(summarize("donor/acceptor IN a repeat", in_rm))
    print(summarize("not in a repeat", out_rm))
    if in_rm and out_rm:
        delta = statistics.median(out_rm) - statistics.median(in_rm)
        print(f"\n    median(non-repeat) - median(repeat) = {delta:+.1f}")
        print("    Expected if the score works: POSITIVE (repeat junctions score lower).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
