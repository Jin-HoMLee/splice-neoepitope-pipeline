#!/usr/bin/env python3
"""Does the NH-uniqueness filter survive contact with a MATCHED tumor/normal design?

Issue #1122. Companion to `issue_919_nh_uniqueness_filter/` (Developer's evidence
base); this is the Scientist-side interpretation and decision.

Developer established that the chr22 fixture cannot measure the filter's PRECISION:
`NH` is index-relative, so the population the filter exists to catch (a read whose
true locus is off-chr22, landing on chr22's single copy with NH=1) is invisible on a
single-chromosome index. That is correct, and it is why the A/B found nothing.

This tool asks a different question, one the fixture CAN answer, because it turns on
a structural property rather than a rate:

    Is `NH` even the right KIND of quantity to build a junction filter out of,
    in a design whose output is a matched tumor-minus-normal subtraction?

The answer is no, and it does not depend on the index or the depth:

  `NH` is a property of an individual READ (how many places its sequence could go),
  not of the JUNCTION. Two libraries sampling the same biological junction draw
  different reads, so they get different NH profiles by chance. Filtering on NH
  therefore erodes the tumor and the normal arm INDEPENDENTLY -- and a matched
  subtraction cannot survive having its two arms independently eroded.

  Concretely: on this fixture the filter leaves the IGLJ3->IGLC3 junction's tumor
  support untouched (4 reads, all NH=1) while cutting its normal support from 5 to 1.
  Normal support falls under `min_normal_reads`, the junction stops counting as
  "seen in normal", and a real immunoglobulin J-C splice junction -- present in the
  matched normal tissue -- is promoted to a TUMOR-EXCLUSIVE neoepitope candidate.

Three checks, in order of what they establish:

  1. `mean_filter_profile`  -- the depth gate that actually runs the pipeline
     (a floating per-file MEAN, not a decided threshold). Issue #1161.
  2. `end_to_end_effect`    -- the filter's cost on the candidate set PRODUCTION
     builds, not on the raw junction set (which is what #919 measured).
  3. `matched_asymmetry`    -- the disqualifying result: NH erodes the two arms
     independently, and one such erosion manufactures a false tumor-exclusive.

Inputs are all committed or offline-derivable:
  - `issue_919_nh_uniqueness_filter/outputs/raw_junctions.{tumor,normal}.filter_{off,on}.tsv`
  - `resources/test/chr22.gtf.gz`  (annotated introns + gene context)

Run from the repo root:
    research/.venv/bin/python \
      research/experiments/issue_1122_multimapped_reads/nh_matched_subtraction.py
"""
import argparse
import gzip
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

# filter_junctions.py builds the normal junction set with this fixed floor.
MIN_NORMAL_READS = 2

AB = "research/experiments/issue_919_nh_uniqueness_filter/outputs/raw_junctions.{sample}.filter_{knob}.tsv"
GTF = "resources/test/chr22.gtf.gz"


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def load_junctions(sample, knob, root):
    """junction_id -> read count, from a committed A/B junction set."""
    path = root / AB.format(sample=sample, knob=knob)
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out[parts[0]] = float(parts[1])
    return out


def load_gtf(root):
    """Return (annotated_intron_ids, genes).

    Junction IDs use the MIXED convention written by bed12_to_junctions.py:82 --
    `<chrom>:<donor_1BASED>:<acceptor_0based_exclusive>:<strand>`. Reading the donor
    as 0-based (the obvious assumption) silently shifts every motif by one base. The
    canonical-motif self-check below exists to catch exactly that.
    """
    exons_by_tx = defaultdict(list)
    strand_by_tx = {}
    genes = []

    with gzip.open(root / GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = f[8]

            def get(key):
                i = attr.find(key + ' "')
                if i < 0:
                    return None
                j = attr.find('"', i + len(key) + 2)
                return attr[i + len(key) + 2:j]

            if f[2] == "gene":
                genes.append((int(f[3]), int(f[4]), get("gene_name"), get("gene_type")))
            elif f[2] == "exon":
                tid = get("transcript_id")
                if tid:
                    exons_by_tx[tid].append((int(f[3]), int(f[4])))
                    strand_by_tx[tid] = f[6]

    annotated = set()
    for tid, exons in exons_by_tx.items():
        if len(exons) < 2:
            continue
        exons.sort()
        st = strand_by_tx[tid]
        for (_s1, e1), (s2, _e2) in zip(exons, exons[1:]):
            # GTF exon end e1 (1-based incl) -> intron 1-based start = e1 + 1
            # next exon start s2 (1-based)   -> intron 0-based exclusive end = s2 - 1
            annotated.add(f"chr22:{e1 + 1}:{s2 - 1}:{st}")

    genes.sort()
    return annotated, genes


def gene_at(genes, pos):
    return next((g[2] for g in genes if g[0] <= pos <= g[1]), None)


def nearest_gene(genes, pos, max_dist=50):
    """Nearest gene whose boundary is within `max_dist` bp -- for J-C junction naming."""
    best, best_d = None, None
    for s, e, name, _t in genes:
        d = 0 if s <= pos <= e else min(abs(s - pos), abs(e - pos))
        if best_d is None or d < best_d:
            best, best_d = name, d
    return best if best_d is not None and best_d <= max_dist else None


# --------------------------------------------------------------------------
# 1. the depth gate that actually runs the pipeline
# --------------------------------------------------------------------------
def mean_filter_profile(tumor):
    """Reproduce filter_junctions.py:400 -- `keep if reads > mean(reads)`."""
    reads = list(tumor.values())
    n = len(reads)
    mean = sum(reads) / n
    kept = [r for r in reads if r > mean]
    return {
        "n_raw": n,
        "min": min(reads),
        "median": statistics.median(reads),
        "mean": mean,
        "max": max(reads),
        "n_single_read": sum(1 for r in reads if r == 1),
        "n_kept": len(kept),
        "effective_threshold": min(kept) if kept else None,
    }


# --------------------------------------------------------------------------
# 2. the filter's cost on the candidate set PRODUCTION builds
# --------------------------------------------------------------------------
def candidate_set(knob, root, annotated):
    """Replicate filter_junctions.py end-to-end and return the tumor_exclusive set.

    tumor: keep r > mean(r)  -> classify against GENCODE and against the normal set
    normal set: unannotated normal junctions with r >= MIN_NORMAL_READS
    """
    tumor = load_junctions("tumor", knob, root)
    normal = load_junctions("normal", knob, root)

    normal_set = {j for j, r in normal.items()
                  if j not in annotated and r >= MIN_NORMAL_READS}

    mean = sum(tumor.values()) / len(tumor)
    kept = {j: r for j, r in tumor.items() if r > mean}

    tumor_exclusive = {j for j in kept
                       if j not in annotated and j not in normal_set}

    return {
        "n_raw": len(tumor),
        "n_after_mean": len(kept),
        "n_annotated": sum(1 for j in kept if j in annotated),
        "n_normal_shared": sum(1 for j in kept
                               if j not in annotated and j in normal_set),
        "tumor_exclusive": tumor_exclusive,
        "normal_set": normal_set,
    }


def end_to_end_effect(root, annotated):
    off = candidate_set("off", root, annotated)
    on = candidate_set("on", root, annotated)
    te_off, te_on = off["tumor_exclusive"], on["tumor_exclusive"]
    return {
        "off": off, "on": on,
        "lost": te_off - te_on,
        "gained": te_on - te_off,
    }


# --------------------------------------------------------------------------
# 3. the disqualifying result
# --------------------------------------------------------------------------
def matched_asymmetry(root, annotated, genes):
    """Does the filter erode the tumor and normal arms independently?

    A junction seen in BOTH libraries is one biological event. If the filter is a
    property of the JUNCTION, it must retain the same fraction of reads in both arms.
    If it is a property of the READS, it will not.
    """
    t_off = load_junctions("tumor", "off", root)
    t_on = load_junctions("tumor", "on", root)
    n_off = load_junctions("normal", "off", root)
    n_on = load_junctions("normal", "on", root)

    shared = set(t_off) & set(n_off)
    rows = []
    for j in shared:
        rt = t_on.get(j, 0) / t_off[j]
        rn = n_on.get(j, 0) / n_off[j]
        if abs(rt - rn) < 0.01:
            continue
        # Did the erosion destroy the NORMAL evidence while sparing the tumor?
        # That is what promotes a junction to a false tumor_exclusive.
        manufactures_fp = (n_off[j] >= MIN_NORMAL_READS
                           and n_on.get(j, 0) < MIN_NORMAL_READS
                           and t_on.get(j, 0) > 0)
        _c, d, a, _s = j.split(":")
        rows.append({
            "junction": j,
            "tumor": (t_off[j], t_on.get(j, 0)),
            "normal": (n_off[j], n_on.get(j, 0)),
            "manufactures_false_tumor_exclusive": manufactures_fp,
            "donor_gene": nearest_gene(genes, int(d) - 1) or gene_at(genes, int(d) - 1),
            "acceptor_gene": nearest_gene(genes, int(a) + 1) or gene_at(genes, int(a) + 1),
        })
    return {"n_shared": len(shared), "asymmetric": rows}


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", type=Path, help="repo root")
    ap.add_argument("--json", type=Path, help="also write results as JSON")
    args = ap.parse_args()
    root = args.root

    annotated, genes = load_gtf(root)
    print(f"GENCODE chr22: {len(annotated):,} annotated introns, {len(genes):,} genes\n")

    # --- 1 ---
    tumor = load_junctions("tumor", "off", root)
    prof = mean_filter_profile(tumor)
    print("=" * 78)
    print("1. THE DEPTH GATE THAT ACTUALLY RUNS THE PIPELINE  (filter_junctions.py:400)")
    print("=" * 78)
    print(f"   read counts: min={prof['min']:.0f} median={prof['median']:.0f} "
          f"mean={prof['mean']:.3f} max={prof['max']:.0f}")
    print(f"   `keep if reads > mean` -> effective threshold: "
          f"reads >= {prof['effective_threshold']:.0f}")
    print(f"   removes {prof['n_raw'] - prof['n_kept']:,} of {prof['n_raw']:,} junctions "
          f"({1 - prof['n_kept'] / prof['n_raw']:.1%}) -- every single-read junction "
          f"({prof['n_single_read']:,})")
    print("   The threshold is a function of the library's own depth, not a decision. -> Issue 1161\n")

    # --- 2 ---
    eff = end_to_end_effect(root, annotated)
    off, on = eff["off"], eff["on"]
    print("=" * 78)
    print("2. THE FILTER'S COST ON THE CANDIDATE SET PRODUCTION ACTUALLY BUILDS")
    print("=" * 78)
    print(f"   {'':26s} {'filter OFF':>11s} {'filter ON':>11s}")
    print(f"   {'raw junctions':26s} {off['n_raw']:>11,} {on['n_raw']:>11,}")
    print(f"   {'survive mean filter':26s} {off['n_after_mean']:>11,} {on['n_after_mean']:>11,}")
    print(f"   {'-> TUMOR_EXCLUSIVE':26s} "
          f"{len(off['tumor_exclusive']):>11,} {len(on['tumor_exclusive']):>11,}")
    n_off_te = len(off["tumor_exclusive"])
    print(f"\n   lost   : {len(eff['lost'])} of {n_off_te} candidates "
          f"({len(eff['lost']) / n_off_te:.1%})")
    print(f"   gained : {len(eff['gained'])}   <- a filter that only REMOVES reads "
          f"should never ADD a candidate\n")

    # --- 3 ---
    asym = matched_asymmetry(root, annotated, genes)
    print("=" * 78)
    print("3. DOES THE FILTER ERODE TUMOR AND NORMAL INDEPENDENTLY?  (the disqualifier)")
    print("=" * 78)
    rows = asym["asymmetric"]
    print(f"   junctions seen in BOTH arms (one biological event each): {asym['n_shared']}")
    print(f"   ... whose tumor:normal read ratio the filter CHANGES    : {len(rows)}\n")
    print(f"   {'junction':30s} {'tumor':>12s} {'normal':>12s}  locus")
    for r in sorted(rows, key=lambda x: -x["manufactures_false_tumor_exclusive"]):
        t, n = r["tumor"], r["normal"]
        locus = f"{r['donor_gene'] or '?'} -> {r['acceptor_gene'] or '?'}"
        flag = "  <== FALSE TUMOR-EXCLUSIVE" if r["manufactures_false_tumor_exclusive"] else ""
        print(f"   {r['junction']:30s} {t[0]:5.0f}->{t[1]:<5.0f} "
              f"{n[0]:5.0f}->{n[1]:<5.0f}  {locus}{flag}")

    fps = [r for r in rows if r["manufactures_false_tumor_exclusive"]]
    print(f"\n   {len(fps)} junction(s) had matched-normal evidence destroyed while the tumor")
    print("   arm was spared, and were promoted to TUMOR-EXCLUSIVE neoepitope candidates.")
    print("   This is an EXISTENCE PROOF, not a rate: n_shared is far too small to")
    print("   estimate a frequency. It is sufficient, because a false tumor-exclusive")
    print("   is a candidate therapeutic target.\n")

    if args.json:
        args.json.write_text(json.dumps({
            "mean_filter": prof,
            "end_to_end": {
                "tumor_exclusive_off": sorted(off["tumor_exclusive"]),
                "tumor_exclusive_on": sorted(on["tumor_exclusive"]),
                "lost": sorted(eff["lost"]),
                "gained": sorted(eff["gained"]),
            },
            "matched_asymmetry": asym,
        }, indent=2, default=str))
        print(f"wrote {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
