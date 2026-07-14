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

  The SAM spec (hts-specs, SAMtags) defines NH as the "Number of REPORTED alignments
  that contain the QUERY in the current record". It is keyed on the query -- the
  individual READ -- and counts what the ALIGNER CHOSE TO REPORT. So it is not a
  property of the JUNCTION, nor even a clean property of the read's sequence: it is a
  property of one read under one aligner's reporting policy. Two libraries sampling
  the same biological junction draw different reads, so they get different NH profiles
  by chance. Filtering on NH therefore erodes the tumor and the normal arm
  INDEPENDENTLY -- and a matched subtraction cannot survive having its two arms
  independently eroded.

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

Inputs are ALL COMMITTED -- this runs on a fresh clone, no network, no BAM, and no
`prepare_test_data.sh`:
  - `issue_919_nh_uniqueness_filter/outputs/raw_junctions.{tumor,normal}.filter_{off,on}.tsv`
    (Developer's committed A/B junction sets)
  - `inputs/chr22_annotated_introns.tsv` + `inputs/chr22_genes.tsv`
    (derived here from GENCODE chr22; committed for exactly this reason, since
    `resources/test/chr22.gtf.gz` is GITIGNORED and absent from a fresh clone)

Regenerate the two annotation fixtures from the GTF (needs `prepare_test_data.sh`
to have fetched it) with `--regenerate-annotation`.

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
HERE = "research/experiments/issue_1122_multimapped_reads"
INTRONS_TSV = f"{HERE}/inputs/chr22_annotated_introns.tsv"
GENES_TSV = f"{HERE}/inputs/chr22_genes.tsv"

# GITIGNORED (.gitignore: `resources/test/chr22*`) and fetched by
# scripts/prepare_test_data.sh -- so it is ABSENT from a fresh clone. Only
# `--regenerate-annotation` touches it; the committed fixtures above are what the
# analysis actually reads.
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


def load_annotation(root):
    """Return (annotated_intron_ids, genes) from the COMMITTED fixtures."""
    introns_path, genes_path = root / INTRONS_TSV, root / GENES_TSV
    if not introns_path.exists() or not genes_path.exists():
        raise SystemExit(
            f"missing committed annotation fixture ({introns_path} / {genes_path}).\n"
            "Regenerate with --regenerate-annotation (needs scripts/prepare_test_data.sh "
            "to have fetched resources/test/chr22.gtf.gz, which is gitignored)."
        )

    annotated = set()
    with open(introns_path) as fh:
        next(fh)  # header
        for line in fh:
            if line.strip():
                annotated.add(line.strip())

    genes = []
    with open(genes_path) as fh:
        next(fh)  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                genes.append((int(f[0]), int(f[1]), f[2] or None, f[3] or None))
    genes.sort()
    return annotated, genes


def regenerate_annotation(root):
    """Rebuild the committed fixtures from GENCODE chr22. Provenance for load_annotation().

    Junction IDs use the MIXED convention written by bed12_to_junctions.py:82 --
    `<chrom>:<donor_1BASED>:<acceptor_0based_exclusive>:<strand>`. Reading the donor
    as 0-based (the obvious assumption) silently shifts every motif by one base; the
    canonical-motif falsifier (GENCODE introns must come back ~98% GT-AG) is what
    catches it. Do not "simplify" the +1 away.
    """
    exons_by_tx = defaultdict(list)
    strand_by_tx = {}
    genes = []

    gtf_path = root / GTF
    if not gtf_path.exists():
        raise SystemExit(
            f"{gtf_path} not found -- it is gitignored; run scripts/prepare_test_data.sh first."
        )

    with gzip.open(gtf_path, "rt") as fh:
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

    introns_path, genes_path = root / INTRONS_TSV, root / GENES_TSV
    introns_path.parent.mkdir(parents=True, exist_ok=True)
    with open(introns_path, "w") as fh:
        fh.write("junction_id\n")
        for j in sorted(annotated):
            fh.write(j + "\n")
    with open(genes_path, "w") as fh:
        fh.write("start\tend\tgene_name\tgene_type\n")
        for s, e, name, gtype in genes:
            fh.write(f"{s}\t{e}\t{name or ''}\t{gtype or ''}\n")

    print(f"regenerated {introns_path} ({len(annotated):,} introns)")
    print(f"regenerated {genes_path} ({len(genes):,} genes)")
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
    if not n:
        raise SystemExit("empty tumor junction file -- nothing to profile")
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

    if not tumor:
        raise SystemExit(f"empty tumor junction file for knob={knob}")

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
def matched_asymmetry(root, annotated, genes, gained):
    """Does the filter erode the tumor and normal arms independently?

    A junction seen in BOTH libraries is one biological event. If the filter is a
    property of the JUNCTION, it must retain the same fraction of reads in both arms.
    If it is a property of the READS, it will not.

    `gained` is the AUTHORITATIVE set of junctions promoted to tumor_exclusive by the
    filter, computed by end_to_end_effect() through the full production path (mean
    filter, GENCODE discard, normal subtraction). We take it as an argument rather than
    re-deriving a local predicate: a hand-rolled "would this become a false positive?"
    heuristic agrees with the real computation today and could silently drift from it on
    the next fixture. The authoritative set is the one that decides; `normal_destroyed`
    below only explains WHY, and is never the thing being claimed.
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
        # WHY it happens: the erosion destroyed the normal evidence while sparing the
        # tumor. Explanatory only -- `promoted_to_tumor_exclusive` below is the claim.
        normal_destroyed = (n_off[j] >= MIN_NORMAL_READS
                            and n_on.get(j, 0) < MIN_NORMAL_READS
                            and t_on.get(j, 0) > 0)
        _c, d, a, _s = j.split(":")
        rows.append({
            "junction": j,
            "tumor": (t_off[j], t_on.get(j, 0)),
            "normal": (n_off[j], n_on.get(j, 0)),
            # THE CLAIM: membership in the authoritative `gained` set.
            "promoted_to_tumor_exclusive": j in gained,
            # THE EXPLANATION: never asserted on its own.
            "normal_evidence_destroyed": normal_destroyed,
            "donor_gene": nearest_gene(genes, int(d) - 1),
            "acceptor_gene": nearest_gene(genes, int(a) + 1),
        })
    return {"n_shared": len(shared), "asymmetric": rows}


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", type=Path, help="repo root")
    ap.add_argument("--json", type=Path, help="also write results as JSON")
    ap.add_argument("--regenerate-annotation", action="store_true",
                    help="rebuild inputs/chr22_{annotated_introns,genes}.tsv from the "
                         "GENCODE GTF (gitignored; needs scripts/prepare_test_data.sh)")
    args = ap.parse_args()
    root = args.root

    if args.regenerate_annotation:
        annotated, genes = regenerate_annotation(root)
    else:
        annotated, genes = load_annotation(root)
    print(f"GENCODE chr22: {len(annotated):,} annotated introns, {len(genes):,} genes "
          f"(committed fixtures)\n")

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
    asym = matched_asymmetry(root, annotated, genes, eff["gained"])
    print("=" * 78)
    print("3. DOES THE FILTER ERODE TUMOR AND NORMAL INDEPENDENTLY?  (the disqualifier)")
    print("=" * 78)
    rows = asym["asymmetric"]
    print(f"   junctions seen in BOTH arms (one biological event each): {asym['n_shared']}")
    print(f"   ... whose tumor:normal read ratio the filter CHANGES    : {len(rows)}\n")
    print(f"   {'junction':30s} {'tumor':>12s} {'normal':>12s}  locus")
    for r in sorted(rows, key=lambda x: -x["promoted_to_tumor_exclusive"]):
        t, n = r["tumor"], r["normal"]
        locus = f"{r['donor_gene'] or '?'} -> {r['acceptor_gene'] or '?'}"
        flag = "  <== FALSE TUMOR-EXCLUSIVE" if r["promoted_to_tumor_exclusive"] else ""
        print(f"   {r['junction']:30s} {t[0]:5.0f}->{t[1]:<5.0f} "
              f"{n[0]:5.0f}->{n[1]:<5.0f}  {locus}{flag}")

    fps = [r for r in rows if r["promoted_to_tumor_exclusive"]]
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
