"""
Issue #299 — identify NEJ_GNAS / NEJ_RPL22 candidates from Snaptron GTEx hg19.

Strategy: filter Snaptron junctions in each gene region for the molecular
signature described in Kwok et al. (Nature 2025):
- NEJ_GNAS  (chr20, + strand): A3 acceptor shift, loss of 2 nt -> frameshift
- NEJ_RPL22 (chr1, - strand) : A3 acceptor shift, loss of 6 nt -> in-frame

For an A3 (alternative 3' splice site) event, the donor (5' splice site) is
shared with a canonical junction; the acceptor (3' splice site) is shifted.

In Snaptron coords (intron interval, 1-based-end), for a + strand gene:
  donor = start, acceptor = end   -> NEJ acceptor offset shifts `end`
For a - strand gene:
  donor = end,   acceptor = start -> NEJ acceptor offset shifts `start`

The shift sign follows transcription direction. "Loss of N nt" means the
mature transcript is N nt shorter than canonical, i.e. the acceptor moved
deeper into the intron (more intron retained near the exon boundary):
  + strand: end shifts +N (acceptor moves to higher coord, intron becomes longer)
  - strand: start shifts -N (acceptor moves to lower coord, intron becomes longer)
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent

CASES = [
    {"gene": "GNAS",  "strand": "+", "loss_nt": 2, "tsv": HERE / "snaptron_GNAS_gtex_hg19.tsv"},
    {"gene": "RPL22", "strand": "-", "loss_nt": 6, "tsv": HERE / "snaptron_RPL22_gtex_hg19.tsv"},
]

# Snaptron GTEx hg19 cohort size (samples in the dataset)
N_SAMPLES_GTEX_HG19 = 9662  # canonical figure from Snaptron docs / Kwok used n=9166

def load_junctions(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def is_donor_annotated(j, strand):
    """For +strand: 5'ss is the 'left' (start) end. For -strand: 5'ss is the 'right' (end) end."""
    return (j["left_annotated"] if strand == "+" else j["right_annotated"]) != "0"

def is_acceptor_annotated(j, strand):
    """For +strand: 3'ss is the 'right' (end) end. For -strand: 3'ss is the 'left' (start) end."""
    return (j["right_annotated"] if strand == "+" else j["left_annotated"]) != "0"

def find_a3_candidates(junctions, strand, loss_nt):
    """
    Find junctions where:
      - donor (5'ss) is annotated in at least one DB
      - acceptor (3'ss) is shifted by exactly `loss_nt` from a canonical acceptor at the same donor
    Non-annotated overall (annotated=0) — i.e., this exact donor-acceptor pairing is novel.
    """
    # Build map: donor position -> set of canonical acceptor positions
    # Canonical = both endpoints annotated (fully annotated junction in some DB OR endpoints both annotated)
    # Use the more permissive criterion: both endpoints annotated.
    canonical_acceptors_by_donor = defaultdict(set)
    for j in junctions:
        if j["strand"] != strand:
            continue
        if is_donor_annotated(j, strand) and is_acceptor_annotated(j, strand):
            donor = int(j["start"]) if strand == "+" else int(j["end"])
            acceptor = int(j["end"]) if strand == "+" else int(j["start"])
            canonical_acceptors_by_donor[donor].add(acceptor)

    candidates = []
    for j in junctions:
        if j["strand"] != strand:
            continue
        if j["annotated"] != "0":
            continue  # exclude already-fully-annotated junctions (they're canonical isoforms)
        if not is_donor_annotated(j, strand):
            continue  # NEJ requires at least the donor to be canonical
        donor = int(j["start"]) if strand == "+" else int(j["end"])
        canonical_acceptors = canonical_acceptors_by_donor.get(donor, set())
        if not canonical_acceptors:
            continue
        sign = +1 if strand == "+" else -1
        expected_offset = sign * loss_nt
        acceptor = int(j["end"]) if strand == "+" else int(j["start"])
        matching = sorted(a for a in canonical_acceptors if (acceptor - a) == expected_offset)
        if matching:
            candidates.append({
                "junc": j,
                "donor": donor,
                "nej_acceptor": acceptor,
                "canonical_acceptors": matching,
            })
    return candidates

def parse_samples(samples_str):
    """Parse Snaptron samples field ',id1:count1,id2:count2,...' -> dict."""
    if not samples_str or samples_str == ",":
        return {}
    parts = [p for p in samples_str.split(",") if p]
    out = {}
    for p in parts:
        sid, _, cnt = p.partition(":")
        if cnt:
            out[sid] = int(cnt)
    return out

def compute_psr_gtex(nej_samples, canonical_samples_per_donor, threshold=0.01):
    """
    Per Kwok et al.: PSR = fraction of cohort samples in which the NJ has
    read frequency >=1% relative to the canonical junction at the same donor.

    For samples where the canonical has 0 reads, the NEJ ratio is undefined;
    Kwok's pipeline counts those samples as 'expressing the NEJ' if NEJ>0
    (lenient interpretation, but exact semantics depend on their script).

    Returns: (n_psr_pass, n_total_samples_with_nej, ratios_per_sample)
    """
    n_psr_pass = 0
    n_with_nej = 0
    ratios = []
    for sid, nej_cnt in nej_samples.items():
        if nej_cnt < 1:
            continue
        n_with_nej += 1
        canon_cnt = canonical_samples_per_donor.get(sid, 0)
        if canon_cnt == 0:
            # Undefined ratio; treat as passing (NEJ visible without canonical)
            ratios.append(("inf", nej_cnt, 0))
            n_psr_pass += 1
        else:
            ratio = nej_cnt / canon_cnt
            ratios.append((ratio, nej_cnt, canon_cnt))
            if ratio >= threshold:
                n_psr_pass += 1
    return n_psr_pass, n_with_nej, ratios

def main():
    print(f"# Issue #299 — NEJ candidate identification + PSR_GTEx (hg19, n~{N_SAMPLES_GTEX_HG19} GTEx samples)")
    print()

    for case in CASES:
        gene, strand, loss_nt, tsv = case["gene"], case["strand"], case["loss_nt"], case["tsv"]
        junctions = load_junctions(tsv)
        print(f"## {gene} ({strand} strand, A3 loss of {loss_nt} nt)")
        print(f"  Total junctions in region: {len(junctions)}")
        fully_annot = sum(1 for j in junctions if j["annotated"] != "0")
        both_ends_annot = sum(1 for j in junctions if is_donor_annotated(j, strand) and is_acceptor_annotated(j, strand))
        print(f"  Fully-annotated junctions (annotated != 0): {fully_annot}")
        print(f"  Both-endpoints-annotated junctions (canonical pool): {both_ends_annot}")
        print(f"  Non-annotated (annotated == 0, i.e. novel pairing): {len(junctions) - fully_annot}")

        cands = find_a3_candidates(junctions, strand, loss_nt)
        print(f"  A3 candidates (non-annotated, donor shared with annotated, acceptor shifted by {loss_nt} nt): {len(cands)}")

        # For each candidate, compute PSR_GTEx
        # Build per-donor canonical sample-counts: aggregate all annotated junctions sharing the donor
        annot_samples_by_donor = defaultdict(lambda: defaultdict(int))
        for j in junctions:
            if j["strand"] != strand:
                continue
            if not (is_donor_annotated(j, strand) and is_acceptor_annotated(j, strand)):
                continue
            donor = int(j["start"]) if strand == "+" else int(j["end"])
            for sid, cnt in parse_samples(j["samples"]).items():
                annot_samples_by_donor[donor][sid] += cnt

        for c in cands:
            j = c["junc"]
            donor = c["donor"]
            nej_samples = parse_samples(j["samples"])
            canon_samples = annot_samples_by_donor.get(donor, {})
            n_pass, n_with, ratios = compute_psr_gtex(nej_samples, canon_samples)
            psr = n_pass / N_SAMPLES_GTEX_HG19 * 100
            our_filter = "REJECT" if n_with >= 1 else "PASS"
            print(f"\n  Candidate: chr={j['chromosome']} strand={j['strand']} {j['start']}-{j['end']} "
                  f"len={j['length']} motif={j['left_motif']}-{j['right_motif']}")
            print(f"    snaptron_id={j['snaptron_id']}  samples_count={j['samples_count']}  coverage_sum={j['coverage_sum']}")
            print(f"    Donor (shared): {donor}")
            print(f"    Canonical acceptor(s) at this donor: {c['canonical_acceptors']}")
            print(f"    NEJ acceptor: {c['nej_acceptor']}  (shifted by {loss_nt} nt)")
            print(f"    Samples with NEJ read >=1: {n_with}")
            print(f"    Samples with NEJ/canonical ratio >=1%: {n_pass}")
            print(f"    PSR_GTEx (Kwok-style, threshold 1%): {psr:.4f}%  (= {n_pass}/{N_SAMPLES_GTEX_HG19})")
            print(f"    Our min_read_count: 1 filter verdict: {our_filter}")

if __name__ == "__main__":
    main()
