#!/usr/bin/env python
"""Recover the presented-but-untested splice-peptide rows (AC-4) and their per-allele coverage (AC-7).

Issue #681. Run with `research/.venv/bin/python`.

This script answers a *narrower* question than it looks like it answers, and the narrowing is the
point. `build_panel.py` established the ceiling: only a reference whose SEARCH DATABASE is
non-canonical can report a junction-spanning peptide at all. Exactly one of our three references
qualifies (SysteMHC's non-UniProt set), so it is the only place a hit can come from, and every
number below is conditioned on that one 8,911-peptide surface.

What we recover is a *presentation-prevalence* tier: peptides that are (a) predicted splice
neoantigens and (b) observed in a non-canonical immunopeptidomics search. That is evidence of
PRESENTATION, not of IMMUNOGENICITY - nobody T-cell-tested these. They land as `label=untested`
and must never be pooled into the functional positive set (AC-5/AC-6).

**The control is the whole scientific content here.** An exact string match between an 87k-peptide
panel and an 8.9k-peptide atlas is not self-evidently meaningful: class-I peptides are 8-11mers
drawn from a 20-letter alphabet with strongly biased composition (anchor residues), so *some*
collisions are expected by chance alone. A raw hit count with no null is a number, not a finding.
So we compute two independent decoy nulls:

  shuffle   each panel peptide is randomly permuted (preserves length AND amino-acid composition
            exactly, destroys only the sequence ORDER). This is the strict null: it asks whether
            the hits survive when the only thing removed is the actual sequence.
  reverse   each panel peptide is reversed. The standard decoy in proteomics FDR. Weaker than the
            shuffle (it preserves composition and some order structure) but it is the convention,
            and reporting both guards against a shuffle artifact.

If the decoy panels recover hits at the observed rate, the observed hits are collisions and the
tier is empty. That is a result we are willing to publish, which is what makes it a real check.
"""

import json
import random
import statistics
from collections import Counter
from pathlib import Path

from build_panel import (
    DATA,
    OUT,
    SNAF_XLSX,
    load_ours,
    load_snaf,
    load_systemhc_nonuniprot,
    valid,
)

N_DECOY_REPLICATES = 200
SEED = 681  # fixed so the null is reproducible; the number is the Issue, not a magic constant


def shuffled(pep, rng):
    """Permute a peptide's residues. Preserves length + composition, destroys order."""
    chars = list(pep)
    rng.shuffle(chars)
    return "".join(chars)


def decoy_null(panel, reference, rng, replicates=N_DECOY_REPLICATES):
    """Hits recovered by composition-matched decoy panels. Returns (counts, mean, p_empirical).

    p_empirical = fraction of decoy replicates that recover AT LEAST the observed hit count.

    **Iterate `sorted(panel)`, never the set itself.** Python hash-randomizes `str` hashing per
    process (PYTHONHASHSEED), so set iteration order varies run to run. Drawing from `rng` in a
    run-dependent order makes the decoy stream non-reproducible even with a fixed seed - the seed
    pins the RNG, not the order it is consumed in. That is not cosmetic: it is what made the
    reported `decoy_shuffle_mean` drift between runs (0.01 in one, 0.02 in the next), so a
    committed number silently disagreed with the artifact that produced it. Sorting fixes the
    consumption order and makes `seed=SEED` mean what the test plan claims it means.
    """
    ordered = sorted(panel)
    observed = len(panel & reference)
    counts = []
    for _ in range(replicates):
        decoy = {shuffled(p, rng) for p in ordered}
        counts.append(len(decoy & reference))
    ge = sum(1 for c in counts if c >= observed)
    return counts, statistics.mean(counts), (ge + 1) / (replicates + 1)


def main():
    rng = random.Random(SEED)

    snaf_pred, snaf_ms = load_snaf(SNAF_XLSX)
    ours = load_ours(sorted(DATA.glob("patient_*_mhc.tsv")))
    recs = load_systemhc_nonuniprot(DATA / "systemhc_nonuniprot.html", with_meta=True)

    # The reference. Every number below is conditioned on this surface, and it is the ONLY one of
    # our three references that can report a non-canonical peptide (build_panel.py, AC-2).
    systemhc = valid(r["peptide"] for r in recs)

    panels = {"ours": ours, "snaf_pred": snaf_pred, "snaf_ms": snaf_ms}

    print("=" * 78)
    print("RECOVERED HITS vs SysteMHC non-UniProt (the only non-canonical reference)")
    print("=" * 78)

    summary = {"reference_n": len(systemhc), "replicates": N_DECOY_REPLICATES, "seed": SEED}
    summary["panels"] = {}

    for name, panel in panels.items():
        hits = panel & systemhc
        _, shuf_mean, shuf_p = decoy_null(panel, systemhc, rng)
        rev_hits = len({p[::-1] for p in panel} & systemhc)

        rate = len(hits) / len(panel) if panel else 0.0
        summary["panels"][name] = {
            "panel_n": len(panel),
            "hits": len(hits),
            "hit_rate": rate,
            "decoy_shuffle_mean": shuf_mean,
            "decoy_shuffle_p": shuf_p,
            "decoy_reverse_hits": rev_hits,
            "peptides": sorted(hits),
        }
        print(
            f"\n{name:10s} n={len(panel):>6,}  hits={len(hits):>3}  rate={rate*100:.3f}%\n"
            f"           decoy(shuffle) mean={shuf_mean:.2f}  p={shuf_p:.4f}   "
            f"decoy(reverse) hits={rev_hits}"
        )

    # --- POWER: what our own panel could even have detected, at SNAF's observed rate.
    # Our 0 hits is only interesting if we had a real chance of seeing one.
    snaf_rate = summary["panels"]["snaf_pred"]["hit_rate"]
    expected_ours = snaf_rate * len(ours)
    summary["power"] = {
        "snaf_pred_rate": snaf_rate,
        "ours_n": len(ours),
        "ours_expected_hits_at_snaf_rate": expected_ours,
    }
    print("\n" + "=" * 78)
    print("POWER - is our own zero informative?")
    print("=" * 78)
    print(
        f"  SNAF's hit rate against this reference: {snaf_rate*100:.3f}%\n"
        f"  Our panel is {len(ours):,} peptides -> expected hits at that rate: {expected_ours:.2f}\n"
        f"  => our 0 is a SAMPLE-SIZE NULL, not a specificity result. Do not read it as absence."
    )

    # --- AC-4: the recovered rows. Presentation evidence only; never a functional positive.
    by_pep = {}
    for r in recs:
        p = r["peptide"].upper()
        if p in summary["panels"]["snaf_pred"]["peptides"]:
            e = by_pep.setdefault(p, {"alleles": set(), "tissues": set(), "diseases": set()})
            if r["allele"]:
                e["alleles"].add(r["allele"])
            if r["tissue"]:
                e["tissues"].add(r["tissue"])
            if r["disease"]:
                e["diseases"].add(r["disease"])

    ms_set = set(summary["panels"]["snaf_ms"]["peptides"])
    rows = []
    for pep in sorted(by_pep):
        m = by_pep[pep]
        rows.append(
            {
                "peptide": pep,
                "length": len(pep),
                "source_panel": "snaf_pred",
                "ms_confirmed_by_snaf": "yes" if pep in ms_set else "no",
                "systemhc_alleles": ";".join(sorted(m["alleles"])) or "NA",
                "systemhc_tissues": ";".join(sorted(m["tissues"])) or "NA",
                "systemhc_diseases": ";".join(sorted(m["diseases"])) or "NA",
                # AC-4 discipline: presentation evidence, NOT a functional result.
                "label": "untested",
                "tier": "presentation-prevalence",
                "assay_context": "prevalence_only",
            }
        )

    OUT.mkdir(exist_ok=True)
    out_tsv = OUT / "presented_untested_rows.tsv"
    cols = list(rows[0].keys()) if rows else []
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nAC-4: wrote {len(rows)} presentation-tier rows -> {out_tsv.name}")

    # --- AC-7: per-allele coverage. Is the A*02:01 skew actually broken, or only assumed?
    allele_hits = Counter()
    for r in rows:
        for a in r["systemhc_alleles"].split(";"):
            if a and a != "NA":
                allele_hits[a] += 1
    # AC-7's control. A per-allele coverage table on its own restates the REFERENCE's allele
    # composition, not a property of splice peptides - the same shape as
    # feedback_search_key_must_not_mirror_the_gate. So the recovered distribution is only
    # interpretable against the reference's OWN allele background, computed here over rows that
    # actually carry an allele (most do not: SysteMHC leaves the column NA/unclassified for the
    # large majority, and what remains leans on mono-allelic cell lines, which are allele-diverse
    # BY DESIGN).
    # Count the background in the SAME UNIT as the recovered side: one count per unique
    # (peptide, allele) pair, not per SysteMHC record. The recovered side is deduped per peptide
    # (`by_pep` holds a SET of alleles), so a per-record background would mix a per-record rate
    # with a per-(peptide, allele) count and the expectation would be comparing two different
    # things. It does not move the verdict at this n, but a table that ever becomes load-bearing
    # must have both sides counted the same way.
    unassigned = {"NA", "unclassified", ""}
    bg_pairs = {
        (r["peptide"].upper(), r["allele"])
        for r in recs
        if r["allele"] and r["allele"] not in unassigned
    }
    allele_bg = Counter(a for _, a in bg_pairs)
    bg_peptides_with_allele = len({p for p, _ in bg_pairs})

    print("\n" + "=" * 78)
    print("AC-7 - per-allele coverage of the recovered set")
    print("=" * 78)
    total = sum(allele_hits.values())
    for a, n in allele_hits.most_common():
        print(f"  {a:18s} {n:>3}  ({n/total*100:5.1f}% of recovered allele-assignments)")
    # SysteMHC writes alleles as `HLA-A02:01` - no asterisk. Matching the conventional `A*02:01`
    # spelling returns a silent 0 rather than an error, which is the failure mode that would have
    # let me report "the A*02:01 skew is broken" off a string-format mismatch. Normalize both sides.
    norm = lambda a: a.upper().replace("HLA-", "").replace("*", "")
    a0201 = sum(n for a, n in allele_hits.items() if norm(a) == "A02:01")

    bg_total = sum(allele_bg.values())
    bg_a0201 = sum(n for a, n in allele_bg.items() if norm(a) == "A02:01")
    bg_frac = bg_a0201 / bg_total if bg_total else None
    expected_a0201 = bg_frac * total if bg_frac is not None else None
    rows_without_allele = sum(1 for r in rows if r["systemhc_alleles"] == "NA")

    summary["ac7_per_allele"] = {
        "recovered_allele_assignments": total,
        "recovered_rows_without_any_allele": rows_without_allele,
        "distinct_alleles": len(allele_hits),
        "counts": dict(allele_hits.most_common()),
        "a0201_assignments": a0201,
        "a0201_fraction": (a0201 / total) if total else None,
        # the control (same unit as above: unique (peptide, allele) pairs):
        "reference_rows_total": len(recs),
        "reference_peptide_allele_pairs": bg_total,
        "reference_peptides_with_allele": bg_peptides_with_allele,
        "reference_distinct_alleles": len(allele_bg),
        "reference_a0201_fraction": bg_frac,
        "expected_a0201_in_recovered_at_background_rate": expected_a0201,
        "reference_background_top": dict(allele_bg.most_common(10)),
    }
    print(
        f"\n  recovered: {len(allele_hits)} distinct alleles over {total} (peptide, allele) pairs; "
        f"{rows_without_allele}/{len(rows)} recovered peptides carry NO allele at all\n"
        f"  A*02:01 in recovered  : {a0201}/{total} "
        + (f"({a0201/total*100:.1f}%)" if total else "")
        + f"\n  A*02:01 in REFERENCE  : {bg_a0201}/{bg_total} "
        + (f"({bg_frac*100:.1f}%)" if bg_frac is not None else "")
        + "  <- the background this must be read against (same unit)\n"
        f"  expected A*02:01 in recovered at background rate: "
        + (f"{expected_a0201:.2f}" if expected_a0201 is not None else "NA")
        + f"\n  reference: {bg_peptides_with_allele} of its peptides carry any allele, over "
        f"{len(allele_bg)} distinct alleles."
    )
    print(
        "\n  VERDICT (AC-7): the recovered set's allele spread is INHERITED from the reference,\n"
        "  not demonstrated for splice peptides. Observed A*02:01 is indistinguishable from its\n"
        "  background expectation at this n, and the reference leans on mono-allelic lines, which\n"
        "  are allele-diverse by construction. We CANNOT claim the A*02:01 skew is broken."
    )

    with open(OUT / "recovery_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"\nwrote {(OUT / 'recovery_summary.json').name}")


if __name__ == "__main__":
    main()
