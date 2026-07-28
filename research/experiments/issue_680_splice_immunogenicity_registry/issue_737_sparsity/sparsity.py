#!/usr/bin/env python3
"""Sparsity analysis of the #680 splice-immunogenicity registry (Issue #1230).

This is the cadence artifact: it must be recomputed every time the registry grows.
It was previously a notebook a human re-executed and then hand-reconciled across four
files, which is how Issue #1069 shipped a stale number that only a bot review caught.

`outputs/sparsity_stats.json` is the single source of truth. Prose readers (the
writeup, the slides, the registry README section) are checked against it by
`check_prose_consistency`, so a number can drift in exactly one place before a test
goes red.

Usage:
    research/.venv/bin/python sparsity.py           # recompute and rewrite stats.json
    research/.venv/bin/python sparsity.py --check   # verify without writing
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
REGISTRY = EXPERIMENT / "registry.tsv"
DECOYS = EXPERIMENT / "decoy_negatives" / "presented_decoys_681.tsv"
OUTPUTS = HERE / "outputs"
STATS_PATH = OUTPUTS / "sparsity_stats.json"
WRITEUP_PATH = HERE / "sparsity_writeup.md"
SLIDES_PATH = HERE / "slides.qmd"
README_PATH = EXPERIMENT / "README.md"

# Rows that cleared both curation gates. "Dual gate" is a tier whitelist; there is
# deliberately no separate `label` test here, matching the published analysis.
GATE2_TIERS = [
    "functional-scorable",
    "functional-nonscorable",
    "candidate-negative",
    "hard-negative-true-splice",
]
SCORABLE_TIER = "functional-scorable"

POWER_AUCS = [0.70, 0.75, 0.80, 0.85]
UNRESOLVED = "(unresolved)"


def canon_study(source):
    """Collapse the two source strings that appear in split forms."""
    source = str(source)
    if source.startswith("SNAF"):
        return "SNAF (Li 2024)"
    if source.startswith("IRIS"):
        return "IRIS (Pan/Xing)"
    return source


def concentration(counts):
    """Concentration record for one facet, from a {category: count} mapping.

    `hhi` is the Herfindahl-Hirschman index over unrounded shares; `effective_n` is
    the inverse-Simpson count (1/HHI) computed from the unrounded index and only then
    rounded, so it is not a rounded value of a rounded value.

    `counts` key order is preserved as given (callers pass descending-count order).
    `top` is resolved with an explicit tie-break on the category name so dict
    iteration order can never decide which contributor is reported as largest.
    """
    if not counts:
        raise ValueError("concentration() needs at least one category")

    total = int(sum(counts.values()))
    if total <= 0:
        raise ValueError(f"concentration() needs a positive total, got {total}")

    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0])))
    top_name, top_count = ordered[0]
    distinct = len(ordered)

    hhi = float(sum((count / total) ** 2 for count in counts.values()))
    top_share = float(top_count / total)
    if distinct >= 2:
        top2_share = float((top_count + ordered[1][1]) / total)
    else:
        top2_share = top_share

    return {
        "n": total,
        "distinct": distinct,
        "top": top_name,
        "top_share": round(top_share, 3),
        "top2_share": round(top2_share, 3),
        "hhi": round(hhi, 3),
        "effective_n": round(1.0 / hhi, 2),
        "counts": dict(counts),
    }


def _ordered_counts(series, facet):
    """Descending-count mapping for a column, with blanks made visible.

    `value_counts` silently drops NaN, which would shrink the denominator below the
    row count and inflate the effective number. Filling first keeps every row in the
    denominator; the assertion is the loud backstop if that ever stops holding.
    """
    filled = series.fillna(UNRESOLVED)
    counts = filled.value_counts()
    if int(counts.sum()) != len(filled):
        raise ValueError(
            f"{facet}: {len(filled) - int(counts.sum())} rows dropped from the denominator"
        )
    return {str(k): int(v) for k, v in counts.to_dict().items()}


def n_per_arm(auc, alpha=0.05, power=0.80, null=0.5):
    """Hanley-McNeil sample size per arm to detect `auc` against the `null` AUC.

    Balanced arms, two-sided alpha. Ceiled once at the end.
    """
    q1 = auc / (2 - auc)
    q2 = 2 * auc * auc / (1 + auc)
    v1 = (q1 - auc**2) + (q2 - auc**2)
    q1_null = null / (2 - null)
    q2_null = 2 * null * null / (1 + null)
    v0 = (q1_null - null**2) + (q2_null - null**2)
    z_alpha = norm.ppf(1 - alpha / 2)
    z_beta = norm.ppf(power)
    numerator = z_alpha * np.sqrt(v0) + z_beta * np.sqrt(v1)
    return int(np.ceil((numerator / (auc - null)) ** 2))


def power_n_per_arm(aucs):
    """Sample-size curve keyed by the AUC's float repr ("0.7", not "0.70")."""
    return {repr(float(auc)): n_per_arm(float(auc)) for auc in aucs}


def compute_stats(registry_path=REGISTRY, decoys_path=DECOYS):
    registry = pd.read_csv(registry_path, sep="\t")
    scorable = registry[registry.tier == SCORABLE_TIER].assign(
        study=lambda d: d.source.map(canon_study)
    )

    stats = {
        "registry_rows_total": int(len(registry)),
        "dual_gate_rows": int(len(registry[registry.tier.isin(GATE2_TIERS)])),
        "scorable_positives": int(len(scorable)),
        "by_study": concentration(_ordered_counts(scorable.study, "study")),
        "by_allele": concentration(_ordered_counts(scorable.hla, "allele")),
        "by_mechanism": concentration(
            _ordered_counts(scorable.splice_mechanism_canonical, "mechanism")
        ),
    }

    negatives = registry[registry.label == "negative"]
    hard = int((negatives.evidence_strength == "hard").sum())
    soft = int((negatives.evidence_strength == "soft").sum())
    # Negatives graded `na` fall in neither bucket, so tested_total is intentionally
    # smaller than the negative row count.
    tier2 = int(len(pd.read_csv(decoys_path, sep="\t")))
    stats["negatives"] = {
        "hard_true_negative": hard,
        "soft_negative": soft,
        "tested_total": hard + soft,
        "tier2_presented_untested": tier2,
        "usable_decoy_ceiling": hard + soft + tier2,
    }

    stats["power_n_per_arm"] = power_n_per_arm(POWER_AUCS)
    return stats


# Every headline number a prose reader may state, and where it comes from.
#
# Each claim carries a LIST of patterns because the three readers word the same fact
# differently ("the effective number of alleles is 1.26" vs "effective ~ 1.26
# alleles"); the first pattern that matches wins. `precision` None means an exact
# integer match; an int rounds the stats value to that many decimals first, so prose
# may legitimately quote 3.9 for a stored 3.93.
def _claim_catalogue(stats):
    def pct(share):
        return round(float(share) * 100)

    return {
        "registry rows total": (
            [r"Of\s+\*{0,2}(\d+)\*{0,2}\s+registry rows"],
            stats["registry_rows_total"], None),
        "dual-gate rows": (
            [r"registry rows,\s+\*{0,2}(\d+)\*{0,2}\s+pass both gates"],
            stats["dual_gate_rows"], None),
        "scorable positives": (
            [r"\*{0,2}(\d+)\s+are scorable", r"(\d+)\s+scorable positives"],
            stats["scorable_positives"], None),
        "effective study count": (
            [r"effective number of independent studies is\s+\*{0,2}(\d+(?:\.\d+)?)",
             r"really worth\s*~\s*\*{0,2}(\d+(?:\.\d+)?)",
             r"[Ee]ffective\s*[~≈]\s*\*{0,2}(\d+(?:\.\d+)?)\*{0,2}\s+independent studies"],
            stats["by_study"]["effective_n"], 1),
        "study HHI": (
            [r"effective number of independent studies is[^(]*\(HHI\s*=\s*(\d+(?:\.\d+)?)\)"],
            stats["by_study"]["hhi"], 3),
        "effective allele count": (
            [r"effective number of alleles is\s+\*{0,2}(\d+(?:\.\d+)?)",
             r"effective\s*[~≈]\s*\*{0,2}(\d+(?:\.\d+)?)\s*alleles"],
            stats["by_allele"]["effective_n"], 2),
        "allele HHI": (
            [r"effective number of alleles is[^(]*\(HHI\s*=\s*(\d+(?:\.\d+)?)\)"],
            stats["by_allele"]["hhi"], 3),
        "effective mechanism count": (
            [r"effective number of mechanisms is\s+\*{0,2}(\d+(?:\.\d+)?)",
             r"[Mm]echanism spread is healthier\s*\(effective\s*[~≈]\s*\*{0,2}(\d+(?:\.\d+)?)"],
            stats["by_mechanism"]["effective_n"], 2),
        "top study share (%)": (
            [r"contributes\s+\d+\s+peptides\s+\((\d+)%\)",
             r"[Tt]op study[^=]*=\s*\*{0,2}(\d+)%"],
            pct(stats["by_study"]["top_share"]), None),
        "top-two study share (%)": (
            [r"top two studies supply\*{0,2}\s+\*{0,2}(\d+)%",
             r"top two[^=]*=\s*\*{0,2}(\d+)%"],
            pct(stats["by_study"]["top2_share"]), None),
        "dominant allele share (%)": (
            [r"\((\d+)%\)\*{0,2}\s+are restricted",
             r"\d+/\d+\s*\((\d+)%\)"],
            pct(stats["by_allele"]["top_share"]), None),
        "hard true-negative count": (
            # The slides spell this one out ("**One** hard true-negative"), so the
            # pattern admits a number word and `_parse_count` normalises it.
            [r"\*{0,2}(\d+|[Oo]ne|[Tt]wo|[Tt]hree|[Zz]ero)\*{0,2}\s+hard true-negative"],
            stats["negatives"]["hard_true_negative"], None),
        "power n at AUC 0.75": (
            [r"\*{0,2}(\d+)\s+negatives to detect AUC\s*=\s*0\.75",
             r"needs\s+\*{0,2}(\d+)-\d+\*{0,2}\s+negatives"],
            stats["power_n_per_arm"]["0.75"], None),
        "power n at AUC 0.70": (
            [r"(\d+)\s+to detect AUC\s*=\s*0\.70",
             r"needs\s+\*{0,2}\d+-(\d+)\*{0,2}\s+negatives"],
            stats["power_n_per_arm"]["0.7"], None),
    }


# Which claims each reader is responsible for. A reader is only held to the numbers
# it actually states, but once listed the claim must be present: deleting the
# sentence must not become a way to make the check pass.
READER_CLAIMS = {
    "sparsity_writeup.md": [
        "registry rows total", "dual-gate rows", "scorable positives",
        "effective study count", "study HHI", "effective allele count", "allele HHI",
        "effective mechanism count", "top study share (%)", "top-two study share (%)",
        "dominant allele share (%)", "power n at AUC 0.75", "power n at AUC 0.70",
    ],
    "slides.qmd": [
        "registry rows total", "dual-gate rows", "scorable positives",
        "effective study count", "effective allele count", "effective mechanism count",
        "top study share (%)", "top-two study share (%)", "dominant allele share (%)",
        "hard true-negative count", "power n at AUC 0.75", "power n at AUC 0.70",
    ],
    "registry README.md": [
        "scorable positives", "effective study count", "effective allele count",
        "effective mechanism count", "top study share (%)",
        "top-two study share (%)", "dominant allele share (%)",
        "hard true-negative count", "power n at AUC 0.75", "power n at AUC 0.70",
    ],
}


NUMBER_WORDS = {"zero": 0, "one": 1, "two": 2, "three": 3}


def _parse_count(text):
    """Integer from prose, which may spell a small number out."""
    key = text.strip().lower()
    if key in NUMBER_WORDS:
        return NUMBER_WORDS[key]
    return int(text)


def check_prose_consistency(text, stats, claim_names=None):
    """Mismatches between numbers asserted in prose and the stats single source.

    Defaults to the writeup's claim set. A claim whose patterns all fail to match is
    reported as missing rather than skipped.
    """
    catalogue = _claim_catalogue(stats)
    if claim_names is None:
        claim_names = READER_CLAIMS["sparsity_writeup.md"]

    mismatches = []
    for name in claim_names:
        patterns, expected, precision = catalogue[name]
        # EVERY occurrence is checked, not just the first. A document that states the
        # same fact twice can go stale in one place only, which is exactly the drift
        # Issue #1069 shipped; stopping at the first match would not see it.
        found = [m.group(1) for p in patterns for m in re.finditer(p, text)]
        if not found:
            mismatches.append(f"{name}: claim not found (expected {expected})")
            continue
        for found_text in found:
            if precision is None:
                ok = _parse_count(found_text) == int(expected)
                shown = expected
            else:
                shown = round(float(expected), precision)
                ok = float(found_text) == shown
            if not ok:
                mismatches.append(
                    f"{name}: prose says {found_text}, stats say {shown}"
                )
    return mismatches


def check_all_readers(stats, writeup=None, slides=None, readme=None):
    """Run the canary over every reader that quotes a headline number."""
    paths = {
        "sparsity_writeup.md": writeup or WRITEUP_PATH,
        "slides.qmd": slides or SLIDES_PATH,
        "registry README.md": readme or README_PATH,
    }
    return {
        reader: check_prose_consistency(
            path.read_text(), stats, READER_CLAIMS[reader]
        )
        for reader, path in paths.items()
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the committed stats and prose without writing anything",
    )
    args = parser.parse_args(argv)

    stats = compute_stats()
    serialised = json.dumps(stats, indent=2, default=str)

    if args.check:
        problems = []
        if STATS_PATH.exists() and json.loads(STATS_PATH.read_text()) != stats:
            problems.append(f"{STATS_PATH.name} is stale; re-run without --check")
        for reader, reader_problems in check_all_readers(stats).items():
            problems.extend(f"{reader}: {p}" for p in reader_problems)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            f"{'FAIL' if problems else 'OK'}: {len(problems)} inconsistency(ies)",
            file=sys.stderr,
        )
        return 1 if problems else 0

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    STATS_PATH.write_text(serialised)
    print(
        f"{stats['scorable_positives']} scorable positives from "
        f"{stats['registry_rows_total']} registry rows -> {STATS_PATH}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
