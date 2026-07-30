#!/usr/bin/env python3
"""Regenerate the #737 sparsity deck figures from sparsity_stats.json.

Run: research/.venv/bin/python research/experiments/issue_680_splice_immunogenicity_registry/issue_737_sparsity/figures/_regenerate_figures.py

`outputs/sparsity_stats.json` is canonical (Issue #1230): every number drawn here is
read from it, never recomputed from registry.tsv. That is deliberate. The figures render
numbers visually, so they are the one reader `sparsity.py`'s prose canary cannot check
(it reads text, not pixels) - deriving them from the same single source is what
substitutes for that check. A second, independent computation path would reintroduce
exactly the drift Issue #1069 shipped.

Deterministic: no randomness, no network. Re-run after any registry.tsv edit (recompute
stats first with sparsity.py) and commit the PNGs.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: this must run in CI and without a display
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
STATS_PATH = EXPERIMENT / "outputs" / "sparsity_stats.json"

ACCENT, MUTE, WARN = "#2c6fbb", "#9bb8d6", "#c4453a"

plt.rcParams.update({
    "font.size": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 130,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
})


def pareto_series(counts):
    """Descending (labels, values, cumulative-share) for a {category: count} mapping."""
    if not counts:
        raise ValueError("pareto_series() needs at least one category")
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0])))
    labels = [str(k) for k, _ in ordered]
    values = [int(v) for _, v in ordered]
    total = sum(values)
    running = 0
    cumulative = []
    for value in values:
        running += value
        cumulative.append(running / total)
    return labels, values, cumulative


def power_at(stats, auc):
    """Sample size for an AUC, keyed the way JSON stores it.

    The notebook indexed this table with float keys (`need[0.75]`). Round-tripping
    through JSON turns them into strings, so a float lookup raises. Normalising here
    keeps callers from silently depending on which side of the file they are on.
    """
    return stats["power_n_per_arm"][repr(float(auc))]


def _figure_study_pareto(stats, outdir):
    facet = stats["by_study"]
    labels, values, cumulative = pareto_series(facet["counts"])

    fig, ax = plt.subplots(figsize=(9, 4.6))
    x = range(len(values))
    ax.bar(x, values, color=[ACCENT if i < 2 else MUTE for i in x])
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_ylabel("scorable positive peptides")
    ax.set_title(
        f"Positive base = {facet['n']} peptides from {facet['distinct']} studies "
        f"(effective ≈ {facet['effective_n']})",
        fontsize=12,
    )

    ax2 = ax.twinx()
    ax2.plot(x, cumulative, color=WARN, marker="o", lw=2)
    ax2.set_ylabel("cumulative share", color=WARN)
    ax2.set_ylim(0, 1.02)
    ax2.spines.top.set_visible(False)
    ax2.axhline(0.8, ls=":", color=WARN, alpha=0.6)
    ax2.annotate(
        f"top 2 studies = {facet['top2_share']:.0%}",
        (1, cumulative[1]),
        xytext=(2.2, 0.55),
        color=WARN,
        arrowprops=dict(arrowstyle="->", color=WARN),
    )

    path = outdir / "fig_study_pareto.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def _figure_allele_mechanism(stats, outdir):
    panels = [
        (stats["by_allele"], "HLA restriction"),
        (stats["by_mechanism"], "Splice mechanism"),
    ]
    # Wide figure + generous wspace so the right panel's long mechanism labels clear
    # the left panel's bars instead of overlapping them.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.75})
    for ax, (facet, title) in zip(axes, panels):
        labels, values, _ = pareto_series(facet["counts"])
        positions = list(range(len(values)))[::-1]
        ax.barh(positions, values, color=[ACCENT] + [MUTE] * (len(values) - 1))
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=10)
        ax.margins(y=0.04)
        ax.set_title(
            f"{title}\ntop {facet['top_share']:.0%} · effective ≈ {facet['effective_n']}",
            fontsize=11,
        )
        ax.set_xlabel("scorable positives")

    fig.suptitle("The positive set is an A*02:01 monoculture", y=1.04, fontsize=13)
    path = outdir / "fig_allele_mechanism.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def _figure_negative_scarcity(stats, outdir):
    negatives = stats["negatives"]
    hard = negatives["hard_true_negative"]
    tested = negatives["tested_total"]
    ceiling = negatives["usable_decoy_ceiling"]
    need = stats["power_n_per_arm"]

    available = {
        "hard\ntrue-neg": hard,
        "tested\n(hard+soft)": tested,
        "pooled\n+ untested decoys": ceiling,
    }

    fig, ax = plt.subplots(figsize=(9, 4.6))
    x = range(len(available))
    ax.bar(x, list(available.values()), color=[WARN, "#d98a84", "#e8c0bc"], width=0.6)
    for i, value in enumerate(available.values()):
        ax.text(i, value + 0.4, str(value), ha="center", fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(available.keys())
    ax.set_ylabel("negative peptides")

    for auc, n in need.items():
        ax.axhline(n, ls="--", color=ACCENT, alpha=0.7)
        ax.text(
            len(available) - 0.45,
            n + 0.3,
            f"need {n}  (AUC {auc})",
            color=ACCENT,
            fontsize=9,
            ha="right",
        )

    ax.set_title(
        "Negative set is the binding constraint:\n"
        f"{hard} hard true-negative vs {power_at(stats, 0.75)}-{power_at(stats, 0.70)} "
        "needed for a powered AUC probe",
        fontsize=12,
    )
    ax.set_ylim(0, max(need.values()) + 6)

    path = outdir / "fig_negative_scarcity.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def regenerate(stats, outdir):
    """Write all three deck figures from `stats`; returns the paths written."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return [
        _figure_study_pareto(stats, outdir),
        _figure_allele_mechanism(stats, outdir),
        _figure_negative_scarcity(stats, outdir),
    ]


def main(argv=None):
    if not STATS_PATH.exists():
        print(
            f"{STATS_PATH} not found; run sparsity.py first (it is canonical)",
            file=sys.stderr,
        )
        return 1
    stats = json.loads(STATS_PATH.read_text())
    written = regenerate(stats, HERE)
    for path in written:
        print(f"wrote {path.relative_to(EXPERIMENT)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
