#!/usr/bin/env python3
"""Regenerate the #1105 registry-arc deck figures from registry.tsv.

Run: research/.venv/bin/python research/experiments/issue_680_splice_immunogenicity_registry/figures/_regenerate_figures.py

The registry is canonical: every registry-derived number on a slide comes from here, so
no such figure or prose number is hand-typed (Sandve et al. 2013, rules 1/7/9). The one
non-registry number is Zhao's 139 candidate antigens (it contributes zero registry rows),
hardcoded below and cited to its Supplementary Table S1. Deterministic - no randomness,
no network. Re-run after any registry.tsv edit and commit the PNGs.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
REGISTRY = HERE.parent / "registry.tsv"

INK = "#1b1b1b"
MUTED = "#9aa0a6"
ACCENT = "#c1443b"
COOL = "#3b6ea5"

plt.rcParams.update({
    "font.size": 13,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "figure.dpi": 160,
    "savefig.bbox": "tight",
})


def _load() -> pd.DataFrame:
    return pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")


def fig_identity_matrix(df: pd.DataFrame) -> None:
    """The 2x2 that explains the whole schema change: which identity columns a row has.

    The bottom-left cell (junction published, sequence withheld) is the shape the old
    peptide-keyed schema could not represent at all. Zhao 2025's Table S1 is the motivating
    case (103 junctions), but it is still unfolded: those rows are predicted-only, so they
    have no legal (label, tier) yet. See PROVENANCE.md and issue 1125.
    """
    has_pep = df.peptide.str.strip() != ""
    has_jun = df.junction_id.str.strip() != ""
    counts = [[int((has_pep & has_jun).sum()), int((has_pep & ~has_jun).sum())],
              [int((~has_pep & has_jun).sum()), int((~has_pep & ~has_jun).sum())]]

    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    for r in range(2):
        for c in range(2):
            n = counts[r][c]
            blocked = (r == 1 and c == 0)
            face = "#fdecea" if blocked else ("#eef2f7" if n else "#f6f6f6")
            edge = ACCENT if blocked else MUTED
            ax.add_patch(plt.Rectangle((c, 1 - r), 1, 1, facecolor=face,
                                       edgecolor=edge, linewidth=2.2 if blocked else 1))
            ax.text(c + .5, 1 - r + .60, str(n), ha="center", va="center",
                    fontsize=30, color=ACCENT if blocked else INK, fontweight="bold")
            if blocked:
                ax.text(c + .5, 1 - r + .27, "used to be impossible to store\n(Zhao 2025: 103 junctions,\npredicted-only, still unfolded)",
                        ha="center", va="center", fontsize=9.5, color=ACCENT)

    ax.set_xticks([.5, 1.5]); ax.set_xticklabels(["junction_id present", "junction_id absent"])
    ax.set_yticks([1.5, .5]); ax.set_yticklabels(["peptide\npresent", "peptide\nabsent"])
    ax.set_xlim(0, 2); ax.set_ylim(0, 2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    fig.savefig(HERE / "fig_identity_matrix.png")
    plt.close(fig)


def fig_allele_skew(df: pd.DataFrame) -> None:
    """Three quarters of the registry is one allele. This is the scoring-skew figure."""
    hla = df.hla.replace("", "(none published)").value_counts()
    fig, ax = plt.subplots(figsize=(7.6, 4.3))
    colors = [ACCENT if a == "HLA-A*02:01" else COOL for a in hla.index]
    ax.barh(range(len(hla)), hla.values, color=colors)
    ax.set_yticks(range(len(hla))); ax.set_yticklabels(hla.index, fontsize=10.5)
    ax.invert_yaxis()
    ax.set_xlabel("registry rows")
    for i, v in enumerate(hla.values):
        ax.text(v + 0.8, i, str(v), va="center", fontsize=10)
    total = int(hla.sum())
    ax.set_title(f"HLA-A*02:01 carries {int(hla.iloc[0])} of {total} rows", loc="left", fontsize=13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(HERE / "fig_allele_skew.png")
    plt.close(fig)


def fig_evidence_balance(df: pd.DataFrame) -> None:
    """Positives are abundant; experimental negatives are the binding constraint."""
    scorable_pos = int(((df.label == "positive") & (df.tier == "functional-scorable")).sum())
    hard_neg = int((df.tier == "hard-negative-true-splice").sum())
    soft_neg = int((df.tier == "candidate-negative").sum())

    fig, ax = plt.subplots(figsize=(8.2, 3.9))
    bars = ["scorable\npositives", "soft negatives\n(donor IVS)", "hard negatives\n(MS-presented)"]
    vals = [scorable_pos, soft_neg, hard_neg]
    ax.bar(bars, vals, color=[COOL, "#e0a458", ACCENT], width=.5)
    ax.tick_params(axis="x", labelsize=11)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.2, str(v), ha="center", fontsize=13, fontweight="bold")
    ax.set_ylabel("registry rows")
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_title(f"{scorable_pos} positives are checked by {soft_neg + hard_neg} experimental negatives",
                 loc="left", fontsize=13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(HERE / "fig_evidence_balance.png")
    plt.close(fig)


def main() -> int:
    if not REGISTRY.exists():
        print(f"registry not found: {REGISTRY}", file=sys.stderr)
        return 1
    df = _load()
    fig_identity_matrix(df)
    fig_allele_skew(df)
    fig_evidence_balance(df)
    print(f"regenerated 3 figures from {len(df)} registry rows -> {HERE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
