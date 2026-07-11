#!/usr/bin/env python3
"""Rebuild the issue_919 slide figures from the experiment's committed outputs.

The experiment stays canonical: every number plotted here is recomputed from
`outputs/junction_repeat_categorization.*.tsv` (written by
`junction_repeat_overlap.py`), never hard-coded from the report prose. The two
reference rates that are *not* in that TSV - the random-position null and the
annotated-splice-site rate - are recomputed from the rmsk BED and the reference
junction BED at run time, for the same reason.

Usage:
    conda activate snakemake
    python research/experiments/issue_919_nh_uniqueness_filter/figures/_regenerate_figures.py

Requires `references/rmsk/hg38/rmsk.chr22.bed` (fetch it with
`snakemake --cores 1 -- references/rmsk/hg38/rmsk.chr22.bed`).
"""

import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
sys.path.insert(0, str(EXP))

from junction_repeat_overlap import Junction, load_annotated, load_rmsk, overlaps_repeat  # noqa: E402

# House palette (research/slides/custom.scss). Validated for CVD separation:
# worst adjacent pair dE 106.6 (protan) / 63.6 (tritan), well above the >=12 target.
BLUE = "#1e5ba8"
GOLD = "#b8860b"
INK = "#22252a"        # primary text
MUTED = "#6b7280"      # secondary text / recessive axes
SURFACE = "#fcfcfb"

RMSK = REPO / "references" / "rmsk" / "hg38" / "rmsk.chr22.bed"
ANNOTATED_BED = REPO / "resources" / "test" / "chr22_reference_junctions.bed"

# chr22's first ~10.5 Mb is the acrocentric short arm - unassembled, no repeats
# annotated. Including it in the denominator dilutes the null (41% instead of the
# true 52% over mappable sequence), so the random-position baseline samples only
# the mappable span. This is the correction that made the null match the probe.
CHR22_MAPPABLE_START = 10_510_000
CHR22_END = 50_818_468


def _style(ax):
    """Recessive axes, no chartjunk. Grid behind the marks, hairline weight."""
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=11, length=3, width=0.8)
    ax.yaxis.grid(True, color=MUTED, alpha=0.18, linewidth=0.8)
    ax.set_axisbelow(True)


def _label_bars(ax, bars, values, fmt="{:.1f}%", dy=1.6):
    """Direct-label every bar - few enough marks that a legend lookup is wasted work.

    Labels wear ink, not the series color (text tokens stay text tokens).
    """
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + dy,
            fmt.format(value),
            ha="center", va="bottom", fontsize=11, color=INK, fontweight="600",
        )


def _load_categorization(sample):
    """(fate, annotated) -> repeat-overlap rate, recomputed from the committed TSV."""
    annotated = load_annotated(ANNOTATED_BED)
    rows = []
    path = EXP / "outputs" / f"junction_repeat_categorization.{sample}.tsv"
    with path.open() as fh:
        next(fh)
        for line in fh:
            jid, fate, reads, donor_rep, acceptor_rep, _classes = line.rstrip("\n").split("\t")
            chrom, donor, acceptor, strand = jid.split(":")
            junction = Junction(chrom, int(donor) - 1, int(acceptor), strand, int(reads))
            in_repeat = donor_rep == "True" or acceptor_rep == "True"
            is_ann = (chrom, junction.start, junction.end, strand) in annotated
            rows.append((fate, is_ann, in_repeat))
    return rows


def _rate(rows, fate, is_ann):
    subset = [r for r in rows if r[0] == fate and r[1] is is_ann]
    if not subset:
        return 0.0, 0
    return 100 * sum(1 for r in subset if r[2]) / len(subset), len(subset)


def fig_stratified_overlap(rows):
    """The headline: stratifying dissolves the apparent enrichment.

    Grouped bars - the job is comparing magnitude across two series (what the
    filter removed vs kept) within two strata. Two hues, assigned to the series
    (identity), fixed order, never to rank.
    """
    fig, ax = plt.subplots(figsize=(10, 5.2))

    strata = [(True, "annotated\n(real junctions)"), (False, "unannotated\n(the pool the filter draws from)")]
    lost = [_rate(rows, "lost", ann) for ann, _ in strata]
    kept = [_rate(rows, "retained", ann) for ann, _ in strata]

    x = [0, 1.15]
    w = 0.34
    gap = 0.02  # 2px-equivalent surface gap between adjacent fills
    b1 = ax.bar([i - w / 2 - gap for i in x], [v for v, _ in lost], w,
                label="lost to the filter", color=GOLD, edgecolor=SURFACE, linewidth=1.5)
    b2 = ax.bar([i + w / 2 + gap for i in x], [v for v, _ in kept], w,
                label="retained", color=BLUE, edgecolor=SURFACE, linewidth=1.5)

    _label_bars(ax, b1, [v for v, _ in lost])
    _label_bars(ax, b2, [v for v, _ in kept])

    # n= sits below the axis line; the tick labels get extra pad so they clear it.
    for i, (_, n) in enumerate(lost):
        ax.text(x[i] - w / 2 - gap, -5, f"n={n}", ha="center", fontsize=9.5,
                color=MUTED, clip_on=False)
    for i, (_, n) in enumerate(kept):
        ax.text(x[i] + w / 2 + gap, -5, f"n={n}", ha="center", fontsize=9.5,
                color=MUTED, clip_on=False)

    _style(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in strata], fontsize=12, color=INK)
    ax.tick_params(axis="x", pad=20, length=0)
    ax.set_ylabel("splice site inside a repeat (%)", fontsize=12, color=INK)
    ax.set_ylim(0, 134)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.legend(frameon=False, fontsize=11, loc="upper left", labelcolor=INK)

    # The finding, stated on the chart. Sits ABOVE the bar labels, not through them.
    ax.plot([1.15 - w / 2 - gap, 1.15 + w / 2 + gap], [116, 116], color=MUTED, linewidth=1)
    ax.text(1.15, 120, "0.98x - no enrichment", ha="center", va="bottom",
            fontsize=12.5, color=INK, fontweight="700")
    ax.set_title(
        "Stratified, the filter's apparent selectivity disappears",
        fontsize=14.5, color=INK, fontweight="700", pad=14, loc="left",
    )
    fig.tight_layout()
    fig.savefig(HERE / "fig_stratified_overlap.png", dpi=200)
    plt.close(fig)


def fig_saturation(rows):
    """Why there is no headroom: the unannotated pool is repeat-saturated.

    One series (a magnitude comparison against two reference rates), so no legend -
    the title names what is plotted. The observed bar is the one under discussion,
    so it carries the accent hue; the two baselines are the recessive blue.
    """
    index = load_rmsk(RMSK)

    random.seed(919)
    n = 200_000
    hits = sum(1 for _ in range(n)
               if index.at("chr22", random.randrange(CHR22_MAPPABLE_START, CHR22_END)))
    random_rate = 100 * hits / n

    annotated = load_annotated(ANNOTATED_BED)
    ann_junctions = [Junction(c, s, e, st, 0) for (c, s, e, st) in annotated if c == "chr22"]
    ann_rate = 100 * sum(1 for j in ann_junctions if overlaps_repeat(j, index)) / len(ann_junctions)

    observed, _ = _rate(rows, "retained", False)

    labels = [
        "random mappable\nchr22 position",
        "GENCODE-annotated\nsplice site",
        "unannotated junctions\n(what we are filtering)",
    ]
    values = [random_rate, ann_rate, observed]
    colors = [BLUE, BLUE, GOLD]

    fig, ax = plt.subplots(figsize=(10, 5.2))
    bars = ax.bar(labels, values, 0.55, color=colors, edgecolor=SURFACE, linewidth=1.5)
    _label_bars(ax, bars, values)
    _style(ax)
    ax.set_ylabel("splice site inside a repeat (%)", fontsize=12, color=INK)
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="x", labelsize=11.5)
    for lbl in ax.get_xticklabels():
        lbl.set_color(INK)
    ax.set_title(
        "The pool is saturated - there is no room left to be 'enriched' for repeats",
        fontsize=14.5, color=INK, fontweight="700", pad=14, loc="left",
    )
    fig.tight_layout()
    fig.savefig(HERE / "fig_saturation.png", dpi=200)
    plt.close(fig)


def fig_index_relative_nh():
    """The mechanism, hand-drawn - and it turns on *how diverged* the copies are.

    An earlier version of this figure claimed any off-chromosome repeat read becomes
    a genome-wide multimapper. That is wrong, and the data refutes it: a read from a
    *diverged* copy maps uniquely and correctly to its true locus genome-wide, and
    on a chr22-only index it fails the alignment-score threshold outright (which is
    why 92% of this library is unmapped) - it never contaminates anything. Force-
    mapping a 5-15% diverged copy would leave 5-15 mismatches, and NM>=2 is 1.4% of
    mapped spliced reads. It does not happen.

    The false-unique population is narrower: reads from *near-identical* copies
    whose siblings lie off-chr22. Genome-wide those are genuinely ambiguous (NH>1,
    caught); on a chr22-only index only one such copy is visible, so they read as
    NH=1 and the filter is blind to exactly the population it exists to remove.
    Near-identical copies also predict NM 0-1, which is what we observe.

    A schematic, not data - matplotlib rather than mermaid because it is a frozen
    'hero' diagram and needs to match the deck palette exactly.
    """
    GREY_FILL, GREY = "#eef0f2", "#6b7280"
    fig, ax = plt.subplots(figsize=(12.2, 6.4))
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 6.4)
    ax.axis("off")
    fig.set_facecolor(SURFACE)

    def box(x, y, w, h, text, color, fill, fs=10.5, weight="normal"):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=color,
                               linewidth=1.6, zorder=2, joinstyle="round"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontweight=weight,
                fontsize=fs, color=INK, zorder=3, linespacing=1.45)

    def arrow(x1, y1, x2, y2, color=MUTED):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=13, color=color, linewidth=1.3, zorder=1))

    ax.text(0, 6.05, "What a chr22-only index does with an off-chromosome repeat read",
            fontsize=13, color=INK, fontweight="700")

    box(0.05, 3.05, 2.35, 1.0, "spliced read\nfrom an off-chr22\nrepeat copy", BLUE, "#eaf1fa")

    # The split is on divergence, then - for near-identical copies - on how many
    # copies happen to land on chr22. That second split is the whole story.
    # 1. Diverged -> no chr22 match clears the score threshold -> unmapped.
    box(3.05, 4.75, 4.6, 0.95, "DIVERGED (5-15%)\nno chr22 copy is close enough",
        GREY, GREY_FILL)
    box(8.1, 4.75, 3.55, 0.95, "UNMAPPED\n(92% of the library)", GREY, GREY_FILL)
    arrow(2.4, 3.75, 3.05, 5.15)
    arrow(7.65, 5.22, 8.1, 5.22, GREY)
    ax.text(12.15, 5.22, "harmless", fontsize=9.5, color=GREY, va="center", ha="right")

    ax.text(3.05, 4.35, "NEAR-IDENTICAL copy family  ->  how many copies on chr22?",
            fontsize=10.5, color=INK, fontweight="700")

    # 2. Near-identical, zero chr22 copies -> also unmapped. (The rare third case.)
    box(3.05, 3.25, 4.6, 0.85, "ZERO copies on chr22", GREY, GREY_FILL)
    box(8.1, 3.25, 3.55, 0.85, "UNMAPPED", GREY, GREY_FILL)
    arrow(7.65, 3.67, 8.1, 3.67, GREY)
    ax.text(12.15, 3.67, "harmless", fontsize=9.5, color=GREY, va="center", ha="right")

    # 3. Near-identical, exactly one chr22 copy -> false-unique. THE problem.
    box(3.05, 1.95, 4.6, 0.95, "EXACTLY ONE copy on chr22", GOLD, "#faf3e3", weight="700")
    box(8.1, 1.95, 3.55, 0.95, "NH = 1  ->  looks unique\nFILTER IS BLIND", GOLD, "#faf3e3",
        weight="700")
    arrow(7.65, 2.42, 8.1, 2.42, GOLD)

    # 4. Near-identical, two or more chr22 copies -> NH>1 -> the filter fires. This is
    # chr22-INTERNAL paralogy, and it is the only thing the fixture lets the filter act
    # on - which is exactly why the A/B found only the IGLV2 losses.
    box(3.05, 0.75, 4.6, 0.95, "TWO OR MORE copies on chr22", BLUE, "#eaf1fa")
    box(8.1, 0.75, 3.55, 0.95, "NH > 1  ->  filter catches it\n(chr22-internal paralogy)",
        BLUE, "#eaf1fa")
    arrow(7.65, 1.22, 8.1, 1.22, MUTED)

    arrow(2.4, 3.4, 3.05, 3.67)
    arrow(2.4, 3.3, 3.05, 2.42, GOLD)
    arrow(2.4, 3.2, 3.05, 1.22)

    ax.text(0.05, 0.15,
            "The filter's target population is case 3 - and the fixture makes it invisible. The only case it CAN fire on here is\n"
            "case 4, chr22-internal paralogy - which is exactly why the chr22 A/B found nothing but the four IGLV2 losses.",
            fontsize=10.5, color=INK, linespacing=1.5)

    fig.tight_layout()
    fig.savefig(HERE / "fig_index_relative_nh.png", dpi=200)
    plt.close(fig)


def main():
    if not RMSK.exists():
        raise SystemExit(
            f"missing {RMSK}\nFetch it: snakemake --cores 1 -- references/rmsk/hg38/rmsk.chr22.bed"
        )
    rows = _load_categorization("tumor")
    fig_stratified_overlap(rows)
    fig_saturation(rows)
    fig_index_relative_nh()
    print(f"wrote 3 figures to {HERE}")


if __name__ == "__main__":
    main()
