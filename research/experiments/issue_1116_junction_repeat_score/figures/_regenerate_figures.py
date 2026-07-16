"""Deck-only figures for the Issue #1116 slide deck.

Per `docs/research_artifact_conventions.md`: data-derived figures come from the
experiment's own script (`make_figures.py`, canonical); *deck-only* explanatory
graphics live here.

This file makes ONE figure: the schematic that explains what the score actually
measures. The four-panel results figure assumes you already know what
"repeat-embedded" means - which is the one thing a reader does not know.

Both examples are REAL rows from the chr22 tumor run, not invented sequence.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT = Path(__file__).resolve().parent

BLUE = "#2563EB"    # exon / trustworthy
ORANGE = "#EA580C"  # the suspicious match
GREY = "#9ca3af"    # intron
INK = "#1f2328"
MUTED = "#6b7280"
LIGHT = "#e5e7eb"


def _seq_box(ax, x, y, w, h, text, color, fill, fontsize=10.5):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=color,
                           linewidth=1.6, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, family="monospace", color=INK, zorder=4)


def panel(ax, title, verdict, verdict_color, left_anchor, intron_3p, hamming, explain):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 30)
    ax.axis("off")

    ax.text(0, 28, title, fontsize=12, fontweight="bold", color=INK)

    # The reference: exon | intron | exon
    y = 17
    ax.add_patch(Rectangle((2, y), 22, 3.2, facecolor="#dbeafe", edgecolor=BLUE,
                           linewidth=1.5, zorder=2))
    ax.text(13, y + 1.6, "EXON (left)", ha="center", va="center", fontsize=9, color=BLUE)

    ax.add_patch(Rectangle((24, y + 1.0), 52, 1.2, facecolor=LIGHT, edgecolor=GREY,
                           linewidth=1.2, zorder=2))
    ax.text(50, y + 3.9, "INTRON  (spliced out)", ha="center", fontsize=9, color=MUTED)

    ax.add_patch(Rectangle((76, y), 22, 3.2, facecolor="#dbeafe", edgecolor=BLUE,
                           linewidth=1.5, zorder=2))
    ax.text(87, y + 1.6, "EXON (right)", ha="center", va="center", fontsize=9, color=BLUE)

    # The two windows we compare: the left anchor, and the FAR end of the intron.
    _seq_box(ax, 12, 9.5, 12, 3.4, left_anchor, BLUE, "#dbeafe")
    ax.text(18, 8.2, "left anchor\n(exonic, 10 bp)", ha="center", va="top",
            fontsize=8.5, color=MUTED)

    _seq_box(ax, 64, 9.5, 12, 3.4, intron_3p, verdict_color, "#f3f4f6")
    ax.text(70, 8.2, "last 10 bp\nOF THE INTRON", ha="center", va="top",
            fontsize=8.5, color=MUTED)

    # Dotted leaders from the reference down to the two windows.
    for x0, x1 in ((18, 18), (70, 70)):
        ax.plot([x0, x1], [y, 12.9], color=MUTED, linewidth=0.8, linestyle=":", zorder=1)

    # The comparison.
    arrow = FancyArrowPatch((24.5, 11.2), (63.5, 11.2), arrowstyle="<->",
                            mutation_scale=13, color=verdict_color, linewidth=1.8, zorder=3)
    ax.add_patch(arrow)
    ax.text(44, 12.4, f"Hamming distance = {hamming}", ha="center", fontsize=11,
            fontweight="bold", color=verdict_color)

    ax.text(44, 4.4, verdict, ha="center", fontsize=11.5, fontweight="bold",
            color=verdict_color)
    ax.text(44, 1.2, explain, ha="center", fontsize=9.5, color=MUTED, wrap=True)


def main():
    fig, axes = plt.subplots(2, 1, figsize=(13, 7.4))
    fig.patch.set_facecolor("white")

    panel(
        axes[0],
        "1.  A GENUINE junction: the anchor looks nothing like the intron's far end",
        "HIGH distance  ->  trustworthy",
        BLUE,
        "AAAAAAAAAA", "TTTTTTTTTT", 10,
        "The read could only have been placed here. There is no other position in this region that fits its sequence.",
    )

    panel(
        axes[1],
        "2.  A REPEAT-EMBEDDED junction: the anchor is IDENTICAL to the intron's far end",
        "LOW distance  ->  suspicious",
        ORANGE,
        "TGGGGCATAG", "TGGGGCATAG", 0,
        "Real row from our chr22 tumor run (chr22:10854248:10854704:+). The aligner had NO sequence basis to place the\n"
        "anchor outside the intron rather than inside it - so this 'junction' may be an artifact of a repeat, not a splice event.",
    )

    fig.suptitle(
        "What the score measures: could the aligner have put this read somewhere else?",
        fontsize=13.5, fontweight="bold", color=INK, x=0.012, ha="left", y=0.985,
    )
    fig.text(0.012, 0.935,
             "Compare each exonic ANCHOR against the intron sequence at the OPPOSITE splice site - "
             "the placement the aligner could have shifted into.",
             fontsize=10, color=MUTED, ha="left")

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    path = OUT / "concept_repeat_embedding.png"
    fig.savefig(path, dpi=160, facecolor="white")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
