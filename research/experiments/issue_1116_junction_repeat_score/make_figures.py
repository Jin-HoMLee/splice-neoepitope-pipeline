"""Figures for the junction repeat-embedding score.

[Issue #1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116).

Four panels, each answering one question:

  A  Does the score separate real from spurious junctions?      (the result)
  B  Is it measuring repeats at all?                            (the FALSIFIER)
  C  Does the NH gate select for this signal?                   (it does not)
  D  What would a threshold actually cost us?                   (the DECISION)

Panel D is the one that matters for the open question, and it is drawn on the
**real production candidate set** (`novel_junctions.tsv`, the `tumor_exclusive`
junctions) rather than a proxy - because "85% of unannotated junctions" is
abstract, and "you would delete N of your 155 actual candidates" is not.

Palette: two categorical hues, fixed order, validated (CVD dE 124.7, all checks
pass). Color carries identity only; every series is also directly labeled, so the
figure survives grayscale printing and colorblind readers.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from evaluate_chr22 import (  # noqa: E402
    CHR22_FA,
    CHR22_GTF,
    EXP919,
    annotated_introns,
    load_chr22,
    parse_junction_id,
    read_junctions,
)
from repeat_score import (  # noqa: E402
    JunctionTooCloseToContigEnd,
    repeat_embedding_score,
)

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "figures"
CANDIDATES = REPO / "results/patient_001_test/junctions/novel_junctions.tsv"

# Validated categorical pair (node scripts/validate_palette.js): all checks PASS.
BLUE = "#2563EB"   # the "real / trustworthy" series
ORANGE = "#EA580C"  # the "suspect / spurious" series
INK = "#1f2328"
MUTED = "#6b7280"
GRID = "#e5e7eb"

MAXH = 10  # min_hamming is bounded by the anchor length


def _style(ax):
    ax.set_facecolor("white")
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9, length=0)


def _dist(values):
    """Percent of the group at each min_hamming value."""
    n = len(values)
    return [100 * sum(1 for v in values if v == h) / n if n else 0 for h in range(MAXH + 1)]


def paired_bars(ax, a_vals, b_vals, a_label, b_label, title, subtitle):
    xs = range(MAXH + 1)
    w = 0.4
    # 2px surface gap between adjacent bars: achieved by the width/offset pair.
    ax.bar([x - w / 2 - 0.01 for x in xs], _dist(a_vals), width=w,
           color=BLUE, zorder=2, label=f"{a_label} (n={len(a_vals):,})")
    ax.bar([x + w / 2 + 0.01 for x in xs], _dist(b_vals), width=w,
           color=ORANGE, zorder=2, label=f"{b_label} (n={len(b_vals):,})")
    # Title sits ABOVE the subtitle: pad clears the subtitle line, which is
    # anchored just above the axes. Rendering the figure is what caught these
    # two colliding - the numbers could never have told me.
    ax.set_title(title, fontsize=11, color=INK, loc="left", fontweight="bold", pad=24)
    ax.text(0, 1.015, subtitle, transform=ax.transAxes, fontsize=8.5, color=MUTED)
    ax.set_xlabel("min_hamming   (LOW = repeat-embedded = suspicious)", fontsize=9, color=MUTED)
    ax.set_ylabel("% of junctions in group", fontsize=9, color=MUTED)
    ax.set_xticks(list(xs))
    # Headroom so the legend never sits on the tallest bar (it did, first render).
    ax.set_ylim(0, max(max(_dist(a_vals)), max(_dist(b_vals))) * 1.32)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper center")
    _style(ax)


def main():
    OUT.mkdir(exist_ok=True)
    ref = load_chr22(CHR22_FA)
    annotated = annotated_introns(CHR22_GTF)
    off = read_junctions(EXP919 / "raw_junctions.tumor.filter_off.tsv")
    on = read_junctions(EXP919 / "raw_junctions.tumor.filter_on.tsv")
    removed = set(off) - set(on)

    scored = {}
    for jid in off:
        _c, s, e, _st = parse_junction_id(jid)
        try:
            scored[jid] = repeat_embedding_score(ref, s, e).min_hamming
        except (ValueError, JunctionTooCloseToContigEnd):
            pass

    ann = [v for j, v in scored.items() if j in annotated]
    unann = [v for j, v in scored.items() if j not in annotated]
    nh_rm = [v for j, v in scored.items() if j in removed]
    nh_keep = [v for j, v in scored.items() if j not in removed]

    # RepeatMasker: the independent method (from the #919 run).
    rm = {}
    with (EXP919 / "junction_repeat_categorization.tumor.tsv").open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        di, ai = hdr.index("donor_in_repeat"), hdr.index("acceptor_in_repeat")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rm[f[0]] = f[di] == "True" or f[ai] == "True"
    in_rm = [v for j, v in scored.items() if rm.get(j) is True]
    out_rm = [v for j, v in scored.items() if rm.get(j) is False]

    # The REAL production candidates.
    cands = []
    with CANDIDATES.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        oi = hdr.index("junction_origin")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[oi] == "tumor_exclusive" and f[0] in scored:
                cands.append(scored[f[0]])
    n_cand = len(cands)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9.5))
    fig.patch.set_facecolor("white")

    paired_bars(
        axes[0][0], ann, unann, "GENCODE-annotated (real)", "unannotated",
        "A. The score separates real junctions from the rest",
        "Annotated junctions are almost never repeat-embedded. Unannotated ones almost always are.",
    )
    paired_bars(
        axes[0][1], out_rm, in_rm, "not in a repeat", "donor/acceptor in a repeat",
        "B. FALSIFIER: it agrees with an independent method",
        "Groups defined by RepeatMasker, not by this score. If these overlapped, the score would be meaningless.",
    )
    paired_bars(
        axes[1][0], nh_keep, nh_rm, "KEPT by the NH gate", "removed by the NH gate",
        "C. The NH gate does not select for this signal",
        "70% of what it KEEPS is repeat-embedded too. A real repeat filter would separate these two.",
    )

    # Panel D - the decision curve, on the real candidate set.
    ax = axes[1][1]
    ths = list(range(0, MAXH + 1))
    real_kept = [100 * sum(1 for v in ann if v > t) / len(ann) for t in ths]
    cand_kept = [100 * sum(1 for v in cands if v > t) / n_cand for t in ths]
    ax.plot(ths, real_kept, color=BLUE, linewidth=2, marker="o", markersize=5,
            zorder=3, label="GENCODE-annotated junctions retained")
    ax.plot(ths, cand_kept, color=ORANGE, linewidth=2, marker="o", markersize=5,
            zorder=3, label=f"production tumor_exclusive candidates retained (n={n_cand})")
    ax.axvline(2, color=MUTED, linewidth=1, linestyle="--", zorder=1)
    ax.text(2.1, 104, "threshold = 2", fontsize=8.5, color=MUTED)
    ax.annotate(
        f"only {cand_kept[2]:.0f}% of candidates survive\n"
        f"({round(cand_kept[2] / 100 * n_cand)} of {n_cand})",
        xy=(2, cand_kept[2]), xytext=(4.0, 30), fontsize=9.5, color=ORANGE, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2),
    )
    ax.annotate(
        f"but {real_kept[2]:.0f}% of known-real\njunctions survive",
        xy=(2, real_kept[2]), xytext=(3.6, 62), fontsize=9.5, color=BLUE, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2),
    )
    ax.set_title("D. What a threshold would actually cost", fontsize=11, color=INK,
                 loc="left", fontweight="bold", pad=24)
    ax.text(0, 1.015, "THE OPEN QUESTION. Candidates are unannotated by definition, so the score cuts "
                      "hardest exactly where we hunt.",
            transform=ax.transAxes, fontsize=8.5, color=MUTED)
    ax.set_xlabel("filter threshold: keep junctions with min_hamming > t", fontsize=9, color=MUTED)
    ax.set_ylabel("% retained", fontsize=9, color=MUTED)
    ax.set_xticks(ths)
    ax.set_ylim(-3, 112)
    # Legend below the axes: every in-plot corner is occupied by a curve or an
    # annotation, and a legend sitting on the data is the collision this figure
    # already made once.
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=1)
    _style(ax)

    fig.suptitle(
        "Junction repeat-embedding score on chr22 (Issue #1116)   |   LOW min_hamming = repeat-embedded = suspicious",
        fontsize=12.5, color=INK, x=0.008, ha="left", fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    path = OUT / "repeat_score_chr22.png"
    fig.savefig(path, dpi=160, facecolor="white")
    print(f"wrote {path}")

    # The numbers behind panel D, so the figure is not the only record.
    print(f"\nDecision table (production tumor_exclusive candidates, n={n_cand}):")
    print("  threshold   candidates kept   known-real kept")
    for t in range(0, 5):
        print(f"  > {t}         {cand_kept[t]:5.1f}%           {real_kept[t]:5.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
