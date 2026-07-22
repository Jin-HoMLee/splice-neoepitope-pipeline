#!/usr/bin/env python3
"""Regenerate the slide figures for the issue #1204 ORF-FASTA emitter deck.

Three figures, all offline-regenerable from committed inputs:
  1. dag_ms_search_db.png  - the REAL Snakemake rulegraph (_rulegraph.dot, cached
     from `snakemake --rulegraph`), recolored so the two new rules stand out.
  2. emitter_schematic.png - the ORF-stretch extraction concept (matplotlib).
  3. recovery_spectrum.png - the annotated b/y spectrum of the recovered junction
     peptide (matplotlib), masses reused from make_spectrum.py.

Run: python research/experiments/issue_1204_orf_fasta/figures/_regenerate_figures.py
Requires: matplotlib, graphviz `dot` on PATH. Reads outputs/orf_run_summary.json.
"""

import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(EXP))
from make_spectrum import _PROTON, _RESIDUE, _WATER  # noqa: E402

# Project palette (research/slides/custom.scss); blue/gold validated as a
# categorical pair (dataviz validate_palette: CVD dE 28, normal dE 32).
BLUE = "#1e5ba8"
GOLD = "#b8860b"
INK = "#222222"
MUTED = "#8a8a8a"
LIGHT = "#e8eef6"

SUMMARY = json.loads((EXP / "outputs" / "orf_run_summary.json").read_text())


def fig_dag():
    """Recolor the cached real rulegraph: new rules gold, existing muted."""
    dot = (HERE / "_rulegraph.dot").read_text()
    new_rules = ("emit_orf_fasta", "assemble_contigs_wide")
    out_lines = []
    for line in dot.splitlines():
        if "label = " in line and "[label" in line:
            is_new = any(r in line for r in new_rules)
            fill = GOLD if is_new else LIGHT
            fontcolor = "white" if is_new else INK
            # replace the trailing color/style attributes with our styling
            head = line.split("[label")[0]
            label = line.split('label = "')[1].split('"')[0]
            out_lines.append(
                f'{head}[label = "{label}", style="rounded,filled", '
                f'fillcolor="{fill}", color="{BLUE if is_new else MUTED}", '
                f'fontcolor="{fontcolor}", penwidth=2];'
            )
        else:
            out_lines.append(line)
    recolored = HERE / "_rulegraph_colored.dot"
    recolored.write_text("\n".join(out_lines))
    subprocess.run(
        ["dot", "-Tpng", "-Gdpi=150", str(recolored),
         "-o", str(HERE / "dag_ms_search_db.png")],
        check=True,
    )
    recolored.unlink()
    print("wrote dag_ms_search_db.png")


def fig_emitter_schematic():
    """The ORF-stretch extraction concept: 60-nt contig, breakpoint, 3 frames,
    keep the single breakpoint-crossing stop-free stretch."""
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # contig bar (upstream / downstream halves)
    ax.add_patch(Rectangle((0, 8.2), 30, 1.1, facecolor=LIGHT, edgecolor=BLUE, lw=1.5))
    ax.add_patch(Rectangle((30, 8.2), 30, 1.1, facecolor="#f6efe0", edgecolor=GOLD, lw=1.5))
    ax.text(15, 8.75, "upstream flank (30 nt)", ha="center", va="center", fontsize=9, color=INK)
    ax.text(45, 8.75, "downstream flank (30 nt)", ha="center", va="center", fontsize=9, color=INK)
    # breakpoint
    ax.plot([30, 30], [7.6, 9.7], color="#c0392b", lw=2.2, zorder=5)
    ax.text(30, 9.95, "junction breakpoint", ha="center", va="bottom", fontsize=9,
            color="#c0392b", fontweight="bold")

    # three reading frames as codon tracks
    frame_y = [6.2, 4.7, 3.2]
    for fi, y in enumerate(frame_y):
        ax.text(-0.5, y + 0.35, f"frame {fi}", ha="right", va="center", fontsize=9, color=MUTED)
        x = fi
        while x + 3 <= 60:
            ax.add_patch(Rectangle((x, y), 3, 0.7, facecolor="white", edgecolor="#cccccc", lw=0.6))
            x += 3

    # Stop codons (red, "*") are what split each frame into stop-free stretches.
    # A stretch is KEPT only if it is stop-free AND spans the breakpoint
    # (nt_start < 30 < nt_end). Every box below is bounded by the stops drawn.
    STOP = "#c0392b"

    def stop_cell(x, y):
        ax.add_patch(Rectangle((x, y), 3, 0.7, facecolor=STOP, edgecolor="white",
                               lw=0.6, zorder=6))
        ax.text(x + 1.5, y + 0.33, "*", ha="center", va="center", fontsize=13,
                color="white", fontweight="bold", zorder=7)

    # frame 0: stops at 15 & 42 bound a stop-free stretch (18-42) that spans 30 -> KEPT
    stop_cell(15, 6.2)
    stop_cell(42, 6.2)
    ax.add_patch(FancyBboxPatch((18, 6.05), 24, 0.8, boxstyle="round,pad=0.02,rounding_size=0.3",
                                facecolor=BLUE, edgecolor="none", alpha=0.85, zorder=4))
    ax.text(30, 6.45, "kept: crosses breakpoint", ha="center", va="center", fontsize=8.5,
            color="white", fontweight="bold", zorder=5)

    # frame 1: stop at 28 ends the stop-free run before the breakpoint -> upstream-only
    stop_cell(10, 4.7)
    stop_cell(28, 4.7)
    ax.add_patch(FancyBboxPatch((13, 4.55), 15, 0.8, boxstyle="round,pad=0.02,rounding_size=0.3",
                                facecolor="#dcdcdc", edgecolor="none", alpha=0.9, zorder=4))
    ax.text(20.5, 4.95, "✕ upstream-only → dropped", ha="center", va="center",
            fontsize=8, color="#555555", zorder=5)

    # frame 2: stop at 29 starts the stop-free run after the breakpoint -> downstream-only
    stop_cell(29, 3.2)
    stop_cell(47, 3.2)
    ax.add_patch(FancyBboxPatch((32, 3.05), 15, 0.8, boxstyle="round,pad=0.02,rounding_size=0.3",
                                facecolor="#dcdcdc", edgecolor="none", alpha=0.9, zorder=4))
    ax.text(39.5, 3.45, "✕ downstream-only → dropped", ha="center", va="center",
            fontsize=8, color="#555555", zorder=5)

    # legend for the stop-codon cell
    ax.add_patch(Rectangle((0, 7.12), 1.5, 0.6, facecolor=STOP, edgecolor="white", lw=0.6))
    ax.text(0.75, 7.42, "*", ha="center", va="center", fontsize=11, color="white",
            fontweight="bold")
    ax.text(1.9, 7.42, "= stop codon (splits the frame into stop-free stretches)",
            ha="left", va="center", fontsize=8, color=MUTED)

    ax.text(30, 1.5, "keep iff  nt_start < 30 < nt_end   (stop-free, "
            "no X, ≥ 8 aa)  →  FASTA mini-protein",
            ha="center", va="center", fontsize=10, color=INK,
            bbox=dict(boxstyle="round,pad=0.4", facecolor=LIGHT, edgecolor=BLUE, lw=1))
    fig.tight_layout()
    fig.savefig(HERE / "emitter_schematic.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote emitter_schematic.png")


def _by_ions(pep):
    b = []
    run = 0.0
    for a in pep[:-1]:
        run += _RESIDUE[a]
        b.append(run + _PROTON)
    y = []
    run = 0.0
    for a in reversed(pep[1:]):
        run += _RESIDUE[a]
        y.append(run + _WATER + _PROTON)
    return b, y


def fig_recovery_spectrum():
    """Annotated b/y fragment spectrum of the recovered junction peptide."""
    pep = SUMMARY["recovered_peptide"]
    b, y = _by_ions(pep)
    fig, ax = plt.subplots(figsize=(10, 4.8))
    for i, mz in enumerate(b):
        ax.vlines(mz, 0, 0.55 + 0.05 * (i % 3), color=BLUE, lw=2)
        ax.text(mz, 0.62 + 0.05 * (i % 3), f"b{i+1}", ha="center", va="bottom",
                fontsize=7.5, color=BLUE)
    for j, mz in enumerate(y):
        ax.vlines(mz, 0, 0.55 + 0.05 * (j % 3), color=GOLD, lw=2)
        ax.text(mz, 0.62 + 0.05 * (j % 3), f"y{j+1}", ha="center", va="bottom",
                fontsize=7.5, color=GOLD)
    ax.set_xlim(0, max(max(b), max(y)) * 1.08)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("fragment m/z", fontsize=10, color=INK)
    ax.set_yticks([])
    for s in ("top", "left", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    # legend via proxy handles (identity not color-alone: also labeled on peaks)
    ax.plot([], [], color=BLUE, lw=2, label="b-ions (N-terminal)")
    ax.plot([], [], color=GOLD, lw=2, label="y-ions (C-terminal)")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title(
        f"Recovered junction peptide  {pep}  (rank 1, hyperscore "
        f"{SUMMARY['recovered_hyperscore']})",
        fontsize=11, color=INK, pad=10,
    )
    fig.tight_layout()
    fig.savefig(HERE / "recovery_spectrum.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote recovery_spectrum.png")


if __name__ == "__main__":
    fig_dag()
    fig_emitter_schematic()
    fig_recovery_spectrum()
    print("all figures regenerated")
