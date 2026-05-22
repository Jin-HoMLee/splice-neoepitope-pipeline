"""Regenerate slide figures for the VDJdb panel feature doc (Issue #204).

Reads the frozen chr22 panel outputs from `../data/` and renders:
    allele_coverage.png  — per-allele n_exact_matches (log) + n_in_panel bars,
                           colored by panel_status (ok / low_coverage / empty)

Run from repo root with a Python env that has pandas + matplotlib:
    research/.venv/bin/python docs/features/issue_204_vdjdb_panel/figures/_regenerate_figures.py

Data sources (frozen at PR #464 merge — chr22 test config, sample SRR9143066):
    docs/features/issue_204_vdjdb_panel/data/panel.tsv
    docs/features/issue_204_vdjdb_panel/data/panel_qc.tsv

Re-run after editing if the chr22 outputs change (re-copy from
`results/SRR9143066/tcr_panel/vdjdb/` into `../data/` first).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUT_DIR = Path(__file__).resolve().parent
DATA_DIR = OUT_DIR.parent / "data"

STATUS_COLORS = {
    "ok": "#2a9d8f",
    "low_coverage": "#e9c46a",
    "empty": "#e76f51",
}


def render_allele_coverage(qc_path: Path, out_path: Path) -> None:
    qc = pd.read_csv(qc_path, sep="\t")
    qc = qc.sort_values("allele").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = range(len(qc))
    matches_bars = ax.bar(
        x,
        qc["n_exact_matches"].clip(lower=0.5),  # avoid log(0)
        color=[STATUS_COLORS[s] for s in qc["panel_status"]],
        edgecolor="black",
        linewidth=0.5,
        label="n_exact_matches",
    )
    ax.bar(
        x,
        qc["n_in_panel"],
        color="white",
        edgecolor="black",
        linewidth=0.8,
        hatch="//",
        label="n_in_panel",
    )
    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(qc["allele"], rotation=0)
    ax.set_ylabel("Count (log scale)")
    ax.set_title(
        "VDJdb per-allele coverage — chr22 test sample (SRR9143066)\n"
        "Filled = total VDJdb hits · hatched = placed in panel (cap 10)"
    )
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(axis="y", linestyle=":", alpha=0.5)

    # Annotate panel_status above each bar
    for i, row in qc.iterrows():
        ax.text(
            i,
            max(row["n_exact_matches"], 1) * 1.3,
            row["panel_status"],
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=STATUS_COLORS[row["panel_status"]],
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    render_allele_coverage(
        qc_path=DATA_DIR / "panel_qc.tsv",
        out_path=OUT_DIR / "allele_coverage.png",
    )
