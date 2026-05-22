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

dag.svg
-------
The accompanying `dag.svg` is the real Snakemake rule-graph for the chr22 test
config, regenerated via `snakevision` (NOT this script — it requires the
`snakemake` conda env, not the research venv). Default `rule all` does not pull
`fetch_vdjdb_panel` into the DAG yet (report wiring is Issue #206), so an
explicit panel.tsv target is required:

    conda activate snakemake
    snakemake --configfile config/test_config.yaml --rulegraph --forceall -- \\
        results/patient_001_test/reports/report.html \\
        results/patient_001_test/tcr_panel/vdjdb/panel.tsv \\
        > /tmp/dag.dot
    snakevision --style scale=15 --style node_radius=10 --style edge_stroke_width=3 \\
        -o docs/features/issue_204_vdjdb_panel/figures/dag.svg /tmp/dag.dot

The `--` terminator is required: argparse `nargs="+"` on `--configfile` will
otherwise swallow the target as a (non-existent) configfile (CLAUDE.md gotcha).

After regenerating, manually highlight the `fetch_vdjdb_panel` node — find:

    <circle cx="225.0" cy="315.0" fill="#8cd9d9" id="Nfetch_vdjdb_panel"
            r="10" stroke="white" stroke-width="2.0" />
    <text ... x="285.0" y="315.0">fetch_vdjdb_panel</text>

and replace with:

    <circle cx="225.0" cy="315.0" fill="#8cd9d9" id="Nfetch_vdjdb_panel"
            r="14" stroke="#ff6b35" stroke-width="4" />
    <text ... fill="#ff6b35" font-weight="bold" x="285.0" y="315.0">fetch_vdjdb_panel  ← new (PR #457)</text>

(Coordinates may shift if the DAG structure changes — re-locate by `id=`.)
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
