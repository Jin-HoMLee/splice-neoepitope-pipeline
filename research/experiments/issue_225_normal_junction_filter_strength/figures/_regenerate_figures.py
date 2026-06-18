"""Regenerate the 2 deck-only figures for Issue #225.

Run from repo root:
    conda activate splice-neoepitope-alphagenome
    python research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py

Outputs (next to this script):
    pr_curve.png      — Precision-Recall curve over universe-restricted F1 sweep
    caught_bar.png    — Tumor caught (% of tumor total) per filter source + union

Logic mirrors research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
sections 2(b) and 3-4 verbatim — same data sources, same universe construction.

Re-run after editing the notebook to keep slide figures in sync. Both the notebook
and this script remain reproducible from the same cached parquet.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve

REPO_ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT_DIR = REPO_ROOT / "research" / "experiments" / "issue_225_normal_junction_filter_strength"
OUTPUTS_DIR = EXPERIMENT_DIR / "outputs"
FIGURES_DIR = Path(__file__).resolve().parent
TARGET_CHROM = "chr22"

MATCHED_NORMAL_TSV = REPO_ROOT / "results" / "patient_001_test" / "alignment" / "SRR9143065_test" / "junctions.tsv"
AG_PARQUET = REPO_ROOT / "research" / "experiments" / "issue_224_alphagenome_exp1" / "outputs" / "chr22_stomach_predicted_junctions.parquet"
GENCODE_GTF = REPO_ROOT / "resources" / "test" / "chr22.gtf.gz"


def load_pipeline_junctions(tsv_path: Path) -> pd.DataFrame:
    raw = pd.read_csv(tsv_path, sep="\t", header=None, names=["key", "count"])
    parts = raw["key"].str.split(":", expand=True)
    parts.columns = ["chrom", "donor_1based", "acceptor_0based_excl", "strand"]
    return pd.DataFrame({
        "chrom": parts["chrom"],
        "donor": parts["donor_1based"].astype(int) - 1,
        "acceptor": parts["acceptor_0based_excl"].astype(int),
        "strand": parts["strand"],
        "count": raw["count"].astype(int),
    })


def gencode_introns_chr22(gtf_path: Path) -> set[tuple]:
    transcripts: dict[str, list] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]) - 1, int(f[4]), f[6], f[8]
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            transcripts.setdefault(m.group(1), []).append((chrom, start, end, strand))
    introns: set[tuple] = set()
    for exons in transcripts.values():
        exons.sort(key=lambda e: e[1])
        for a, b in zip(exons, exons[1:]):
            if a[0] != b[0] or a[3] != b[3]:
                continue
            donor, acceptor = a[2], b[1]
            if donor >= acceptor:
                continue
            introns.add((a[0], donor, acceptor, a[3]))
    return introns


def to_key_set(df: pd.DataFrame) -> set[tuple]:
    return set(map(tuple, df[["chrom", "donor", "acceptor", "strand"]].itertuples(index=False, name=None)))


def make_pr_curve() -> None:
    """Plot PR curve from universe-restricted F1 sweep (matches notebook §2(b))."""
    mn = load_pipeline_junctions(MATCHED_NORMAL_TSV)
    mn_set = to_key_set(mn)
    gencode_introns = gencode_introns_chr22(GENCODE_GTF)
    ag = pd.read_parquet(AG_PARQUET)
    score_col = "score" if "score" in ag.columns else next(
        (c for c in ag.columns if "score" in c.lower() or "prob" in c.lower()), None
    )
    assert score_col is not None, f"No score column in AG parquet; cols: {ag.columns.tolist()}"
    ag_dedup = ag.groupby(["chrom", "donor", "acceptor", "strand"], as_index=False)[score_col].max()
    ag_lookup = dict(zip(
        map(tuple, ag_dedup[["chrom", "donor", "acceptor", "strand"]].itertuples(index=False, name=None)),
        ag_dedup[score_col].to_numpy(),
    ))
    universe_keys = list(gencode_introns)
    scores = np.array([ag_lookup.get(k, 0.0) for k in universe_keys], dtype=float)
    labels = np.array([1 if k in mn_set else 0 for k in universe_keys], dtype=np.int8)

    p, r, t = precision_recall_curve(labels, scores)
    denom = p + r
    f1 = np.where(denom > 0, 2 * p * r / np.where(denom > 0, denom, 1.0), 0.0)
    best_idx = int(np.argmax(f1[:-1]))
    ap = float(average_precision_score(labels, scores))
    baseline = labels.sum() / len(labels)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(r, p, lw=2, color="#1e5ba8", label=f"AG vs MN ∩ GENCODE  (AP = {ap:.3f})")
    ax.axhline(baseline, ls="--", color="grey", alpha=0.6, label=f"Baseline (prevalence = {baseline:.3f})")
    ax.plot(r[best_idx], p[best_idx], "o", color="crimson", markersize=10,
            label=f"F1-max = {f1[best_idx]:.3f} at τ = {t[best_idx]:.3f}\n(P = {p[best_idx]:.3f}, R = {r[best_idx]:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("AlphaGenome PR curve — chr22 universe-restricted F1 sweep\n(universe = GENCODE chr22 introns; positives = MN ∩ GENCODE)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    out = FIGURES_DIR / "pr_curve.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}  (F1={f1[best_idx]:.4f}, τ={t[best_idx]:.4f}, AP={ap:.4f})")


def make_caught_bar() -> None:
    """Bar chart of tumor caught per filter source + union (matches notebook §3-§4)."""
    df = pd.read_csv(OUTPUTS_DIR / "filter_overlap_table.tsv", sep="\t")
    labels = ["MN", "GTEx", "AG", "Caught by any (union)"]
    rows = df[df["region"].isin(labels)].set_index("region").loc[labels].reset_index()

    colors = ["#1e5ba8", "#2e8b57", "#b8860b", "#666666"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(rows["region"], rows["pct_of_tumor"], color=colors)
    for bar, n, pct in zip(bars, rows["n"], rows["pct_of_tumor"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                f"{int(n):,}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=11)
    ax.set_ylabel("% of tumor junctions caught")
    ax.set_title("Tumor junctions caught by each normal-filter source\n(patient_001 chr22; tumor n = 1,872)")
    ax.set_ylim(0, max(rows["pct_of_tumor"]) * 1.18)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "caught_bar.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def main() -> None:
    print(f"Repo root:      {REPO_ROOT}")
    print(f"Experiment dir: {EXPERIMENT_DIR.relative_to(REPO_ROOT)}")
    print(f"Figures dir:    {FIGURES_DIR.relative_to(REPO_ROOT)}")
    make_pr_curve()
    make_caught_bar()
    print("Done.")


if __name__ == "__main__":
    main()
