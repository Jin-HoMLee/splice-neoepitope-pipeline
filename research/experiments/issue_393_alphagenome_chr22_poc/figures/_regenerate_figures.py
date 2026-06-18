"""Regenerate the 3 pilot-deck figures from the cached AG parquet + raw GTF/TSV.

Run from repo root:
    conda activate splice-neoepitope-alphagenome
    python research/slides/issue_393_alphagenome_chr22_poc/figures/_regenerate_figures.py

Outputs (next to this script):
    pr_curve.png            — Precision-Recall curve + AP + baseline
    prf1_vs_threshold.png   — Precision / Recall / F1 vs threshold τ
    bootstrap_f1.png        — Bootstrap F1 histogram + 95% CI band

Logic mirrors `research/experiments/issue_224_alphagenome_exp1/notebook.ipynb`
§1–§5 verbatim — same data sources, same universe construction, same AP/bootstrap
semantics (sklearn AP convention; F1-at-fixed-τ positives-only bootstrap).

Re-run after editing the notebook to keep slide figures in sync. Both the notebook
and this script remain reproducible from the same cached parquet.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = Path(__file__).resolve().parent
MATCHED_NORMAL_TSV = REPO_ROOT / "results" / "patient_001_test" / "alignment" / "SRR9143065_test" / "junctions.tsv"
GENCODE_GTF = REPO_ROOT / "resources" / "test" / "chr22.gtf.gz"
AG_PARQUET = REPO_ROOT / "research" / "experiments" / "issue_224_alphagenome_exp1" / "outputs" / "chr22_stomach_predicted_junctions.parquet"
TARGET_CHROM = "chr22"


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


def parse_exons(gtf_path: Path, target_chrom: str = TARGET_CHROM) -> pd.DataFrame:
    rows = []
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon" or fields[0] != target_chrom:
                continue
            start_1, end_1 = int(fields[3]), int(fields[4])
            strand = fields[6]
            tx = None
            for kv in fields[8].split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id "):
                    tx = kv.split('"')[1]
                    break
            if tx is None:
                continue
            rows.append({"chrom": fields[0], "start_0": start_1 - 1, "end_0": end_1, "strand": strand, "transcript_id": tx})
    return pd.DataFrame(rows)


def exons_to_introns(exons: pd.DataFrame) -> pd.DataFrame:
    introns = []
    for (tx_id, chrom, strand), grp in exons.groupby(["transcript_id", "chrom", "strand"]):
        sorted_exons = grp.sort_values("start_0").reset_index(drop=True)
        for i in range(len(sorted_exons) - 1):
            donor = int(sorted_exons.loc[i, "end_0"])
            acceptor = int(sorted_exons.loc[i + 1, "start_0"])
            introns.append({"chrom": chrom, "donor": donor, "acceptor": acceptor, "strand": strand})
    return pd.DataFrame(introns).drop_duplicates().reset_index(drop=True)


def junction_keys(df: pd.DataFrame) -> set[tuple]:
    return set(zip(df["chrom"], df["donor"].astype(int), df["acceptor"].astype(int), df["strand"]))


def pr_f1_universe(uni_df: pd.DataFrame, threshold: float) -> dict:
    pred_pos = uni_df["ag_score"].values >= threshold
    actual_pos = uni_df["label"].values == 1
    tp = int(np.sum(pred_pos & actual_pos))
    fp = int(np.sum(pred_pos & ~actual_pos))
    fn = int(np.sum(~pred_pos & actual_pos))
    tn = int(np.sum(~pred_pos & ~actual_pos))
    p = tp / (tp + fp) if (tp + fp) else 0.0
    r = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"threshold": threshold, "tp": tp, "fp": fp, "fn": fn, "tn": tn, "precision": p, "recall": r, "f1": f1}


def average_precision_numpy(y_true: np.ndarray, scores: np.ndarray):
    order = np.argsort(-scores, kind="stable")
    y_sorted = y_true[order]
    s_sorted = scores[order]
    tp_cum = np.cumsum(y_sorted == 1)
    fp_cum = np.cumsum(y_sorted == 0)
    total_pos = int((y_true == 1).sum())
    last_of_tied = np.r_[np.diff(s_sorted) != 0, True]
    precision = tp_cum[last_of_tied] / (tp_cum[last_of_tied] + fp_cum[last_of_tied])
    recall = tp_cum[last_of_tied] / total_pos
    recall_with_anchor = np.r_[0.0, recall]
    ap = float(np.sum(np.diff(recall_with_anchor) * precision))
    return ap, precision, recall


def main() -> None:
    print(f"Repo root: {REPO_ROOT}")
    print(f"Out dir:   {OUT_DIR}")

    matched_normal = load_pipeline_junctions(MATCHED_NORMAL_TSV)
    exons = parse_exons(GENCODE_GTF)
    annotated = exons_to_introns(exons)
    predicted = pd.read_parquet(AG_PARQUET)

    mn_keys = junction_keys(matched_normal)
    ann_keys_set = junction_keys(annotated)
    ground_truth_keys = mn_keys & ann_keys_set

    ag_score_map = {
        (r.chrom, int(r.donor), int(r.acceptor), r.strand): float(r.score)
        for r in predicted.itertuples(index=False)
    }
    ann_keys_list = list(zip(
        annotated["chrom"], annotated["donor"].astype(int),
        annotated["acceptor"].astype(int), annotated["strand"],
    ))
    universe = annotated.copy()
    universe["ag_score"] = [ag_score_map.get(k, 0.0) for k in ann_keys_list]
    universe["label"] = [int(k in ground_truth_keys) for k in ann_keys_list]

    print(f"Universe={len(universe):,} positives={int(universe['label'].sum())} negatives={int((universe['label']==0).sum())}")

    thresholds = np.unique(universe["ag_score"].values)
    metrics = pd.DataFrame([pr_f1_universe(universe, t) for t in thresholds])
    ap, precision_curve, recall_curve = average_precision_numpy(
        universe["label"].values.astype(int), universe["ag_score"].values,
    )
    best_idx = int(metrics["f1"].idxmax())
    best = metrics.loc[best_idx]
    best_threshold = float(best["threshold"])
    print(f"AP={ap:.4f} best_F1={best['f1']:.4f} at τ={best_threshold:.4f} P={best['precision']:.4f} R={best['recall']:.4f}")

    # Bootstrap CI (1000 iter, F1-at-fixed-τ, positives-only)
    rng = np.random.default_rng(seed=42)
    n_boot = 1000
    pos = universe[universe["label"] == 1].reset_index(drop=True)
    neg = universe[universe["label"] == 0].reset_index(drop=True)
    f1_boot = np.empty(n_boot)
    for b in range(n_boot):
        ps = pos.iloc[rng.integers(0, len(pos), size=len(pos))]
        boot = pd.concat([ps, neg], ignore_index=True)
        f1_boot[b] = pr_f1_universe(boot, best_threshold)["f1"]
    f1_lo, f1_hi = np.percentile(f1_boot, [2.5, 97.5])
    print(f"F1 95% bootstrap CI: [{f1_lo:.4f}, {f1_hi:.4f}]")

    baseline = float(universe["label"].mean())

    # Figure 1 — P-R curve
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(recall_curve, precision_curve, marker=".", linestyle="-", alpha=0.75, color="steelblue")
    ax.axhline(y=baseline, color="grey", linestyle="--", alpha=0.6, label=f"baseline (random) = {baseline:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Precision–Recall curve (AP = {ap:.3f})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "pr_curve.png", dpi=160)
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'pr_curve.png'}")

    # Figure 2 — P/R/F1 vs threshold
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(metrics["threshold"], metrics["precision"], label="Precision", alpha=0.8)
    ax.plot(metrics["threshold"], metrics["recall"], label="Recall", alpha=0.8)
    ax.plot(metrics["threshold"], metrics["f1"], label="F1", linewidth=2.2, color="black")
    ax.axvline(x=best_threshold, color="red", linestyle=":", alpha=0.7, label=f"best τ = {best_threshold:.3f}")
    ax.set_xlabel("AlphaGenome threshold τ")
    ax.set_ylabel("Score")
    ax.set_title("Precision / Recall / F1 vs threshold")
    ax.set_xscale("symlog", linthresh=0.01)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "prf1_vs_threshold.png", dpi=160)
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'prf1_vs_threshold.png'}")

    # Figure 3 — Bootstrap F1
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(f1_boot, bins=40, color="steelblue", alpha=0.75, edgecolor="white")
    ax.axvline(x=best["f1"], color="red", linestyle="-", linewidth=2.2, label=f"point F1 = {best['f1']:.3f}")
    ax.axvline(x=f1_lo, color="black", linestyle="--", alpha=0.7, label=f"95% CI = [{f1_lo:.3f}, {f1_hi:.3f}]")
    ax.axvline(x=f1_hi, color="black", linestyle="--", alpha=0.7)
    ax.set_xlabel("F1 (bootstrap samples)")
    ax.set_ylabel("Count")
    ax.set_title(f"Bootstrap F1 distribution (n={n_boot}, F1-at-fixed-τ, positives-only)")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "bootstrap_f1.png", dpi=160)
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'bootstrap_f1.png'}")


if __name__ == "__main__":
    main()
