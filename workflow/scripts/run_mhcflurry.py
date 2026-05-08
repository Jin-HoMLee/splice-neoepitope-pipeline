#!/usr/bin/env python3
"""run_mhcflurry.py — Wrapper to run MHCflurry 2.x on junction-spanning peptides.

MHCflurry is an open-source MHC-I binding predictor that achieves
state-of-the-art performance.  Unlike NetMHCPan, it does not require
academic registration and can be installed via pip.

Reference:
  O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
  of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
  Cell Systems, 11(1), 42-48.e7.

The script:
  1. Reads junction-spanning peptides from the TSV produced by translate_peptides.py.
  2. Runs MHCflurry Class1PresentationPredictor.predict() with the patient's full
     HLA genotype (all alleles at once), yielding one best-allele prediction per
     peptide: ic50_nM, processing_score, presentation_score, presentation_percentile.
  3. Assigns a presentation_class label (strong/weak/non) from presentation_percentile
     using configurable thresholds (default strong ≤ 0.5%, weak ≤ 2.0%).

Allele sources (in priority order):
  1. ``--alleles-tsv`` / ``snakemake.input.alleles_tsv`` — alleles.tsv from
     aggregate_hla_alleles (patient-specific HLA typing via OptiType).
  2. ``--alleles`` / ``snakemake.params.fallback_alleles`` — explicit allele list,
     drawn from config.mhcflurry.fallback_alleles when HLA typing is disabled.

Output TSV columns:
  contig_key  start_nt  peptide  best_allele  ic50_nM
  processing_score  presentation_score  presentation_percentile  presentation_class
  {allele}_presentation_score  {allele}_presentation_percentile  (one pair per allele)
  genotype_presentation_score  n_strong_alleles  best_presentation_percentile

Usage (standalone, explicit alleles):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/mhc_presentation.tsv \\
      --alleles HLA-A*02:01 HLA-B*07:02 HLA-C*07:02

Usage (standalone, alleles from HLA typing):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/mhc_presentation.tsv \\
      --alleles-tsv results/hla_typing/patient_001/alleles.tsv

Usage (Snakemake):
  Called automatically by the ``run_mhcflurry`` rule.
"""

import argparse
import csv
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Allele loading
# ---------------------------------------------------------------------------

def _load_alleles_from_tsv(tsv_path: str) -> list[str]:
    """Load unique alleles from an alleles.tsv produced by aggregate_hla_alleles.

    Columns: locus, allele1, allele2. Returns a deduplicated list in input
    order (homozygous patients have the same allele in allele1 and allele2 —
    deduplication ensures it is predicted only once).
    """
    alleles: list[str] = []
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for col in ("allele1", "allele2"):
                a = (row.get(col) or "").strip()
                if a:
                    alleles.append(a)
    seen: set[str] = set()
    unique: list[str] = []
    for a in alleles:
        if a not in seen:
            seen.add(a)
            unique.append(a)
    return unique


# ---------------------------------------------------------------------------
# MHCflurry prediction
# ---------------------------------------------------------------------------

def _load_mhcflurry_predictor():
    """Load and return the MHCflurry Class1PresentationPredictor (once per run)."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        log.error(
            "MHCflurry is not installed. Install it with: "
            "pip install mhcflurry && mhcflurry-downloads fetch"
        )
        raise
    try:
        predictor = Class1PresentationPredictor.load()
    except Exception as exc:
        log.error(
            "Failed to load MHCflurry models. Run: mhcflurry-downloads fetch\n%s", exc
        )
        raise
    return predictor


def _normalise_allele(predictor, allele: str) -> str:
    """Normalise an allele string to the format MHCflurry expects.

    Falls back gracefully when the allele is not in the supported list.
    """
    normalised = allele.replace("HLA-", "").replace("*", "").replace(":", "")
    if len(normalised) >= 4:
        mhcflurry_allele = f"HLA-{normalised[0]}*{normalised[1:3]}:{normalised[3:5]}"
    else:
        mhcflurry_allele = allele

    supported = predictor.supported_alleles
    if mhcflurry_allele not in supported:
        alt = allele.replace("-", "").replace("*", "").replace(":", "")
        log.warning(
            "Allele %s not directly supported, trying alternatives. "
            "Supported alleles: %d total",
            mhcflurry_allele, len(supported),
        )
        matching = [a for a in supported if alt[:5] in a.replace("*", "").replace(":", "")]
        if matching:
            mhcflurry_allele = matching[0]
            log.info("Using matched allele: %s", mhcflurry_allele)
        else:
            log.warning("No matching allele found, using original: %s", mhcflurry_allele)

    log.info("Using MHCflurry with allele: %s (normalised: %s)", allele, mhcflurry_allele)
    return mhcflurry_allele


def _run_mhcflurry_predictions(
    peptides: list[str],
    alleles: list[str],
    predictor=None,
) -> pd.DataFrame:
    """Run MHCflurry predictions on a list of peptides for a patient HLA genotype.

    Class1PresentationPredictor.predict() takes the full genotype (≤6 alleles)
    and returns one best-allele prediction per peptide.

    Args:
        peptides:  List of peptide sequences.
        alleles:   Patient HLA genotype (≤6 alleles, e.g. A/B/C loci × 2 each).
        predictor: Pre-loaded Class1PresentationPredictor. Loaded fresh if None.

    Returns:
        DataFrame with columns: peptide, affinity, best_allele, processing_score,
        presentation_score, presentation_percentile (plus metadata columns
        peptide_num/sample_name, dropped downstream).
    """
    if predictor is None:
        predictor = _load_mhcflurry_predictor()

    normalised = [_normalise_allele(predictor, a) for a in alleles]

    if len(normalised) > 6:
        raise ValueError(
            f"Class1PresentationPredictor.predict() supports at most 6 alleles "
            f"(got {len(normalised)}). Check alleles.tsv for unexpected loci."
        )

    log.info(
        "Running MHCflurry predictions for %d peptides against %d-allele genotype...",
        len(peptides), len(normalised),
    )

    return predictor.predict(
        peptides=peptides,
        alleles=normalised,
    )


def _compute_per_allele_features(
    peptides: list[str],
    alleles: list[str],
    hla_c_weight: float,
    strong_threshold: float,
) -> pd.DataFrame:
    """Compute per-allele presentation scores and genotype-level features for each peptide.

    Makes one predict([allele]) call per allele using the module-level predictor
    cache, then computes:
      - {allele}_presentation_score       — absolute presentation probability (0–1)
      - {allele}_presentation_percentile  — rank percentile (lower = better)
      - genotype_presentation_score       — 1 − ∏(1 − wᵢ·pᵢ), pᵢ = presentation_score
      - n_strong_alleles                  — allele count with percentile ≤ strong_threshold
      - best_presentation_percentile      — min percentile across all alleles

    HLA-C alleles are weighted by hla_c_weight; HLA-A/B alleles by 1.0.
    """
    per_allele: dict[str, dict[str, float]] = {p: {} for p in peptides}

    for allele in alleles:
        result = _run_mhcflurry_predictions(peptides, [allele], predictor=_worker_predictor)
        pep_to_score = result.set_index("peptide")["presentation_score"].to_dict()
        pep_to_pct = result.set_index("peptide")["presentation_percentile"].to_dict()
        for p in peptides:
            score = pep_to_score.get(p, 0.0)
            if not 0.0 <= score <= 1.0:
                log.warning(
                    "MHCflurry returned presentation_score=%.6f for peptide %s / %s "
                    "— outside the defined [0, 1] range; value will be used as-is.",
                    score, p, allele,
                )
            per_allele[p][f"{allele}_presentation_score"] = score
            per_allele[p][f"{allele}_presentation_percentile"] = pep_to_pct.get(p, 100.0)

    rows = []
    for pep in peptides:
        rec = per_allele[pep]

        product = 1.0
        n_strong = 0
        best_pct = 100.0
        for allele in alleles:
            w = hla_c_weight if allele.startswith("HLA-C") else 1.0
            p = rec.get(f"{allele}_presentation_score", 0.0)
            pct = rec.get(f"{allele}_presentation_percentile", 100.0)
            product *= (1.0 - w * p)
            if pct <= strong_threshold:
                n_strong += 1
            if pct < best_pct:
                best_pct = pct

        row: dict = {"peptide": pep}
        row.update(rec)
        row["genotype_presentation_score"] = round(1.0 - product, 6)
        row["n_strong_alleles"] = n_strong
        row["best_presentation_percentile"] = best_pct
        rows.append(row)

    df = pd.DataFrame(rows)
    log.info(
        "Per-allele features: %d peptides, %d allele(s); "
        "genotype_presentation_score range [%.4f, %.4f]; %d with n_strong_alleles >= 1",
        len(peptides), len(alleles),
        df["genotype_presentation_score"].min(), df["genotype_presentation_score"].max(),
        (df["n_strong_alleles"] >= 1).sum(),
    )
    return df


# ---------------------------------------------------------------------------
# Binder classification
# ---------------------------------------------------------------------------

def classify_by_percentile(
    percentile: float,
    strong_threshold: float = 0.5,
    weak_threshold: float = 2.0,
) -> str:
    """Classify a peptide by percentile rank (lower = better, per allele).

    Boundaries are inclusive: percentile <= strong_threshold → strong,
    percentile <= weak_threshold → weak, otherwise non.
    """
    if percentile <= strong_threshold:
        return "strong"
    if percentile <= weak_threshold:
        return "weak"
    return "non"


# ---------------------------------------------------------------------------
# GPU detection
# ---------------------------------------------------------------------------

def _has_gpu() -> bool:
    """Return True if a CUDA GPU is available and can execute PyTorch kernels.

    Uses PyTorch (MHCflurry's actual inference backend). torch.cuda.is_available()
    returns True even when the GPU's SM version isn't in the PyTorch build's compiled
    architectures, so we run a minimal smoke-test kernel to catch that case early.
    """
    try:
        import torch
        if not torch.cuda.is_available():
            return False
        t = torch.zeros(2, device="cuda")
        torch.nn.functional.relu(t)  # real kernel dispatch — fails on SM mismatch
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Predictor cache and worker helpers
# ---------------------------------------------------------------------------

# Module-level predictor cache: populated by _load_predictor() before predict().
_worker_predictor = None


def _load_predictor() -> None:
    """Load predictor into module-level cache."""
    global _worker_predictor
    _worker_predictor = _load_mhcflurry_predictor()


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------

def _write_zero_stats(stats_output_path: str | Path | None) -> None:
    """Emit a zero-count mhc-affinity stats TSV for empty-input cases.

    Required so the cross-step aggregator (``aggregate_filtering_stats``)
    finds the file even when zero peptides reach the prediction step.
    """
    if stats_output_path is None:
        return
    stats_output_path = Path(stats_output_path)
    stats_output_path.parent.mkdir(parents=True, exist_ok=True)
    with stats_output_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["category", "count"])
        writer.writerow(["strong_presenters", 0])
        writer.writerow(["weak_presenters", 0])
    log.info("MHC-affinity stats (zero-count) written to %s", stats_output_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_tsv: str | Path,
    output_tsv: str | Path,
    alleles: list[str] | None = None,
    alleles_tsv: str | None = None,
    presentation_percentile_strong: float = 0.5,
    presentation_percentile_weak: float = 2.0,
    hla_c_weight: float = 0.5,
    stats_output_path: str | Path | None = None,
) -> None:
    """Run MHCflurry on junction-spanning peptides and write results to TSV.

    Runs Class1PresentationPredictor.predict() with the full patient HLA genotype
    (all alleles at once) for best-allele columns, then makes one per-allele call to
    compute genotype_presentation_score and related breadth features.

    Args:
        peptides_tsv:                  TSV of junction-spanning peptides.
        output_tsv:                    Destination TSV file.
        alleles:                       HLA alleles to predict. Ignored when
                                       alleles_tsv is provided.
        alleles_tsv:                   Path to alleles.tsv from aggregate_hla_alleles.
                                       Takes precedence over alleles.
        presentation_percentile_strong: Percentile threshold for strong presentation_class
                                       and n_strong_alleles counting.
        presentation_percentile_weak:  Percentile threshold for weak presentation_class.
        hla_c_weight:                  Weight for HLA-C alleles in genotype_presentation_score
                                       formula (default 0.5, reflecting ~50% lower surface
                                       density vs HLA-A/B).
        stats_output_path:             Optional destination TSV for the
                                       mhc-affinity funnel slice (Issue #215).
                                       Two columns — ``category, count`` —
                                       feeding the cross-step aggregator.
                                       Categories: ``strong_presenters``,
                                       ``weak_presenters`` (presenter
                                       vocabulary per CLAUDE.md, distinct from
                                       affinity-only "binders").

    Raises:
        ValueError: If neither alleles nor alleles_tsv is provided, or if
                    hla_c_weight is outside [0, 1].
    """
    if not 0.0 <= hla_c_weight <= 1.0:
        raise ValueError(
            f"hla_c_weight must be in [0, 1]; got {hla_c_weight}. "
            "Values outside this range cannot be interpreted as locus weights."
        )

    # Resolve alleles: alleles_tsv > alleles (caller must supply one)
    if alleles_tsv:
        resolved_alleles = _load_alleles_from_tsv(alleles_tsv)
        log.info(
            "Loaded %d alleles from %s: %s",
            len(resolved_alleles), alleles_tsv, resolved_alleles,
        )
        if not resolved_alleles:
            raise ValueError(
                f"alleles_tsv {alleles_tsv!r} contained no valid alleles. "
                "Check that the file has allele1/allele2 columns with non-empty values."
            )
    elif alleles:
        resolved_alleles = alleles
    else:
        raise ValueError("Either alleles or alleles_tsv must be provided.")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Read junction-spanning peptides
    peptides_df = pd.read_csv(peptides_tsv, sep="\t")

    if peptides_df.empty:
        log.warning("No peptides found in %s", peptides_tsv)
        per_allele_cols: list[str] = []
        for a in resolved_alleles:
            per_allele_cols += [f"{a}_presentation_score", f"{a}_presentation_percentile"]
        empty_df = pd.DataFrame(columns=[
            "contig_key", "start_nt", "peptide", "best_allele",
            "ic50_nM", "processing_score", "presentation_score",
            "presentation_percentile", "presentation_class",
        ] + per_allele_cols + [
            "genotype_presentation_score", "n_strong_alleles", "best_presentation_percentile",
        ])
        empty_df.to_csv(output_tsv, sep="\t", index=False)
        _write_zero_stats(stats_output_path)
        return

    unique_peptides = peptides_df["peptide"].unique().tolist()
    log.info(
        "Extracted %d peptides (%d unique) from %s",
        len(peptides_df), len(unique_peptides), peptides_tsv,
    )

    # Step 2: Run genotype-level prediction (single call, all alleles at once).
    if _has_gpu():
        log.info("GPU detected — loading predictor in main process")
    else:
        log.info(
            "No GPU — running %d allele(s) on CPU (all cores available)",
            len(resolved_alleles),
        )
    _load_predictor()

    pred_df = _run_mhcflurry_predictions(
        unique_peptides, resolved_alleles, predictor=_worker_predictor
    )

    # Rename and keep only needed columns; drop metadata (peptide_num, sample_name)
    pred_df = pred_df.rename(columns={"affinity": "ic50_nM"})[
        ["peptide", "best_allele", "ic50_nM", "processing_score",
         "presentation_score", "presentation_percentile"]
    ]
    pred_df["ic50_nM"] = pred_df["ic50_nM"].fillna(float("inf"))
    pred_df["presentation_class"] = pred_df["presentation_percentile"].apply(
        lambda v: classify_by_percentile(v, presentation_percentile_strong, presentation_percentile_weak)
    )

    # Step 3: Per-allele calls to compute genotype_presentation_score and breadth features
    log.info(
        "Running per-allele predictions for genotype features (%d allele(s))...",
        len(resolved_alleles),
    )
    per_allele_df = _compute_per_allele_features(
        unique_peptides, resolved_alleles, hla_c_weight, presentation_percentile_strong
    )
    pred_df = pred_df.merge(per_allele_df, on="peptide", how="left")

    # Step 4: Merge with peptides_df to restore contig_key/start_nt, then write output
    per_allele_extra_cols = [c for c in per_allele_df.columns if c != "peptide"]
    output_cols = (
        ["contig_key", "start_nt", "peptide", "best_allele",
         "ic50_nM", "processing_score", "presentation_score",
         "presentation_percentile", "presentation_class"]
        + per_allele_extra_cols
    )
    df = peptides_df.merge(pred_df, on="peptide", how="left")[output_cols]

    df.to_csv(output_tsv, sep="\t", index=False)
    n_strong = int((df["presentation_class"] == "strong").sum())
    n_weak = int((df["presentation_class"] == "weak").sum())
    n_non = int((df["presentation_class"] == "non").sum())
    log.info(
        "Predictions: %d rows (%d unique peptides, %d allele(s)); "
        "genotype_presentation_score range [%.4f, %.4f]; "
        "%d strong / %d weak / %d non → %s",
        len(df), len(unique_peptides), len(resolved_alleles),
        df["genotype_presentation_score"].min(), df["genotype_presentation_score"].max(),
        n_strong, n_weak, n_non, output_tsv,
    )

    if stats_output_path is not None:
        stats_output_path = Path(stats_output_path)
        stats_output_path.parent.mkdir(parents=True, exist_ok=True)
        with stats_output_path.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["category", "count"])
            writer.writerow(["strong_presenters", n_strong])
            writer.writerow(["weak_presenters", n_weak])
        log.info("MHC-affinity stats written to %s", stats_output_path)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    # Use patient-specific alleles when HLA typing is enabled (alleles_tsv input
    # is only wired in mhcflurry.smk when config.hla.enabled is true).
    alleles_tsv = getattr(snakemake.input, "alleles_tsv", None)  # type: ignore[name-defined]  # noqa: F821
    run_prediction(
        peptides_tsv=snakemake.input.peptides_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.mhc_presentation_tsv,  # type: ignore[name-defined]  # noqa: F821
        alleles=None if alleles_tsv else list(snakemake.params.fallback_alleles),  # type: ignore[name-defined]  # noqa: F821
        alleles_tsv=alleles_tsv,
        presentation_percentile_strong=float(snakemake.params.presentation_percentile_strong),  # type: ignore[name-defined]  # noqa: F821
        presentation_percentile_weak=float(snakemake.params.presentation_percentile_weak),  # type: ignore[name-defined]  # noqa: F821
        hla_c_weight=float(snakemake.params.hla_c_weight),  # type: ignore[name-defined]  # noqa: F821
        stats_output_path=getattr(snakemake.output, "stats", None),  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run MHCflurry epitope prediction on junction-spanning peptides."
    )
    parser.add_argument("--peptides-tsv", required=True, help="Input peptides TSV")
    parser.add_argument("--output", required=True, help="Output predictions TSV")
    parser.add_argument(
        "--alleles", nargs="+",
        help="HLA alleles to predict (one or more). Ignored when --alleles-tsv is given.",
    )
    parser.add_argument(
        "--alleles-tsv", default=None,
        help="Path to alleles.tsv from aggregate_hla_alleles (overrides --alleles).",
    )
    parser.add_argument("--presentation-percentile-strong", type=float, default=0.5)
    parser.add_argument("--presentation-percentile-weak", type=float, default=2.0)
    parser.add_argument(
        "--hla-c-weight", type=float, default=0.5,
        help="Weight for HLA-C alleles in genotype_presentation_score formula (default 0.5).",
    )
    parser.add_argument(
        "--stats-output", default=None,
        help="Optional MHC-affinity stats TSV (Issue #215)",
    )
    args = parser.parse_args()

    if not args.alleles and not args.alleles_tsv:
        parser.error("Provide --alleles or --alleles-tsv.")

    run_prediction(
        peptides_tsv=args.peptides_tsv,
        output_tsv=args.output,
        alleles=args.alleles,
        alleles_tsv=args.alleles_tsv,
        presentation_percentile_strong=args.presentation_percentile_strong,
        presentation_percentile_weak=args.presentation_percentile_weak,
        hla_c_weight=args.hla_c_weight,
        stats_output_path=args.stats_output,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
