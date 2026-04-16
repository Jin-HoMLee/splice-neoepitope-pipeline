#!/usr/bin/env python3
"""run_mhcflurry.py — Wrapper to run MHCflurry 2.x on junction-spanning 9-mers.

MHCflurry is an open-source MHC-I binding predictor that achieves
state-of-the-art performance.  Unlike NetMHCPan, it does not require
academic registration and can be installed via pip.

Reference:
  O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
  of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
  Cell Systems, 11(1), 42-48.e7.

The script:
  1. Reads junction-spanning 9-mers from the TSV produced by translate_peptides.py.
  2. Runs MHCflurry affinity prediction for each HLA allele.
  3. Classifies each 9-mer as a strong binder (IC50 <= 50 nM), weak binder
     (IC50 <= 500 nM), or non-binder.

Allele sources (in priority order):
  1. ``--alleles-tsv`` / ``snakemake.input.alleles_tsv`` — alleles.tsv from
     aggregate_hla_alleles (patient-specific HLA typing via OptiType).
  2. ``--alleles`` / ``snakemake.params.fallback_alleles`` — explicit allele list,
     drawn from config.mhcflurry.fallback_alleles when HLA typing is disabled.

Output TSV columns:
  contig_key  start_nt  peptide  allele  ic50_nM  percentile_rank  binder_class

Usage (standalone, explicit alleles):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/predictions.tsv \\
      --alleles HLA-A*02:01 HLA-B*07:02 HLA-C*07:02 \\
      --ic50-strong 50 \\
      --ic50-weak 500

Usage (standalone, alleles from HLA typing):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/predictions.tsv \\
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
    """Load and return the MHCflurry Class1AffinityPredictor (once per run)."""
    try:
        from mhcflurry import Class1AffinityPredictor
    except ImportError:
        log.error(
            "MHCflurry is not installed. Install it with: "
            "pip install mhcflurry && mhcflurry-downloads fetch"
        )
        raise
    try:
        predictor = Class1AffinityPredictor.load()
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
    allele: str,
    predictor=None,
) -> pd.DataFrame:
    """Run MHCflurry predictions on a list of peptides.

    Args:
        peptides:  List of peptide sequences.
        allele:    HLA allele in MHCflurry format (e.g., 'HLA-A*02:01').
        predictor: Pre-loaded Class1AffinityPredictor. Loaded fresh if None.

    Returns:
        DataFrame with columns: peptide, allele, prediction, prediction_percentile.
    """
    if predictor is None:
        predictor = _load_mhcflurry_predictor()

    mhcflurry_allele = _normalise_allele(predictor, allele)

    log.info("Running MHCflurry predictions for %d peptides...", len(peptides))

    # predict_to_dataframe() returns a DataFrame with columns:
    # peptide, allele, prediction, prediction_low, prediction_high, prediction_percentile
    # (mhcflurry 2.2.x renamed mhcflurry_affinity → prediction)
    return predictor.predict_to_dataframe(
        peptides=peptides,
        alleles=[mhcflurry_allele] * len(peptides),
    )


# ---------------------------------------------------------------------------
# Binder classification
# ---------------------------------------------------------------------------

def classify(ic50: float, ic50_strong: float = 50.0, ic50_weak: float = 500.0) -> str:
    """Classify a peptide as a strong, weak, or non-binder by IC50 (nM).

    Boundaries are inclusive: IC50 <= ic50_strong → strong,
    IC50 <= ic50_weak → weak, otherwise non.
    """
    if ic50 <= ic50_strong:
        return "strong"
    if ic50 <= ic50_weak:
        return "weak"
    return "non"


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_tsv: str | Path,
    output_tsv: str | Path,
    alleles: list[str] | None = None,
    alleles_tsv: str | None = None,
    ic50_strong: float = 50.0,
    ic50_weak: float = 500.0,
) -> None:
    """Run MHCflurry on junction-spanning 9-mers and write results to TSV.

    Predictions are run for every resolved allele and concatenated into a
    single output TSV (one row per peptide × allele combination).

    Args:
        peptides_tsv:   TSV of junction-spanning 9-mers (contig_key, start_nt, peptide).
        output_tsv:     Destination TSV file.
        alleles:        HLA alleles to predict (MHCflurry format, e.g. ['HLA-A*02:01']).
                        Ignored when alleles_tsv is provided.
        alleles_tsv:    Path to alleles.tsv from aggregate_hla_alleles. When provided,
                        alleles are loaded from the file (takes precedence over alleles).
        ic50_strong:    Strong-binder IC50 threshold (nM).
        ic50_weak:      Weak-binder IC50 threshold (nM).

    Raises:
        ValueError: If neither alleles nor alleles_tsv is provided.
    """
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

    # Step 1: Read junction-spanning 9-mers
    peptides_df = pd.read_csv(peptides_tsv, sep="\t")

    if peptides_df.empty:
        log.warning("No peptides found in %s", peptides_tsv)
        empty_df = pd.DataFrame(
            columns=["contig_key", "start_nt", "peptide", "allele",
                     "ic50_nM", "percentile_rank", "binder_class"]
        )
        empty_df.to_csv(output_tsv, sep="\t", index=False)
        return

    unique_peptides = peptides_df["peptide"].unique().tolist()
    log.info(
        "Extracted %d 9-mers (%d unique) from %s",
        len(peptides_df), len(unique_peptides), peptides_tsv,
    )

    # Step 2: Load the predictor once, then run for each allele
    predictor = _load_mhcflurry_predictor()

    allele_dfs = []
    for allele in resolved_alleles:
        pred_df = _run_mhcflurry_predictions(unique_peptides, allele, predictor=predictor)

        # Rename mhcflurry output columns and merge with the full peptides table
        pred_df = pred_df.rename(columns={
            "prediction": "ic50_nM",
            "prediction_percentile": "percentile_rank",
        })[["peptide", "ic50_nM", "percentile_rank"]]

        merged = peptides_df.merge(pred_df, on="peptide", how="left")
        merged["allele"] = allele
        merged["ic50_nM"] = merged["ic50_nM"].fillna(float("inf"))
        merged["percentile_rank"] = merged["percentile_rank"].fillna(float("nan"))
        merged["binder_class"] = merged["ic50_nM"].apply(
            lambda v: classify(v, ic50_strong, ic50_weak)
        )
        allele_dfs.append(merged)

    # Step 3: Write combined output
    df = pd.concat(allele_dfs, ignore_index=True)[
        ["contig_key", "start_nt", "peptide", "allele",
         "ic50_nM", "percentile_rank", "binder_class"]
    ]
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Predictions: %d total (%d allele(s)), %d strong, %d weak binders → %s",
        len(df),
        len(resolved_alleles),
        (df["binder_class"] == "strong").sum(),
        (df["binder_class"] == "weak").sum(),
        output_tsv,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    # Use patient-specific alleles when HLA typing is enabled (alleles_tsv input
    # is only wired in predict.smk when config.hla.enabled is true).
    alleles_tsv = getattr(snakemake.input, "alleles_tsv", None)  # type: ignore[name-defined]  # noqa: F821
    run_prediction(
        peptides_tsv=snakemake.input.peptides_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        alleles=None if alleles_tsv else list(snakemake.params.fallback_alleles),  # type: ignore[name-defined]  # noqa: F821
        alleles_tsv=alleles_tsv,
        ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
        ic50_weak=float(snakemake.params.ic50_weak),  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run MHCflurry epitope prediction on junction-spanning 9-mers."
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
    parser.add_argument("--ic50-strong", type=float, default=50.0)
    parser.add_argument("--ic50-weak", type=float, default=500.0)
    args = parser.parse_args()

    if not args.alleles and not args.alleles_tsv:
        parser.error("Provide --alleles or --alleles-tsv.")

    run_prediction(
        peptides_tsv=args.peptides_tsv,
        output_tsv=args.output,
        alleles=args.alleles,
        alleles_tsv=args.alleles_tsv,
        ic50_strong=args.ic50_strong,
        ic50_weak=args.ic50_weak,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
