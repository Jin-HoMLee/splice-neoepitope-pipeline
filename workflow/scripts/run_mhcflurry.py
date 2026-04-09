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
  2. Runs MHCflurry affinity prediction for the specified HLA allele.
  3. Classifies each 9-mer as a strong binder (IC50 < 50 nM), weak binder
     (IC50 < 500 nM), or non-binder.

Output TSV columns:
  contig_key  start_nt  peptide  allele  ic50_nM  percentile_rank  binder_class

Usage (standalone):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/TCGA-BRCA/peptides.tsv \\
      --output results/predictions/TCGA-BRCA/predictions.tsv \\
      --allele HLA-A*02:01 \\
      --ic50-strong 50 \\
      --ic50-weak 500

Usage (Snakemake):
  Called automatically by the ``run_mhcflurry`` rule.
"""

import argparse
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
# MHCflurry prediction
# ---------------------------------------------------------------------------

def _run_mhcflurry_predictions(
    peptides: list[str],
    allele: str,
) -> pd.DataFrame:
    """Run MHCflurry predictions on a list of peptides.

    Args:
        peptides: List of peptide sequences.
        allele:   HLA allele in MHCflurry format (e.g., 'HLA-A*02:01').

    Returns:
        DataFrame with columns: peptide, allele, prediction, prediction_percentile.
    """
    try:
        from mhcflurry import Class1AffinityPredictor
    except ImportError:
        log.error(
            "MHCflurry is not installed. Install it with: "
            "pip install mhcflurry && mhcflurry-downloads fetch"
        )
        raise

    # Normalise allele format (MHCflurry uses HLA-A*02:01 format)
    normalised_allele = allele.replace("HLA-", "").replace("*", "").replace(":", "")
    if len(normalised_allele) >= 4:
        mhcflurry_allele = f"HLA-{normalised_allele[0]}*{normalised_allele[1:3]}:{normalised_allele[3:5]}"
    else:
        mhcflurry_allele = allele

    log.info("Using MHCflurry with allele: %s (normalised: %s)", allele, mhcflurry_allele)

    try:
        predictor = Class1AffinityPredictor.load()
    except Exception as exc:
        log.error(
            "Failed to load MHCflurry models. Run: mhcflurry-downloads fetch\n%s", exc
        )
        raise

    supported_alleles = predictor.supported_alleles
    if mhcflurry_allele not in supported_alleles:
        alt_allele = allele.replace("-", "").replace("*", "").replace(":", "")
        log.warning(
            "Allele %s not directly supported, trying alternatives. "
            "Supported alleles: %d total",
            mhcflurry_allele, len(supported_alleles),
        )
        matching = [a for a in supported_alleles if alt_allele[:5] in a.replace("*", "").replace(":", "")]
        if matching:
            mhcflurry_allele = matching[0]
            log.info("Using matched allele: %s", mhcflurry_allele)
        else:
            log.warning("No matching allele found, using original: %s", mhcflurry_allele)

    log.info("Running MHCflurry predictions for %d peptides...", len(peptides))

    # predict_to_dataframe() returns a DataFrame with columns:
    # peptide, allele, prediction, prediction_low, prediction_high, prediction_percentile
    # (mhcflurry 2.2.x renamed mhcflurry_affinity → prediction)
    return predictor.predict_to_dataframe(
        peptides=peptides,
        alleles=[mhcflurry_allele] * len(peptides),
    )


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_tsv: str | Path,
    output_tsv: str | Path,
    allele: str = "HLA-A*02:01",
    ic50_strong: float = 50.0,
    ic50_weak: float = 500.0,
) -> None:
    """Run MHCflurry on junction-spanning 9-mers and write results to TSV.

    Args:
        peptides_tsv:   TSV of junction-spanning 9-mers (contig_key, start_nt, peptide).
        output_tsv:     Destination TSV file.
        allele:         HLA allele (MHCflurry format, e.g., HLA-A*02:01).
        ic50_strong:    Strong-binder IC50 threshold (nM).
        ic50_weak:      Weak-binder IC50 threshold (nM).
    """
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

    # Step 2: Run MHCflurry predictions on unique peptides
    predictions_df = _run_mhcflurry_predictions(unique_peptides, allele)

    peptide_to_prediction = {
        row["peptide"]: {
            "ic50_nM": row["prediction"],
            "percentile_rank": row.get("prediction_percentile", float("nan")),
        }
        for _, row in predictions_df.iterrows()
    }

    # Step 3: Build output — join predictions back to all rows
    def classify(ic50: float) -> str:
        if ic50 <= ic50_strong:
            return "strong"
        if ic50 <= ic50_weak:
            return "weak"
        return "non"

    output_records = []
    for _, row in peptides_df.iterrows():
        pred = peptide_to_prediction.get(
            row["peptide"], {"ic50_nM": float("inf"), "percentile_rank": float("nan")}
        )
        output_records.append({
            "contig_key":      row["contig_key"],
            "start_nt":        row["start_nt"],
            "peptide":         row["peptide"],
            "allele":          allele,
            "ic50_nM":         pred["ic50_nM"],
            "percentile_rank": pred["percentile_rank"],
            "binder_class":    classify(pred["ic50_nM"]),
        })

    df = pd.DataFrame(
        output_records,
        columns=["contig_key", "start_nt", "peptide", "allele",
                 "ic50_nM", "percentile_rank", "binder_class"],
    )
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Predictions: %d total, %d strong, %d weak binders → %s",
        len(df),
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

    run_prediction(
        peptides_tsv=snakemake.input.peptides_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        allele=snakemake.params.hla_allele,  # type: ignore[name-defined]  # noqa: F821
        ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
        ic50_weak=float(snakemake.params.ic50_weak),  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run MHCflurry epitope prediction on junction-spanning 9-mers."
    )
    parser.add_argument("--peptides-tsv", required=True, help="Input peptides TSV")
    parser.add_argument("--output", required=True, help="Output predictions TSV")
    parser.add_argument("--allele", default="HLA-A*02:01", help="HLA allele")
    parser.add_argument("--ic50-strong", type=float, default=50.0)
    parser.add_argument("--ic50-weak", type=float, default=500.0)
    args = parser.parse_args()

    run_prediction(
        peptides_tsv=args.peptides_tsv,
        output_tsv=args.output,
        allele=args.allele,
        ic50_strong=args.ic50_strong,
        ic50_weak=args.ic50_weak,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
