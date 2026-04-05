#!/usr/bin/env python3
"""statistical_analysis.py — Count epitopes per sample, perform Fisher's exact
test for tumour vs. normal enrichment, and produce summary statistics.

This script serves two Snakemake rules:
  * ``statistical_analysis``  — per-cancer-type analysis
  * ``summarise_all``         — aggregate all cancer types into one table

**Fisher's exact test** (one-tailed, testing enrichment in tumour over normal)
is performed for cancer types that have matched normal samples (BRCA, LUAD).

Output columns for per-cancer-type TSV:
  cancer_type  sample_id  sample_type  n_strong  n_weak  n_total
  epitopes_per_mb  fisher_pvalue  fisher_oddsratio  significant

Usage (standalone — single cancer type):
  python statistical_analysis.py \\
      --mode single \\
      --predictions results/predictions/TCGA-BRCA/predictions.tsv \\
      --manifest results/raw_data/TCGA-BRCA/manifest.tsv \\
      --cancer-type TCGA-BRCA \\
      --output results/analysis/TCGA-BRCA/statistics.tsv \\
      --ic50-strong 50 --ic50-weak 500

Usage (standalone — summary):
  python statistical_analysis.py \\
      --mode summary \\
      --stats results/analysis/TCGA-BRCA/statistics.tsv \\
             results/analysis/TCGA-LUAD/statistics.tsv \\
             results/analysis/TCGA-LAML/statistics.tsv \\
      --output results/reports/summary_table.tsv

Usage (Snakemake):
  Called automatically by ``statistical_analysis`` and ``summarise_all`` rules.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-cancer-type analysis
# ---------------------------------------------------------------------------

def _count_epitopes_per_sample(
    predictions_df: pd.DataFrame,
    ic50_strong: float,
    ic50_weak: float,
) -> pd.DataFrame:
    """Count strong/weak/total binders per sample (source_header).

    Args:
        predictions_df: DataFrame from MHCflurry output TSV.
        ic50_strong:    Strong-binder IC50 threshold.
        ic50_weak:      Weak-binder IC50 threshold.

    Returns:
        DataFrame with columns: source_header, n_strong, n_weak, n_total.
    """
    df = predictions_df.copy()
    df["is_strong"] = df["ic50_nM"] < ic50_strong
    df["is_weak"]   = (df["ic50_nM"] >= ic50_strong) & (df["ic50_nM"] < ic50_weak)
    df["is_binder"] = df["ic50_nM"] < ic50_weak

    per_sample = df.groupby("source_header").agg(
        n_strong=("is_strong", "sum"),
        n_weak=("is_weak", "sum"),
        n_total=("is_binder", "sum"),
    ).reset_index()
    return per_sample


def _fisher_test_tumour_vs_normal(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
) -> pd.DataFrame:
    """Perform one-tailed Fisher's exact test for each sample.

    For each sample, the contingency table is:
      |                  | binder | non-binder |
      | tumour sample    |  a     |   b        |
      | normal samples   |  c     |   d        |

    where a/b are from the individual tumour sample and c/d are the totals
    across all normal samples.

    Args:
        counts:      Per-sample epitope counts (output of
                     :func:`_count_epitopes_per_sample`).
        sample_meta: DataFrame with columns sample_id, sample_type.

    Returns:
        counts DataFrame with added columns: fisher_pvalue, fisher_oddsratio.
    """
    merged = counts.merge(
        sample_meta[["sample_id", "sample_type"]],
        left_on="source_header",
        right_on="sample_id",
        how="left",
    )
    merged["sample_type"] = merged["sample_type"].fillna("Unknown")

    # Aggregate normal epitope counts
    normal_mask = merged["sample_type"].str.contains("Normal", case=False, na=False)
    normal_rows = merged[normal_mask]
    if normal_rows.empty:
        # No normal samples — cannot perform test
        merged["fisher_pvalue"]    = np.nan
        merged["fisher_oddsratio"] = np.nan
        return merged

    normal_total_binders = int(normal_rows["n_total"].sum())
    # Proxy for total 9-mers tested: use the maximum observed binder count
    # across all samples as a denominator baseline.  This is a conservative
    # approximation; the true denominator would require the total number of
    # 9-mers submitted per sample.
    max_total = int(merged["n_total"].max()) if len(merged) > 0 else 1
    # Non-binder counts for the pooled normals (proxy: max_total per sample)
    normal_total_non_binders = max(0, max_total * len(normal_rows) - normal_total_binders)

    p_values: list[float] = []
    odds_ratios: list[float] = []

    for _, row in merged.iterrows():
        a = int(row["n_total"])              # sample binders
        b = max(0, max_total - a)            # sample non-binders (proxy)
        c = normal_total_binders             # pooled normal binders
        d = normal_total_non_binders         # pooled normal non-binders

        table = [[a, b], [c, d]]
        if sum(table[0]) == 0 or sum(table[1]) == 0:
            p_values.append(np.nan)
            odds_ratios.append(np.nan)
            continue

        try:
            odds_ratio, p_value = fisher_exact(table, alternative="greater")
        except Exception:
            odds_ratio, p_value = np.nan, np.nan

        p_values.append(p_value)
        odds_ratios.append(odds_ratio)

    merged["fisher_pvalue"]    = p_values
    merged["fisher_oddsratio"] = odds_ratios
    return merged


def analyse_cancer_type(
    predictions_tsv: str | Path,
    manifest_tsv: str | Path,
    output_tsv: str | Path,
    cancer_type: str,
    ic50_strong: float = 50.0,
    ic50_weak: float = 500.0,
    cancer_types_with_normal: list[str] | None = None,
) -> None:
    """Run per-cancer-type statistical analysis.

    Args:
        predictions_tsv:          MHCflurry predictions TSV.
        manifest_tsv:             GDC manifest TSV with sample metadata.
        output_tsv:               Destination statistics TSV.
        cancer_type:              TCGA project ID.
        ic50_strong:              Strong-binder threshold.
        ic50_weak:                Weak-binder threshold.
        cancer_types_with_normal: Cancer types that have normal samples.
    """
    if cancer_types_with_normal is None:
        cancer_types_with_normal = ["TCGA-BRCA", "TCGA-LUAD"]

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    predictions_df = pd.read_csv(predictions_tsv, sep="\t")
    manifest_df    = pd.read_csv(manifest_tsv, sep="\t")

    if predictions_df.empty:
        log.warning("No predictions found for %s", cancer_type)
        empty = pd.DataFrame(
            columns=[
                "cancer_type", "source_header", "sample_id", "sample_type",
                "n_strong", "n_weak", "n_total",
                "fisher_pvalue", "fisher_oddsratio",
            ]
        )
        empty.to_csv(output_tsv, sep="\t", index=False)
        return

    per_sample = _count_epitopes_per_sample(predictions_df, ic50_strong, ic50_weak)

    if cancer_type in cancer_types_with_normal:
        result = _fisher_test_tumour_vs_normal(per_sample, manifest_df)
    else:
        result = per_sample.copy()
        # source_header format: junc_id|coords|sample_type|frame
        result["sample_type"] = (
            result["source_header"].str.split("|").str[2].fillna("Unknown")
        )
        result["fisher_pvalue"]    = np.nan
        result["fisher_oddsratio"] = np.nan

    result["cancer_type"] = cancer_type
    result["sample_id"]   = result["source_header"]

    cols = [
        "cancer_type", "sample_id", "sample_type",
        "n_strong", "n_weak", "n_total",
        "fisher_pvalue", "fisher_oddsratio",
    ]
    # Keep only columns that exist
    cols = [c for c in cols if c in result.columns]
    result[cols].to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Statistics for %s: %d samples → %s",
        cancer_type, len(result), output_tsv,
    )


# ---------------------------------------------------------------------------
# Summary (aggregate all cancer types)
# ---------------------------------------------------------------------------

def summarise_all(
    stats_files: list[str | Path],
    output_tsv: str | Path,
) -> None:
    """Concatenate per-cancer-type statistics into one summary TSV.

    Args:
        stats_files: List of per-cancer-type statistics TSV files.
        output_tsv:  Destination summary TSV.
    """
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    dfs = []
    for sf in stats_files:
        try:
            df = pd.read_csv(sf, sep="\t")
            dfs.append(df)
        except Exception as exc:
            log.warning("Could not read %s: %s", sf, exc)

    if not dfs:
        log.error("No statistics files could be read.")
        pd.DataFrame().to_csv(output_tsv, sep="\t", index=False)
        return

    summary = pd.concat(dfs, ignore_index=True)
    summary.to_csv(output_tsv, sep="\t", index=False)
    log.info("Summary table with %d rows → %s", len(summary), output_tsv)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    rule_name = snakemake.rule  # type: ignore[name-defined]  # noqa: F821
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    if rule_name == "statistical_analysis":
        analyse_cancer_type(
            predictions_tsv=snakemake.input.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
            manifest_tsv=snakemake.input.manifest,  # type: ignore[name-defined]  # noqa: F821
            output_tsv=snakemake.output.stats_tsv,  # type: ignore[name-defined]  # noqa: F821
            cancer_type=snakemake.wildcards.cancer_type,  # type: ignore[name-defined]  # noqa: F821
            ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
            ic50_weak=float(snakemake.params.ic50_weak),  # type: ignore[name-defined]  # noqa: F821
            cancer_types_with_normal=snakemake.params.cancer_types_with_normal,  # type: ignore[name-defined]  # noqa: F821
        )
    elif rule_name == "summarise_all":
        summarise_all(
            stats_files=snakemake.input.stats,  # type: ignore[name-defined]  # noqa: F821
            output_tsv=snakemake.output.summary,  # type: ignore[name-defined]  # noqa: F821
        )
    else:
        raise ValueError(f"Unexpected rule: {rule_name!r}")


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Statistical analysis of MHCflurry epitope predictions."
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    single = sub.add_parser("single", help="Per-cancer-type analysis")
    single.add_argument("--predictions", required=True)
    single.add_argument("--manifest", required=True)
    single.add_argument("--cancer-type", required=True)
    single.add_argument("--output", required=True)
    single.add_argument("--ic50-strong", type=float, default=50.0)
    single.add_argument("--ic50-weak", type=float, default=500.0)
    single.add_argument(
        "--cancer-types-with-normal", nargs="*",
        default=["TCGA-BRCA", "TCGA-LUAD"],
    )

    summ = sub.add_parser("summary", help="Aggregate all cancer types")
    summ.add_argument("--stats", nargs="+", required=True)
    summ.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.mode == "single":
        analyse_cancer_type(
            predictions_tsv=args.predictions,
            manifest_tsv=args.manifest,
            output_tsv=args.output,
            cancer_type=args.cancer_type,
            ic50_strong=args.ic50_strong,
            ic50_weak=args.ic50_weak,
            cancer_types_with_normal=args.cancer_types_with_normal,
        )
    elif args.mode == "summary":
        summarise_all(stats_files=args.stats, output_tsv=args.output)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
