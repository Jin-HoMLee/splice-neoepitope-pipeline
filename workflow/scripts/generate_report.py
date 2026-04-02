#!/usr/bin/env python3
"""generate_report.py — Generate an HTML report with epitope distribution plots,
Fisher's test visualisations, and annotated summary tables.

The report mirrors the key figures from the original Lee et al. (2015) paper:
  * Figure 2: epitope counts per sample (tumour vs. normal box plots)
  * Figure 3: -log10(p) volcano-style plots from Fisher's exact test

Output: self-contained HTML file (plots embedded as base64 PNG).

Usage (standalone):
  python generate_report.py \\
      --stats results/analysis/TCGA-BRCA/statistics.tsv \\
      --predictions results/predictions/TCGA-BRCA/predictions.tsv \\
      --novel-junctions results/junctions/TCGA-BRCA/novel_junctions.tsv \\
      --cancer-type TCGA-BRCA \\
      --output results/reports/TCGA-BRCA/report.html

Usage (Snakemake):
  Called automatically by the ``generate_report`` rule.
"""

import argparse
import base64
import io
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Use a clean seaborn style
sns.set_theme(style="whitegrid", palette="muted")


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _fig_to_base64(fig: plt.Figure) -> str:
    """Encode a matplotlib Figure as a base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def plot_epitope_counts(stats_df: pd.DataFrame, cancer_type: str) -> plt.Figure:
    """Box plot of epitope counts per sample, coloured by sample type.

    Args:
        stats_df:    Statistics DataFrame with columns n_total, sample_type.
        cancer_type: Cancer type label for the plot title.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    # Map sample_type to tumour/normal labels
    def _label(st: str) -> str:
        st_lower = str(st).lower()
        if "normal" in st_lower:
            return "Normal"
        return "Tumour"

    plot_df = stats_df.copy()
    plot_df["group"] = plot_df["sample_type"].apply(_label)

    order = ["Tumour", "Normal"] if "Normal" in plot_df["group"].values else ["Tumour"]

    sns.boxplot(
        data=plot_df,
        x="group",
        y="n_total",
        order=order,
        ax=ax,
        width=0.4,
    )
    ax.set_title(f"{cancer_type} — Epitope counts per sample")
    ax.set_xlabel("Sample type")
    ax.set_ylabel("Number of epitopes (IC50 < 500 nM)")
    plt.tight_layout()
    return fig


def plot_fisher_pvalues(stats_df: pd.DataFrame, cancer_type: str) -> plt.Figure:
    """Bar chart of -log10(Fisher's p-value) for each tumour sample.

    Args:
        stats_df:    Statistics DataFrame with fisher_pvalue column.
        cancer_type: Cancer type label for the plot title.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    tumour_mask = ~stats_df["sample_type"].str.contains(
        "Normal", case=False, na=False
    )
    plot_df = stats_df[tumour_mask].copy()

    if plot_df["fisher_pvalue"].isna().all():
        ax.text(0.5, 0.5, "No normal samples — Fisher's test not applicable",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{cancer_type} — Fisher's exact test")
        plt.tight_layout()
        return fig

    plot_df = plot_df.dropna(subset=["fisher_pvalue"])
    plot_df["-log10_p"] = -np.log10(plot_df["fisher_pvalue"].clip(lower=1e-300))
    plot_df = plot_df.sort_values("-log10_p", ascending=False).reset_index(drop=True)

    ax.bar(range(len(plot_df)), plot_df["-log10_p"], color="steelblue", alpha=0.7)
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", label="p = 0.05")
    ax.set_title(f"{cancer_type} — Fisher's exact test (tumour enrichment)")
    ax.set_xlabel("Tumour sample (ranked)")
    ax.set_ylabel("-log10(p-value)")
    ax.legend()
    plt.tight_layout()
    return fig


def plot_binder_breakdown(predictions_df: pd.DataFrame, cancer_type: str) -> plt.Figure:
    """Stacked bar chart showing strong vs. weak binders per sample type.

    Args:
        predictions_df: MHCflurry predictions DataFrame.
        cancer_type:    Cancer type label.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(5, 4))

    if predictions_df.empty:
        ax.text(0.5, 0.5, "No predictions available",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{cancer_type} — Binder breakdown")
        return fig

    counts = predictions_df["binder_class"].value_counts()
    labels  = [c for c in ["strong", "weak", "non"] if c in counts.index]
    values  = [counts[c] for c in labels]
    colors  = {"strong": "#d62728", "weak": "#ff7f0e", "non": "#aec7e8"}

    ax.bar(labels, values, color=[colors[l] for l in labels])
    ax.set_title(f"{cancer_type} — Binder classification")
    ax.set_xlabel("Binder class")
    ax.set_ylabel("Count of predicted 9-mers")
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# HTML report builder
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Splice Neoepitope Report — {cancer_type}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2em; color: #333; }}
    h1   {{ color: #2c3e50; }}
    h2   {{ color: #34495e; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
    th {{ background-color: #f2f2f2; }}
    tr:nth-child(even) {{ background-color: #f9f9f9; }}
    img {{ max-width: 100%; height: auto; margin: 1em 0; }}
    .figure-row {{ display: flex; flex-wrap: wrap; gap: 1em; }}
    .figure-block {{ flex: 1; min-width: 280px; }}
  </style>
</head>
<body>
  <h1>Splice Neoepitope Pipeline Report</h1>
  <h2>Cancer type: {cancer_type}</h2>
  <p>Generated by the modernised splice neoepitope pipeline
     (reimplementation of Lee et al., 2015).</p>

  <h2>Summary statistics</h2>
  {summary_table}

  <h2>Figures</h2>
  <div class="figure-row">
    <div class="figure-block">
      <h3>Epitope counts per sample</h3>
      <img src="data:image/png;base64,{fig_counts}" alt="Epitope counts"/>
    </div>
    <div class="figure-block">
      <h3>Fisher's exact test (tumour enrichment)</h3>
      <img src="data:image/png;base64,{fig_fisher}" alt="Fisher p-values"/>
    </div>
    <div class="figure-block">
      <h3>Binder classification</h3>
      <img src="data:image/png;base64,{fig_binders}" alt="Binder breakdown"/>
    </div>
  </div>

  <h2>Top strong binders (IC50 &lt; 50 nM)</h2>
  {strong_table}

  <hr/>
  <p><small>Pipeline source: <a href="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline">
  splice-neoepitope-pipeline</a>.
  Original work: Jin-Ho Lee, Seoul National University, 2015.</small></p>
</body>
</html>
"""


def _df_to_html_table(df: pd.DataFrame, max_rows: int = 100) -> str:
    """Convert a DataFrame to an HTML table string (truncated if large)."""
    if df.empty:
        return "<p><em>No data.</em></p>"
    truncated = df.head(max_rows)
    html = truncated.to_html(index=False, border=0, classes="", escape=True)
    if len(df) > max_rows:
        html += f"<p><em>Showing {max_rows} of {len(df)} rows.</em></p>"
    return html


def generate_report(
    stats_tsv: str | Path,
    predictions_tsv: str | Path,
    novel_junctions_tsv: str | Path,
    output_html: str | Path,
    cancer_type: str,
) -> None:
    """Generate the HTML report.

    Args:
        stats_tsv:           Per-sample statistics TSV.
        predictions_tsv:     MHCflurry predictions TSV.
        novel_junctions_tsv: Novel junctions TSV.
        output_html:         Destination HTML report.
        cancer_type:         TCGA project ID.
    """
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    stats_df       = pd.read_csv(stats_tsv, sep="\t")
    predictions_df = pd.read_csv(predictions_tsv, sep="\t")

    # Figures
    fig_counts  = plot_epitope_counts(stats_df, cancer_type)
    fig_fisher  = plot_fisher_pvalues(stats_df, cancer_type)
    fig_binders = plot_binder_breakdown(predictions_df, cancer_type)

    b64_counts  = _fig_to_base64(fig_counts)
    b64_fisher  = _fig_to_base64(fig_fisher)
    b64_binders = _fig_to_base64(fig_binders)

    plt.close("all")

    # Summary table
    summary_cols = [c for c in [
        "sample_id", "sample_type", "n_strong", "n_weak", "n_total",
        "fisher_pvalue", "fisher_oddsratio",
    ] if c in stats_df.columns]
    summary_table_html = _df_to_html_table(stats_df[summary_cols])

    # Top strong binders
    strong_df = predictions_df[predictions_df["binder_class"] == "strong"].sort_values(
        "ic50_nM"
    )
    strong_cols = [c for c in [
        "source_header", "peptide_9mer", "allele", "ic50_nM", "rank",
    ] if c in strong_df.columns]
    strong_table_html = _df_to_html_table(strong_df[strong_cols], max_rows=50)

    html = _HTML_TEMPLATE.format(
        cancer_type=cancer_type,
        summary_table=summary_table_html,
        fig_counts=b64_counts,
        fig_fisher=b64_fisher,
        fig_binders=b64_binders,
        strong_table=strong_table_html,
    )

    output_html.write_text(html, encoding="utf-8")
    log.info("Report written to %s", output_html)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    generate_report(
        stats_tsv=snakemake.input.stats_tsv,  # type: ignore[name-defined]  # noqa: F821
        predictions_tsv=snakemake.input.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        novel_junctions_tsv=snakemake.input.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        output_html=snakemake.output.report_html,  # type: ignore[name-defined]  # noqa: F821
        cancer_type=snakemake.params.cancer_type,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate HTML report from splice neoepitope pipeline results."
    )
    parser.add_argument("--stats", required=True, help="Statistics TSV")
    parser.add_argument("--predictions", required=True, help="Predictions TSV")
    parser.add_argument("--novel-junctions", required=True, help="Novel junctions TSV")
    parser.add_argument("--cancer-type", required=True)
    parser.add_argument("--output", required=True, help="Output HTML report")
    args = parser.parse_args()

    generate_report(
        stats_tsv=args.stats,
        predictions_tsv=args.predictions,
        novel_junctions_tsv=args.novel_junctions,
        output_html=args.output,
        cancer_type=args.cancer_type,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
