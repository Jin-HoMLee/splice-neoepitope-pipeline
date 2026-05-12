#!/usr/bin/env python3
"""Cross-patient cohort aggregator for per-patient ``report.tsv`` files.

Pipeline runs are single-patient by design (see `workflow/rules/common.smk`),
so combining patients into a cohort table is a post-pipeline research-time
step rather than a Snakemake rule.

Each per-patient ``report.tsv`` carries the unified long-format schema
``patient_id | stage | metric | value | notes`` (see
``workflow/scripts/generate_report.py::_build_report_tsv``), so cohort
aggregation is a plain concatenation. Output is sorted by
(``patient_id``, ``stage``, ``metric``) for stable diffs.

Usage:
    # explicit input paths
    python research/scripts/aggregate_cohort.py \\
        --inputs results/patient_001/reports/report.tsv \\
                 results/patient_002/reports/report.tsv \\
        --output research/cohort_summary.tsv

    # patient IDs + results root (paths resolved as {root}/{id}/reports/report.tsv)
    python research/scripts/aggregate_cohort.py \\
        --patients patient_001 patient_002 \\
        --results-root results \\
        --output research/cohort_summary.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ("patient_id", "stage", "metric", "value", "notes")
SORT_KEYS = ["patient_id", "stage", "metric"]


def resolve_inputs(
    inputs: list[str] | None,
    patients: list[str] | None,
    results_root: str,
) -> list[Path]:
    if inputs:
        return [Path(p) for p in inputs]
    return [
        Path(results_root) / pid / "reports" / "report.tsv" for pid in patients
    ]


def load_report(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"report.tsv not found: {path}")
    df = pd.read_csv(path, sep="\t")
    missing = set(REQUIRED_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: report.tsv is missing required columns "
            f"{sorted(missing)}. Got: {sorted(df.columns)}"
        )
    return df


def aggregate(paths: list[Path]) -> pd.DataFrame:
    frames = [load_report(p) for p in paths]
    cohort = pd.concat(frames, ignore_index=True)
    return cohort.sort_values(SORT_KEYS, kind="stable").reset_index(drop=True)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--inputs", nargs="+", metavar="PATH",
        help="Explicit list of per-patient report.tsv paths.",
    )
    src.add_argument(
        "--patients", nargs="+", metavar="ID",
        help="Patient IDs; paths resolved as {results-root}/{ID}/reports/report.tsv.",
    )
    p.add_argument(
        "--results-root", default="results",
        help="Root directory for patient results (default: results/). Only used with --patients.",
    )
    p.add_argument(
        "--output", "-o", required=True, type=Path,
        help="Output cohort TSV path.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    paths = resolve_inputs(args.inputs, args.patients, args.results_root)
    cohort = aggregate(paths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(args.output, sep="\t", index=False)
    print(
        f"Wrote {len(cohort)} rows from {len(paths)} patient(s) "
        f"({cohort['patient_id'].nunique()} unique) to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
