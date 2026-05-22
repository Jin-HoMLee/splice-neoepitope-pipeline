#!/usr/bin/env python3
"""aggregate_filtering_stats.py — Combine per-step funnel stats into a
single per-patient ``filtering_stats.tsv`` (Issue #215).

Reads the per-step stats TSVs emitted by each pipeline stage and writes a
unified long-format TSV with one row per (sample_id, step, category):

    patient_id  sample_id  sample_type  step  category  count

Per-sample steps (``junction-filter``) preserve their existing
``sample_id`` / ``sample_type`` columns. Per-patient steps
(``contig-assemble``, ``peptide-translate``, ``proteome-filter``,
``mhc-affinity``) leave those columns empty since they operate on the
patient-level junction/peptide aggregate.

The funnel categories are listed in Issue #215. Numeric values from the
``junction-filter`` distribution rows (``min_reads``, ``mean_reads``,
``median_reads``, ``max_reads``) are descriptive and not part of the
funnel sum — they coexist in the same long-format TSV under the
``junction-filter`` step.

Usage (Snakemake):
  Called automatically by the ``aggregate_filtering_stats`` rule.
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

UNIFIED_COLUMNS = [
    "patient_id", "sample_id", "sample_type", "step", "category", "count",
]


_PER_SAMPLE_REQUIRED = {"sample_id", "sample_type", "category", "count"}
_PER_PATIENT_REQUIRED = {"category", "count"}


def _read_per_sample_stats(path: Path, step: str) -> "pd.DataFrame":
    """Load a stats file that already has sample_id / sample_type columns."""
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    missing = _PER_SAMPLE_REQUIRED - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: per-sample stats file missing required columns "
            f"{sorted(missing)}. Got: {sorted(df.columns)}"
        )
    df["step"] = step
    return df


def _read_per_patient_stats(path: Path, step: str) -> "pd.DataFrame":
    """Load a stats file with just (category, count). Sample columns left blank."""
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    missing = _PER_PATIENT_REQUIRED - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: per-patient stats file missing required columns "
            f"{sorted(missing)}. Got: {sorted(df.columns)}"
        )
    df["sample_id"] = ""
    df["sample_type"] = ""
    df["step"] = step
    return df


def aggregate(
    patient_id: str,
    junction_filter_tsv: str | Path,
    contig_assemble_tsv: str | Path,
    translate_tsv: str | Path,
    mhc_tsv: str | Path,
    output_tsv: str | Path,
    proteome_tsv: str | Path | None = None,
) -> None:
    """Concatenate per-step stats files into a unified filtering_stats.tsv.

    ``proteome_tsv`` is optional: when ``proteome_filter.enabled`` is false in
    the pipeline config, the proteome step does not run, no stats file is
    produced, and the proteome-filter rows are simply omitted from the audit
    trail.
    """
    import pandas as pd

    frames: list["pd.DataFrame"] = [
        _read_per_sample_stats(Path(junction_filter_tsv), "junction-filter"),
        _read_per_patient_stats(Path(contig_assemble_tsv), "contig-assemble"),
        _read_per_patient_stats(Path(translate_tsv), "peptide-translate"),
    ]
    if proteome_tsv is not None:
        frames.append(_read_per_patient_stats(Path(proteome_tsv), "proteome-filter"))
    frames.append(_read_per_patient_stats(Path(mhc_tsv), "mhc-affinity"))

    df = pd.concat(frames, ignore_index=True, sort=False)
    df["patient_id"] = patient_id
    df = df[UNIFIED_COLUMNS]

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info("Filtering audit trail written to %s (%d rows)", output_tsv, len(df))


def _snakemake_main() -> None:
    sm = snakemake  # type: ignore[name-defined]  # noqa: F821
    log_file = sm.log[0]
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    aggregate(
        patient_id=sm.wildcards.patient_id,
        junction_filter_tsv=sm.input.junction_filter,
        contig_assemble_tsv=sm.input.contig_assemble,
        translate_tsv=sm.input.translate,
        mhc_tsv=sm.input.mhc,
        output_tsv=sm.output.filtering_stats,
        proteome_tsv=getattr(sm.input, "proteome", None),
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate per-step pipeline stats into filtering_stats.tsv (Issue #215).",
    )
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--junction-filter", required=True)
    parser.add_argument("--contig-assemble", required=True)
    parser.add_argument("--translate", required=True)
    parser.add_argument(
        "--proteome",
        default=None,
        help="Optional. Omit when proteome_filter is disabled in the pipeline config.",
    )
    parser.add_argument("--mhc", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    aggregate(
        patient_id=args.patient_id,
        junction_filter_tsv=args.junction_filter,
        contig_assemble_tsv=args.contig_assemble,
        translate_tsv=args.translate,
        mhc_tsv=args.mhc,
        output_tsv=args.output,
        proteome_tsv=args.proteome,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
