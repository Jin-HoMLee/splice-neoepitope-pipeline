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
  2. Runs MHCflurry Class1PresentationPredictor for each HLA allele, producing five
     scores per peptide × allele: ic50_nM, affinity_percentile, processing_score,
     presentation_score, presentation_percentile.
  3. Assigns two classification labels using per-allele percentile thresholds:
       binder_class        — based on affinity_percentile
       presentation_class  — based on presentation_percentile
     Both use the same strong/weak/non threshold pair (default 0.5% / 2.0%).

Allele sources (in priority order):
  1. ``--alleles-tsv`` / ``snakemake.input.alleles_tsv`` — alleles.tsv from
     aggregate_hla_alleles (patient-specific HLA typing via OptiType).
  2. ``--alleles`` / ``snakemake.params.fallback_alleles`` — explicit allele list,
     drawn from config.mhcflurry.fallback_alleles when HLA typing is disabled.

Output TSV columns:
  contig_key  start_nt  peptide  allele  ic50_nM  affinity_percentile
  processing_score  presentation_score  presentation_percentile
  binder_class  presentation_class

Usage (standalone, explicit alleles):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/mhc_affinity.tsv \\
      --alleles HLA-A*02:01 HLA-B*07:02 HLA-C*07:02

Usage (standalone, alleles from HLA typing):
  python run_mhcflurry.py \\
      --peptides-tsv results/peptides/patient_001/peptides.tsv \\
      --output results/predictions/patient_001/mhc_affinity.tsv \\
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
    allele: str,
    predictor=None,
) -> pd.DataFrame:
    """Run MHCflurry predictions on a list of peptides.

    Args:
        peptides:  List of peptide sequences.
        allele:    HLA allele in MHCflurry format (e.g., 'HLA-A*02:01').
        predictor: Pre-loaded Class1PresentationPredictor. Loaded fresh if None.

    Returns:
        DataFrame with columns from Class1PresentationPredictor.predict_to_dataframe:
        peptide, allele, affinity, affinity_percentile, processing_score,
        presentation_score, presentation_percentile.
    """
    if predictor is None:
        predictor = _load_mhcflurry_predictor()

    mhcflurry_allele = _normalise_allele(predictor, allele)

    log.info("Running MHCflurry predictions for %d peptides...", len(peptides))

    return predictor.predict_to_dataframe(
        peptides=peptides,
        alleles=[mhcflurry_allele] * len(peptides),
    )


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

# Module-level predictor cache: populated by _load_predictor_for_gpu() or
# _load_predictor_for_cpu() before the per-allele prediction loop.
_worker_predictor = None


def _load_predictor_for_gpu() -> None:
    """Load predictor into module-level cache for GPU execution.

    Does NOT restrict TF/BLAS thread counts — GPU path offloads matrix ops to
    the device; TF's CPU threads handle batching/preprocessing and benefit from
    all available vCPUs.
    """
    global _worker_predictor
    _worker_predictor = _load_mhcflurry_predictor()


def _load_predictor_for_cpu() -> None:
    """Load predictor into module-level cache for CPU execution.

    Sets TF/BLAS thread env vars to allow TF to use all available cores for
    the single sequential allele loop (no parallel workers competing for cores).
    """
    global _worker_predictor
    _worker_predictor = _load_mhcflurry_predictor()


def _predict_allele_worker(
    allele: str,
    unique_peptides: list[str],
    affinity_percentile_strong: float,
    affinity_percentile_weak: float,
    presentation_percentile_strong: float,
    presentation_percentile_weak: float,
) -> "pd.DataFrame":
    """Predict binding for one allele; returns a lean per-allele DataFrame.

    Columns: peptide, ic50_nM, affinity_percentile, processing_score,
             presentation_score, presentation_percentile, allele,
             binder_class, presentation_class.

    peptides_df (with contig_key/start_nt) is kept in the parent process to
    avoid pickling it into every worker.
    """
    pred_df = _run_mhcflurry_predictions(unique_peptides, allele, predictor=_worker_predictor)
    # Class1PresentationPredictor returns:
    #   affinity (IC50 nM), affinity_percentile, processing_score,
    #   presentation_score, presentation_percentile
    pred_df = pred_df.rename(columns={"affinity": "ic50_nM"})[
        ["peptide", "ic50_nM", "affinity_percentile",
         "processing_score", "presentation_score", "presentation_percentile"]
    ]
    pred_df["allele"] = allele
    pred_df["ic50_nM"] = pred_df["ic50_nM"].fillna(float("inf"))
    pred_df["binder_class"] = pred_df["affinity_percentile"].apply(
        lambda v: classify_by_percentile(v, affinity_percentile_strong, affinity_percentile_weak)
    )
    pred_df["presentation_class"] = pred_df["presentation_percentile"].apply(
        lambda v: classify_by_percentile(v, presentation_percentile_strong, presentation_percentile_weak)
    )
    return pred_df


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_tsv: str | Path,
    output_tsv: str | Path,
    alleles: list[str] | None = None,
    alleles_tsv: str | None = None,
    affinity_percentile_strong: float = 0.5,
    affinity_percentile_weak: float = 2.0,
    presentation_percentile_strong: float = 0.5,
    presentation_percentile_weak: float = 2.0,
) -> None:
    """Run MHCflurry on junction-spanning peptides and write results to TSV.

    Predictions are run for every resolved allele and concatenated into a
    single output TSV (one row per peptide × allele combination).

    Args:
        peptides_tsv:                  TSV of junction-spanning peptides.
        output_tsv:                    Destination TSV file.
        alleles:                       HLA alleles to predict. Ignored when
                                       alleles_tsv is provided.
        alleles_tsv:                   Path to alleles.tsv from aggregate_hla_alleles.
                                       Takes precedence over alleles.
        affinity_percentile_strong:    Affinity percentile threshold for strong binder_class.
        affinity_percentile_weak:      Affinity percentile threshold for weak binder_class.
        presentation_percentile_strong: Presentation percentile threshold for strong presentation_class.
        presentation_percentile_weak:  Presentation percentile threshold for weak presentation_class.

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

    # Step 1: Read junction-spanning peptides
    peptides_df = pd.read_csv(peptides_tsv, sep="\t")

    if peptides_df.empty:
        log.warning("No peptides found in %s", peptides_tsv)
        empty_df = pd.DataFrame(
            columns=["contig_key", "start_nt", "peptide", "allele",
                     "ic50_nM", "affinity_percentile", "processing_score",
                     "presentation_score", "presentation_percentile",
                     "binder_class", "presentation_class"]
        )
        empty_df.to_csv(output_tsv, sep="\t", index=False)
        return

    unique_peptides = peptides_df["peptide"].unique().tolist()
    log.info(
        "Extracted %d peptides (%d unique) from %s",
        len(peptides_df), len(unique_peptides), peptides_tsv,
    )

    # Step 2: Run predictions for each allele sequentially in the main process.
    #
    # GPU: multiple CUDA contexts on a single GPU conflict across subprocesses.
    # CPU: spawning N worker processes each loading a full TF model copy causes
    #      OOM (6 alleles × ~8 GB model ≈ 48 GB on a 52 GB VM). Sequential also
    #      lets TF use all available cores per call, which is faster than N
    #      single-threaded workers (OMP_NUM_THREADS=1) competing for the same cores.
    if _has_gpu():
        log.info(
            "GPU detected — running %d allele(s) sequentially in main process",
            len(resolved_alleles),
        )
        _load_predictor_for_gpu()
    else:
        log.info(
            "No GPU — running %d allele(s) sequentially on CPU (all cores available per allele)",
            len(resolved_alleles),
        )
        _load_predictor_for_cpu()

    pred_dfs = [
        _predict_allele_worker(
            allele, unique_peptides,
            affinity_percentile_strong, affinity_percentile_weak,
            presentation_percentile_strong, presentation_percentile_weak,
        )
        for allele in resolved_alleles
    ]

    # Step 3: Merge per-allele predictions with peptides_df (held in parent to
    # avoid pickling the full table into every worker), then concatenate.
    allele_dfs = [peptides_df.merge(pred_df, on="peptide", how="left") for pred_df in pred_dfs]

    # Step 4: Write combined output
    df = pd.concat(allele_dfs, ignore_index=True)[
        ["contig_key", "start_nt", "peptide", "allele",
         "ic50_nM", "affinity_percentile", "processing_score",
         "presentation_score", "presentation_percentile",
         "binder_class", "presentation_class"]
    ]
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Predictions: %d total (%d allele(s)), "
        "%d strong/%d weak binder_class, "
        "%d strong/%d weak presentation_class → %s",
        len(df),
        len(resolved_alleles),
        (df["binder_class"] == "strong").sum(),
        (df["binder_class"] == "weak").sum(),
        (df["presentation_class"] == "strong").sum(),
        (df["presentation_class"] == "weak").sum(),
        output_tsv,
    )


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
        output_tsv=snakemake.output.mhc_affinity_tsv,  # type: ignore[name-defined]  # noqa: F821
        alleles=None if alleles_tsv else list(snakemake.params.fallback_alleles),  # type: ignore[name-defined]  # noqa: F821
        alleles_tsv=alleles_tsv,
        affinity_percentile_strong=float(snakemake.params.affinity_percentile_strong),  # type: ignore[name-defined]  # noqa: F821
        affinity_percentile_weak=float(snakemake.params.affinity_percentile_weak),  # type: ignore[name-defined]  # noqa: F821
        presentation_percentile_strong=float(snakemake.params.presentation_percentile_strong),  # type: ignore[name-defined]  # noqa: F821
        presentation_percentile_weak=float(snakemake.params.presentation_percentile_weak),  # type: ignore[name-defined]  # noqa: F821
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
    parser.add_argument("--affinity-percentile-strong", type=float, default=0.5)
    parser.add_argument("--affinity-percentile-weak", type=float, default=2.0)
    parser.add_argument("--presentation-percentile-strong", type=float, default=0.5)
    parser.add_argument("--presentation-percentile-weak", type=float, default=2.0)
    args = parser.parse_args()

    if not args.alleles and not args.alleles_tsv:
        parser.error("Provide --alleles or --alleles-tsv.")

    run_prediction(
        peptides_tsv=args.peptides_tsv,
        output_tsv=args.output,
        alleles=args.alleles,
        alleles_tsv=args.alleles_tsv,
        affinity_percentile_strong=args.affinity_percentile_strong,
        affinity_percentile_weak=args.affinity_percentile_weak,
        presentation_percentile_strong=args.presentation_percentile_strong,
        presentation_percentile_weak=args.presentation_percentile_weak,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
