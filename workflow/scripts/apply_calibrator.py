#!/usr/bin/env python3
"""apply_calibrator.py — wire the fitted immunogenicity calibrator into the pipeline.

Slots between MHCflurry (`presentation.tsv`, carries `genotype_presentation_score`)
and TCRdock. Emits two new columns onto the presentation TSV:

  * ``calibrated_immunogenicity_log_odds`` — ``np.interp(score, cx, cy)`` over the
    serialized centered-isotonic knots in ``calibrator_v1.joblib``. This is a
    **provisional, SECONDARY** ranking signal: ``genotype_presentation_score`` stays
    the primary ranker, so this script only ADDS columns and never re-sorts rows.
    The "provisional" status is discharged by Issue #870 (validate the SNV->splice
    transfer on measured labels); accuracy for splice stays open at Issue #680.
  * ``out_of_calibration_support`` (bool) — True where the score is outside the
    interpolation support ``[cx[0], cx[-1]]`` (both the floor-clip and ceil-clip),
    where ``np.interp`` returns a flat constant rather than a fitted gradient. NaN
    scores are also flagged (a missing score cannot be certified in-support).

The support thresholds are read FROM the artifact (``cx[0]`` / ``cx[-1]``), never
hard-coded, so a future calibrator refit cannot silently desync the flag from the
curve.

Calibrator fitting + training lives in the research experiment
(``research/experiments/issue_547_immunogenicity_calibration/``); production only
applies the serialized knots — no research-dir import, no fitting code here.

Issue #709 (https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709).
"""

import argparse
import logging
from pathlib import Path
from typing import Optional, Tuple, Union

import joblib
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_SCORE_COL = "genotype_presentation_score"
KNOWN_CALIBRATOR_VERSIONS = {"v1"}


def calibrate(scores, cx, cy) -> Tuple[np.ndarray, np.ndarray]:
    """Apply the centered-isotonic calibrator knots to an array of scores.

    Args:
        scores: genotype_presentation_score values (array-like; NaN allowed).
        cx, cy: calibrator knot x/y arrays (``cx`` strictly increasing).

    Returns:
        ``(log_odds, out_of_support)`` — ``log_odds`` is ``np.interp(scores, cx, cy)``
        (flat-clipped outside ``[cx[0], cx[-1]]``, NaN preserved); ``out_of_support``
        is a boolean array, True outside ``[cx[0], cx[-1]]`` or where the score is NaN.
    """
    scores = np.asarray(scores, dtype=float)
    cx = np.asarray(cx, dtype=float)
    cy = np.asarray(cy, dtype=float)

    log_odds = np.interp(scores, cx, cy)  # np.interp(nan) -> nan; clips to endpoints
    nan_mask = np.isnan(scores)
    # np.interp already returns nan for nan input, but be explicit (defensive).
    log_odds = np.where(nan_mask, np.nan, log_odds)
    out_of_support = (scores < cx[0]) | (scores > cx[-1]) | nan_mask
    return log_odds, out_of_support


def load_calibrator_knots(calibrator_path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """Load ``(cx, cy)`` knot arrays from a serialized calibrator artifact.

    The artifact is a joblib dict (``{cx, cy, version, ...}``) produced by
    ``PresentationCalibrator.save``. Only the knots are needed to apply it. An
    unrecognised ``version`` warns (the cx/cy format is stable) rather than failing.
    """
    d = joblib.load(calibrator_path)
    version = d.get("version", "unknown")
    if version not in KNOWN_CALIBRATOR_VERSIONS:
        log.warning(
            "calibrator version %r not in known set %s; applying on cx/cy knots anyway",
            version,
            sorted(KNOWN_CALIBRATOR_VERSIONS),
        )
    cx = np.asarray(d["cx"], dtype=float)
    cy = np.asarray(d["cy"], dtype=float)
    if cx.ndim != 1 or cx.shape != cy.shape or cx.size < 2:
        raise ValueError(
            f"malformed calibrator knots in {calibrator_path}: "
            f"cx{cx.shape} cy{cy.shape}"
        )
    return cx, cy


def apply_calibrator(
    input_tsv: Union[str, Path],
    calibrator_path: Union[str, Path],
    output_tsv: Union[str, Path],
    score_col: str = DEFAULT_SCORE_COL,
) -> int:
    """Read ``input_tsv``, append calibrated columns, write ``output_tsv``.

    Returns the number of rows written. Row order and all original columns are
    preserved unchanged (secondary-signal guardrail).
    """
    input_tsv = Path(input_tsv)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_tsv, sep="\t")
    if score_col not in df.columns:
        raise KeyError(
            f"input TSV {input_tsv} lacks required score column {score_col!r}"
        )

    cx, cy = load_calibrator_knots(calibrator_path)
    log_odds, out_of_support = calibrate(df[score_col].to_numpy(dtype=float), cx, cy)
    df["calibrated_immunogenicity_log_odds"] = log_odds
    df["out_of_calibration_support"] = out_of_support

    df.to_csv(output_tsv, sep="\t", index=False)
    n_oos = int(np.asarray(out_of_support).sum())
    log.info(
        "Calibrated %d rows (%d out-of-support, support=[%.4f, %.4f]) -> %s",
        len(df),
        n_oos,
        cx[0],
        cx[-1],
        output_tsv,
    )
    return len(df)


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Append calibrated_immunogenicity_log_odds + "
        "out_of_calibration_support to an MHCflurry presentation TSV."
    )
    parser.add_argument("--input", required=True, help="MHCflurry presentation.tsv")
    parser.add_argument(
        "--calibrator", required=True, help="serialized calibrator (.joblib)"
    )
    parser.add_argument("--output", required=True, help="augmented TSV output")
    parser.add_argument(
        "--score-col",
        default=DEFAULT_SCORE_COL,
        help=f"score column to calibrate (default: {DEFAULT_SCORE_COL})",
    )
    args = parser.parse_args()
    apply_calibrator(args.input, args.calibrator, args.output, args.score_col)


if __name__ == "__main__":
    _cli_main()
