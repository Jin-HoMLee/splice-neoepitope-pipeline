#!/usr/bin/env python3
"""apply_calibrator.py — wire the fitted immunogenicity calibrator into the pipeline.

Slots between MHCflurry (`presentation.tsv`, carries `genotype_presentation_score`)
and TCRdock. Emits two new columns onto the presentation TSV:

  * ``calibrated_immunogenicity_log_odds`` — the centered-isotonic knots in
    ``calibrator_v1.joblib`` applied via ``interp_monotone_extrapolate`` (in-range
    interpolation, monotone linear extrapolation out of range; Issue #805). This is
    a **provisional, SECONDARY** ranking signal: ``genotype_presentation_score``
    stays the primary ranker, so this script only ADDS columns and never re-sorts
    rows. The "provisional" status is discharged by Issue #870 (validate the
    SNV->splice transfer on measured labels); accuracy for splice stays open at
    Issue #680.
  * ``out_of_calibration_support`` (bool) — True where the score is outside the
    interpolation support ``[cx[0], cx[-1]]``. Out-of-support scores are now
    extrapolated (ordered) rather than clipped-flat (Issue #805), but they are
    still flagged here, so the value is ranking-usable while remaining marked as not
    calibration-accurate. NaN scores are also flagged (a missing score cannot be
    certified in-support).

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


def interp_monotone_extrapolate(scores, cx, cy) -> np.ndarray:
    """``np.interp`` but with monotone linear extrapolation beyond ``[cx[0], cx[-1]]``.

    ``np.interp`` clips every out-of-range score to the nearest endpoint value, so
    all scores past a knot collapse to one log-odds - ties that shed ranking
    resolution exactly at the top of the candidate list, where neoepitope ranking
    matters most (Issue #805, empirically: calibrated AUPRC fell *below* raw). Here
    each out-of-range score is instead extrapolated off the terminal knot segment's
    slope, so a strictly higher score always maps to a strictly higher log-odds
    (unless that terminal segment is flat). The extrapolated value is NOT
    calibration-accurate out of support - it exists only to preserve ordering; the
    ``out_of_support`` flag (unchanged) is what marks it untrustworthy as a
    probability. NaN in -> NaN out. Falls back to plain ``np.interp`` with < 2 knots
    (no segment to take a slope from).

    Keep in sync with ``PresentationCalibrator.transform`` in
    ``research/experiments/issue_547_immunogenicity_calibration/calibrator.py`` -
    the two are deliberately decoupled (no research-dir import here) but must apply
    the same map so the offline eval matches production.
    """
    scores = np.atleast_1d(np.asarray(scores, dtype=float))
    cx = np.asarray(cx, dtype=float)
    cy = np.asarray(cy, dtype=float)
    out = np.interp(scores, cx, cy)  # correct in-range; clipped out-of-range; nan->nan
    if cx.size >= 2:
        left_slope = (cy[1] - cy[0]) / (cx[1] - cx[0])
        right_slope = (cy[-1] - cy[-2]) / (cx[-1] - cx[-2])
        out = np.where(scores < cx[0], cy[0] + (scores - cx[0]) * left_slope, out)
        out = np.where(scores > cx[-1], cy[-1] + (scores - cx[-1]) * right_slope, out)
    return out


def calibrate(scores, cx, cy) -> Tuple[np.ndarray, np.ndarray]:
    """Apply the centered-isotonic calibrator knots to an array of scores.

    Args:
        scores: genotype_presentation_score values (array-like; NaN allowed).
        cx, cy: calibrator knot x/y arrays (``cx`` strictly increasing).

    Returns:
        ``(log_odds, out_of_support)`` — ``log_odds`` is
        ``interp_monotone_extrapolate(scores, cx, cy)`` (linearly extrapolated
        outside ``[cx[0], cx[-1]]`` so ordering is preserved, NaN preserved);
        ``out_of_support`` is a boolean array, True outside ``[cx[0], cx[-1]]`` or
        where the score is NaN. Extrapolation changes only the out-of-range
        *value* (ordered instead of clipped-flat); the support flag is unchanged,
        so the #826 applicability gate behaves identically.
    """
    scores = np.asarray(scores, dtype=float)
    cx = np.asarray(cx, dtype=float)
    cy = np.asarray(cy, dtype=float)

    log_odds = interp_monotone_extrapolate(scores, cx, cy)
    nan_mask = np.isnan(scores)
    # interp_monotone_extrapolate preserves nan, but be explicit (defensive).
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
    # np.interp requires strictly-increasing xp; on unsorted xp it returns
    # silently-wrong values (no error). centered_isotonic() produces sorted cx
    # today, but guard against a malformed/resorted future artifact.
    if not np.all(np.diff(cx) > 0):
        raise ValueError(
            f"calibrator cx knots must be strictly increasing "
            f"(np.interp requires sorted xp): {calibrator_path}"
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
