"""Tests for apply_calibrator.py — wire the fitted immunogenicity calibrator into
the pipeline as a post-MHCflurry / pre-TCRdock column.

Issue #709 (https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709).

The calibrator artifact (`calibrator_v1.joblib`) is a serialized dict of centered-
isotonic knots {cx, cy}. Applying it is `np.interp(score, cx, cy)` — flat-clipped
outside the support `[cx[0], cx[-1]]`. Per the Scientist's #826 contract:

  * `calibrated_immunogenicity_log_odds` is a **provisional secondary** signal —
    `genotype_presentation_score` stays primary, so the rule only ADDS columns and
    must NOT re-sort rows.
  * `out_of_calibration_support` flags rows outside `[cx[0], cx[-1]]` (both the
    floor-clip and ceil-clip), where `np.interp` returns a flat constant.
"""

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest

from apply_calibrator import apply_calibrator, calibrate, load_calibrator_knots

# A tiny synthetic calibrator: knots chosen so np.interp is hand-checkable.
#   cx = [0.1, 0.5, 0.9]   cy = [-5.0, -1.0, 0.0]
_CX = [0.1, 0.5, 0.9]
_CY = [-5.0, -1.0, 0.0]


# --------------------------------------------------------------------------- #
# Pure-logic tests (numpy only) — calibrate(scores, cx, cy)
# --------------------------------------------------------------------------- #
def test_calibrate_in_support_interpolates():
    log_odds, oos = calibrate([0.5, 0.3], _CX, _CY)
    # 0.5 -> exact knot -1.0; 0.3 -> halfway on [0.1,0.5] segment -> -3.0
    assert log_odds[0] == pytest.approx(-1.0)
    assert log_odds[1] == pytest.approx(-3.0)
    assert not oos[0] and not oos[1]


def test_calibrate_floor_extrapolates_and_is_flagged():
    # Issue #805: below cx[0] extrapolates off the first segment slope
    # (= (cy[1]-cy[0])/(cx[1]-cx[0]) = 10.0), not clipped flat at cy[0].
    log_odds, oos = calibrate([0.05], _CX, _CY)
    assert log_odds[0] == pytest.approx(-5.0 + (0.05 - 0.1) * 10.0)  # -5.5
    assert bool(oos[0]) is True  # still flagged out of support


def test_calibrate_ceil_extrapolates_and_is_flagged():
    # Above cx[-1] extrapolates off the last segment slope
    # (= (cy[-1]-cy[-2])/(cx[-1]-cx[-2]) = 2.5), not clipped flat at cy[-1].
    log_odds, oos = calibrate([0.95], _CX, _CY)
    assert log_odds[0] == pytest.approx(0.0 + (0.95 - 0.9) * 2.5)  # 0.125
    assert bool(oos[0]) is True


def test_calibrate_extrapolation_preserves_order_and_breaks_ties():
    # The whole point of #805: two distinct out-of-range scores must map to two
    # distinct, correctly-ordered log-odds (np.interp clipping tied them).
    log_odds, oos = calibrate([0.02, 0.05, 0.95, 1.00], _CX, _CY)
    below_lo, below_hi, above_lo, above_hi = log_odds
    assert below_lo < below_hi          # ordered below the floor (no tie)
    assert above_lo < above_hi          # ordered above the ceiling (no tie)
    assert below_hi < above_lo          # global monotonicity across the range
    assert len(set(log_odds)) == 4      # all distinct - clipping would give ties
    assert all(bool(f) for f in oos)    # all four still flagged out of support


def test_calibrate_boundary_is_in_support():
    # Exactly cx[0]/cx[-1] is IN support (flag is strict < / >).
    _, oos = calibrate([0.1, 0.9], _CX, _CY)
    assert not oos[0] and not oos[1]


def test_calibrate_nan_score_is_flagged_out_of_support():
    # NaN presentation score must not silently read as in-support (NaN < x is False).
    log_odds, oos = calibrate([np.nan], _CX, _CY)
    assert np.isnan(log_odds[0])
    assert bool(oos[0]) is True


# --------------------------------------------------------------------------- #
# IO round-trip — apply_calibrator(input_tsv, calibrator_path, output_tsv)
# --------------------------------------------------------------------------- #
@pytest.fixture
def synthetic_calibrator(tmp_path):
    path = tmp_path / "calibrator_v1.joblib"
    joblib.dump({"cx": np.array(_CX), "cy": np.array(_CY), "version": "v1"}, path)
    return path


@pytest.fixture
def presentation_tsv(tmp_path):
    df = pd.DataFrame(
        {
            "peptide": ["AAA", "BBB", "CCC"],
            "genotype_presentation_score": [0.5, 0.05, 0.95],
            "presentation_class": ["weak", "non", "strong"],
        }
    )
    path = tmp_path / "presentation.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_apply_calibrator_adds_two_columns(synthetic_calibrator, presentation_tsv, tmp_path):
    out = tmp_path / "presentation_calibrated.tsv"
    apply_calibrator(presentation_tsv, synthetic_calibrator, out)
    res = pd.read_csv(out, sep="\t")
    assert "calibrated_immunogenicity_log_odds" in res.columns
    assert "out_of_calibration_support" in res.columns
    # values: 0.5->-1.0 (in); 0.05->-5.5 (floor, extrapolated off slope 10.0, oos);
    # 0.95->0.125 (ceil, extrapolated off slope 2.5, oos). Out-of-support flags
    # unchanged by the #805 extrapolation - only the out-of-range *values* differ.
    assert res["calibrated_immunogenicity_log_odds"].tolist() == pytest.approx([-1.0, -5.5, 0.125])
    assert res["out_of_calibration_support"].tolist() == [False, True, True]


def test_apply_calibrator_preserves_original_columns_and_order(
    synthetic_calibrator, presentation_tsv, tmp_path
):
    # Secondary-signal guardrail: original columns + ROW ORDER unchanged (no re-rank).
    out = tmp_path / "presentation_calibrated.tsv"
    apply_calibrator(presentation_tsv, synthetic_calibrator, out)
    res = pd.read_csv(out, sep="\t")
    assert res["peptide"].tolist() == ["AAA", "BBB", "CCC"]
    assert res["presentation_class"].tolist() == ["weak", "non", "strong"]


def test_apply_calibrator_empty_input_writes_header(synthetic_calibrator, tmp_path):
    empty = tmp_path / "empty.tsv"
    pd.DataFrame({"genotype_presentation_score": pd.Series([], dtype=float)}).to_csv(
        empty, sep="\t", index=False
    )
    out = tmp_path / "out.tsv"
    apply_calibrator(empty, synthetic_calibrator, out)
    res = pd.read_csv(out, sep="\t")
    assert "calibrated_immunogenicity_log_odds" in res.columns
    assert "out_of_calibration_support" in res.columns
    assert len(res) == 0


# --------------------------------------------------------------------------- #
# Fail-loud guards
# --------------------------------------------------------------------------- #
def test_load_calibrator_knots_rejects_shape_mismatch(tmp_path):
    bad = tmp_path / "bad.joblib"
    joblib.dump({"cx": np.array([0.1, 0.5]), "cy": np.array([-5.0]), "version": "v1"}, bad)
    with pytest.raises(ValueError):
        load_calibrator_knots(bad)


def test_load_calibrator_knots_rejects_non_monotonic_cx(tmp_path):
    # np.interp silently returns wrong values on unsorted xp — guard must reject it.
    bad = tmp_path / "nonmono.joblib"
    joblib.dump(
        {"cx": np.array([0.5, 0.1, 0.9]), "cy": np.array(_CY), "version": "v1"}, bad
    )
    with pytest.raises(ValueError):
        load_calibrator_knots(bad)


def test_apply_calibrator_missing_score_col_raises(synthetic_calibrator, tmp_path):
    bad_tsv = tmp_path / "noscore.tsv"
    pd.DataFrame({"peptide": ["AAA"]}).to_csv(bad_tsv, sep="\t", index=False)
    with pytest.raises(KeyError):
        apply_calibrator(bad_tsv, synthetic_calibrator, tmp_path / "out.tsv")
