import json
from pathlib import Path
import numpy as np
import pytest
from sklearn.isotonic import IsotonicRegression
from calibrator import centered_isotonic

FIX = json.loads((Path(__file__).parent / "fixtures" / "cir_fixtures.json").read_text())

@pytest.mark.parametrize("name", list(FIX["cases"].keys()))
def test_matches_r_cir(name):
    """Port reproduces R `cir`'s centred-isotonic *function* (live oracle).

    We compare the calibration function the calibrator actually consumes —
    ``np.interp(score, cx, cy)`` with endpoint clipping (calibrator.py:109) —
    NOT the exact knot encoding. R's `cir` (2.5.1) keeps boundary-anchor knots
    (first/last x) inside a collapsed level-set that the compact port drops;
    both define the identical interpolant once np.interp clips to endpoints
    (verified max|Δ| = 0 over a dense grid). Asserting exact knot arrays would
    test an arbitrary representation choice, not the curve — and would re-create
    the self-validating trap the old hand-computed fixtures fell into. Fixtures
    are now live R `cir` output (see fixtures/cir_fixtures.json r_version/cir_version).
    """
    c = FIX["cases"][name]
    x, y, w = np.array(c["x"], float), np.array(c["y"], float), np.array(c["w"], float)
    y_iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit_transform(x, y, sample_weight=w)
    cx, cy = centered_isotonic(x, y_iso, w)
    grid = np.linspace(float(x.min()), float(x.max()), 1001)
    got = np.interp(grid, cx, cy)
    want = np.interp(grid, np.array(c["cir_x"], float), np.array(c["cir_y"], float))
    np.testing.assert_allclose(got, want, atol=1e-6)

def test_centred_knots_are_monotone():
    x = np.array([1, 2, 3, 4, 5, 6.])
    y_iso = np.array([0.1, 0.2, 0.2, 0.2, 0.5, 0.9])
    # Asymmetric weights: plateau at y=0.2 spans x=2,3,4 with weights 2,5,1
    w = np.array([1., 2., 5., 1., 1., 1.])
    cx, cy = centered_isotonic(x, y_iso, w)
    assert np.all(np.diff(cy) >= 0)
    assert np.all(np.diff(cx) > 0)
    # Plateau centroid: (2*2 + 3*5 + 4*1)/(2+5+1) = (4+15+4)/8 = 23/8 = 2.875
    # knots: x=1 (singleton), x=2.875 (plateau), x=5 (singleton), x=6 (singleton)
    np.testing.assert_allclose(cx[1], 23.0 / 8.0, atol=1e-9)


def test_asymmetric_weighted_plateau():
    """Plateau with asymmetric weights: exact weighted centroid must be returned."""
    # Plateau at y=0.3 spans x=1,2,3 with weights 1,2,5
    # centroid = (1*1 + 2*2 + 3*5)/(1+2+5) = (1+4+15)/8 = 20/8 = 2.5
    # x=4 with y=0.9 is a singleton, stays at 4.0
    cx, cy = centered_isotonic(
        x=[1, 2, 3, 4],
        y_iso=[0.3, 0.3, 0.3, 0.9],
        w=[1, 2, 5, 1],
    )
    np.testing.assert_allclose(cx, np.array([2.5, 4.0]), atol=1e-9)
    np.testing.assert_allclose(cy, np.array([0.3, 0.9]), atol=1e-9)
