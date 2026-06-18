import json
from pathlib import Path
import numpy as np
import pytest
from sklearn.isotonic import IsotonicRegression
from calibrator import centered_isotonic

FIX = json.loads((Path(__file__).parent / "fixtures" / "cir_fixtures.json").read_text())

@pytest.mark.parametrize("name", list(FIX["cases"].keys()))
def test_matches_r_cir(name):
    c = FIX["cases"][name]
    x, y, w = np.array(c["x"], float), np.array(c["y"], float), np.array(c["w"], float)
    y_iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit_transform(x, y, sample_weight=w)
    cx, cy = centered_isotonic(x, y_iso, w)
    # the port must reproduce R cir's centred knots within tolerance
    np.testing.assert_allclose(cx, np.array(c["cir_x"], float), atol=1e-6)
    np.testing.assert_allclose(cy, np.array(c["cir_y"], float), atol=1e-6)

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
