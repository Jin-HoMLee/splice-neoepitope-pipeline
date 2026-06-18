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
    w = np.ones_like(x)
    cx, cy = centered_isotonic(x, y_iso, w)
    assert np.all(np.diff(cy) >= 0)
    assert np.all(np.diff(cx) > 0)
