import numpy as np
import pytest
from calibrator import PresentationCalibrator

def _synthetic(seed=0, n_pos=200, n_neg=4000):
    rng = np.random.default_rng(seed)
    pos = np.clip(rng.beta(6.0, 2.2, n_pos), 0, 1)   # immunogenic skew high
    neg = np.clip(rng.beta(2.0, 4.0, n_neg), 0, 1)   # non-immunogenic skew low
    scores = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(int)
    return scores, labels

@pytest.mark.parametrize("mode", ["adaptive", "fixed"])
def test_transform_is_monotone_nondecreasing(mode):
    scores, labels = _synthetic()
    cal = PresentationCalibrator(kde_mode=mode).fit(scores, labels)
    grid = np.linspace(0, 1, 200)
    out = cal.transform(grid)
    assert np.all(np.diff(out) >= -1e-9)

def test_higher_score_higher_logodds():
    scores, labels = _synthetic()
    cal = PresentationCalibrator().fit(scores, labels)
    assert cal.transform([0.9])[0] > cal.transform([0.1])[0]

def test_out_of_range_clips_to_boundary():
    scores, labels = _synthetic()
    cal = PresentationCalibrator().fit(scores, labels)
    lo, hi = cal.score_range_
    assert cal.transform([hi + 5.0])[0] == pytest.approx(cal.transform([hi])[0])
    assert cal.transform([lo - 5.0])[0] == pytest.approx(cal.transform([lo])[0])

def test_prior_uses_true_counts_not_subsample():
    scores, labels = _synthetic(n_pos=200, n_neg=4000)
    # same data, but declare the TRUE base rate as much rarer
    cal_sub = PresentationCalibrator().fit(scores, labels)
    cal_true = PresentationCalibrator().fit(scores, labels, n_pos_true=200, n_neg_true=400000)
    # rarer true prior => uniformly lower log-odds intercept
    g = np.linspace(0.2, 0.8, 50)
    assert np.all(cal_true.transform(g) < cal_sub.transform(g))
