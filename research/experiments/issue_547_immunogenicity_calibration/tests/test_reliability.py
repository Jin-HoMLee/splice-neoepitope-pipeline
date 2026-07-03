"""Unit tests for reliability.jeffreys_ci (Jeffreys / boundary-corrected binomial CI).

Motivation (Issue #804): the reliability-plot error bars used a Wilson score
interval, which *under-covers* as the observed rate p -> 0. Our lowest-prevalence
held-out cohort (NCI, prevalence ~ 2.45e-4) sits exactly in that regime. The
Jeffreys interval (Beta(1/2, 1/2) prior; Brown, Cai & DasGupta, *Statist. Sci.*
2001) has coverage far closer to nominal near the boundary.

Correctness anchors:
  * Jeffreys is *defined* by Beta quantiles, so exact numeric literals (computed
    once from scipy.stats.beta) pin the formula against silent drift.
  * The Brown-Cai-DasGupta boundary rule: lower limit is exactly 0 when there are
    no successes, upper limit is exactly 1 when all trials are successes.
  * CI width is driven by the Kish *effective* n, not the weight-inflated raw
    count (the reason the helper takes n_eff separately).
  * A seeded Monte-Carlo coverage check reproduces the whole point of the change:
    at low p, Jeffreys reaches nominal coverage where Wilson under-covers.
"""
import numpy as np
import pytest

from reliability import jeffreys_ci


# ---- boundary behavior (Brown-Cai-DasGupta) --------------------------------

def test_no_successes_gives_zero_lower_bound():
    # x_eff = 0 -> lower limit pinned to exactly 0, finite positive upper.
    lo, hi = jeffreys_ci(0.0, 10.0, 10.0)
    assert lo == 0.0
    assert 0.0 < hi < 1.0
    assert hi == pytest.approx(0.21719626750921053, abs=1e-12)


def test_all_successes_gives_unit_upper_bound():
    # x_eff = n_eff -> upper limit pinned to exactly 1, positive lower.
    lo, hi = jeffreys_ci(5.0, 5.0, 5.0)
    assert hi == 1.0
    assert 0.0 < lo < 1.0
    assert lo == pytest.approx(0.6206228577009606, abs=1e-12)


# ---- known-answer literals (pin the Beta-quantile formula) ------------------

@pytest.mark.parametrize(
    "x, n, expected_lo, expected_hi",
    [
        (1.0, 10.0, 0.011011673763161135, 0.38131477106661626),
        (3.0, 50.0, 0.017186649071151135, 0.15153256302766024),
        (250.0, 1000.0, 0.22391109474754972, 0.27753363123662056),
    ],
)
def test_matches_beta_quantile_reference(x, n, expected_lo, expected_hi):
    # Unweighted case: n_eff == n, so x_eff == x.
    lo, hi = jeffreys_ci(x, n, n)
    assert lo == pytest.approx(expected_lo, abs=1e-12)
    assert hi == pytest.approx(expected_hi, abs=1e-12)


# ---- interval validity ------------------------------------------------------

@pytest.mark.parametrize("x, n", [(1, 20), (5, 20), (10, 20), (19, 20), (7, 100)])
def test_interval_is_ordered_and_within_unit(x, n):
    lo, hi = jeffreys_ci(float(x), float(n), float(n))
    p = x / n
    assert 0.0 <= lo <= p <= hi <= 1.0


def test_degenerate_inputs_return_uninformative_interval():
    assert jeffreys_ci(0.0, 0.0, 0.0) == (0.0, 1.0)
    assert jeffreys_ci(3.0, 10.0, 0.0) == (0.0, 1.0)
    assert jeffreys_ci(3.0, -1.0, 10.0) == (0.0, 1.0)


# ---- Kish effective-n drives the width, not the raw weighted count ----------

def test_effective_n_controls_width():
    # Same observed rate p = 0.1, but a larger effective n must give a strictly
    # narrower interval. This exercises the n_eff plumbing: the raw weighted
    # count (n_pos_w / n_total_w) fixes only p; n_eff sets the variance.
    lo_small, hi_small = jeffreys_ci(1.0, 10.0, n_eff=10.0)
    lo_large, hi_large = jeffreys_ci(1.0, 10.0, n_eff=1000.0)
    assert (hi_large - lo_large) < (hi_small - lo_small)


# ---- the reason for the change: better boundary coverage than Wilson --------

def _wilson_ci(x, n, z=1.96):
    """Reference Wilson score interval (the method being replaced)."""
    if n <= 0:
        return 0.0, 1.0
    p = x / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return centre - margin, centre + margin


def test_coverage_beats_wilson_near_zero():
    # NCI-like low-prevalence regime. Seeded so the assertion is deterministic.
    rng = np.random.default_rng(20260703)
    p_true, n, alpha, n_sim = 0.002, 500, 0.05, 20000
    xs = rng.binomial(n, p_true, n_sim)

    jeff_hits = wilson_hits = 0
    for x in xs:
        lo, hi = jeffreys_ci(float(x), float(n), float(n), alpha=alpha)
        jeff_hits += lo <= p_true <= hi
        lo, hi = _wilson_ci(x, n)
        wilson_hits += lo <= p_true <= hi

    jeff_cov = jeff_hits / n_sim
    wilson_cov = wilson_hits / n_sim

    # Jeffreys reaches nominal coverage; Wilson under-covers at the boundary.
    assert jeff_cov >= 0.95, f"Jeffreys under-covered: {jeff_cov:.4f}"
    assert wilson_cov < 0.95, f"Wilson unexpectedly covered: {wilson_cov:.4f}"
    assert jeff_cov > wilson_cov
