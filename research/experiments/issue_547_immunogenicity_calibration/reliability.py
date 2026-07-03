"""Reliability-plot confidence intervals for the immunogenicity calibrator.

The reliability diagram (predicted log-odds vs empirical logit of the observed
rate) carries per-bin error bars. Those CIs are a *visual diagnostic* on the
held-out SNV cohorts, not load-bearing for any deployment decision (#592).

We use the **Jeffreys** interval (equal-tailed Beta(1/2, 1/2) credible interval;
Brown, Cai & DasGupta, "Interval Estimation for a Binomial Proportion",
*Statist. Sci.* 16(2), 2001) rather than the Wilson score interval. Wilson
under-covers as the observed rate p -> 0, and our lowest-prevalence held-out
cohort (NCI, prevalence ~ 2.45e-4) sits exactly in that boundary regime where
Jeffreys has coverage far closer to nominal (Issue #804).

Width uses the Kish (1965) effective sample size n_eff = (Sum w)^2 / Sum(w^2),
passed in separately: the large negative-reweighting weights inflate Sum(w)
thousands-fold while n_eff stays near the *sampled* n, so a naive weighted count
would badly understate the variance.
"""
from scipy.stats import beta


def jeffreys_ci(n_pos_w, n_total_w, n_eff, alpha=0.05):
    """Jeffreys credible interval for a binomial proportion on the effective-n scale.

    Parameters
    ----------
    n_pos_w, n_total_w
        Weighted positive and total counts in the bin; their ratio is the observed
        rate p. (Only the ratio matters - the raw weighted magnitudes, which the
        negative-reweighting inflates, do not set the interval width.)
    n_eff
        Kish effective sample size for the bin; this alone sets the CI width.
    alpha
        Two-sided miscoverage (default 0.05 -> a 95% interval).

    Returns
    -------
    (lo, hi) : tuple of float
        The Jeffreys interval, with the Brown-Cai-DasGupta boundary rule applied:
        the lower limit is exactly 0 when there are no effective successes and the
        upper limit is exactly 1 when every effective trial is a success.

    Degenerate inputs (non-positive total weight or effective n) return the
    uninformative interval (0.0, 1.0).
    """
    if n_total_w <= 0 or n_eff <= 0:
        return 0.0, 1.0

    p = n_pos_w / n_total_w
    p = min(max(p, 0.0), 1.0)

    # Effective successes / failures on the Kish scale, plus the Jeffreys Beta(1/2)
    # prior pseudo-counts.
    x_eff = p * n_eff
    a = x_eff + 0.5
    b = (n_eff - x_eff) + 0.5

    lo = 0.0 if x_eff <= 0.0 else float(beta.ppf(alpha / 2, a, b))
    hi = 1.0 if x_eff >= n_eff else float(beta.ppf(1 - alpha / 2, a, b))
    return lo, hi
