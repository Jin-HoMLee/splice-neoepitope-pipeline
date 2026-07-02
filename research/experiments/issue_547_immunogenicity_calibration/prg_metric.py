"""Area Under the Precision-Recall-Gain curve (AUPRG), weight-aware.

Method: Flach, P.A. & Kull, M. (2015), "Precision-Recall-Gain Curves: PR Analysis
Done Right", Advances in Neural Information Processing Systems 28 (NeurIPS).
Reference implementation: https://github.com/meeliskull/prg (clean-room reimplemented
here from the published definitions; validated against its canonical example, which
returns AUPRG = 0.432462103574 -- see tests/test_prg_metric.py).

Why we report it (Issue #803): the LOCO discrimination metric was AUPRC / prevalence
("lift"), a crude baseline correction with the interpolation + baseline artifacts PRG
was designed to fix. Precision-gain and recall-gain baseline precision and recall
against the prevalence on a harmonic (F-score) scale, so AUPRG = 0 for a random
(non-discriminating) ranker and 1 for a perfect ranker *independent of prevalence* --
values are directly comparable and meaningfully averaged across cohorts.

Weighting: the published implementation is unweighted. We generalise the contingency
counts (TP/FP/FN/TN) to sums of sample weights so the metric stays consistent with the
negative-reweighted AUPRC used in the same evaluation (compute_sample_weights).  With
unit weights it reduces to the reference exactly; integer weights equal physically
replicating the corresponding rows (tested).
"""

import numpy as np


def _precision_gain(tp, fp, ratio):
    """1 - (n_pos/n_neg) * (fp/tp).  ratio = n_pos/n_neg = pi/(1-pi)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 - ratio * (fp / tp)


def _recall_gain(tp, fn, ratio):
    """1 - (n_pos/n_neg) * (fn/tp)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 - ratio * (fn / tp)


def _weighted_segments(scores, labels, weights):
    """Group tied scores (descending) into segments, returning the weighted
    positive and negative mass falling on each distinct score."""
    order = np.argsort(-scores, kind="mergesort")
    s, lab, w = scores[order], labels[order], weights[order]
    seg_pos, seg_neg = [], []
    i, n = 0, len(s)
    while i < n:
        j, p, q = i, 0.0, 0.0
        while j < n and s[j] == s[i]:
            if lab[j] == 1:
                p += w[j]
            else:
                q += w[j]
            j += 1
        seg_pos.append(p)
        seg_neg.append(q)
        i = j
    return np.asarray(seg_pos), np.asarray(seg_neg)


def _prg_curve(labels, scores, weights):
    """Recall-gain / precision-gain at every threshold (origin prepended), with
    crossing points inserted where the curve crosses recall_gain = 0 and where it
    crosses precision_gain = 0 within the recall_gain >= 0 region.  Both crossings
    are placed by interpolating the contingency table linearly and recomputing the
    (hyperbolic) gain there -- matching the reference so the trapezoidal area is
    exact for the piecewise curve."""
    n_pos = weights[labels == 1].sum()
    n_neg = weights[labels == 0].sum()
    n = n_pos + n_neg
    ratio = n_pos / n_neg

    seg_pos, seg_neg = _weighted_segments(scores, labels, weights)
    tp = np.insert(np.cumsum(seg_pos), 0, 0.0)
    fp = np.insert(np.cumsum(seg_neg), 0, 0.0)
    fn = n_pos - tp

    rg = _recall_gain(tp, fn, ratio)
    pg = _precision_gain(tp, fp, ratio)

    tp, fp, rg, pg = _insert_recall_gain_crossing(tp, fp, rg, pg, n_pos, n_neg, n, ratio)
    tp, fp, rg, pg = _insert_precision_gain_crossings(tp, fp, rg, pg, n_pos, n_neg, n, ratio)
    return rg, pg


def _insert_recall_gain_crossing(tp, fp, rg, pg, n_pos, n_neg, n, ratio):
    """Insert the point where recall_gain = 0 (i.e. tp = n_pos^2 / n) so the area
    from recall_gain 0 up to the first non-negative operating point is captured."""
    nonneg = np.where(rg >= 0)[0]
    if len(nonneg) == 0 or rg[nonneg[0]] == 0:
        return tp, fp, rg, pg
    j = nonneg[0]
    tp_star = n_pos * n_pos / n
    dtp = tp[j] - tp[j - 1]
    if dtp <= 0:
        return tp, fp, rg, pg
    alpha = (tp_star - tp[j - 1]) / dtp
    fp_star = fp[j - 1] + alpha * (fp[j] - fp[j - 1])
    tp = np.insert(tp, j, tp_star)
    fp = np.insert(fp, j, fp_star)
    rg = np.insert(rg, j, 0.0)
    pg = np.insert(pg, j, _precision_gain(tp_star, fp_star, ratio))
    return tp, fp, rg, pg


def _insert_precision_gain_crossings(tp, fp, rg, pg, n_pos, n_neg, n, ratio):
    """Insert points where precision_gain crosses 0 within recall_gain >= 0, so the
    piecewise-linear trapezoid follows the true (hyperbolic) curve through the axis."""
    i = 1
    while i < len(pg):
        if rg[i - 1] >= 0 and np.isfinite(pg[i - 1]) and np.isfinite(pg[i]) and pg[i - 1] * pg[i] < 0:
            # precision_gain = 0  <=>  fp = tp * n_neg / n_pos.  Interpolate the
            # contingency point linearly to that condition.
            dtp = tp[i] - tp[i - 1]
            dfp = fp[i] - fp[i - 1]
            denom = dfp - (n_neg / n_pos) * dtp
            if denom != 0:
                alpha = ((n_neg / n_pos) * tp[i - 1] - fp[i - 1]) / denom
                if 0.0 < alpha < 1.0:
                    tp_c = tp[i - 1] + alpha * dtp
                    fp_c = fp[i - 1] + alpha * dfp
                    tp = np.insert(tp, i, tp_c)
                    fp = np.insert(fp, i, fp_c)
                    rg = np.insert(rg, i, _recall_gain(tp_c, n_pos - tp_c, ratio))
                    pg = np.insert(pg, i, 0.0)
                    i += 1
        i += 1
    return tp, fp, rg, pg


def auprg(labels, scores, weights=None):
    """Area under the Precision-Recall-Gain curve.

    Parameters
    ----------
    labels : array of {0, 1}
    scores : array; higher score = ranked more positive.
    weights : optional per-sample weights (default all 1.0).

    Returns
    -------
    float
        AUPRG.  1.0 = perfect ranker, 0.0 = random / non-discriminating,
        < 0 = worse than random.  ``nan`` if ``labels`` is single-class.
    """
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    n_lab = labels.sum()
    if n_lab == 0 or n_lab == len(labels):
        return float("nan")
    if weights is None:
        weights = np.ones(len(labels), dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)

    rg, pg = _prg_curve(labels, scores, weights)
    area = 0.0
    for i in range(1, len(rg)):
        if np.isnan(rg[i - 1]) or rg[i - 1] < 0:
            continue
        width = rg[i] - rg[i - 1]
        height = (pg[i] + pg[i - 1]) / 2.0
        area += width * height
    return float(area)
