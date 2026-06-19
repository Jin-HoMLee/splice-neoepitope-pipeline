"""Immunogenicity calibrator: KDE -> density-ratio + prior -> isotonic ->
centered isotonic (Oron & Flournoy 2017). Reimplemented from the NeoGuider
paper (Wei et al. 2026); no NeoGuider code vendored."""
import joblib
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.isotonic import IsotonicRegression


def centered_isotonic(x, y_iso, w):
    """Centered isotonic regression (Oron & Flournoy 2017).

    Collapse each flat level-set of the isotonic fit `y_iso` to a single knot
    whose x is the weight-centroid of the pooled level-set points; keep the
    isotonic y. Returns sorted (centre_x, centre_y).
    """
    x = np.asarray(x, float)
    y_iso = np.asarray(y_iso, float)
    w = np.asarray(w, float)
    order = np.argsort(x)
    x, y_iso, w = x[order], y_iso[order], w[order]
    cx, cy = [], []
    i, n = 0, len(x)
    while i < n:
        j = i
        # atol-only comparison: PAVA level-sets share an *identical* value, so match on
        # absolute equality. The np.isclose default (rtol=1e-5) could merge two genuinely
        # distinct adjacent levels that happen to be within 1e-5 relative.
        while j + 1 < n and np.isclose(y_iso[j + 1], y_iso[i], rtol=0, atol=1e-10):
            j += 1
        sl = slice(i, j + 1)
        wsum = w[sl].sum()
        cx.append((x[sl] * w[sl]).sum() / wsum if wsum > 0 else x[sl].mean())
        cy.append(y_iso[i])
        i = j + 1
    return np.array(cx), np.array(cy)


_EPS = 1e-12


def fixed_kde(samples):
    """Fixed-bandwidth Gaussian KDE (Scott's rule) — the baseline."""
    kde = gaussian_kde(np.asarray(samples, float))
    return lambda grid: kde(np.atleast_1d(np.asarray(grid, float)))


def adaptive_kde(samples):
    """Abramson sample-point variable-bandwidth Gaussian KDE.

    Local bandwidth ∝ 1/sqrt(pilot density), normalised to the geometric mean.
    """
    s = np.asarray(samples, float)
    pilot = gaussian_kde(s)
    f_pilot = np.clip(pilot(s), _EPS, None)
    g = np.exp(np.mean(np.log(f_pilot)))
    lam = np.sqrt(g / f_pilot)                 # per-sample local factor
    std = np.std(s) if len(s) > 1 else 1.0
    bw = pilot.factor * std * lam              # per-sample bandwidth
    bw = np.clip(bw, _EPS, None)

    def density(grid):
        x = np.atleast_1d(np.asarray(grid, float))
        u = (x[:, None] - s[None, :]) / bw[None, :]
        k = np.exp(-0.5 * u ** 2) / (np.sqrt(2 * np.pi) * bw[None, :])
        return k.mean(axis=1)
    return density


class PresentationCalibrator:
    """Calibrate genotype_presentation_score -> calibrated_immunogenicity_log_odds."""

    CALIBRATOR_VERSION = "v1"

    def __init__(self, kde_mode="adaptive", n_grid=512):
        if kde_mode not in ("adaptive", "fixed"):
            raise ValueError("kde_mode must be 'adaptive' or 'fixed'")
        self.kde_mode = kde_mode
        self.n_grid = n_grid

    def fit(self, scores, labels, n_pos_true=None, n_neg_true=None, fit_cohorts=None):
        scores = np.asarray(scores, float)
        labels = np.asarray(labels, int)
        pos, neg = scores[labels == 1], scores[labels == 0]
        if len(pos) < 2 or len(neg) < 2:
            raise ValueError("need >=2 samples per class")
        n_pos_true = len(pos) if n_pos_true is None else n_pos_true
        n_neg_true = len(neg) if n_neg_true is None else n_neg_true
        self.prior_ = float(np.log(n_pos_true / n_neg_true))

        make = adaptive_kde if self.kde_mode == "adaptive" else fixed_kde
        kde_pos, kde_neg = make(pos), make(neg)

        lo, hi = float(scores.min()), float(scores.max())
        self.score_range_ = (lo, hi)
        grid = np.linspace(lo, hi, self.n_grid)
        d_pos = np.clip(kde_pos(grid), _EPS, None)
        d_neg = np.clip(kde_neg(grid), _EPS, None)
        raw = np.log(d_pos) - np.log(d_neg) + self.prior_         # log-odds(grid)
        # local data support as isotonic weights (true-count-scaled densities)
        w = n_pos_true * d_pos + n_neg_true * d_neg
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
        y_iso = iso.fit_transform(grid, raw, sample_weight=w)
        self.cx_, self.cy_ = centered_isotonic(grid, y_iso, w)
        self.fit_cohorts_ = list(fit_cohorts) if fit_cohorts is not None else None
        return self

    def transform(self, scores):
        if not hasattr(self, "cx_"):
            raise RuntimeError("PresentationCalibrator must be fit before transform()")
        scores = np.atleast_1d(np.asarray(scores, float))
        return np.interp(scores, self.cx_, self.cy_)  # np.interp clips to endpoints

    def save(self, path):
        joblib.dump({
            "cx": self.cx_, "cy": self.cy_, "prior": self.prior_,
            "score_range": self.score_range_, "kde_mode": self.kde_mode,
            "n_grid": self.n_grid,
            "fit_cohorts": self.fit_cohorts_, "version": self.CALIBRATOR_VERSION,
        }, path)

    @classmethod
    def load(cls, path):
        d = joblib.load(path)
        # n_grid: default 512 for artifacts saved before it was serialised (transform()
        # uses only cx_/cy_, so older artifacts round-trip identically).
        obj = cls(kde_mode=d["kde_mode"], n_grid=d.get("n_grid", 512))
        obj.cx_, obj.cy_ = d["cx"], d["cy"]
        obj.prior_ = d["prior"]
        obj.score_range_ = tuple(d["score_range"])
        obj.fit_cohorts_ = d["fit_cohorts"]
        return obj
