"""Immunogenicity calibrator: KDE -> density-ratio + prior -> isotonic ->
centered isotonic (Oron & Flournoy 2017). Reimplemented from the NeoGuider
paper (Wei et al. 2026); no NeoGuider code vendored."""
import numpy as np


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
        while j + 1 < n and np.isclose(y_iso[j + 1], y_iso[i]):
            j += 1
        sl = slice(i, j + 1)
        wsum = w[sl].sum()
        cx.append((x[sl] * w[sl]).sum() / wsum if wsum > 0 else x[sl].mean())
        cy.append(y_iso[i])
        i = j + 1
    return np.array(cx), np.array(cy)
