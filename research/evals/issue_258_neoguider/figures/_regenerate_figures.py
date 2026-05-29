"""Regenerate the deck-only figures for the NeoGuider KDE + centered-isotonic
backup (appendix) slides — Issue #258.

Run from repo root with any interpreter carrying numpy/scipy/sklearn/matplotlib
(the research venv, or the pyenv global that has the notebook stack):

    ~/.pyenv/shims/python research/evals/issue_258_neoguider/figures/_regenerate_figures.py

Outputs (next to this script):
    kde_densities.png     — two class-conditional KDEs (immunogenic vs not) over presentation_score
    density_ratio.png     — log density-ratio (= per-feature log-odds shape) vs score
    imbalance_prior.png    — class imbalance enters ONLY as an additive vertical shift (prior log-odds)
    isotonic_vs_cir.png   — the core: PAVA step function (plateau) vs centered-isotonic smooth curve
    calibrated_curve.png  — final calibrated log-odds curve + how a NEW peptide's score maps
    combine_schematic.png — block diagram: per-feature CIR log-odds → logistic regression → rank
    integration_map.png   — pipeline integration map (NeoGuider calibration layer position)

These are PEDAGOGICAL figures (synthetic data, fixed seed; no pipeline data is read).
The isotonic + centered-isotonic transforms are faithful — sklearn IsotonicRegression and a
from-scratch centered-isotonic-regression (Oron & Flournoy 2017) implementation. The KDE is the
one deliberate simplification: NeoGuider uses ONE shared *adaptive*-bandwidth kernel to jointly
estimate both class densities (bandwidth tracks the spacing of consecutive positive examples,
then is iteratively refined); here we draw two independent fixed-bandwidth scipy.gaussian_kde
curves for visual clarity, and the slide says so. Method source: Wei et al., Genome Medicine
2025, doi:10.1186/s13073-025-01592-9.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.isotonic import IsotonicRegression

FIGURES_DIR = Path(__file__).resolve().parent

# --- house style (legible on a 1280x720 reveal.js slide) -------------------
plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 15,
    "axes.titlesize": 17,
    "axes.labelsize": 15,
    "legend.fontsize": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
C_POS = "#c44e52"   # immunogenic (positive) — warm red
C_NEG = "#4c72b0"   # non-immunogenic (negative) — steel blue
C_IR = "#8c8c8c"    # isotonic step — grey
C_CIR = "#b8860b"   # centered isotonic — gold (matches the deck's calibration node)

RNG = np.random.default_rng(258)


# --- shared synthetic model -------------------------------------------------
def make_data():
    """Synthetic labelled cohort: presentation_score in [0,1] for two classes.

    Positives skew high, negatives skew low-to-mid, with a deliberate overlap zone.
    Mild 1:6 ratio for legible KDEs (the real regime is ~1:100; the imbalance
    figure handles the base-rate point analytically).
    """
    pos = np.clip(RNG.beta(6.0, 2.2, size=70), 0, 1)        # immunogenic, skew high
    neg = np.clip(RNG.beta(2.0, 4.0, size=420), 0, 1)        # non-immunogenic, skew low
    return pos, neg


def class_kdes(pos, neg):
    return gaussian_kde(pos), gaussian_kde(neg)


def binned_empirical_logodds(pos, neg, n_bins=12):
    """Per-bin empirical log-odds with a continuity correction → naturally noisy,
    so PAVA has a real monotonicity violation to pool into a plateau."""
    edges = np.linspace(0, 1, n_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    pos_c, _ = np.histogram(pos, bins=edges)
    neg_c, _ = np.histogram(neg, bins=edges)
    # keep bins with any data
    keep = (pos_c + neg_c) > 0
    mids, pos_c, neg_c = mids[keep], pos_c[keep], neg_c[keep]
    logodds = np.log((pos_c + 0.5) / (neg_c + 0.5))
    weight = pos_c + neg_c
    return mids, logodds, weight


def centered_isotonic(x, y_iso, w):
    """Centered isotonic regression (Oron & Flournoy 2017).

    Collapse each flat level-set of the isotonic fit to a single point whose
    x is the sample-size-weighted centroid of the pooled design points; keep its
    isotonic y; linearly interpolate between those centred points.
    Returns (centre_x, centre_y).
    """
    x, y_iso, w = np.asarray(x), np.asarray(y_iso), np.asarray(w)
    cx, cy = [], []
    i, n = 0, len(x)
    while i < n:
        j = i
        while j + 1 < n and np.isclose(y_iso[j + 1], y_iso[i]):
            j += 1
        sl = slice(i, j + 1)
        wsum = w[sl].sum()
        cx.append((x[sl] * w[sl]).sum() / wsum)
        cy.append(y_iso[i])
        i = j + 1
    return np.array(cx), np.array(cy)


# --- figure 1: two class-conditional KDEs -----------------------------------
def fig_kde_densities(pos, neg, kpos, kneg):
    grid = np.linspace(0, 1, 400)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.fill_between(grid, kneg(grid), color=C_NEG, alpha=0.30)
    ax.fill_between(grid, kpos(grid), color=C_POS, alpha=0.30)
    ax.plot(grid, kneg(grid), color=C_NEG, lw=2.4, label=r"non-immunogenic  $p(s\,|\,\mathrm{neg})$")
    ax.plot(grid, kpos(grid), color=C_POS, lw=2.4, label=r"immunogenic  $p(s\,|\,\mathrm{pos})$")
    # rug
    ax.plot(neg, np.full_like(neg, -0.06), "|", color=C_NEG, alpha=0.5, markersize=8)
    ax.plot(pos, np.full_like(pos, -0.16), "|", color=C_POS, alpha=0.6, markersize=8)
    ax.axvspan(0.45, 0.85, color="0.5", alpha=0.10)
    ax.text(0.65, ax.get_ylim()[1] * 0.62, "overlap zone", ha="center", color="0.4", fontsize=12, style="italic")
    ax.set_xlabel("MHCflurry presentation_score")
    ax.set_ylabel("class-conditional density")
    ax.set_title("Step 1a — class-conditional densities (illustrative)")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig(FIGURES_DIR / "kde_densities.png")
    plt.close(fig)


# --- figure 2: log density ratio = per-feature log-odds shape ----------------
def support_grid(kpos, kneg, thresh=1e-3):
    """Grid restricted to where BOTH class densities have real support.

    The density ratio is only statistically meaningful where both KDEs carry
    mass. Outside that region one density underflows and the log-ratio explodes
    (the old eps-on-both-terms guard masked this by bending the curve, which
    silently broke monotonicity). Masking to the joint-support region keeps the
    plotted log-odds well-behaved and monotone, as the method intends.
    """
    grid = np.linspace(0.0, 1.0, 400)
    dp, dn = kpos(grid), kneg(grid)
    keep = (dp > thresh) & (dn > thresh)
    return grid[keep], dp[keep], dn[keep]


def fig_density_ratio(kpos, kneg):
    grid, dp, dn = support_grid(kpos, kneg)
    log_lr = np.log(dp / dn)
    full = np.linspace(0, 1, 400)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))

    # --- faint original class densities (secondary axis) — the visual bridge to
    #     Step 1a: the gold curve below is literally log(red / blue). ---------------
    ax2 = ax.twinx()
    ax2.fill_between(full, kpos(full), color=C_POS, alpha=0.06, zorder=0)
    ax2.fill_between(full, kneg(full), color=C_NEG, alpha=0.06, zorder=0)
    ax2.plot(full, kpos(full), color=C_POS, lw=1.3, ls="--", alpha=0.32, zorder=1)
    ax2.plot(full, kneg(full), color=C_NEG, lw=1.3, ls="--", alpha=0.32, zorder=1)
    ax2.set_ylim(0, max(kpos(full).max(), kneg(full).max()) * 1.05)
    ax2.set_ylabel("class density  (faint — from Step 1a)", color="0.6", fontsize=10.5)
    ax2.tick_params(axis="y", labelcolor="0.65", labelsize=9.5)
    ax2.text(0.31, kneg(0.31) * 0.45, r"$p(s\,|\,\mathrm{neg})$", color=C_NEG, alpha=0.65, fontsize=11)
    ax2.text(0.82, kpos(0.82) * 0.5, r"$p(s\,|\,\mathrm{pos})$", color=C_POS, alpha=0.65, fontsize=11)

    # --- the log-odds curve (primary axis, raised above the faint densities) ------
    ax.axhline(0, color="0.7", lw=1, zorder=2)
    ax.fill_between(grid, 0, log_lr, where=(log_lr >= 0), color=C_POS, alpha=0.08, zorder=2)
    ax.fill_between(grid, 0, log_lr, where=(log_lr < 0), color=C_NEG, alpha=0.08, zorder=2)
    ax.plot(grid, log_lr, color=C_CIR, lw=2.8, zorder=4)
    ax.text(grid.max() - 0.05, max(log_lr) * 0.5, "evidence FOR\nimmunogenicity", color=C_POS, fontsize=11.5, ha="center", zorder=5)
    ax.text(grid.min() + 0.05, min(log_lr) * 0.5, "evidence\nAGAINST", color=C_NEG, fontsize=11.5, ha="center", zorder=5)
    ax.set_xlabel("MHCflurry presentation_score")
    ax.set_ylabel(r"$\log\,\dfrac{p(s\,|\,\mathrm{pos})}{p(s\,|\,\mathrm{neg})}$  (log-odds, pre-prior)")
    ax.set_title("Step 1b — log-odds = log of the two densities' ratio")
    ax.set_xlim(0, 1)
    ax.set_zorder(ax2.get_zorder() + 1)   # gold curve on top of the faint densities
    ax.patch.set_visible(False)           # let ax2 show through ax's background
    fig.savefig(FIGURES_DIR / "density_ratio.png")
    plt.close(fig)


# --- figure 3: class imbalance is just an additive vertical shift ------------
def fig_imbalance_prior(kpos, kneg):
    grid, dp, dn = support_grid(kpos, kneg)
    log_lr = np.log(dp / dn)
    prior_balanced = 0.0
    prior_imbal = np.log(1 / 100)  # ~ -4.6
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.axhline(0, color="0.8", lw=1)
    ax.plot(grid, log_lr + prior_balanced, color="#55a868", lw=2.6, label="balanced cohort  (prior log-odds = 0)")
    ax.plot(grid, log_lr + prior_imbal, color=C_CIR, lw=2.6, label="1:100 imbalance  (prior log-odds ≈ −4.6)")
    # annotate the constant gap
    xa = 0.7
    ya0 = np.interp(xa, grid, log_lr)
    ax.annotate("", xy=(xa, ya0), xytext=(xa, ya0 + prior_imbal),
                arrowprops=dict(arrowstyle="<->", color="0.4"))
    ax.text(xa + 0.015, ya0 + prior_imbal / 2, "same constant\nshift everywhere", fontsize=11.5, color="0.3", va="center")
    ax.set_xlabel("MHCflurry presentation_score")
    ax.set_ylabel("posterior log-odds of immunogenicity")
    ax.set_title("Imbalance shifts the curve — it never bends it")
    ax.set_xlim(0, 1)
    ax.legend(loc="lower right", frameon=False)
    fig.savefig(FIGURES_DIR / "imbalance_prior.png")
    plt.close(fig)


# --- figure 4: isotonic step vs centered-isotonic curve (THE core) -----------
def fig_isotonic_vs_cir(pos, neg):
    x, y, w = binned_empirical_logodds(pos, neg)
    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    y_iso = ir.fit_transform(x, y, sample_weight=w)
    cx, cy = centered_isotonic(x, y_iso, w)

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    # raw noisy empirical points
    ax.scatter(x, y, s=18 + w, color="0.45", alpha=0.7, zorder=3,
               label="raw empirical log-odds (noisy, non-monotone)")
    # isotonic step function (where="post": the level-set value holds until the
    # next design point, i.e. a constant run over each pooled interval)
    ax.step(x, y_iso, where="post", color=C_IR, lw=2.4, zorder=4,
            label="isotonic regression (PAVA) — step function")
    # highlight every plateau (consecutive equal y_iso)
    for i in range(len(x) - 1):
        if np.isclose(y_iso[i], y_iso[i + 1]):
            ax.axvspan(x[i] - 0.005, x[i + 1] + 0.005, color=C_IR, alpha=0.10, zorder=1)
    # CIR curve
    ax.plot(cx, cy, "-o", color=C_CIR, lw=2.8, ms=7, zorder=5,
            label="centered isotonic (CIR) — centroid + interpolate")

    # annotate the DEFINING CIR move on the first (widest) plateau:
    # the whole flat run collapses to ONE node at its sample-weighted x-centroid.
    plateau_start = next((i for i in range(len(x) - 1) if np.isclose(y_iso[i], y_iso[i + 1])), None)
    if plateau_start is not None:
        plateau_end = plateau_start
        while plateau_end + 1 < len(x) and np.isclose(y_iso[plateau_end + 1], y_iso[plateau_start]):
            plateau_end += 1
        sl = slice(plateau_start, plateau_end + 1)
        yp = y_iso[plateau_start]
        cxp = (x[sl] * w[sl]).sum() / w[sl].sum()   # the centroid this plateau collapses to
        xmid = 0.5 * (x[plateau_start] + x[plateau_end])
        # label the flat run (above the span)
        ax.annotate("flat plateau (ties → unrankable)", xy=(xmid, yp),
                    xytext=(xmid, yp + 2.0), fontsize=10.5, color=C_IR, ha="center",
                    arrowprops=dict(arrowstyle="->", color=C_IR))
        # ring + arrow to the single gold centroid node it collapses to
        ax.plot(cxp, yp, "o", color=C_CIR, ms=15, mfc="none", mew=2.4, zorder=6)
        ax.annotate("collapses to ONE node\nat the sample-weighted centroid",
                    xy=(cxp, yp), xytext=(cxp + 0.07, yp + 0.9),
                    fontsize=10.5, color=C_CIR, ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=C_CIR, lw=1.9))
    ax.set_xlabel("MHCflurry presentation_score")
    ax.set_ylabel("per-feature immunogenicity log-odds")
    ax.set_title("Step 2 — isotonic (staircase) → centered isotonic (smooth)")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", frameon=False, fontsize=11)
    fig.savefig(FIGURES_DIR / "isotonic_vs_cir.png")
    plt.close(fig)


# --- figure 5: final calibrated curve + new-peptide lookup -------------------
def fig_calibrated_curve(pos, neg):
    x, y, w = binned_empirical_logodds(pos, neg)
    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    y_iso = ir.fit_transform(x, y, sample_weight=w)
    cx, cy = centered_isotonic(x, y_iso, w)

    def calibrate(s):
        return np.interp(s, cx, cy)

    grid = np.linspace(cx.min(), cx.max(), 300)
    fig, ax = plt.subplots(figsize=(8.4, 4.7))
    ax.axhline(0, color="0.8", lw=1)
    ax.plot(grid, calibrate(grid), color=C_CIR, lw=3, zorder=3)
    ax.plot(cx, cy, "o", color=C_CIR, ms=6, zorder=4)
    # two new peptides — labels parked in the empty upper-left, clear of the title
    for s_new, lab, xytext in [(0.66, "peptide A", (0.31, 0.7)), (0.83, "peptide B", (0.31, 2.0))]:
        lo = float(calibrate(s_new))
        ax.plot([s_new, s_new], [ax.get_ylim()[0], lo], "--", color="0.45", lw=1.4, zorder=2)
        ax.plot([grid.min(), s_new], [lo, lo], "--", color="0.45", lw=1.4, zorder=2)
        ax.plot(s_new, lo, "o", color=C_POS, ms=9, zorder=5)
        ax.annotate(f"{lab}: score {s_new}\n→ log-odds {lo:+.2f}", xy=(s_new, lo),
                    xytext=xytext, fontsize=11.5, color=C_POS,
                    arrowprops=dict(arrowstyle="->", color=C_POS))
    ax.set_xlabel("MHCflurry presentation_score (new peptide)")
    ax.set_ylabel("calibrated immunogenicity log-odds")
    ax.set_title("Step 3a — look up a new peptide on the calibrated curve")
    fig.savefig(FIGURES_DIR / "calibrated_curve.png")
    plt.close(fig)


# --- figure 6: combine per-feature log-odds via logistic regression ----------
def fig_combine_schematic():
    """Block diagram: per-feature CIR log-odds → logistic regression → rank.

    No data — makes the headline correction visible: the per-feature log-odds are
    INPUTS to a logistic regression that LEARNS per-feature weights, NOT summed.
    (A naive-Bayes sum is just the degenerate all-weights-equal special case.)
    """
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    feats = [
        ("presentation_score", "+1.3"),
        ("agretopicity", "-0.4"),
        ("foreignness", "+0.9"),
        ("… more features", ""),
    ]
    weights = [r"$w_1$", r"$w_2$", r"$w_3$", r"$w_k$"]

    fig, ax = plt.subplots(figsize=(9.4, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # per-feature CIR log-odds boxes (left)
    box_x, box_w, box_h = 0.3, 3.05, 0.95
    ys = [5.6, 4.2, 2.8, 1.4]
    for (name, val), y in zip(feats, ys):
        fc = "#f3ede0" if val else "#efefef"
        ax.add_patch(FancyBboxPatch((box_x, y - box_h / 2), box_w, box_h,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=fc, ec=C_CIR if val else "0.6", lw=1.6))
        label = f"{name}\nCIR log-odds {val}" if val else name
        ax.text(box_x + box_w / 2, y, label, ha="center", va="center", fontsize=11, color="0.15")

    # combiner node (center)
    comb_x, comb_y, comb_w, comb_h = 4.8, 3.5, 2.5, 2.1
    ax.add_patch(FancyBboxPatch((comb_x, comb_y - comb_h / 2), comb_w, comb_h,
                                boxstyle="round,pad=0.03,rounding_size=0.15",
                                fc="#ffe6a3", ec="#b8860b", lw=2.4))
    ax.text(comb_x + comb_w / 2, comb_y + 0.4, "logistic\nregression", ha="center", va="center",
            fontsize=14, fontweight="bold", color="#7a5c00")
    ax.text(comb_x + comb_w / 2, comb_y - 0.55, r"learns weights $w_1 \dots w_k$",
            ha="center", va="center", fontsize=10.5, color="#7a5c00")

    # arrows feature → combiner, labelled with the learned weights
    for wlbl, y in zip(weights, ys):
        ytip = comb_y + (y - comb_y) * 0.30
        ax.add_patch(FancyArrowPatch((box_x + box_w, y), (comb_x, ytip),
                                     arrowstyle="-|>", mutation_scale=15, color="0.45", lw=1.5))
        ax.text((box_x + box_w + comb_x) / 2, (y + ytip) / 2 + 0.16, wlbl,
                fontsize=12, color="#b8860b", ha="center")

    # output arrow → rank
    ax.add_patch(FancyArrowPatch((comb_x + comb_w, comb_y), (8.55, comb_y),
                                 arrowstyle="-|>", mutation_scale=18, color="0.3", lw=1.8))
    ax.add_patch(FancyBboxPatch((8.6, comb_y - 0.75), 1.25, 1.5,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#d7ecd9", ec="#55a868", lw=1.8))
    ax.text(8.6 + 0.625, comb_y, "immuno-\ngenicity\nrank", ha="center", va="center",
            fontsize=11, color="#2f5d3a")

    ax.text(5.0, 0.3,
            r"Not a naive-Bayes sum (that is just the special case $w_i = 1$) — the weights are learned.",
            ha="center", va="center", fontsize=11.5, style="italic", color="0.3")
    ax.set_title("Step 3b — combine per-feature log-odds into one ranking", fontsize=17)
    fig.savefig(FIGURES_DIR / "combine_schematic.png")
    plt.close(fig)


# --- figure 7: integration map (replaces the main-deck mermaid flowchart) ----
def fig_integration_map():
    """Where NeoGuider's calibration layer plugs into our scoring stack.

    Matplotlib reimplementation of the 'Where it would plug in' mermaid flowchart,
    in the deck's house style. Process stages are boxes; the data passed between
    stages rides on the arrows (so the two NeoGuider interfaces stay visible
    without two extra boxes). No data — a block diagram.
    """
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    NF, NE = "#eef1f5", "#9aa6b2"   # neutral fill / edge — our existing stages
    boxes = [
        ("Tumor-exclusive\njunctions", NF, NE),
        ("Peptide assembly\n+ flanking", NF, NE),
        ("MHCflurry 2.x\npresentation", NF, NE),
        ("NeoGuider\nrank-calibration\n(KDE + CIR)", "#ffe6a3", "#b8860b"),
        ("TCRdock\n(top candidate)", "#a3d8ff", "#1e5ba8"),
    ]
    pitch, bw = 3.3, 2.18
    centers = [1.35 + i * pitch for i in range(len(boxes))]
    y = 2.0

    fig, ax = plt.subplots(figsize=(12.8, 2.95))
    ax.set_xlim(0, centers[-1] + 1.35)
    ax.set_ylim(0, 3.55)
    ax.axis("off")

    for (label, fc, ec), cx in zip(boxes, centers):
        h = 1.5 if label.count("\n") >= 2 else 1.12
        lw = 2.9 if fc == "#ffe6a3" else 1.6
        ax.add_patch(FancyBboxPatch((cx - bw / 2, y - h / 2), bw, h,
                                    boxstyle="round,pad=0.03,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=lw, zorder=3))
        ax.text(cx, y, label, ha="center", va="center", fontsize=13, color="0.12", zorder=4)

    # the data interface NeoGuider consumes / emits, labelled on its two arrows
    arrow_labels = {2: "presentation_score\ngenotype_presentation_score",
                    3: "calibrated\nlog-odds"}
    for i, (cx0, cx1) in enumerate(zip(centers[:-1], centers[1:])):
        ax.add_patch(FancyArrowPatch((cx0 + bw / 2, y), (cx1 - bw / 2, y),
                                     arrowstyle="-|>", mutation_scale=17, color="0.4", lw=1.7, zorder=2))
        if i in arrow_labels:
            ax.text((cx0 + cx1) / 2, y + 0.92, arrow_labels[i], ha="center", va="bottom",
                    fontsize=9.5, color="#5a6b7a", style="italic", zorder=4)

    # emphasis: the new layer
    ax.annotate("the only new layer\n(slots in post-MHCflurry)", xy=(centers[3], y - 0.8),
                xytext=(centers[3], 0.4), ha="center", va="center", fontsize=11.5,
                color="#9c7a14", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#b8860b", lw=1.9), zorder=5)

    # no figure title — the slide header ("Where it would plug in") already names it
    fig.savefig(FIGURES_DIR / "integration_map.png")
    plt.close(fig)


def main():
    pos, neg = make_data()
    kpos, kneg = class_kdes(pos, neg)
    fig_kde_densities(pos, neg, kpos, kneg)
    fig_density_ratio(kpos, kneg)
    fig_imbalance_prior(kpos, kneg)
    fig_isotonic_vs_cir(pos, neg)
    fig_calibrated_curve(pos, neg)
    fig_combine_schematic()
    fig_integration_map()
    print("Wrote 7 figures to", FIGURES_DIR)


if __name__ == "__main__":
    main()
