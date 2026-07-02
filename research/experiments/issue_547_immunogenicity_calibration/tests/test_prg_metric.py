"""Unit tests for prg_metric.auprg (weight-aware Area Under Precision-Recall-Gain).

Correctness anchors:
  * The NeurIPS canonical example (Flach & Kull 2015) has a published AUPRG of
    0.432462103574 (meeliskull/prg reference notebook). Reproducing it exactly
    validates the curve construction + crossing-point handling.
  * Perfect ranker -> 1.0; a non-discriminating (all-tied) ranker -> 0.0;
    an inverted ranker -> < 0.  These are prevalence-independent by design.
  * Integer sample weights must equal physically replicating those rows
    (the property that lets us reuse it with the negative-reweighting weights).
"""
import numpy as np
import pytest

from prg_metric import auprg


# ---- published known-answer -------------------------------------------------

def test_canonical_neurips_example_matches_published_value():
    labels = np.array([1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0,
                       1, 0, 0, 1, 0, 0, 0, 1, 0, 1])
    scores = np.arange(1, 26)[::-1].astype(float)
    assert auprg(labels, scores) == pytest.approx(0.432462103574, abs=1e-9)


# ---- analytical anchors -----------------------------------------------------

def test_perfect_ranker_is_one():
    labels = np.array([0, 0, 0, 0, 1, 1, 1])
    scores = np.arange(len(labels), dtype=float)  # positives score highest
    assert auprg(labels, scores) == pytest.approx(1.0, abs=1e-9)


def test_non_discriminating_all_tied_is_zero():
    # every item shares one score -> a single operating point -> AUPRG 0.
    labels = np.array([1, 1, 0, 0, 0, 0])
    scores = np.zeros(len(labels))
    assert auprg(labels, scores) == pytest.approx(0.0, abs=1e-12)


def test_inverted_ranker_is_negative():
    labels = np.array([0, 0, 0, 0, 1, 1, 1])
    scores = np.arange(len(labels), dtype=float)
    assert auprg(labels, -scores) < 0.0


def test_single_class_returns_nan():
    assert np.isnan(auprg(np.ones(5, dtype=int), np.arange(5.0)))
    assert np.isnan(auprg(np.zeros(5, dtype=int), np.arange(5.0)))


# ---- prevalence-invariance --------------------------------------------------

def test_perfect_ranker_is_one_regardless_of_prevalence():
    for n_pos, n_neg in [(1, 99), (50, 50), (95, 5)]:
        labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(int)
        scores = np.concatenate([np.arange(n_neg, n_neg + n_pos),
                                 np.arange(n_neg)]).astype(float)
        assert auprg(labels, scores) == pytest.approx(1.0, abs=1e-9)


# ---- weight semantics -------------------------------------------------------

def test_integer_weights_equal_physical_replication():
    rng = np.random.default_rng(7)
    n = 40
    labels = rng.integers(0, 2, n)
    if labels.sum() in (0, n):            # guard the degenerate draw
        labels[0], labels[1] = 1, 0
    scores = rng.normal(size=n)
    weights = rng.integers(1, 5, n).astype(float)

    rep_labels = np.repeat(labels, weights.astype(int))
    rep_scores = np.repeat(scores, weights.astype(int))

    assert auprg(labels, scores, weights) == pytest.approx(
        auprg(rep_labels, rep_scores), abs=1e-9)


def test_uniform_weights_equal_unweighted():
    rng = np.random.default_rng(11)
    labels = rng.integers(0, 2, 50)
    scores = rng.normal(size=50)
    w = np.full(50, 3.0)
    assert auprg(labels, scores, w) == pytest.approx(auprg(labels, scores), abs=1e-9)


# ---- robustness -------------------------------------------------------------

def test_tie_heavy_scores_run_and_stay_bounded():
    # isotonic-calibrated log-odds carry many tied values; must not error.
    labels = np.array([1, 0, 1, 0, 1, 0, 1, 0, 0, 0])
    scores = np.array([3, 3, 3, 2, 2, 2, 1, 1, 1, 1], dtype=float)
    val = auprg(labels, scores)
    assert -1.0 <= val <= 1.0


def test_better_ranking_scores_higher():
    rng = np.random.default_rng(3)
    n_pos, n_neg = 30, 90
    labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(int)
    signal = np.concatenate([rng.normal(1.5, 1, n_pos), rng.normal(0, 1, n_neg)])
    noise = rng.normal(size=n_pos + n_neg)
    assert auprg(labels, signal) > auprg(labels, noise)
