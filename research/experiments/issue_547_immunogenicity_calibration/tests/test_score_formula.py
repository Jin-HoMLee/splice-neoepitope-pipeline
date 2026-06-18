"""Pure-Python unit tests for genotype_score and normalize_hla.

These tests MUST NOT import mhcflurry or pandas — only stdlib + score_cohort.
They run under research/.venv (Python 3.14, no mhcflurry).
"""
import math
import pytest
from score_cohort import genotype_score, normalize_hla


# ---------------------------------------------------------------------------
# normalize_hla tests
# ---------------------------------------------------------------------------

def test_normalize_hla_bare_supertype():
    """A*01:01 -> HLA-A*01:01"""
    assert normalize_hla("A*01:01") == "HLA-A*01:01"


def test_normalize_hla_with_hla_prefix_no_star():
    """HLA-A01:01 -> HLA-A*01:01 (insert missing star)"""
    assert normalize_hla("HLA-A01:01") == "HLA-A*01:01"


def test_normalize_hla_already_canonical():
    """HLA-A*01:01 -> HLA-A*01:01 (no-op)"""
    assert normalize_hla("HLA-A*01:01") == "HLA-A*01:01"


def test_normalize_hla_c_locus():
    """C*07:02 -> HLA-C*07:02"""
    assert normalize_hla("C*07:02") == "HLA-C*07:02"


# ---------------------------------------------------------------------------
# genotype_score tests
# ---------------------------------------------------------------------------

def test_genotype_score_formula_match():
    """genotype_score matches the analytic formula: log(1 + sum(weights * scores)).

    Weights: A=1.0, B=1.0, C=0.5. For two alleles with known scores, verify
    the formula is applied exactly.
    """
    per_allele = {
        "HLA-A*01:01": 0.8,
        "HLA-B*07:02": 0.6,
    }
    expected = math.log1p(1.0 * 0.8 + 1.0 * 0.6)
    result = genotype_score(per_allele)
    assert abs(result - expected) < 1e-9, f"got {result}, expected {expected}"


def test_genotype_score_single_allele():
    """Single allele with C locus uses weight 0.5."""
    per_allele = {"HLA-C*07:02": 0.4}
    expected = math.log1p(0.5 * 0.4)
    result = genotype_score(per_allele)
    assert abs(result - expected) < 1e-9, f"got {result}, expected {expected}"


def test_genotype_score_empty_returns_zero():
    """Empty per_allele dict -> genotype_score 0.0."""
    assert genotype_score({}) == 0.0


def test_genotype_score_c_weight_lower_than_ab():
    """Swapping a C-locus score with an A-locus score yields a lower result (weight 0.5 vs 1.0)."""
    score_a = genotype_score({"HLA-A*01:01": 0.5})
    score_c = genotype_score({"HLA-C*01:02": 0.5})
    assert score_a > score_c
