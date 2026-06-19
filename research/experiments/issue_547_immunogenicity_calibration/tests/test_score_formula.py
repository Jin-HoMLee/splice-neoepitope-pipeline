"""Pure-Python unit tests for genotype_score and normalize_hla.

These tests MUST NOT import mhcflurry or pandas — only stdlib + score_cohort.
They run under research/.venv (Python 3.14, no mhcflurry).
"""
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
# genotype_score tests — product formula: 1 - prod(1 - w_i * p_i)
# ---------------------------------------------------------------------------

def test_genotype_score_formula_match():
    """genotype_score matches the product formula: 1 - prod(1 - w_i * p_i).

    Weights: A/B = 1.0, C = 0.5 (hla_c_weight default).
    Two alleles: HLA-A*02:01 (w=1.0, p=0.8) and HLA-C*07:02 (w=0.5, p=0.6).
    Expected: 1 - (1 - 1.0*0.8) * (1 - 0.5*0.6) = 1 - 0.2 * 0.7 = 0.86
    """
    per_allele = {
        "HLA-A*02:01": 0.8,
        "HLA-C*07:02": 0.6,
    }
    expected = 1.0 - (1.0 - 1.0 * 0.8) * (1.0 - 0.5 * 0.6)  # = 0.86
    result = genotype_score(per_allele)
    assert abs(result - expected) < 1e-9, f"got {result}, expected {expected}"


def test_genotype_score_single_c_allele():
    """Single C-locus allele with p=0.4: 1 - (1 - 0.5*0.4) = 0.2."""
    per_allele = {"HLA-C*07:02": 0.4}
    expected = 1.0 - (1.0 - 0.5 * 0.4)  # = 0.2
    result = genotype_score(per_allele)
    assert abs(result - expected) < 1e-9, f"got {result}, expected {expected}"


def test_genotype_score_empty_returns_zero():
    """Empty per_allele dict -> genotype_score 0.0."""
    assert genotype_score({}) == 0.0


def test_genotype_score_c_weight_lower_than_ab():
    """C-locus contributes less than A-locus for the same presentation score (product form)."""
    score_a = genotype_score({"HLA-A*01:01": 0.5})
    score_c = genotype_score({"HLA-C*01:02": 0.5})
    assert score_a > score_c
