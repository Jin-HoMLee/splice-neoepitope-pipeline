"""Dedup + merge on the coalesced key (#1086 AC6).

Row identity is COALESCE(junction_id, peptide), disambiguated by hla. Dedup keys on
the (junction_id, peptide, hla) triple with either identity column allowed null, so
the two shapes - coordinate-first and sequence-first - share one rule.
"""
import pandas as pd
import pytest

from registry_dedup import duplicate_keys, junction_view, row_identity

# The real one-junction-multiple-peptides case, already in registry.tsv: one Kim 2025
# constitutive-intron event yields three distinct A*02:01 peptides (two of them a
# nested 9/10-mer pair).
CI_JUNCTION = "ci@16:719606:720123:+|16:719606:719607:+"


def _rows(*triples):
    return pd.DataFrame(
        [{"junction_id": j, "peptide": p, "hla": h} for j, p, h in triples], dtype=str
    ).fillna("")


# --- identity ---------------------------------------------------------------------


def test_identity_coalesces_junction_then_peptide():
    assert row_identity({"junction_id": "chr1:1-2:+", "peptide": "AAA"}) == "chr1:1-2:+"
    assert row_identity({"junction_id": "", "peptide": "AAA"}) == "AAA"
    assert row_identity({"junction_id": " ", "peptide": "AAA"}) == "AAA"


def test_identity_rejects_a_row_with_neither():
    with pytest.raises(ValueError, match="at-least-one-non-null"):
        row_identity({"junction_id": "", "peptide": ""})


# --- dedup ------------------------------------------------------------------------


def test_exact_triple_duplicate_is_flagged():
    df = _rows((CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01"),
               (CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01"))
    assert duplicate_keys(df) == [(CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01")]


def test_one_junction_with_several_peptides_is_not_a_duplicate():
    """The merge rule must NOT collapse distinct peptides onto their shared junction:
    immunogenicity is a property of the peptide-HLA pair, so these are three rows."""
    df = _rows((CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01"),
               (CI_JUNCTION, "FLWPGLGPSV", "HLA-A*02:01"),
               (CI_JUNCTION, "ILGSLTWSC", "HLA-A*02:01"))
    assert duplicate_keys(df) == []


def test_same_peptide_on_two_alleles_is_not_a_duplicate():
    df = _rows(("", "AAAAAAAAA", "HLA-A*02:01"), ("", "AAAAAAAAA", "HLA-B*07:02"))
    assert duplicate_keys(df) == []


def test_two_peptide_null_rows_on_different_junctions_are_distinct():
    df = _rows(("chr1:1-2:+", "", ""), ("chr1:3-4:+", "", ""))
    assert duplicate_keys(df) == []


def test_duplicate_peptide_null_junction_rows_are_flagged():
    df = _rows(("chr1:1-2:+", "", ""), ("chr1:1-2:+", "", ""))
    assert duplicate_keys(df) == [("chr1:1-2:+", "", "")]


def test_a_junction_null_row_never_collides_with_a_peptide_null_row():
    """Both have one empty identity column; the coalesced key must keep them apart
    rather than matching on the shared blank."""
    df = _rows(("", "AAAAAAAAA", ""), ("chr1:1-2:+", "", ""))
    assert duplicate_keys(df) == []


# --- junction-resolution view -----------------------------------------------------


def test_junction_view_groups_peptides_under_their_junction():
    df = _rows((CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01"),
               (CI_JUNCTION, "ILGSLTWSC", "HLA-A*02:01"),
               ("chr1:1-2:+", "", ""))
    view = junction_view(df)
    assert view[CI_JUNCTION] == ["FLWPGLGPS", "ILGSLTWSC"]
    assert view["chr1:1-2:+"] == []


def test_junction_view_dedupes_a_peptide_presented_on_two_alleles():
    """Two registry rows, one peptide: immunogenicity is per peptide-HLA pair, but at
    junction resolution the allele is gone and the peptide must not be double-listed."""
    df = _rows((CI_JUNCTION, "FLWPGLGPS", "HLA-A*02:01"),
               (CI_JUNCTION, "FLWPGLGPS", "HLA-B*07:02"))
    assert junction_view(df)[CI_JUNCTION] == ["FLWPGLGPS"]


def test_junction_view_omits_rows_with_no_junction():
    df = _rows(("", "AAAAAAAAA", "HLA-A*02:01"))
    assert junction_view(df) == {}


def test_live_registry_has_no_duplicate_keys_and_one_shared_junction():
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    assert duplicate_keys(df) == []
    shared = {j: p for j, p in junction_view(df).items() if len(p) > 1}
    assert shared == {CI_JUNCTION: ["FLWPGLGPS", "FLWPGLGPSV", "ILGSLTWSC"]}
