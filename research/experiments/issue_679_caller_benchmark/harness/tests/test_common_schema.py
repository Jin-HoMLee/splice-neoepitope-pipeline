"""Tests for the benchmark harness common output schema (#966, AC-2)."""

import pytest

from common_schema import (
    CommonRecord,
    normalize_junction_id,
    parse_junction_string,
    validate,
)


def _peptide_record(**overrides):
    base = dict(caller="asneo", junction_id=None, genome_build="hg19",
                strand=None, peptide="LEQGTHPKFQ", record_level="peptide")
    base.update(overrides)
    return CommonRecord(**base)


def test_validate_peptide_level_allows_null_junction():
    validate(_peptide_record())  # must not raise


def test_validate_junction_level_requires_junction_id():
    # Matched pair: flip only record_level; opposite outcomes.
    bad = _peptide_record(record_level="junction")   # junction_id/strand still None
    with pytest.raises(ValueError):
        validate(bad)
    good = _peptide_record(record_level="junction",
                           junction_id="chr22:100-200:+", strand="+")
    validate(good)  # must not raise


def test_validate_peptide_level_rejects_a_junction_id():
    with pytest.raises(ValueError):
        validate(_peptide_record(junction_id="chr22:100-200:+", strand="+"))


def test_validate_unknown_record_level_raises():
    with pytest.raises(ValueError):
        validate(_peptide_record(record_level="isoform"))


class TestNormalizeJunctionId:
    def test_adds_ucsc_chr_prefix_to_ensembl_contig(self):
        # Ensembl-style bare contig "2" -> UCSC "chr2"; canonical genomic-order key.
        assert (
            normalize_junction_id("2", 152389996, 152392205, "-")
            == "chr2:152389996-152392205:-"
        )

    def test_keeps_existing_chr_prefix(self):
        assert normalize_junction_id("chr2", 100, 200, "+") == "chr2:100-200:+"

    def test_orders_coordinates_genomically_regardless_of_input_order(self):
        # A junction is identified by the same key no matter which endpoint is
        # passed first, so two tools that disagree on argument order still collide.
        assert normalize_junction_id("chr7", 200, 100, "-") == "chr7:100-200:-"

    def test_ensembl_mt_maps_to_ucsc_chrm_not_chrmt(self):
        # UCSC mitochondrial contig is "chrM"; Ensembl names it "MT". Bare
        # prefixing would give "chrMT", so a tool emitting "MT" and one emitting
        # "chrM" would not collide on one key. Normalize both to "chrM".
        assert normalize_junction_id("MT", 100, 200, "+") == "chrM:100-200:+"
        assert normalize_junction_id("chrM", 100, 200, "+") == "chrM:100-200:+"


class TestParseJunctionString:
    def test_parses_splice2neo_style_string(self):
        # splice2neo emits "chr2:152389996-152392205:-".
        assert parse_junction_string("chr2:152389996-152392205:-") == (
            "chr2",
            152389996,
            152392205,
            "-",
        )
