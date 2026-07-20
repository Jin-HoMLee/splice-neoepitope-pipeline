"""Tests for the benchmark harness common output schema (#966, AC-2)."""

from common_schema import normalize_junction_id, parse_junction_string


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


class TestParseJunctionString:
    def test_parses_splice2neo_style_string(self):
        # splice2neo emits "chr2:152389996-152392205:-".
        assert parse_junction_string("chr2:152389996-152392205:-") == (
            "chr2",
            152389996,
            152392205,
            "-",
        )
