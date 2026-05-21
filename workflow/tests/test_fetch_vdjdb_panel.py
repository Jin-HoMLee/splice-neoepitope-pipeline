"""Tests for fetch_vdjdb_panel.py — VDJdb panel construction."""

from pathlib import Path

import pytest

from fetch_vdjdb_panel import load_and_filter_vdjdb, normalize_allele_to_4digit


FIXTURE_PATH = Path(__file__).parent / "fixtures" / "vdjdb_mini.tsv"


class TestNormalizeAlleleTo4Digit:
    @pytest.mark.parametrize("inp,expected", [
        ("HLA-A*02:01", "HLA-A*02:01"),       # already 4-digit
        ("HLA-A*02:01:110", "HLA-A*02:01"),   # 6-digit truncates
        ("HLA-B*15:63", "HLA-B*15:63"),       # rare allele, already 4-digit
        ("HLA-C*07:01:01:03", "HLA-C*07:01"), # 8-digit truncates
    ])
    def test_4digit_or_longer_normalizes(self, inp, expected):
        assert normalize_allele_to_4digit(inp) == expected

    @pytest.mark.parametrize("inp", [
        "HLA-A*02",      # 2-digit only — excluded
        "HLA-A",         # no allele subtype at all
        "",              # empty
    ])
    def test_2digit_or_invalid_returns_none(self, inp):
        assert normalize_allele_to_4digit(inp) is None


class TestLoadAndFilterVdjdb:
    def test_filters_to_homosapiens_mhci_minscore(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        # Expected pass rows: 1, 2, 3, 5, 6, 9, 10 — i.e. 7 rows
        # (row 4 fails score=0; row 7 fails species=Mouse; row 8 fails MHCII)
        assert len(df) == 7
        assert (df["species"] == "HomoSapiens").all()
        assert (df["mhc.class"] == "MHCI").all()
        assert (df["vdjdb.score"] >= 2).all()

    def test_min_score_threshold_is_inclusive(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=3)
        # Drops the score=2 rows (rows 3 and 6 in fixture)
        assert (df["vdjdb.score"] >= 3).all()
        assert len(df) == 5

    def test_normalizes_mhc_a_inplace(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        # Row 2's HLA-A*02:01:110 must be normalized to HLA-A*02:01
        a0201_rows = df[df["mhc.a_4digit"] == "HLA-A*02:01"]
        assert len(a0201_rows) == 5  # rows 1, 2, 3, 9, 10

    def test_drops_rows_with_unnormalizable_allele(self, tmp_path):
        # Synthetic fixture with a 2-digit-only allele — should be dropped
        tsv = tmp_path / "vdjdb_partial.tsv"
        header = open(FIXTURE_PATH).readline()
        # 35 cols total: 7 empty (1-7) + 4 set (8-11: species, mhc.a, mhc.b, mhc.class)
        # + 22 empty (12-33) + vdjdb.score (34) + TCR_hash (35).
        bad_row = "\t".join([""] * 7 + ["HomoSapiens", "HLA-A*02", "B2M", "MHCI"] + [""] * 22 + ["3", ""]) + "\n"
        tsv.write_text(header + bad_row)
        df = load_and_filter_vdjdb(tsv, min_score=2)
        assert len(df) == 0
