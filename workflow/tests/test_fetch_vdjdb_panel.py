"""Tests for fetch_vdjdb_panel.py — VDJdb panel construction."""

import pytest

from fetch_vdjdb_panel import normalize_allele_to_4digit


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
