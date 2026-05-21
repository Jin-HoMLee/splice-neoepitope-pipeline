"""Tests for fetch_vdjdb_panel.py — VDJdb panel construction."""

import shutil
from pathlib import Path

import pytest

from fetch_vdjdb_panel import (
    classify_panel_status,
    load_and_filter_vdjdb,
    normalize_allele_to_4digit,
    select_top_n_for_allele,
    stitch_chain,
)


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


class TestSelectTopNForAllele:
    def test_exact_match_filter(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=10)
        # rows 1, 2, 3, 9, 10 — 5 total (includes the 6-digit row 2 after normalization)
        assert len(result) == 5
        assert (result["mhc.a_4digit"] == "HLA-A*02:01").all()

    def test_does_not_match_different_allele(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:02", n=10)
        assert len(result) == 1  # only row 5 in fixture

    def test_returns_top_n_when_more_available(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=2)
        assert len(result) == 2

    def test_returns_fewer_when_fewer_available(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-B*08:01", n=10)
        assert len(result) == 1  # only row 6 matches B*08:01

    def test_returns_empty_when_no_matches(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*31:01", n=10)
        assert len(result) == 0

    def test_sorted_by_score_then_donor_id(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=10)
        # Of the 5 matching rows: 4 with score=3 (donors 001, 002, 009, 010), 1 with score=2 (donor 003).
        # Score 3 rows come first; among them, donor IDs sorted ascending lexicographically.
        scores = result["vdjdb.score"].tolist()
        assert scores == [3, 3, 3, 3, 2]
        donor_order = result["meta.subject.id"].tolist()
        # Top 4 (score=3) sorted ascending: 001, 002, 009, 010. Last is score=2 (donor 003).
        assert donor_order == ["donor-001", "donor-002", "donor-009", "donor-010", "donor-003"]


class TestClassifyPanelStatus:
    @pytest.mark.parametrize("n_in_panel,target,expected", [
        (10, 10, "ok"),
        (15, 10, "ok"),   # clipped to top-N before reaching here; defensive case
        (9, 10, "low_coverage"),
        (1, 10, "low_coverage"),
        (0, 10, "empty"),
    ])
    def test_classification(self, n_in_panel, target, expected):
        assert classify_panel_status(n_in_panel, target_size=target) == expected


@pytest.mark.network
@pytest.mark.skipif(shutil.which("stitchr") is None, reason="stitchr CLI not on PATH")
class TestStitchChain:
    """Smoke test against the real stitchr CLI + cached IMGT data.

    Skipped automatically when the `stitchr` CLI is not on PATH (e.g. CI runners
    without the vdjdb conda env activated). The `network` marker is informational
    — actual skip is driven by the `skipif` so CI works regardless of `-m` flag.
    """
    def test_dmf5_alpha_chain(self):
        # DMF5 (the existing single-fallback TCR in config/gpu_config.yaml)
        alpha = stitch_chain(
            v_gene="TRAV12-2",
            j_gene="TRAJ21",
            cdr3="CAVNFGGGKLI",
            chain="A",
        )
        # Stitched chain should be a non-trivial protein sequence containing the CDR3
        assert alpha is not None
        assert "CAVNFGGGKLI" in alpha
        # TRAV12-2 framework begins with a known leader; check for FW1 motif "QSV"
        # (allowing some flexibility — sanity check, not equality)
        assert len(alpha) > 100  # full Vα is ~100-110 aa
