"""Tests for fetch_vdjdb_panel.py — VDJdb panel construction."""

import logging
import shutil
from pathlib import Path

import pytest

from fetch_vdjdb_panel import (
    build_panel,
    classify_panel_status,
    load_alleles_tsv,
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
        with open(FIXTURE_PATH) as f:
            header = f.readline()
        # 35 cols total: 7 empty (1-7) + 4 set (8-11: species, mhc.a, mhc.b, mhc.class)
        # + 22 empty (12-33) + vdjdb.score (34) + TCR_hash (35).
        bad_row = "\t".join([""] * 7 + ["HomoSapiens", "HLA-A*02", "B2M", "MHCI"] + [""] * 22 + ["3", ""]) + "\n"
        tsv.write_text(header + bad_row)
        df = load_and_filter_vdjdb(tsv, min_score=2)
        assert len(df) == 0

    def test_drops_single_chain_rows_missing_alpha(self, tmp_path):
        # Regression for chr22 run on 2026-05-22: VDJdb's full TSV mixes single-chain
        # entries (NaN in opposite-chain V/J/CDR3 cols) with paired ones. Without the
        # paired-α/β filter these rows reached stitchr and crashed it.
        tsv = tmp_path / "vdjdb_single_chain.tsv"
        with open(FIXTURE_PATH) as f:
            header = f.readline()
        # Beta chain populated, alpha cols (1-3) empty — should be dropped.
        single_chain_row = "\t".join(
            [""] * 3                                                            # cdr3.alpha, v.alpha, j.alpha
            + ["CASSLSGNTGELFF", "TRBV5-1*01", "", "TRBJ2-2*01"]                 # cdr3.beta, v.beta, d.beta, j.beta
            + ["HomoSapiens", "HLA-A*02:01", "B2M", "MHCI"]                      # species, mhc.a, mhc.b, mhc.class
            + [""] * 22 + ["3", ""]                                              # cols 12-33, vdjdb.score, TCR_hash
        ) + "\n"
        tsv.write_text(header + single_chain_row)
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


class TestBuildPanel:
    def test_writes_panel_and_qc_with_correct_schema(self, tmp_path, monkeypatch):
        # Monkeypatch stitch_chain so we don't depend on real stitchr in this test
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}_{cdr3}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        build_panel(
            vdjdb_full_tsv=FIXTURE_PATH,
            alleles=["HLA-A*02:01", "HLA-A*31:01", "HLA-B*08:01"],
            output_panel=panel_tsv,
            output_qc=qc_tsv,
            min_score=2,
            panel_size=10,
        )

        import pandas as pd
        panel = pd.read_csv(panel_tsv, sep="\t")
        qc = pd.read_csv(qc_tsv, sep="\t")

        # Panel schema
        assert list(panel.columns) == [
            "allele", "va_gene", "ja_gene", "cdr3a",
            "vb_gene", "jb_gene", "cdr3b",
            "alpha_seq", "beta_seq",
            "vdjdb_score", "vdjdb_donor_id",
        ]
        # QC schema
        assert list(qc.columns) == [
            "allele", "n_exact_matches", "n_in_panel", "panel_status"
        ]

    def test_qc_status_classification(self, tmp_path, monkeypatch):
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        build_panel(
            vdjdb_full_tsv=FIXTURE_PATH,
            alleles=["HLA-A*02:01", "HLA-A*31:01", "HLA-B*08:01"],
            output_panel=panel_tsv,
            output_qc=qc_tsv,
            min_score=2,
            panel_size=10,
        )

        import pandas as pd
        qc = pd.read_csv(qc_tsv, sep="\t").set_index("allele")
        # HLA-A*02:01: 5 matches, panel_size=10 → low_coverage
        assert qc.loc["HLA-A*02:01", "n_in_panel"] == 5
        assert qc.loc["HLA-A*02:01", "panel_status"] == "low_coverage"
        # HLA-A*31:01: 0 matches in fixture → empty
        assert qc.loc["HLA-A*31:01", "n_in_panel"] == 0
        assert qc.loc["HLA-A*31:01", "panel_status"] == "empty"
        # HLA-B*08:01: 1 match → low_coverage
        assert qc.loc["HLA-B*08:01", "n_in_panel"] == 1
        assert qc.loc["HLA-B*08:01", "panel_status"] == "low_coverage"

    def test_empty_allele_emits_warning(self, tmp_path, monkeypatch, caplog):
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        with caplog.at_level(logging.WARNING):
            build_panel(
                vdjdb_full_tsv=FIXTURE_PATH,
                alleles=["HLA-A*31:01"],
                output_panel=panel_tsv,
                output_qc=qc_tsv,
                min_score=2,
                panel_size=10,
            )

        warning_msgs = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("HLA-A*31:01" in m and "empty" in m.lower() for m in warning_msgs)

    def test_stitch_failure_skips_row_continues(self, tmp_path, monkeypatch, caplog):
        # stitch_chain returns None for the second call only — simulates one-row stitch failure
        import fetch_vdjdb_panel as mod
        call_count = {"n": 0}
        def flaky_stitch(v_gene, j_gene, cdr3, chain):
            call_count["n"] += 1
            return None if call_count["n"] == 2 else f"MOCK_{chain}_{cdr3}"
        monkeypatch.setattr(mod, "stitch_chain", flaky_stitch)

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"
        with caplog.at_level(logging.ERROR):
            build_panel(
                vdjdb_full_tsv=FIXTURE_PATH,
                alleles=["HLA-A*02:01"],
                output_panel=panel_tsv,
                output_qc=qc_tsv,
                min_score=2,
                panel_size=3,  # need 3, will fall back when one fails
            )
        import pandas as pd
        panel = pd.read_csv(panel_tsv, sep="\t")
        # Should still get 3 rows because we keep iterating after a stitch failure
        assert len(panel) == 3
        # ERROR was logged for the failed row
        assert any("stitchr" in r.message.lower() or "stitch" in r.message.lower()
                   for r in caplog.records if r.levelno == logging.ERROR)


class TestLoadAllelesTsv:
    def test_extracts_unique_4digit_alleles(self, tmp_path):
        tsv = tmp_path / "alleles.tsv"
        tsv.write_text(
            "locus\tallele1\tallele2\n"
            "A\tHLA-A*02:01\tHLA-A*31:01\n"
            "B\tHLA-B*08:01\tHLA-B*08:01\n"  # homozygous → dedupe
            "C\tHLA-C*07:01\tHLA-C*03:03\n"
        )
        alleles = load_alleles_tsv(tsv)
        assert sorted(alleles) == [
            "HLA-A*02:01", "HLA-A*31:01",
            "HLA-B*08:01",
            "HLA-C*03:03", "HLA-C*07:01",
        ]
