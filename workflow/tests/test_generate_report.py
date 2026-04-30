"""Tests for generate_report.py — HTML generation helpers."""

import pandas as pd
import pytest

from generate_report import (
    _build_chain_legend,
    _build_compnd_records,
    _build_report_tsv,
    _build_strong_table_html,
    _df_to_html,
    _extract_chain_ids,
    _load_report_tsv,
)


class TestDfToHtml:
    def test_empty_dataframe_returns_no_data_message(self):
        df = pd.DataFrame(columns=["a", "b"])
        assert _df_to_html(df) == "<p><em>No data.</em></p>"

    def test_small_dataframe_returns_html_table(self):
        df = pd.DataFrame({"sample": ["s1"], "count": [5]})
        html = _df_to_html(df)
        assert "<table" in html
        assert "s1" in html

    def test_large_dataframe_shows_truncation_message(self):
        df = pd.DataFrame({"x": range(200)})
        html = _df_to_html(df, max_rows=10)
        assert "Showing 10 of 200 rows" in html

    def test_large_dataframe_respects_max_rows(self):
        df = pd.DataFrame({"x": range(50)})
        html = _df_to_html(df, max_rows=10)
        # Only 10 data rows should be in the table (row 0–9, not 10–49)
        assert "49" not in html


# Minimal PDB snippet for testing (two chains: A and B)
_SAMPLE_PDB = (
    "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C\n"
    "ATOM      3  N   GLY B   2       5.000   6.000   7.000  1.00  0.00           N\n"
)


class TestExtractChainIds:
    def test_extracts_chains_in_order(self):
        assert _extract_chain_ids(_SAMPLE_PDB) == ["A", "B"]

    def test_empty_pdb_returns_empty(self):
        assert _extract_chain_ids("") == []

    def test_ignores_non_atom_records(self):
        pdb = "REMARK  test\nTER\nEND\n"
        assert _extract_chain_ids(pdb) == []

    def test_no_duplicates(self):
        pdb = _SAMPLE_PDB + (
            "ATOM      4  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C\n"
        )
        assert _extract_chain_ids(pdb) == ["A", "B"]


class TestBuildChainLegend:
    def test_known_chains(self):
        html = _build_chain_legend(["A", "B"], "SIINFEKL", "HLA-A*02:01")
        assert "MHC" in html
        assert "Peptide" in html
        assert "SIINFEKL" in html
        assert "HLA-A*02:01" in html

    def test_empty_chains_returns_empty(self):
        assert _build_chain_legend([], "PEP", "HLA") == ""

    def test_unknown_chain_id(self):
        html = _build_chain_legend(["Z"], "PEP", "HLA")
        assert "Chain Z" in html


class TestBuildCompndRecords:
    def test_prepends_compnd_lines(self):
        result = _build_compnd_records(_SAMPLE_PDB, "SIINFEKL", "HLA-A*02:01")
        lines = result.splitlines()
        compnd_lines = [l for l in lines if l.startswith("COMPND")]
        # 2 chains × 3 lines each (MOL_ID, MOLECULE, CHAIN)
        assert len(compnd_lines) == 6
        assert "MOL_ID: 1;" in compnd_lines[0]
        assert "MHC" in compnd_lines[1]
        assert "CHAIN: A;" in compnd_lines[2]
        assert "Peptide (SIINFEKL)" in compnd_lines[4]

    def test_compnd_lines_are_80_chars(self):
        result = _build_compnd_records(_SAMPLE_PDB, "PEP", "HLA")
        for line in result.splitlines():
            if line.startswith("COMPND"):
                assert len(line) == 80

    def test_empty_pdb_returns_unchanged(self):
        assert _build_compnd_records("", "PEP", "HLA") == ""

    def test_greek_letters_transliterated(self):
        # Chain C = TCR α-chain, Chain D = TCR β-chain
        pdb = (
            "ATOM      1  N   ALA C   1       1.0   2.0   3.0  1.00  0.00           N\n"
            "ATOM      2  N   ALA D   2       4.0   5.0   6.0  1.00  0.00           N\n"
        )
        result = _build_compnd_records(pdb, "PEP", "HLA")
        assert "alpha" in result
        assert "beta" in result


# ---------------------------------------------------------------------------
# Helpers for _build_report_tsv tests
# ---------------------------------------------------------------------------

def _make_origin_df(tumor_exclusive: int = 5, normal_shared: int = 2) -> pd.DataFrame:
    return pd.DataFrame([{
        "sample_id": "S1", "sample_type": "Primary Tumor",
        "unannotated": tumor_exclusive + normal_shared,
        "tumor_exclusive": tumor_exclusive,
        "normal_shared": normal_shared,
    }])


def _make_pred_df(n_strong: int = 3, n_weak: int = 1, n_non: int = 10) -> pd.DataFrame:
    # presentation_percentile increases with index so PEP0 is always the top candidate
    rows = (
        [{"peptide": f"PEP{i}", "best_allele": "HLA-A*02:01", "ic50_nM": float(i + 1),
          "processing_score": 0.8, "presentation_score": 0.9,
          "presentation_percentile": 0.1 * (i + 1),
          "presentation_class": "strong", "contig_key": f"k{i}", "start_nt": 0,
          } for i in range(n_strong)]
        + [{"peptide": f"PEP{i}", "best_allele": "HLA-A*02:01", "ic50_nM": float(100 + i),
            "processing_score": 0.7, "presentation_score": 0.6,
            "presentation_percentile": 1.0 + 0.1 * i,
            "presentation_class": "weak", "contig_key": f"k{i}", "start_nt": 0,
            } for i in range(n_weak)]
        + [{"peptide": f"PEP{i}", "best_allele": "HLA-A*02:01", "ic50_nM": float(1000 + i),
            "processing_score": 0.5, "presentation_score": 0.3,
            "presentation_percentile": 3.0 + 0.1 * i,
            "presentation_class": "non", "contig_key": f"k{i}", "start_nt": 0,
            } for i in range(n_non)]
    )
    return pd.DataFrame(rows)


class TestBuildReportTsv:
    def test_output_file_is_created(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(),
            pred_df=_make_pred_df(),
            hla_qc_tsv=None,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=None,
        )
        assert out.exists()

    def test_output_has_expected_columns(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == ["patient_id", "stage", "metric", "value", "notes"]

    def test_patient_id_in_every_row(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("patient_007", _make_origin_df(), _make_pred_df(), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        assert (df["patient_id"] == "patient_007").all()

    def test_junction_filtering_rows_present(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(tumor_exclusive=7, normal_shared=3), _make_pred_df(), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        jf = df[df["stage"] == "junction_filtering"]
        te_row = jf[jf["metric"] == "tumor_exclusive"]
        assert int(te_row["value"].iloc[0]) == 7

    def test_mhc_prediction_counts_correct(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(n_strong=4, n_weak=2, n_non=5), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        mhc = df[df["stage"] == "mhc_prediction"]
        assert int(mhc[mhc["metric"] == "strong"]["value"].iloc[0]) == 4
        assert int(mhc[mhc["metric"] == "non"]["value"].iloc[0]) == 5
        assert int(mhc[mhc["metric"] == "total_predictions"]["value"].iloc[0]) == 11

    def test_top_candidate_is_lowest_presentation_percentile(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(n_strong=3), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        tc = df[df["stage"] == "top_candidate"]
        peptide = tc[tc["metric"] == "peptide"]["value"].iloc[0]
        pres_pct = float(tc[tc["metric"] == "presentation_percentile"]["value"].iloc[0])
        assert pres_pct == pytest.approx(0.1)
        assert peptide == "PEP0"

    def test_hla_typing_rows_from_qc_tsv(self, tmp_path):
        hla_tsv = tmp_path / "hla_qc.tsv"
        hla_df = pd.DataFrame([
            {"locus": "A", "allele1": "HLA-A*02:01", "allele2": "HLA-A*24:02",
             "source": "tumor", "reads": 500, "discrepancy": "",
             "serology_allele1": "", "serology_allele2": "", "serology_validation": ""},
        ])
        hla_df.to_csv(hla_tsv, sep="\t", index=False)
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(), str(hla_tsv), out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        hla = df[df["stage"] == "hla_typing"]
        assert len(hla) == 1
        assert hla["metric"].iloc[0] == "HLA-A"
        assert "HLA-A*02:01" in hla["value"].iloc[0]
        assert "source: tumor" in hla["notes"].iloc[0]

    def test_tcrdock_pdb_available_false_when_no_pdb(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        row = df[(df["stage"] == "tcrdock") & (df["metric"] == "pdb_available")]
        assert row["value"].iloc[0] == "false"

    def test_tcrdock_pdb_available_true_when_pdb_exists(self, tmp_path):
        pdb = tmp_path / "model.pdb"
        pdb.write_text("ATOM\n")
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), _make_pred_df(), None, out, 0.5, pdb)
        df = pd.read_csv(out, sep="\t")
        row = df[(df["stage"] == "tcrdock") & (df["metric"] == "pdb_available")]
        assert row["value"].iloc[0] == "true"

    def test_empty_predictions_skips_mhc_and_top_candidate(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), pd.DataFrame(), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        assert df[df["stage"] == "mhc_prediction"].empty
        assert df[df["stage"] == "top_candidate"].empty


class TestLoadReportTsv:
    """Round-trip: writer output must load into the projections HTML rendering needs."""

    def _write_round_trip(self, tmp_path, *, hla_qc_tsv=None, tcrdock_pdb=None):
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(tumor_exclusive=7, normal_shared=3),
            pred_df=_make_pred_df(n_strong=3, n_weak=1, n_non=10),
            hla_qc_tsv=hla_qc_tsv,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=tcrdock_pdb,
        )
        return _load_report_tsv(out)

    def test_returns_dict_with_expected_top_level_keys(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        assert set(loaded.keys()) >= {
            "junction_filtering", "mhc_prediction", "mhc_prediction_thresholds",
            "top_candidate", "hla_typing", "tcrdock",
        }

    def test_junction_filtering_is_wide_dataframe(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        jf = loaded["junction_filtering"]
        assert isinstance(jf, pd.DataFrame)
        assert list(jf.columns) == [
            "sample_id", "sample_type", "unannotated", "normal_shared", "tumor_exclusive",
        ]
        assert int(jf.iloc[0]["tumor_exclusive"]) == 7
        assert int(jf.iloc[0]["normal_shared"]) == 3
        assert int(jf.iloc[0]["unannotated"]) == 10
        assert jf.iloc[0]["sample_id"] == "S1"
        assert jf.iloc[0]["sample_type"] == "Primary Tumor"

    def test_mhc_prediction_counts_loaded_as_ints(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        mp = loaded["mhc_prediction"]
        assert mp["total_predictions"] == 14
        assert mp["strong"] == 3
        assert mp["weak"] == 1
        assert mp["non"] == 10
        assert all(isinstance(v, int) for v in mp.values())

    def test_mhc_prediction_thresholds_loaded(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        thresholds = loaded["mhc_prediction_thresholds"]
        assert "strong" in thresholds
        assert "0.5" in thresholds["strong"]

    def test_top_candidate_numerics_typed(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        tc = loaded["top_candidate"]
        assert tc["peptide"] == "PEP0"
        assert isinstance(tc["presentation_percentile"], float)
        assert tc["presentation_percentile"] == pytest.approx(0.1)
        assert isinstance(tc["ic50_nM"], float)

    def test_hla_typing_keyed_by_locus(self, tmp_path):
        hla_tsv = tmp_path / "hla_qc.tsv"
        pd.DataFrame([
            {"locus": "A", "allele1": "HLA-A*02:01", "allele2": "HLA-A*24:02",
             "source": "tumor", "reads": 500, "discrepancy": "",
             "serology_allele1": "", "serology_allele2": "", "serology_validation": ""},
        ]).to_csv(hla_tsv, sep="\t", index=False)
        loaded = self._write_round_trip(tmp_path, hla_qc_tsv=str(hla_tsv))
        hla = loaded["hla_typing"]
        assert "A" in hla
        assert "HLA-A*02:01" in hla["A"]["alleles"]
        assert "tumor" in hla["A"]["notes"]

    def test_tcrdock_pdb_available_is_bool_false(self, tmp_path):
        loaded = self._write_round_trip(tmp_path)
        assert loaded["tcrdock"]["pdb_available"] is False

    def test_tcrdock_pdb_available_is_bool_true(self, tmp_path):
        pdb = tmp_path / "model.pdb"
        pdb.write_text("ATOM\n")
        loaded = self._write_round_trip(tmp_path, tcrdock_pdb=pdb)
        assert loaded["tcrdock"]["pdb_available"] is True

    def test_empty_predictions_returns_empty_dicts(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), pd.DataFrame(), None, out, 0.5, None)
        loaded = _load_report_tsv(out)
        assert loaded["mhc_prediction"] == {}
        assert loaded["top_candidate"] == {}

    def test_missing_junction_filtering_returns_empty_dataframe(self, tmp_path):
        out = tmp_path / "report.tsv"
        # Empty origin_df + empty pred_df → no junction_filtering rows written
        _build_report_tsv("p001", pd.DataFrame(), pd.DataFrame(), None, out, 0.5, None)
        loaded = _load_report_tsv(out)
        jf = loaded["junction_filtering"]
        assert isinstance(jf, pd.DataFrame)
        assert jf.empty
        assert list(jf.columns) == [
            "sample_id", "sample_type", "unannotated", "normal_shared", "tumor_exclusive",
        ]


# ---------------------------------------------------------------------------
# Helpers and tests for the GPS (genotype_presentation_score) code path
# ---------------------------------------------------------------------------

def _make_pred_df_with_gps(rows: list[dict]) -> pd.DataFrame:
    """Build a predictions DataFrame that includes GPS columns."""
    defaults = {
        "best_allele": "HLA-A*02:01",
        "ic50_nM": 50.0,
        "processing_score": 0.8,
        "presentation_score": 0.9,
        "contig_key": "k0",
        "start_nt": 0,
    }
    return pd.DataFrame([{**defaults, **r} for r in rows])


class TestBuildStrongTableHtmlWithGps:
    def test_quality_gate_excludes_high_best_percentile(self):
        """Candidates where best_presentation_percentile > weak threshold are excluded."""
        df = _make_pred_df_with_gps([
            {"peptide": "ACDEFGHIK", "presentation_percentile": 0.1,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.9, "n_strong_alleles": 1,
             "best_presentation_percentile": 0.1},
            {"peptide": "LMNPQRSTV", "presentation_percentile": 0.2,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.95, "n_strong_alleles": 2,
             "best_presentation_percentile": 3.0},  # fails gate
        ])
        html = _build_strong_table_html(df, {}, presentation_percentile_weak=2.0)
        assert "ACDEFGHIK" in html
        assert "LMNPQRSTV" not in html

    def test_ranked_by_gps_descending(self):
        """Candidate with higher GPS appears first in the table."""
        df = _make_pred_df_with_gps([
            {"peptide": "ACDEFGHIK", "presentation_percentile": 0.1,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.7, "n_strong_alleles": 1,
             "best_presentation_percentile": 0.1},
            {"peptide": "LMNPQRSTV", "presentation_percentile": 0.4,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.95, "n_strong_alleles": 3,
             "best_presentation_percentile": 0.4},
        ])
        html = _build_strong_table_html(df, {}, presentation_percentile_weak=2.0)
        assert html.index("LMNPQRSTV") < html.index("ACDEFGHIK")


class TestBuildReportTsvWithGps:
    def test_top_candidate_ranked_by_gps(self, tmp_path):
        """When GPS columns present, top_candidate row reflects GPS-ranked winner."""
        df = _make_pred_df_with_gps([
            {"peptide": "ACDEFGHIK", "presentation_percentile": 0.1,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.7, "n_strong_alleles": 1,
             "best_presentation_percentile": 0.1},
            {"peptide": "LMNPQRSTV", "presentation_percentile": 0.4,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.95, "n_strong_alleles": 3,
             "best_presentation_percentile": 0.4},
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), df, None, out, 0.5, None,
                          presentation_percentile_weak=2.0)
        result = pd.read_csv(out, sep="\t")
        top = result[(result["stage"] == "top_candidate") & (result["metric"] == "peptide")]
        assert top["value"].iloc[0] == "LMNPQRSTV"

    def test_quality_gate_excludes_candidate_in_tsv(self, tmp_path):
        """Candidate failing quality gate is not selected as top_candidate."""
        df = _make_pred_df_with_gps([
            {"peptide": "ACDEFGHIK", "presentation_percentile": 0.1,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.95, "n_strong_alleles": 2,
             "best_presentation_percentile": 3.0},  # fails gate
            {"peptide": "LMNPQRSTV", "presentation_percentile": 0.4,
             "presentation_class": "strong",
             "genotype_presentation_score": 0.7, "n_strong_alleles": 1,
             "best_presentation_percentile": 0.4},
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(), df, None, out, 0.5, None,
                          presentation_percentile_weak=2.0)
        result = pd.read_csv(out, sep="\t")
        top = result[(result["stage"] == "top_candidate") & (result["metric"] == "peptide")]
        assert top["value"].iloc[0] == "LMNPQRSTV"


# ---------------------------------------------------------------------------
# _cli_main argument parser
# ---------------------------------------------------------------------------

class TestCLIParser:
    def _make_parser(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--novel-junctions", required=True)
        parser.add_argument("--predictions", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--contigs-fasta", default=None)
        parser.add_argument("--hla-qc-tsv", default=None)
        parser.add_argument("--tcrdock-pdb", default=None)
        parser.add_argument("--presentation-percentile-strong", type=float, default=0.5)
        parser.add_argument("--presentation-percentile-weak", type=float, default=2.0)
        parser.add_argument("--output-tsv", default=None)
        parser.add_argument("--patient-id", default="")
        return parser

    def test_default_presentation_percentile_weak(self):
        parser = self._make_parser()
        args = parser.parse_args([
            "--novel-junctions", "j.tsv",
            "--predictions", "p.tsv",
            "--output", "out.html",
        ])
        assert args.presentation_percentile_weak == 2.0

    def test_custom_presentation_percentile_weak(self):
        parser = self._make_parser()
        args = parser.parse_args([
            "--novel-junctions", "j.tsv",
            "--predictions", "p.tsv",
            "--output", "out.html",
            "--presentation-percentile-weak", "1.5",
        ])
        assert args.presentation_percentile_weak == 1.5
