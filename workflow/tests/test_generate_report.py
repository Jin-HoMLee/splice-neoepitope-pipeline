"""Tests for generate_report.py — HTML generation helpers."""

import pandas as pd
import pytest

from generate_report import (
    _build_chain_legend,
    _build_compnd_records,
    _build_contig_peek,
    _build_docking_metrics_html,
    _build_filtering_funnel_html,
    _build_report_3d_structure_tsv,
    _build_report_tsv,
    _build_strong_table_html,
    _build_strong_table_html_from_top_candidates,
    _build_structure_section,
    _build_tcr_provenance_html,
    _build_vdjdb_panel_section,
    _df_to_html,
    _extract_chain_ids,
    _load_report_tsv,
    _presenter_counts_html,
    _render_contig_peek,
    generate_report,
    TOP_CANDIDATES_LIMIT,
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

def _make_origin_df(tumor_exclusive: int = 5, normal_shared: int = 2,
                    gtex_pantissue_shared: int = 0) -> pd.DataFrame:
    return pd.DataFrame([{
        "sample_id": "S1", "sample_type": "Primary Tumor",
        "unannotated": tumor_exclusive + normal_shared + gtex_pantissue_shared,
        "tumor_exclusive": tumor_exclusive,
        "normal_shared": normal_shared,
        "gtex_pantissue_shared": gtex_pantissue_shared,
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

    def test_top_candidates_total_and_capped_recorded(self, tmp_path):
        # 12 strong + 1 weak = 13 quality-gated presenters, capped to TOP_CANDIDATES_LIMIT
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(),
                          _make_pred_df(n_strong=12, n_weak=1, n_non=2), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        mhc = df[df["stage"] == "mhc_prediction"]
        total = int(mhc[mhc["metric"] == "top_candidates_total"]["value"].iloc[0])
        capped = int(mhc[mhc["metric"] == "top_candidates_capped"]["value"].iloc[0])
        assert total == 13
        assert capped == TOP_CANDIDATES_LIMIT

    def test_top_candidates_capped_equals_total_when_under_limit(self, tmp_path):
        out = tmp_path / "report.tsv"
        _build_report_tsv("p001", _make_origin_df(),
                          _make_pred_df(n_strong=3, n_weak=1, n_non=2), None, out, 0.5, None)
        df = pd.read_csv(out, sep="\t")
        mhc = df[df["stage"] == "mhc_prediction"]
        total = int(mhc[mhc["metric"] == "top_candidates_total"]["value"].iloc[0])
        capped = int(mhc[mhc["metric"] == "top_candidates_capped"]["value"].iloc[0])
        assert total == 4
        assert capped == 4

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

    # ---------- Junction-funnel patient-level totals (Issue #214) ----------

    def _write_stats_tsv(self, path, rows):
        """rows: list of (sample_id, sample_type, category, count) tuples."""
        df = pd.DataFrame(rows, columns=["sample_id", "sample_type", "category", "count"])
        df.to_csv(path, sep="\t", index=False)

    def test_junction_funnel_totals_emitted_when_stats_provided(self, tmp_path):
        stats = tmp_path / "junction_filter_stats.tsv"
        # Each sample has mean_reads_filtered too — the funnel must reconcile.
        self._write_stats_tsv(stats, [
            ("S1", "Primary Tumor", "junctions_raw", 1000),
            ("S1", "Primary Tumor", "mean_reads_filtered", 100),
            ("S1", "Primary Tumor", "annotated_discarded", 750),
            ("S1", "Primary Tumor", "normal_shared", 40),
            ("S1", "Primary Tumor", "tumor_exclusive", 110),
            ("S2", "Primary Tumor", "junctions_raw", 500),
            ("S2", "Primary Tumor", "mean_reads_filtered", 50),
            ("S2", "Primary Tumor", "annotated_discarded", 370),
            ("S2", "Primary Tumor", "normal_shared", 20),
            ("S2", "Primary Tumor", "tumor_exclusive", 60),
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(),
            pred_df=_make_pred_df(),
            hla_qc_tsv=None,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=None,
            junction_filter_stats_tsv=stats,
        )
        df = pd.read_csv(out, sep="\t")
        jf = df[df["stage"] == "junction_filtering"].set_index("metric")["value"]
        # Patient-level totals are sums across the 2 tumor samples
        assert int(jf["junctions_extracted_total"]) == 1500
        assert int(jf["junctions_mean_reads_filtered"]) == 150
        assert int(jf["junctions_annotated_discarded"]) == 1120
        # unannotated_total = normal_shared + tumor_exclusive across samples
        assert int(jf["junctions_unannotated_total"]) == 40 + 110 + 20 + 60

    def test_junction_funnel_totals_reconcile_in_report_tsv(self, tmp_path):
        """report.tsv consumers should see a closed funnel (PR #240 review observation)."""
        stats = tmp_path / "junction_filter_stats.tsv"
        self._write_stats_tsv(stats, [
            ("S1", "Primary Tumor", "junctions_raw", 1000),
            ("S1", "Primary Tumor", "mean_reads_filtered", 100),
            ("S1", "Primary Tumor", "annotated_discarded", 750),
            ("S1", "Primary Tumor", "normal_shared", 40),
            ("S1", "Primary Tumor", "tumor_exclusive", 110),
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(),
            pred_df=_make_pred_df(),
            hla_qc_tsv=None,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=None,
            junction_filter_stats_tsv=stats,
        )
        df = pd.read_csv(out, sep="\t")
        jf = df[df["stage"] == "junction_filtering"].set_index("metric")["value"]
        # extracted_total == mean_reads_filtered + annotated_discarded + unannotated_total
        assert int(jf["junctions_extracted_total"]) == (
            int(jf["junctions_mean_reads_filtered"])
            + int(jf["junctions_annotated_discarded"])
            + int(jf["junctions_unannotated_total"])
        )

    def test_junction_funnel_unannotated_total_includes_gtex(self, tmp_path):
        """Issue #212: gtex_pantissue_shared is a 3rd unannotated bucket and must be
        folded into junctions_unannotated_total so the report.tsv funnel stays closed."""
        stats = tmp_path / "junction_filter_stats.tsv"
        self._write_stats_tsv(stats, [
            ("S1", "Primary Tumor", "junctions_raw", 1000),
            ("S1", "Primary Tumor", "mean_reads_filtered", 100),
            ("S1", "Primary Tumor", "annotated_discarded", 750),
            ("S1", "Primary Tumor", "normal_shared", 40),
            ("S1", "Primary Tumor", "gtex_pantissue_shared", 30),
            ("S1", "Primary Tumor", "tumor_exclusive", 80),
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(gtex_pantissue_shared=30),
            pred_df=_make_pred_df(),
            hla_qc_tsv=None,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=None,
            junction_filter_stats_tsv=stats,
        )
        df = pd.read_csv(out, sep="\t")
        jf = df[df["stage"] == "junction_filtering"].set_index("metric")["value"]
        # unannotated_total folds in the GTEx bucket: 40 + 30 + 80
        assert int(jf["junctions_unannotated_total"]) == 40 + 30 + 80
        # funnel still closed with the 3rd bucket present
        assert int(jf["junctions_extracted_total"]) == (
            int(jf["junctions_mean_reads_filtered"])
            + int(jf["junctions_annotated_discarded"])
            + int(jf["junctions_unannotated_total"])
        )
        # the per-sample gtex_pantissue_shared metric row is surfaced in report.tsv
        assert int(jf["gtex_pantissue_shared"]) == 30

    def test_junction_funnel_totals_omitted_when_no_stats_path(self, tmp_path):
        """Without stats path, the 4 new totals must not appear (back-compat)."""
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
        df = pd.read_csv(out, sep="\t")
        metrics = set(df[df["stage"] == "junction_filtering"]["metric"])
        assert "junctions_extracted_total" not in metrics
        assert "junctions_mean_reads_filtered" not in metrics
        assert "junctions_annotated_discarded" not in metrics
        assert "junctions_unannotated_total" not in metrics

    def test_junction_funnel_totals_marked_as_patient_level(self, tmp_path):
        """The 4 new rows should have a notes value distinguishing them from per-sample rows."""
        stats = tmp_path / "junction_filter_stats.tsv"
        self._write_stats_tsv(stats, [
            ("S1", "Primary Tumor", "junctions_raw", 100),
            ("S1", "Primary Tumor", "mean_reads_filtered", 10),
            ("S1", "Primary Tumor", "annotated_discarded", 70),
            ("S1", "Primary Tumor", "normal_shared", 5),
            ("S1", "Primary Tumor", "tumor_exclusive", 15),
        ])
        out = tmp_path / "report.tsv"
        _build_report_tsv(
            patient_id="p001",
            origin_df=_make_origin_df(),
            pred_df=_make_pred_df(),
            hla_qc_tsv=None,
            output_tsv=out,
            presentation_percentile_strong=0.5,
            tcrdock_pdb=None,
            junction_filter_stats_tsv=stats,
        )
        df = pd.read_csv(out, sep="\t")
        totals = df[df["metric"].isin([
            "junctions_extracted_total",
            "junctions_mean_reads_filtered",
            "junctions_annotated_discarded",
            "junctions_unannotated_total",
        ])]
        # All 4 rows present
        assert len(totals) == 4
        # Notes field distinguishes them from per-sample rows (which carry sample_id)
        assert (totals["notes"] == "all tumor samples").all()


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
            "sample_id", "sample_type", "unannotated",
            "normal_shared", "gtex_pantissue_shared", "tumor_exclusive",
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

    def test_top_candidates_total_and_capped_exposed(self, tmp_path):
        # _write_round_trip uses 3 strong + 1 weak = 4 gated presenters (under the cap)
        loaded = self._write_round_trip(tmp_path)
        mp = loaded["mhc_prediction"]
        assert mp["top_candidates_total"] == 4
        assert mp["top_candidates_capped"] == 4
        assert isinstance(mp["top_candidates_total"], int)
        assert isinstance(mp["top_candidates_capped"], int)

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
            "sample_id", "sample_type", "unannotated",
            "normal_shared", "gtex_pantissue_shared", "tumor_exclusive",
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


# ---------------------------------------------------------------------------
# Tests for the artefact-driven HTML helpers (Phase 2)
# ---------------------------------------------------------------------------

class TestRenderContigPeek:
    """Round-trip: bracketed plain-text peek → styled HTML with same span classes."""

    def test_empty_input_returns_empty(self):
        assert _render_contig_peek("") == ""

    def test_emits_junction_marker(self):
        # Junction-spanning peptide: | inside [...]
        html = _render_contig_peek("AAA[CC|GG]TTT")
        assert '<span class="junction-mark">|</span>' in html

    def test_classifies_upstream_downstream_outside_peptide(self):
        html = _render_contig_peek("AA[CC|GG]TT")
        # Outside [..], before |: nt-up
        assert html.count('class="nt-up">A</span>') == 2
        # Outside [..], after |: nt-down
        assert html.count('class="nt-down">T</span>') == 2

    def test_classifies_peptide_upstream_downstream_inside_peptide(self):
        html = _render_contig_peek("AA[CC|GG]TT")
        # Inside [..], before |: nt-pep-up
        assert html.count('class="nt-pep-up">C</span>') == 2
        # Inside [..], after |: nt-pep-down
        assert html.count('class="nt-pep-down">G</span>') == 2

    def test_round_trip_with_build_contig_peek(self):
        """The plain-text writer + HTML reader should preserve all nucleotides."""
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"  # 50 nt
        peek = _build_contig_peek(seq, start_nt=20, end_nt_incl=34, upstream_nt=26)
        html = _render_contig_peek(peek)
        for nt in seq:
            assert f'>{nt}</span>' in html
        assert '<span class="junction-mark">|</span>' in html


class TestBuildStrongTableHtmlFromTopCandidates:
    def _make_top_df(self, n_rows: int = 3, with_gps: bool = True) -> pd.DataFrame:
        rows = []
        for i in range(n_rows):
            row = {
                "patient_id": "p001", "rank": i + 1,
                "peptide": f"PEP{i}", "best_allele": "HLA-A*02:01",
                "best_presentation_percentile": 0.1 * (i + 1),
                "ic50_nM": float(50 + i),
                "n_strong_alleles": 1, "presentation_class": "strong",
                "contig_key": f"k{i}", "contig_start_nt": 20,
                "contig_peek": "AAAAAAAAAAAAAAAAAAAAAAAAAA[CCCC|CCCC]TTTTTTTTTTTTTTTTTT",
            }
            if with_gps:
                row["genotype_presentation_score"] = 0.99 - 0.01 * i
            rows.append(row)
        return pd.DataFrame(rows)

    def test_empty_returns_no_results_message(self):
        assert "No strong" in _build_strong_table_html_from_top_candidates(pd.DataFrame())

    def test_renders_one_row_per_candidate(self):
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(n_rows=3))
        assert html.count("<tr>") == 4  # 1 header row + 3 data rows

    def test_includes_peptide_and_allele(self):
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(n_rows=2))
        assert "PEP0" in html and "PEP1" in html
        assert "HLA-A*02:01" in html

    def test_includes_gps_column_when_present(self):
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(with_gps=True))
        assert "GPS" in html

    def test_omits_gps_column_when_absent(self):
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(with_gps=False))
        assert "GPS" not in html

    def test_renders_styled_contig_peek(self):
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(n_rows=1))
        assert '<span class="junction-mark">' in html
        assert '<span class="nt-pep-up">' in html

    def test_truncation_notice_when_total_exceeds_capped(self):
        html = _build_strong_table_html_from_top_candidates(
            self._make_top_df(n_rows=10), total=47, capped=10)
        assert "Showing 10 of 47" in html

    def test_no_truncation_notice_when_total_within_cap(self):
        html = _build_strong_table_html_from_top_candidates(
            self._make_top_df(n_rows=5), total=5, capped=5)
        assert "Showing" not in html

    def test_no_truncation_notice_when_total_unknown(self):
        # Backward-compatible default: callers without report.tsv pass nothing → no notice.
        html = _build_strong_table_html_from_top_candidates(self._make_top_df(n_rows=3))
        assert "Showing" not in html

    def test_no_truncation_notice_when_only_capped_given(self):
        # Guard requires BOTH total and capped; one-of-two must suppress the notice.
        html = _build_strong_table_html_from_top_candidates(
            self._make_top_df(n_rows=5), total=None, capped=10)
        assert "Showing" not in html

    def test_no_truncation_notice_when_only_total_given(self):
        html = _build_strong_table_html_from_top_candidates(
            self._make_top_df(n_rows=5), total=47, capped=None)
        assert "Showing" not in html


class TestBuildFilteringFunnelHtml:
    """Issue #215 — funnel renders for the unified filtering_stats.tsv."""

    def _write_unified_stats(self, tmp_path):
        path = tmp_path / "filtering_stats.tsv"
        pd.DataFrame([
            # Junction-filter — per-sample funnel + distribution rows
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "junctions_raw", "count": 100},
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "tumor_exclusive", "count": 30},
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "mean_reads", "count": 12.5},
            # Downstream patient-level steps
            {"patient_id": "p1", "sample_id": "", "sample_type": "",
             "step": "contig-assemble", "category": "contigs_written", "count": 25},
            {"patient_id": "p1", "sample_id": "", "sample_type": "",
             "step": "mhc-affinity", "category": "strong_presenters", "count": 18},
        ]).to_csv(path, sep="\t", index=False)
        return path

    def test_returns_empty_when_path_missing(self):
        assert _build_filtering_funnel_html(None) == ""

    def test_returns_empty_when_file_does_not_exist(self, tmp_path):
        assert _build_filtering_funnel_html(tmp_path / "nope.tsv") == ""

    def test_renders_section_header(self, tmp_path):
        path = self._write_unified_stats(tmp_path)
        html = _build_filtering_funnel_html(path)
        assert "Filtering funnel" in html
        assert "filtering_stats.tsv" in html

    def test_includes_per_sample_table(self, tmp_path):
        path = self._write_unified_stats(tmp_path)
        html = _build_filtering_funnel_html(path)
        assert "Junction-level (per sample)" in html
        assert "T1" in html
        assert "tumor_exclusive" in html

    def test_includes_distribution_table(self, tmp_path):
        path = self._write_unified_stats(tmp_path)
        html = _build_filtering_funnel_html(path)
        assert "Raw read-count distribution" in html
        assert "mean_reads" in html

    def test_includes_downstream_steps(self, tmp_path):
        path = self._write_unified_stats(tmp_path)
        html = _build_filtering_funnel_html(path)
        assert "Pipeline funnel" in html
        assert "contig-assemble" in html
        assert "strong_presenters" in html

    def test_no_nan_when_funnel_categories_partially_missing(self, tmp_path):
        """Issue #215 follow-up: ``reindex`` introduces NaN columns when a
        funnel/distribution category is absent from the input — fillna runs
        AFTER the reindex now so the output table never renders 'NaN'."""
        path = tmp_path / "filtering_stats.tsv"
        # Only 2 of the 5 funnel categories present (junctions_raw +
        # tumor_exclusive); the other 3 must reindex-fill to 0.
        # Only 1 of the 4 distribution categories present (mean_reads); the
        # other 3 must reindex-fill to 0.
        pd.DataFrame([
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "junctions_raw", "count": 100},
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "tumor_exclusive", "count": 30},
            {"patient_id": "p1", "sample_id": "T1", "sample_type": "Primary Tumor",
             "step": "junction-filter", "category": "mean_reads", "count": 12.5},
        ]).to_csv(path, sep="\t", index=False)

        html = _build_filtering_funnel_html(path)
        assert "NaN" not in html
        assert "T1" in html  # still rendered


class TestPresenterCountsHtml:
    def test_empty_dict_returns_no_predictions(self):
        assert "No predictions" in _presenter_counts_html({})

    def test_only_total_predictions_returns_no_predictions(self):
        assert "No predictions" in _presenter_counts_html({"total_predictions": 100})

    def test_excludes_total_predictions_row(self):
        html = _presenter_counts_html({"total_predictions": 14, "strong": 3, "weak": 1, "non": 10})
        assert "<table" in html
        assert "total_predictions" not in html
        assert "strong" in html and "weak" in html and "non" in html

    def test_excludes_top_candidates_metrics(self):
        html = _presenter_counts_html({
            "total_predictions": 14, "strong": 3, "weak": 1, "non": 10,
            "top_candidates_total": 47, "top_candidates_capped": 10,
        })
        assert "top_candidates_total" not in html
        assert "top_candidates_capped" not in html
        assert "strong" in html and "weak" in html and "non" in html

    def test_only_non_class_metrics_returns_no_predictions(self):
        html = _presenter_counts_html(
            {"total_predictions": 14, "top_candidates_total": 5, "top_candidates_capped": 5})
        assert "No predictions" in html


class TestGenerateReportEndToEnd:
    """Integration test: drives generate_report() through the artefact-driven path."""

    def _write_inputs(self, tmp_path, n_strong=2, n_weak=1, n_non=3):
        junc_tsv = tmp_path / "novel_junctions.tsv"
        pd.DataFrame([
            {"sample_id": "S1", "sample_type": "Primary Tumor",
             "junction_origin": "tumor_exclusive", "contig_key": "k0", "start_nt": 20},
            {"sample_id": "S1", "sample_type": "Primary Tumor",
             "junction_origin": "tumor_exclusive", "contig_key": "k1", "start_nt": 22},
            {"sample_id": "S1", "sample_type": "Primary Tumor",
             "junction_origin": "normal_shared", "contig_key": "k2", "start_nt": 18},
        ]).to_csv(junc_tsv, sep="\t", index=False)

        pred_df = _make_pred_df(n_strong=n_strong, n_weak=n_weak, n_non=n_non)
        pred_tsv = tmp_path / "predictions.tsv"
        pred_df.to_csv(pred_tsv, sep="\t", index=False)

        # One contig per distinct prediction key so the artefact writer finds a
        # sequence for every candidate (default keys k0/k1/k2 stay unchanged).
        contigs_fa = tmp_path / "contigs.fa"
        keys = sorted(set(pred_df["contig_key"]), key=lambda k: (len(k), k))
        contigs_fa.write_text("".join(f">{k}\n" + "ACGTACGT" * 7 + "\n" for k in keys))
        return junc_tsv, pred_tsv, contigs_fa

    def test_full_pipeline_writes_html_and_artefacts(self, tmp_path):
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        html_out = out_dir / "report.html"
        tsv_out = out_dir / "report.tsv"
        top_out = out_dir / "report_top_candidates.tsv"

        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(tsv_out),
            output_top_candidates_tsv=str(top_out),
            patient_id="p001",
        )

        assert html_out.exists()
        assert tsv_out.exists()
        assert top_out.exists()

        html = html_out.read_text()
        # Surface signal: vocabulary and key sections rendered
        assert "Top strong presentations" in html
        assert "presenters" in html  # pipeline-diagram label, presenter vocabulary
        assert "binder" not in html.lower() or "Binding affinity" in html  # IC50 tooltip allowed
        # Junction filtering rendered
        assert "Primary Tumor" in html
        # Top presenter peptide visible (PEP0 has lowest presentation_percentile)
        assert "PEP0" in html

    def test_html_reflects_artefact_contents_not_raw_pred_df(self, tmp_path, monkeypatch):
        """Prove the strong-presenters HTML section reads from
        report_top_candidates.tsv, not the raw pred_df.

        Wrap the artefact writer so that AFTER it produces the file we inject a
        sentinel peptide name into it. The HTML rendering happens after the
        writer and (if the artefact path is genuinely active) reads from the
        modified file. The sentinel can only appear in the HTML if the renderer
        used the artefact — the raw pred_df never contains it.
        """
        import generate_report as gr_mod

        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        html_out = out_dir / "report.html"
        tsv_out = out_dir / "report.tsv"
        top_out = out_dir / "report_top_candidates.tsv"

        sentinel = "SENTINELPEPTIDE"
        original_writer = gr_mod._build_report_top_candidates_tsv

        def writer_with_sentinel(*args, **kwargs):
            original_writer(*args, **kwargs)
            out_path = kwargs["output_tsv"]
            df = pd.read_csv(out_path, sep="\t")
            assert not df.empty, (
                "test data must produce ≥1 presenter for sentinel injection — "
                "if the artefact is empty, the sentinel is never written and "
                "the assertion below would fail for the wrong reason"
            )
            df.loc[0, "peptide"] = sentinel
            df.to_csv(out_path, sep="\t", index=False)

        monkeypatch.setattr(
            gr_mod, "_build_report_top_candidates_tsv", writer_with_sentinel
        )

        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(tsv_out),
            output_top_candidates_tsv=str(top_out),
            patient_id="p001",
        )

        html = html_out.read_text()
        assert sentinel in html, (
            "HTML strong-presenters section must reflect artefact contents — "
            "sentinel peptide was injected into report_top_candidates.tsv "
            "after the writer ran, so its absence proves HTML is rendering "
            "from raw pred_df instead of the artefact."
        )

    def test_full_pipeline_renders_truncation_notice_when_over_cap(self, tmp_path):
        """End-to-end seam: with >TOP_CANDIDATES_LIMIT candidates, the notice must
        flow from _build_report_tsv → report.tsv → _load_report_tsv → call site →
        the artefact-driven renderer's HTML (Issue #226)."""
        # 12 strong + 1 weak = 13 quality-gated candidates, capped to 10.
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(
            tmp_path, n_strong=12, n_weak=1, n_non=2)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        html_out = out_dir / "report.html"
        tsv_out = out_dir / "report.tsv"
        top_out = out_dir / "report_top_candidates.tsv"

        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(tsv_out),
            output_top_candidates_tsv=str(top_out),
            patient_id="p001",
        )

        html = html_out.read_text()
        assert f"Showing {TOP_CANDIDATES_LIMIT} of 13 candidates." in html

    def test_full_pipeline_omits_notice_when_no_report_tsv(self, tmp_path):
        """Backward-compat seam: over-cap data but no report.tsv (output_tsv=None)
        means the renderer gets no total/capped and must omit the notice."""
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(
            tmp_path, n_strong=12, n_weak=1, n_non=2)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        html_out = out_dir / "report.html"
        top_out = out_dir / "report_top_candidates.tsv"

        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=None,
            output_top_candidates_tsv=str(top_out),
            patient_id="p001",
        )

        html = html_out.read_text()
        assert "Showing" not in html


# ---------------------------------------------------------------------------
# Issue #206 — VDJdb panel reference table + matched-TCR provenance in report
# ---------------------------------------------------------------------------

# Minimal 4-chain PDB (A=MHC, B=peptide, C=TCRα, D=TCRβ) for structure tests.
_MINI_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      2  CA  GLY B   1       1.000   1.000   1.000  1.00  0.00           C\n"
    "ATOM      3  CA  SER C   1       2.000   2.000   2.000  1.00  0.00           C\n"
    "ATOM      4  CA  LEU D   1       3.000   3.000   3.000  1.00  0.00           C\n"
)


def _make_panel_df() -> pd.DataFrame:
    """Two-allele VDJdb panel mirroring fetch_vdjdb_panel.PANEL_COLUMNS."""
    return pd.DataFrame([
        {"allele": "HLA-C*07:01", "va_gene": "TRAV29/DV5*01", "ja_gene": "TRAJ29*01",
         "cdr3a": "CAANSGNTPLVF", "vb_gene": "TRBV2*01", "jb_gene": "TRBJ1-1*01",
         "cdr3b": "CASSLMGGGNTIYF", "alpha_seq": "AAAA", "beta_seq": "BBBB",
         "vdjdb_score": 2, "vdjdb_donor_id": "donor_42"},
        {"allele": "HLA-A*02:01", "va_gene": "TRAV12-2*01", "ja_gene": "TRAJ33*01",
         "cdr3a": "CAVNDSWGKLQF", "vb_gene": "TRBV6-5*01", "jb_gene": "TRBJ2-7*01",
         "cdr3b": "CASSYQETQYF", "alpha_seq": "CCCC", "beta_seq": "DDDD",
         "vdjdb_score": 1, "vdjdb_donor_id": "donor_7"},
    ])


def _make_qc_df() -> pd.DataFrame:
    """Per-allele QC mirroring fetch_vdjdb_panel.QC_COLUMNS, with all 3 statuses."""
    return pd.DataFrame([
        {"allele": "HLA-C*07:01", "n_exact_matches": 1, "n_in_panel": 1,
         "panel_status": "low_coverage"},
        {"allele": "HLA-A*02:01", "n_exact_matches": 50, "n_in_panel": 10,
         "panel_status": "ok"},
        {"allele": "HLA-B*07:02", "n_exact_matches": 0, "n_in_panel": 0,
         "panel_status": "empty"},
    ])


def _make_docking_scores_df(panel_status: str = "vdjdb_matched") -> pd.DataFrame:
    """docking_scores.tsv top row mirroring run_tcrdock.collect_outputs output."""
    matched = panel_status != "dmf5_fallback"
    return pd.DataFrame([{
        "peptide": "SQIPRTHSY", "mhc": "HLA-C*07:01",
        "model_2_ptm_plddt": 91.3, "model_2_ptm_pae": 5.2, "pmhc_tcr_pae": 6.1,
        "tcr_va": "TRAV29/DV5*01" if matched else "TRAV12-2",
        "tcr_ja": "TRAJ29*01" if matched else "TRAJ23",
        "tcr_vb": "TRBV2*01" if matched else "TRBV6-4",
        "tcr_jb": "TRBJ1-1*01" if matched else "TRBJ1-1",
        "tcr_cdr3a": "CAANSGNTPLVF" if matched else "CAVRPGGAGPFFVVF",
        "tcr_cdr3b": "CASSLMGGGNTIYF" if matched else "CASSLSFGTEAFF",
        "vdjdb_donor_id": "donor_42" if matched else "",
        "vdjdb_score": 2 if matched else pd.NA,
        "panel_status": panel_status,
    }])


class TestBuildVdjdbPanelSection:
    def test_both_none_returns_empty(self):
        assert _build_vdjdb_panel_section(None, None) == ""

    def test_missing_files_return_empty(self, tmp_path):
        assert _build_vdjdb_panel_section(
            str(tmp_path / "nope.tsv"), str(tmp_path / "nope_qc.tsv")
        ) == ""

    def test_renders_coverage_table_with_status_badges(self, tmp_path):
        qc = tmp_path / "panel_qc.tsv"
        _make_qc_df().to_csv(qc, sep="\t", index=False)
        html = _build_vdjdb_panel_section(None, str(qc))
        assert "Per-allele coverage" in html
        assert "low_coverage" in html and "empty" in html and "ok" in html
        # Status colours present (ok=green, low=orange, empty=red)
        assert "#27ae60" in html and "#e67e22" in html and "#c0392b" in html

    def test_renders_reference_panel_table(self, tmp_path):
        panel = tmp_path / "panel.tsv"
        _make_panel_df().to_csv(panel, sep="\t", index=False)
        html = _build_vdjdb_panel_section(str(panel), None)
        assert "Reference panel" in html
        assert "All 2 TCRs" in html
        assert "CAANSGNTPLVF" in html  # CDR3α
        assert "CASSLMGGGNTIYF" in html  # CDR3β
        assert "donor_42" in html
        assert "TRAV29/DV5*01" in html

    def test_renders_both_tables_and_section_header(self, tmp_path):
        panel = tmp_path / "panel.tsv"
        qc = tmp_path / "panel_qc.tsv"
        _make_panel_df().to_csv(panel, sep="\t", index=False)
        _make_qc_df().to_csv(qc, sep="\t", index=False)
        html = _build_vdjdb_panel_section(str(panel), str(qc))
        assert "VDJdb TCR panel (reference)" in html
        assert "Per-allele coverage" in html
        assert "Reference panel" in html

    def test_empty_dataframes_return_empty(self, tmp_path):
        panel = tmp_path / "panel.tsv"
        qc = tmp_path / "panel_qc.tsv"
        _make_panel_df().head(0).to_csv(panel, sep="\t", index=False)
        _make_qc_df().head(0).to_csv(qc, sep="\t", index=False)
        assert _build_vdjdb_panel_section(str(panel), str(qc)) == ""


class TestBuildTcrProvenanceHtml:
    def test_none_meta_returns_legacy_note(self):
        html = _build_tcr_provenance_html("HLA-C*07:01", None)
        assert "DMF5 fallback (see config)" in html

    def test_meta_without_panel_status_returns_legacy_note(self):
        html = _build_tcr_provenance_html("HLA-C*07:01", {"tcr_va": "X"})
        assert "DMF5 fallback (see config)" in html

    def test_dmf5_fallback_flags_warning_with_allele(self):
        html = _build_tcr_provenance_html(
            "HLA-C*07:01", {"panel_status": "dmf5_fallback"}
        )
        assert "tcr-fallback" in html
        assert "DMF5 fallback" in html
        assert "HLA-C*07:01" in html
        assert "sanity check" in html

    def test_vdjdb_matched_shows_full_provenance(self):
        meta = _make_docking_scores_df("vdjdb_matched").iloc[0].to_dict()
        html = _build_tcr_provenance_html("HLA-C*07:01", meta)
        assert "tcr-matched" in html
        assert "HLA-matched" in html
        assert "TRAV29/DV5*01" in html and "TRAJ29*01" in html
        assert "CAANSGNTPLVF" in html and "CASSLMGGGNTIYF" in html
        assert "donor_42" in html
        assert "confidence score 2" in html

    def test_matched_with_na_score_omits_score_phrase(self):
        meta = {"panel_status": "vdjdb_matched", "tcr_va": "TRAV1",
                "vdjdb_donor_id": "d1", "vdjdb_score": pd.NA}
        html = _build_tcr_provenance_html("HLA-A*02:01", meta)
        assert "HLA-matched" in html
        assert "confidence score" not in html  # NA score suppressed
        assert "donor" in html


class TestBuildDockingMetricsHtml:
    def test_none_returns_empty(self):
        assert _build_docking_metrics_html(None) == ""

    def test_no_metric_columns_returns_empty(self):
        assert _build_docking_metrics_html({"peptide": "X", "panel_status": "ok"}) == ""

    def test_renders_plddt_pae_ptm_rows(self):
        meta = _make_docking_scores_df().iloc[0].to_dict()
        html = _build_docking_metrics_html(meta)
        assert "TCRdock confidence metrics" in html
        assert "model_2_ptm_plddt" in html and "91.3" in html
        assert "pmhc_tcr_pae" in html

    def test_skips_nan_and_blank_metric_values(self):
        html = _build_docking_metrics_html(
            {"model_2_ptm_plddt": float("nan"), "model_2_ptm_pae": ""}
        )
        assert html == ""

    def test_skips_pd_na_metric_value(self):
        # pd.NA is not a float, so the guard must use pd.isna (not isinstance).
        html = _build_docking_metrics_html({"model_2_ptm_plddt": pd.NA})
        assert html == ""


class TestBuildReport3dStructureProvenance:
    def test_manifest_carries_provenance_columns(self, tmp_path):
        scores = tmp_path / "docking_scores.tsv"
        _make_docking_scores_df("vdjdb_matched").to_csv(scores, sep="\t", index=False)
        manifest = tmp_path / "report_3d_structure.tsv"
        _build_report_3d_structure_tsv(
            patient_id="p001",
            docking_scores_tsv=str(scores),
            pdb_relative_path="predictions/tcrdock/top_candidate.pdb",
            output_tsv=str(manifest),
        )
        df = pd.read_csv(manifest, sep="\t")
        assert "panel_status" in df.columns
        assert df.iloc[0]["panel_status"] == "vdjdb_matched"
        assert df.iloc[0]["vdjdb_donor_id"] == "donor_42"
        assert "model_2_ptm_plddt" in df.columns  # metrics still pass through

    def test_manifest_without_provenance_columns_is_backward_compatible(self, tmp_path):
        # Pre-#205 docking_scores.tsv: no provenance columns at all.
        scores = tmp_path / "docking_scores.tsv"
        pd.DataFrame([{"peptide": "SQIPRTHSY", "mhc": "HLA-C*07:01",
                       "model_2_ptm_plddt": 88.0}]).to_csv(scores, sep="\t", index=False)
        manifest = tmp_path / "report_3d_structure.tsv"
        _build_report_3d_structure_tsv(
            patient_id="p001",
            docking_scores_tsv=str(scores),
            pdb_relative_path="x.pdb",
            output_tsv=str(manifest),
        )
        df = pd.read_csv(manifest, sep="\t")
        assert "panel_status" not in df.columns
        assert "model_2_ptm_plddt" in df.columns

    def test_dmf5_fallback_na_score_blanked_in_manifest(self, tmp_path):
        scores = tmp_path / "docking_scores.tsv"
        _make_docking_scores_df("dmf5_fallback").to_csv(scores, sep="\t", index=False)
        manifest = tmp_path / "report_3d_structure.tsv"
        _build_report_3d_structure_tsv(
            patient_id="p001",
            docking_scores_tsv=str(scores),
            pdb_relative_path="x.pdb",
            output_tsv=str(manifest),
        )
        df = pd.read_csv(manifest, sep="\t")
        assert df.iloc[0]["panel_status"] == "dmf5_fallback"
        # NA score round-trips as blank (NaN), not a stray literal.
        assert pd.isna(df.iloc[0]["vdjdb_score"]) or df.iloc[0]["vdjdb_score"] == ""


class TestBuildStructureSectionWithMeta:
    def test_matched_meta_renders_provenance_and_metrics(self, tmp_path):
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        meta = _make_docking_scores_df("vdjdb_matched").iloc[0].to_dict()
        html = _build_structure_section(pdb, "SQIPRTHSY", "HLA-C*07:01", meta=meta)
        assert "HLA-matched" in html
        assert "donor_42" in html
        assert "TCRdock confidence metrics" in html
        assert "91.3" in html
        assert "DMF5 fallback (see config)" not in html

    def test_fallback_meta_flags_dmf5(self, tmp_path):
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        meta = _make_docking_scores_df("dmf5_fallback").iloc[0].to_dict()
        html = _build_structure_section(pdb, "SQIPRTHSY", "HLA-C*07:01", meta=meta)
        assert "tcr-fallback" in html
        assert "DMF5 fallback" in html

    def test_no_meta_keeps_legacy_note(self, tmp_path):
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        html = _build_structure_section(pdb, "SQIPRTHSY", "HLA-C*07:01")
        assert "DMF5 fallback (see config)" in html


class TestGenerateReportWithPanelAndStructure:
    """End-to-end: panel section + matched-TCR provenance flow into report.html."""

    def _write_inputs(self, tmp_path):
        junc_tsv = tmp_path / "novel_junctions.tsv"
        pd.DataFrame([
            {"sample_id": "S1", "sample_type": "Primary Tumor",
             "junction_origin": "tumor_exclusive", "contig_key": "k0", "start_nt": 20},
        ]).to_csv(junc_tsv, sep="\t", index=False)
        pred_df = _make_pred_df(n_strong=2, n_weak=1, n_non=2)
        pred_tsv = tmp_path / "predictions.tsv"
        pred_df.to_csv(pred_tsv, sep="\t", index=False)
        contigs_fa = tmp_path / "contigs.fa"
        keys = sorted(set(pred_df["contig_key"]), key=lambda k: (len(k), k))
        contigs_fa.write_text("".join(f">{k}\n" + "ACGTACGT" * 7 + "\n" for k in keys))
        return junc_tsv, pred_tsv, contigs_fa

    def test_panel_section_appears_in_html(self, tmp_path):
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        panel = tmp_path / "panel.tsv"
        qc = tmp_path / "panel_qc.tsv"
        _make_panel_df().to_csv(panel, sep="\t", index=False)
        _make_qc_df().to_csv(qc, sep="\t", index=False)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            vdjdb_panel_tsv=str(panel),
            vdjdb_panel_qc_tsv=str(qc),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert "VDJdb TCR panel (reference)" in html
        assert "CAANSGNTPLVF" in html
        assert "low_coverage" in html

    def test_structure_section_surfaces_matched_tcr_from_manifest(self, tmp_path):
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        scores = tmp_path / "docking_scores.tsv"
        _make_docking_scores_df("vdjdb_matched").to_csv(scores, sep="\t", index=False)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            output_3d_structure_tsv=str(out_dir / "report_3d_structure.tsv"),
            tcrdock_pdb=str(pdb),
            docking_scores_tsv=str(scores),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert "HLA-matched" in html
        assert "donor_42" in html
        assert "TCRdock confidence metrics" in html
        # Manifest carries the provenance the HTML rendered from.
        manifest = pd.read_csv(out_dir / "report_3d_structure.tsv", sep="\t")
        assert manifest.iloc[0]["panel_status"] == "vdjdb_matched"

    def test_structure_section_flags_dmf5_fallback_end_to_end(self, tmp_path):
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        scores = tmp_path / "docking_scores.tsv"
        _make_docking_scores_df("dmf5_fallback").to_csv(scores, sep="\t", index=False)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            output_3d_structure_tsv=str(out_dir / "report_3d_structure.tsv"),
            tcrdock_pdb=str(pdb),
            docking_scores_tsv=str(scores),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert "tcr-fallback" in html
        assert "DMF5 fallback" in html

    def test_panel_section_reflects_artefact_contents(self, tmp_path):
        """Prove the panel section renders the *file* content, not hardcoded
        values: a unique sentinel CDR3 + donor written into panel.tsv must
        surface in the HTML (the panel section reads panel.tsv directly)."""
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        panel = tmp_path / "panel.tsv"
        qc = tmp_path / "panel_qc.tsv"
        panel_df = _make_panel_df()
        panel_df.loc[0, "cdr3a"] = "SENTINELCDR3AXYZ"
        panel_df.loc[0, "vdjdb_donor_id"] = "SENTINELDONOR999"
        panel_df.to_csv(panel, sep="\t", index=False)
        _make_qc_df().to_csv(qc, sep="\t", index=False)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            vdjdb_panel_tsv=str(panel),
            vdjdb_panel_qc_tsv=str(qc),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert "SENTINELCDR3AXYZ" in html
        assert "SENTINELDONOR999" in html

    def test_generate_report_omits_panel_when_none(self, tmp_path):
        """Backward-compat / HLA-disabled path: with no panel artefacts the
        report still renders, just without the VDJdb panel section."""
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            vdjdb_panel_tsv=None,
            vdjdb_panel_qc_tsv=None,
            patient_id="p001",
        )
        html = html_out.read_text()
        assert html_out.exists()
        assert "VDJdb TCR panel" not in html
        assert "Top strong presentations" in html  # other sections intact

    def test_structure_section_with_legacy_docking_scores(self, tmp_path):
        """Backward-compat end-to-end: a pre-#205 docking_scores.tsv with NO
        provenance columns must flow through generate_report() to the legacy
        DMF5 note, with no HLA-matched provenance leaking in (AC4 fallback)."""
        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        scores = tmp_path / "docking_scores.tsv"
        # Pre-#205: only peptide, mhc, and a metric column — no provenance.
        pd.DataFrame([{"peptide": "SQIPRTHSY", "mhc": "HLA-C*07:01",
                       "model_2_ptm_plddt": 88.0}]).to_csv(scores, sep="\t", index=False)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            output_3d_structure_tsv=str(out_dir / "report_3d_structure.tsv"),
            tcrdock_pdb=str(pdb),
            docking_scores_tsv=str(scores),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert "DMF5 fallback (see config)" in html  # legacy note
        assert "HLA-matched" not in html
        assert "VDJdb donor" not in html
        # Metrics still pass through even on a legacy scores file.
        assert "TCRdock confidence metrics" in html

    def test_structure_section_reads_provenance_from_manifest_artefact(self, tmp_path, monkeypatch):
        """Prove the structure section's provenance reads from the reloaded
        report_3d_structure.tsv, not from docking_scores.tsv directly: inject a
        sentinel donor into the manifest AFTER the writer runs; it can only reach
        the HTML if the renderer reads the manifest artefact."""
        import generate_report as gr_mod

        junc_tsv, pred_tsv, contigs_fa = self._write_inputs(tmp_path)
        out_dir = tmp_path / "results" / "p001" / "reports"
        out_dir.mkdir(parents=True)
        pdb = tmp_path / "top_candidate.pdb"
        pdb.write_text(_MINI_PDB)
        scores = tmp_path / "docking_scores.tsv"
        _make_docking_scores_df("vdjdb_matched").to_csv(scores, sep="\t", index=False)

        sentinel = "SENTINELDONOR777"
        original_writer = gr_mod._build_report_3d_structure_tsv

        def writer_with_sentinel(*args, **kwargs):
            original_writer(*args, **kwargs)
            out_path = kwargs["output_tsv"]
            df = pd.read_csv(out_path, sep="\t")
            assert not df.empty and "vdjdb_donor_id" in df.columns
            df.loc[0, "vdjdb_donor_id"] = sentinel
            df.to_csv(out_path, sep="\t", index=False)

        monkeypatch.setattr(gr_mod, "_build_report_3d_structure_tsv", writer_with_sentinel)

        html_out = out_dir / "report.html"
        generate_report(
            novel_junctions_tsv=str(junc_tsv),
            predictions_tsv=str(pred_tsv),
            output_html=str(html_out),
            contigs_fasta=str(contigs_fa),
            output_tsv=str(out_dir / "report.tsv"),
            output_top_candidates_tsv=str(out_dir / "report_top_candidates.tsv"),
            output_3d_structure_tsv=str(out_dir / "report_3d_structure.tsv"),
            tcrdock_pdb=str(pdb),
            docking_scores_tsv=str(scores),
            patient_id="p001",
        )
        html = html_out.read_text()
        assert sentinel in html, (
            "structure section must render TCR provenance from the reloaded "
            "report_3d_structure.tsv artefact — sentinel injected post-write "
            "is absent, so the HTML is not reading from the manifest."
        )
