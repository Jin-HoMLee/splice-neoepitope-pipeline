"""Tests for generate_report.py — HTML generation helpers."""

import pandas as pd

from generate_report import (
    _build_chain_legend,
    _build_compnd_records,
    _df_to_html,
    _extract_chain_ids,
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
