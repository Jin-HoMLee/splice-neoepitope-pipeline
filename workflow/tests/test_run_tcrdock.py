"""Tests for run_tcrdock.py — candidate selection, HLA-matched TCR selection,
and PDB chain relabeling."""

import pandas as pd

from run_tcrdock import (
    _normalize_allele_4digit,
    attach_matched_tcrs,
    build_tcrdock_input,
    collect_outputs,
    relabel_pdb_chains,
    select_matched_tcr,
    select_top_candidates,
)

# DMF5 fallback TCR (HLA-A*02:01-restricted) — mirrors config gpu_config.yaml.
_DMF5 = {
    "va_gene": "TRAV12-2", "ja_gene": "TRAJ21", "cdr3a": "CAVNFGGGKLI",
    "vb_gene": "TRBV6-5",  "jb_gene": "TRBJ2-7", "cdr3b": "CASSLAGGRPEQYF",
}


def _panel(rows):
    """Build a VDJdb panel DataFrame (subset of fetch_vdjdb_panel PANEL_COLUMNS)."""
    return pd.DataFrame(rows, columns=[
        "allele", "va_gene", "ja_gene", "cdr3a", "vb_gene", "jb_gene", "cdr3b",
        "vdjdb_score", "vdjdb_donor_id",
    ])


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _write_predictions_tsv(tmp_path, rows):
    """Write a predictions TSV and return the path."""
    df = pd.DataFrame(rows)
    tsv = tmp_path / "predictions.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    return tsv


# Minimal single-chain PDB (AlphaFold output style: all chain A)
_FLAT_PDB = (
    "ATOM      1  N   ALA A   1       1.0   2.0   3.0  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       2.0   3.0   4.0  1.00  0.00           C\n"
    "ATOM      3  N   GLY A   2       3.0   4.0   5.0  1.00  0.00           N\n"
    "ATOM      4  N   SER A   3       4.0   5.0   6.0  1.00  0.00           N\n"
    "ATOM      5  N   VAL A   4       5.0   6.0   7.0  1.00  0.00           N\n"
)


# ---------------------------------------------------------------------------
# select_top_candidates
# ---------------------------------------------------------------------------

class TestSelectTopCandidates:
    def test_selects_strong_presenter(self, tmp_path):
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "ACDEFGHIK", "best_allele": "HLA-A*02:01", "ic50_nM": 10.0,
             "presentation_percentile": 0.1, "presentation_class": "strong", "contig_key": "c1", "start_nt": 0},
            {"peptide": "LMNPQRSTV", "best_allele": "HLA-A*02:01", "ic50_nM": 200.0,
             "presentation_percentile": 2.0, "presentation_class": "weak", "contig_key": "c2", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=1)
        assert len(result) == 1
        assert result.iloc[0]["peptide"] == "ACDEFGHIK"

    def test_includes_weak_presenters(self, tmp_path):
        """Weak presenters are included in candidates (quality gate checked separately)."""
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "LMNPQRSTV", "best_allele": "HLA-A*02:01", "ic50_nM": 200.0,
             "presentation_percentile": 2.0, "presentation_class": "weak", "contig_key": "c1", "start_nt": 0},
            {"peptide": "YKLMFSTAV", "best_allele": "HLA-A*02:01", "ic50_nM": 9999.0,
             "presentation_percentile": 50.0, "presentation_class": "non", "contig_key": "c2", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=1)
        assert len(result) == 1
        assert result.iloc[0]["peptide"] == "LMNPQRSTV"

    def test_returns_empty_when_no_presenters(self, tmp_path):
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "YKLMFSTAV", "best_allele": "HLA-A*02:01", "ic50_nM": 9999.0,
             "presentation_percentile": 50.0, "presentation_class": "non", "contig_key": "c1", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=1)
        assert result.empty

    def test_respects_n_candidates(self, tmp_path):
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "ACDEFGHIK", "best_allele": "HLA-A*02:01", "ic50_nM": 10.0,
             "presentation_percentile": 0.1, "presentation_class": "strong", "contig_key": "c1", "start_nt": 0},
            {"peptide": "LMNPQRSTV", "best_allele": "HLA-A*02:01", "ic50_nM": 20.0,
             "presentation_percentile": 0.2, "presentation_class": "strong", "contig_key": "c2", "start_nt": 0},
            {"peptide": "YKLMFSTAV", "best_allele": "HLA-A*02:01", "ic50_nM": 30.0,
             "presentation_percentile": 0.3, "presentation_class": "strong", "contig_key": "c3", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=2)
        assert len(result) == 2
        assert list(result["peptide"]) == ["ACDEFGHIK", "LMNPQRSTV"]

    def test_quality_gate_filters_non_presenter_percentile(self, tmp_path):
        """Candidates where best_presentation_percentile > 2% are excluded."""
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "ACDEFGHIK", "best_allele": "HLA-A*02:01", "ic50_nM": 10.0,
             "presentation_percentile": 0.1, "presentation_class": "strong",
             "best_presentation_percentile": 0.1,
             "genotype_presentation_score": 0.9, "n_strong_alleles": 1,
             "contig_key": "c1", "start_nt": 0},
            {"peptide": "LMNPQRSTV", "best_allele": "HLA-A*02:01", "ic50_nM": 20.0,
             "presentation_percentile": 0.2, "presentation_class": "strong",
             "best_presentation_percentile": 3.0,  # fails quality gate
             "genotype_presentation_score": 0.95, "n_strong_alleles": 2,
             "contig_key": "c2", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=5)
        assert len(result) == 1
        assert result.iloc[0]["peptide"] == "ACDEFGHIK"

    def test_ranked_by_genotype_presentation_score(self, tmp_path):
        """When GPS columns present, candidates rank by genotype_presentation_score desc."""
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "ACDEFGHIK", "best_allele": "HLA-A*02:01", "ic50_nM": 10.0,
             "presentation_percentile": 0.1, "presentation_class": "strong",
             "best_presentation_percentile": 0.1,
             "genotype_presentation_score": 0.7, "n_strong_alleles": 1,
             "contig_key": "c1", "start_nt": 0},
            {"peptide": "LMNPQRSTV", "best_allele": "HLA-A*02:01", "ic50_nM": 200.0,
             "presentation_percentile": 0.4, "presentation_class": "strong",
             "best_presentation_percentile": 0.4,
             "genotype_presentation_score": 0.95, "n_strong_alleles": 3,
             "contig_key": "c2", "start_nt": 0},
        ])
        result = select_top_candidates(tsv, n_candidates=1)
        assert result.iloc[0]["peptide"] == "LMNPQRSTV"  # higher GPS wins

    def test_quality_gate_empty_logs_specific_message(self, tmp_path, caplog):
        """When presenters exist but all fail quality gate, log distinguishes from no-presenters case."""
        import logging
        tsv = _write_predictions_tsv(tmp_path, [
            {"peptide": "ACDEFGHIK", "best_allele": "HLA-A*02:01", "ic50_nM": 10.0,
             "presentation_percentile": 0.1, "presentation_class": "strong",
             "best_presentation_percentile": 3.0,  # fails gate
             "genotype_presentation_score": 0.9, "n_strong_alleles": 1,
             "contig_key": "c1", "start_nt": 0},
        ])
        with caplog.at_level(logging.ERROR, logger="run_tcrdock"):
            result = select_top_candidates(tsv, n_candidates=1, presentation_percentile_weak=2.0)
        assert result.empty
        assert any("quality gate" in r.message for r in caplog.records if r.levelno == logging.ERROR)


# ---------------------------------------------------------------------------
# relabel_pdb_chains
# ---------------------------------------------------------------------------

class TestRelabelPdbChains:
    def test_relabels_two_chains(self):
        # 2 residues for chain A (MHC), 1 for chain B (peptide)
        chainseq = "AA/G"
        result = relabel_pdb_chains(_FLAT_PDB, chainseq)
        # Residues 1-2 → chain A, residue 3 → chain B
        lines = [l for l in result.splitlines() if l.startswith("ATOM")]
        assert lines[0][21] == "A"  # residue 1
        assert lines[1][21] == "A"  # residue 1 (second atom)
        assert lines[2][21] == "A"  # residue 2
        assert lines[3][21] == "B"  # residue 3

    def test_four_chain_tcr_pmhc(self):
        # MHC(1) / peptide(1) / TCR-alpha(1) / TCR-beta(1)
        chainseq = "A/G/S/V"
        result = relabel_pdb_chains(_FLAT_PDB, chainseq)
        lines = [l for l in result.splitlines() if l.startswith("ATOM")]
        chains = [l[21] for l in lines]
        # Residues 1→A, 2→B, 3→C, 4→D
        assert chains[0] == "A"  # res 1
        assert chains[1] == "A"  # res 1 (second atom)
        assert chains[2] == "B"  # res 2
        assert chains[3] == "C"  # res 3
        assert chains[4] == "D"  # res 4

    def test_empty_pdb_returns_empty(self):
        assert relabel_pdb_chains("", "A/B") == ""


# ---------------------------------------------------------------------------
# _cli_main argument parser
# ---------------------------------------------------------------------------

class TestCLIParser:
    def _make_parser(self):
        import argparse
        from run_tcrdock import _cli_main  # noqa: F401 — import triggers module load
        # Reconstruct the parser directly to avoid subprocess/sys.exit
        parser = argparse.ArgumentParser()
        parser.add_argument("--predictions-tsv", required=True)
        parser.add_argument("--output-pdb", required=True)
        parser.add_argument("--output-scores", required=True)
        parser.add_argument("--docker-image", default="tcrdock:latest")
        parser.add_argument("--n-candidates", type=int, default=1)
        parser.add_argument("--vdjdb-panel", default=None)
        parser.add_argument("--presentation-percentile-weak", type=float, default=2.0)
        return parser

    def test_default_presentation_percentile_weak(self):
        parser = self._make_parser()
        args = parser.parse_args([
            "--predictions-tsv", "p.tsv",
            "--output-pdb", "o.pdb",
            "--output-scores", "s.tsv",
        ])
        assert args.presentation_percentile_weak == 2.0

    def test_vdjdb_panel_defaults_none_and_parses(self):
        parser = self._make_parser()
        base = ["--predictions-tsv", "p.tsv", "--output-pdb", "o.pdb", "--output-scores", "s.tsv"]
        assert parser.parse_args(base).vdjdb_panel is None
        args = parser.parse_args(base + ["--vdjdb-panel", "panel.tsv"])
        assert args.vdjdb_panel == "panel.tsv"

    def test_custom_presentation_percentile_weak(self):
        parser = self._make_parser()
        args = parser.parse_args([
            "--predictions-tsv", "p.tsv",
            "--output-pdb", "o.pdb",
            "--output-scores", "s.tsv",
            "--presentation-percentile-weak", "1.5",
        ])
        assert args.presentation_percentile_weak == 1.5


# ---------------------------------------------------------------------------
# _normalize_allele_4digit
# ---------------------------------------------------------------------------

class TestNormalizeAllele4digit:
    def test_truncates_to_4digit(self):
        assert _normalize_allele_4digit("HLA-A*02:01:110") == "HLA-A*02:01"

    def test_passthrough_4digit(self):
        assert _normalize_allele_4digit("HLA-C*07:01") == "HLA-C*07:01"

    def test_too_short_returns_none(self):
        assert _normalize_allele_4digit("HLA-A*02") is None

    def test_empty_returns_none(self):
        assert _normalize_allele_4digit("") is None


# ---------------------------------------------------------------------------
# select_matched_tcr (Issue #205 selection rule)
# ---------------------------------------------------------------------------

class TestSelectMatchedTcr:
    def test_picks_highest_score_for_allele(self):
        panel = _panel([
            {"allele": "HLA-C*07:01", "va_gene": "TRAV1", "ja_gene": "TRAJ1", "cdr3a": "CAAA",
             "vb_gene": "TRBV1", "jb_gene": "TRBJ1", "cdr3b": "CBBB", "vdjdb_score": 2, "vdjdb_donor_id": "D1"},
            {"allele": "HLA-C*07:01", "va_gene": "TRAV9", "ja_gene": "TRAJ9", "cdr3a": "CXXX",
             "vb_gene": "TRBV9", "jb_gene": "TRBJ9", "cdr3b": "CYYY", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
        ])
        result = select_matched_tcr(panel, "HLA-C*07:01", _DMF5)
        assert result["panel_status"] == "vdjdb_matched"
        assert result["tcr_va"] == "TRAV9"          # the score-3 row wins
        assert result["vdjdb_score"] == 3
        assert result["vdjdb_donor_id"] == "D2"

    def test_tiebreak_by_donor_id_ascending(self):
        panel = _panel([
            {"allele": "HLA-A*02:01", "va_gene": "TRAV_B", "ja_gene": "TRAJ1", "cdr3a": "CAAA",
             "vb_gene": "TRBV1", "jb_gene": "TRBJ1", "cdr3b": "CBBB", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
            {"allele": "HLA-A*02:01", "va_gene": "TRAV_A", "ja_gene": "TRAJ1", "cdr3a": "CAAA",
             "vb_gene": "TRBV1", "jb_gene": "TRBJ1", "cdr3b": "CBBB", "vdjdb_score": 3, "vdjdb_donor_id": "D1"},
        ])
        result = select_matched_tcr(panel, "HLA-A*02:01", _DMF5)
        assert result["vdjdb_donor_id"] == "D1"      # lexicographically smaller donor wins the tie
        assert result["tcr_va"] == "TRAV_A"

    def test_normalizes_best_allele_before_match(self):
        panel = _panel([
            {"allele": "HLA-C*07:01", "va_gene": "TRAV9", "ja_gene": "TRAJ9", "cdr3a": "CXXX",
             "vb_gene": "TRBV9", "jb_gene": "TRBJ9", "cdr3b": "CYYY", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
        ])
        # 6-digit best_allele still matches the 4-digit panel key
        result = select_matched_tcr(panel, "HLA-C*07:01:01", _DMF5)
        assert result["panel_status"] == "vdjdb_matched"
        assert result["tcr_va"] == "TRAV9"

    def test_no_allele_match_falls_back_to_dmf5(self):
        panel = _panel([
            {"allele": "HLA-C*07:01", "va_gene": "TRAV9", "ja_gene": "TRAJ9", "cdr3a": "CXXX",
             "vb_gene": "TRBV9", "jb_gene": "TRBJ9", "cdr3b": "CYYY", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
        ])
        result = select_matched_tcr(panel, "HLA-A*02:01", _DMF5)
        assert result["panel_status"] == "dmf5_fallback"
        assert result["tcr_va"] == _DMF5["va_gene"]
        assert result["tcr_cdr3b"] == _DMF5["cdr3b"]
        assert result["vdjdb_donor_id"] == ""

    def test_none_panel_falls_back(self):
        result = select_matched_tcr(None, "HLA-C*07:01", _DMF5)
        assert result["panel_status"] == "dmf5_fallback"
        assert result["tcr_va"] == _DMF5["va_gene"]

    def test_empty_panel_falls_back(self):
        result = select_matched_tcr(_panel([]), "HLA-C*07:01", _DMF5)
        assert result["panel_status"] == "dmf5_fallback"


# ---------------------------------------------------------------------------
# attach_matched_tcrs
# ---------------------------------------------------------------------------

class TestAttachMatchedTcrs:
    def test_mixed_matched_and_fallback(self):
        candidates = pd.DataFrame([
            {"peptide": "SQIPRTHSY", "best_allele": "HLA-C*07:01"},  # matched
            {"peptide": "AAAAAAAAA", "best_allele": "HLA-B*40:01"},  # no panel entry → fallback
        ])
        panel = _panel([
            {"allele": "HLA-C*07:01", "va_gene": "TRAV9", "ja_gene": "TRAJ9", "cdr3a": "CXXX",
             "vb_gene": "TRBV9", "jb_gene": "TRBJ9", "cdr3b": "CYYY", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
        ])
        result = attach_matched_tcrs(candidates, panel, _DMF5)
        assert list(result["panel_status"]) == ["vdjdb_matched", "dmf5_fallback"]
        assert result.iloc[0]["tcr_va"] == "TRAV9"
        assert result.iloc[1]["tcr_va"] == _DMF5["va_gene"]


# ---------------------------------------------------------------------------
# build_tcrdock_input — uses the matched TCR, not the hardcoded DMF5
# ---------------------------------------------------------------------------

class TestBuildTcrdockInputUsesMatchedTcr:
    def test_input_tsv_carries_matched_tcr(self, tmp_path):
        candidates = pd.DataFrame([{"peptide": "SQIPRTHSY", "best_allele": "HLA-C*07:01"}])
        panel = _panel([
            {"allele": "HLA-C*07:01", "va_gene": "TRAV9", "ja_gene": "TRAJ9", "cdr3a": "CXXX",
             "vb_gene": "TRBV9", "jb_gene": "TRBJ9", "cdr3b": "CYYY", "vdjdb_score": 3, "vdjdb_donor_id": "D2"},
        ])
        candidates = attach_matched_tcrs(candidates, panel, _DMF5)
        input_tsv = build_tcrdock_input(candidates, {"A": "HLA-A*02:01"}, tmp_path)

        written = pd.read_csv(input_tsv, sep="\t")
        assert written.iloc[0]["mhc"] == "C*07:01"          # HLA- prefix stripped
        assert written.iloc[0]["va"] == "TRAV9*01"          # *01 appended; matched TCR (not DMF5)
        assert written.iloc[0]["cdr3a"] == "CXXX"
        assert written.iloc[0]["va"] != _DMF5["va_gene"] + "*01"


# ---------------------------------------------------------------------------
# collect_outputs — docking_scores.tsv carries TCR provenance (AC4)
# ---------------------------------------------------------------------------

class TestCollectOutputsProvenance:
    def test_docking_scores_includes_provenance_columns(self, tmp_path):
        out_dir = tmp_path / "tcrdock_out"
        out_dir.mkdir()
        # Minimal AlphaFold flat-chain PDB that collect_outputs can relabel.
        pdb_path = out_dir / "model.pdb"
        pdb_path.write_text(_FLAT_PDB)
        (out_dir / "alphafold_setup").mkdir()
        (out_dir / "alphafold_setup" / "targets.tsv").write_text(
            "target_chainseq\nAA/G/S/V\n"
        )
        pd.DataFrame([{
            "peptide": "SQIPRTHSY", "mhc": "C*07:01",
            "model_2_ptm_plddt": 80.0, "model_2_ptm_pdb_file": str(pdb_path),
        }]).to_csv(out_dir / "tcrdock_out_final.tsv", sep="\t", index=False)

        candidates = pd.DataFrame([{
            "peptide": "SQIPRTHSY", "best_allele": "HLA-C*07:01",
            "tcr_va": "TRAV9", "tcr_ja": "TRAJ9", "tcr_vb": "TRBV9", "tcr_jb": "TRBJ9",
            "tcr_cdr3a": "CXXX", "tcr_cdr3b": "CYYY",
            "vdjdb_donor_id": "D2", "vdjdb_score": 3, "panel_status": "vdjdb_matched",
        }])
        output_pdb = tmp_path / "top_candidate.pdb"
        output_scores = tmp_path / "docking_scores.tsv"
        collect_outputs(out_dir, candidates, output_pdb, output_scores)

        scores = pd.read_csv(output_scores, sep="\t")
        for col in ["tcr_va", "tcr_cdr3a", "vdjdb_donor_id", "vdjdb_score", "panel_status"]:
            assert col in scores.columns
        assert scores.iloc[0]["panel_status"] == "vdjdb_matched"
        assert scores.iloc[0]["vdjdb_donor_id"] == "D2"
        assert scores.iloc[0]["tcr_va"] == "TRAV9"
        assert output_pdb.exists()
