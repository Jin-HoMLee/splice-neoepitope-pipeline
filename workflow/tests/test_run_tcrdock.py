"""Tests for run_tcrdock.py — candidate selection and PDB chain relabeling."""

import pandas as pd

from run_tcrdock import relabel_pdb_chains, select_top_candidates


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
