"""Tests for run_mhcflurry.py — TSV I/O and binder classification."""

import pandas as pd
import pytest

from run_mhcflurry import classify, run_prediction


class TestRunPredictionEmpty:
    """run_prediction should write an empty TSV when input has no rows."""

    def test_empty_input_writes_empty_output(self, tmp_path):
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")

        output_tsv = tmp_path / "predictions.tsv"
        # Skip actual MHCflurry model loading — empty input short-circuits before it
        run_prediction(
            peptides_tsv=peptides_tsv,
            output_tsv=output_tsv,
            alleles=["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"],
            ic50_strong=50.0,
            ic50_weak=500.0,
        )

        df = pd.read_csv(output_tsv, sep="\t")
        assert df.empty
        assert list(df.columns) == [
            "contig_key", "start_nt", "peptide", "allele",
            "ic50_nM", "percentile_rank", "binder_class",
        ]

    def test_output_file_is_created(self, tmp_path):
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")
        output_tsv = tmp_path / "subdir" / "predictions.tsv"

        run_prediction(
            peptides_tsv=peptides_tsv,
            output_tsv=output_tsv,
            alleles=["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"],
        )
        assert output_tsv.exists()


class TestBinderClassification:
    """Classification thresholds: strong <= 50, weak <= 500, non otherwise."""

    @pytest.mark.parametrize("ic50,expected", [
        (10.0, "strong"),
        (49.9, "strong"),
        (50.0, "strong"),   # boundary: ic50 <= ic50_strong → strong
        (50.1, "weak"),
        (499.9, "weak"),
        (500.0, "weak"),    # boundary: ic50 <= ic50_weak → weak
        (500.1, "non"),
        (9999.0, "non"),
    ])
    def test_classify_boundaries(self, ic50, expected):
        """Test the production classify() function boundary behaviour."""
        assert classify(ic50) == expected
