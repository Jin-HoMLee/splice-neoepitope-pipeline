"""Tests for run_mhcflurry.py — TSV I/O, allele loading, and binder classification."""

import csv

import pandas as pd
import pytest

from run_mhcflurry import _load_alleles_from_tsv, classify, run_prediction


# ---------------------------------------------------------------------------
# _load_alleles_from_tsv
# ---------------------------------------------------------------------------

class TestLoadAllelesFromTsv:
    def _write_alleles_tsv(self, path, rows):
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["locus", "allele1", "allele2"], delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_heterozygous_patient(self, tmp_path):
        tsv = self._write_alleles_tsv(tmp_path / "alleles.tsv", [
            {"locus": "A", "allele1": "HLA-A*02:01", "allele2": "HLA-A*24:02"},
            {"locus": "B", "allele1": "HLA-B*07:02", "allele2": "HLA-B*35:01"},
            {"locus": "C", "allele1": "HLA-C*07:01", "allele2": "HLA-C*07:02"},
        ])
        alleles = _load_alleles_from_tsv(str(tsv))
        assert alleles == [
            "HLA-A*02:01", "HLA-A*24:02",
            "HLA-B*07:02", "HLA-B*35:01",
            "HLA-C*07:01", "HLA-C*07:02",
        ]

    def test_homozygous_deduplication(self, tmp_path):
        """Homozygous loci (allele1 == allele2) must appear only once."""
        tsv = self._write_alleles_tsv(tmp_path / "alleles.tsv", [
            {"locus": "A", "allele1": "HLA-A*02:01", "allele2": "HLA-A*02:01"},
            {"locus": "B", "allele1": "HLA-B*07:02", "allele2": "HLA-B*35:01"},
            {"locus": "C", "allele1": "HLA-C*07:02", "allele2": "HLA-C*07:02"},
        ])
        alleles = _load_alleles_from_tsv(str(tsv))
        assert alleles == ["HLA-A*02:01", "HLA-B*07:02", "HLA-B*35:01", "HLA-C*07:02"]

    def test_preserves_order(self, tmp_path):
        """Alleles come out in locus row order, allele1 before allele2."""
        tsv = self._write_alleles_tsv(tmp_path / "alleles.tsv", [
            {"locus": "A", "allele1": "HLA-A*03:01", "allele2": "HLA-A*02:01"},
        ])
        alleles = _load_alleles_from_tsv(str(tsv))
        assert alleles[0] == "HLA-A*03:01"
        assert alleles[1] == "HLA-A*02:01"


# ---------------------------------------------------------------------------
# run_prediction — empty input / file creation
# ---------------------------------------------------------------------------

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
            alleles=["HLA-A*02:01"],
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
            alleles=["HLA-A*02:01"],
        )
        assert output_tsv.exists()

    def test_empty_input_with_alleles_tsv(self, tmp_path):
        """alleles_tsv path is honoured even when peptides TSV is empty."""
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")

        alleles_tsv = tmp_path / "alleles.tsv"
        with alleles_tsv.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["locus", "allele1", "allele2"], delimiter="\t"
            )
            writer.writeheader()
            writer.writerow({"locus": "A", "allele1": "HLA-A*02:01", "allele2": "HLA-A*02:01"})

        output_tsv = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=peptides_tsv,
            output_tsv=output_tsv,
            alleles_tsv=str(alleles_tsv),
        )
        assert output_tsv.exists()
        df = pd.read_csv(output_tsv, sep="\t")
        assert df.empty

    def test_no_alleles_raises(self, tmp_path):
        """run_prediction must raise ValueError when neither alleles nor alleles_tsv given."""
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")
        with pytest.raises(ValueError, match="alleles"):
            run_prediction(
                peptides_tsv=peptides_tsv,
                output_tsv=tmp_path / "out.tsv",
            )


# ---------------------------------------------------------------------------
# Binder classification
# ---------------------------------------------------------------------------

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
