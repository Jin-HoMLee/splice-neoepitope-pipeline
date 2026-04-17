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

    def test_empty_alleles_tsv_raises(self, tmp_path):
        """alleles_tsv with no valid allele values must raise ValueError."""
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")

        alleles_tsv = tmp_path / "alleles.tsv"
        with alleles_tsv.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["locus", "allele1", "allele2"], delimiter="\t"
            )
            writer.writeheader()
            # Blank allele values — simulates misconfigured fallback
            writer.writerow({"locus": "A", "allele1": "", "allele2": ""})

        with pytest.raises(ValueError, match="no valid alleles"):
            run_prediction(
                peptides_tsv=peptides_tsv,
                output_tsv=tmp_path / "out.tsv",
                alleles_tsv=str(alleles_tsv),
            )


# ---------------------------------------------------------------------------
# _worker_init and _predict_allele_worker — isolated unit tests
# ---------------------------------------------------------------------------

class TestWorker:
    """Unit tests for the per-process worker functions.

    Tests run entirely in-process: _worker_predictor is set manually and
    _run_mhcflurry_predictions is monkeypatched, so no MHCflurry models or
    subprocess spawning are needed.
    """

    PEPTIDES = ["ACDEFGHIK", "LMNPQRSTV", "YKLMFSTAV"]

    @pytest.fixture(autouse=True)
    def _patch_mhcflurry(self, monkeypatch):
        import run_mhcflurry

        def fake_predict(peptides, allele, predictor=None):
            ic50 = {"HLA-A*02:01": 30.0, "HLA-B*07:02": 300.0}.get(allele, 9999.0)
            return pd.DataFrame({
                "peptide": peptides,
                "prediction": [ic50] * len(peptides),
                "prediction_percentile": [5.0] * len(peptides),
            })

        monkeypatch.setattr(run_mhcflurry, "_load_mhcflurry_predictor", lambda: object())
        monkeypatch.setattr(run_mhcflurry, "_run_mhcflurry_predictions", fake_predict)
        monkeypatch.setattr(run_mhcflurry, "_worker_predictor", object())

    def test_worker_returns_correct_columns(self):
        from run_mhcflurry import _predict_allele_worker
        df = _predict_allele_worker("HLA-A*02:01", self.PEPTIDES, 50.0, 500.0)
        assert list(df.columns) == ["peptide", "ic50_nM", "percentile_rank", "allele", "binder_class"]

    def test_worker_allele_assigned(self):
        from run_mhcflurry import _predict_allele_worker
        df = _predict_allele_worker("HLA-B*07:02", self.PEPTIDES, 50.0, 500.0)
        assert (df["allele"] == "HLA-B*07:02").all()

    def test_worker_binder_class_strong(self):
        """IC50=30 nM → strong binder (threshold 50 nM)."""
        from run_mhcflurry import _predict_allele_worker
        df = _predict_allele_worker("HLA-A*02:01", self.PEPTIDES, 50.0, 500.0)
        assert (df["binder_class"] == "strong").all()

    def test_worker_binder_class_weak(self):
        """IC50=300 nM → weak binder (threshold 500 nM)."""
        from run_mhcflurry import _predict_allele_worker
        df = _predict_allele_worker("HLA-B*07:02", self.PEPTIDES, 50.0, 500.0)
        assert (df["binder_class"] == "weak").all()

    def test_worker_binder_class_non(self):
        """Unknown allele gets IC50=9999 nM → non-binder."""
        from run_mhcflurry import _predict_allele_worker
        df = _predict_allele_worker("HLA-C*07:02", self.PEPTIDES, 50.0, 500.0)
        assert (df["binder_class"] == "non").all()

    def test_worker_init_sets_env_vars(self, monkeypatch):
        """_worker_init must set TF/BLAS thread-count env vars before model load."""
        import os
        import run_mhcflurry
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS"):
            monkeypatch.delenv(var, raising=False)

        run_mhcflurry._worker_init()

        assert os.environ["OMP_NUM_THREADS"] == "1"
        assert os.environ["TF_NUM_INTRAOP_THREADS"] == "1"
        assert os.environ["TF_NUM_INTEROP_THREADS"] == "1"

    def test_invalid_threads_raises(self, tmp_path):
        """threads=0 must raise ValueError before any workers are spawned."""
        peptides_tsv = tmp_path / "peptides.tsv"
        peptides_tsv.write_text("contig_key\tstart_nt\tpeptide\n")
        with pytest.raises(ValueError, match="threads"):
            run_prediction(
                peptides_tsv=peptides_tsv,
                output_tsv=tmp_path / "out.tsv",
                alleles=["HLA-A*02:01"],
                threads=0,
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
