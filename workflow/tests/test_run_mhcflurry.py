"""Tests for run_mhcflurry.py — TSV I/O, allele loading, and binder classification."""

import csv

import pandas as pd
import pytest

from run_mhcflurry import _load_alleles_from_tsv, classify_by_percentile, run_prediction


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
        )

        df = pd.read_csv(output_tsv, sep="\t")
        assert df.empty
        assert list(df.columns) == [
            "contig_key", "start_nt", "peptide", "best_allele",
            "ic50_nM", "processing_score", "presentation_score",
            "presentation_percentile", "presentation_class",
            "HLA-A*02:01_presentation_score", "HLA-A*02:01_presentation_percentile",
            "genotype_presentation_score", "n_strong_alleles", "best_presentation_percentile",
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
# run_prediction with non-empty input — isolated unit tests
# ---------------------------------------------------------------------------

class TestRunPredictionNonEmpty:
    """Unit tests for run_prediction with real peptide rows.

    _run_mhcflurry_predictions is monkeypatched so no MHCflurry models are
    loaded. Tests verify the full run_prediction orchestration: column renaming,
    classification, contig_key/start_nt merge, and one-row-per-peptide output.
    """

    PEPTIDES = ["ACDEFGHIK", "LMNPQRSTV", "YKLMFSTAV"]

    @pytest.fixture(autouse=True)
    def _patch(self, monkeypatch):
        import run_mhcflurry

        def fake_run_predictions(peptides, alleles, predictor=None):
            # Percentile/IC50 keyed on the first allele in the genotype list
            first = alleles[0] if alleles else "HLA-A*02:01"
            percentile = {"HLA-A*02:01": 0.3, "HLA-B*07:02": 1.5}.get(first, 5.0)
            ic50 = {"HLA-A*02:01": 30.0, "HLA-B*07:02": 300.0}.get(first, 9999.0)
            n = len(peptides)
            return pd.DataFrame({
                "peptide": peptides,
                "affinity": [ic50] * n,
                "best_allele": [first] * n,
                "processing_score": [0.8] * n,
                "presentation_score": [0.9] * n,
                "presentation_percentile": [percentile] * n,
            })

        monkeypatch.setattr(run_mhcflurry, "_load_mhcflurry_predictor", lambda: object())
        monkeypatch.setattr(run_mhcflurry, "_run_mhcflurry_predictions", fake_run_predictions)
        monkeypatch.setattr(run_mhcflurry, "_worker_predictor", object())

    def _make_peptides_tsv(self, tmp_path):
        tsv = tmp_path / "peptides.tsv"
        with tsv.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["contig_key", "start_nt", "peptide"], delimiter="\t"
            )
            writer.writeheader()
            for i, p in enumerate(self.PEPTIDES):
                writer.writerow({"contig_key": f"k{i}", "start_nt": i * 3, "peptide": p})
        return tsv

    def test_output_has_expected_columns(self, tmp_path):
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01"],
        )
        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == [
            "contig_key", "start_nt", "peptide", "best_allele",
            "ic50_nM", "processing_score", "presentation_score",
            "presentation_percentile", "presentation_class",
            "HLA-A*02:01_presentation_score", "HLA-A*02:01_presentation_percentile",
            "genotype_presentation_score", "n_strong_alleles", "best_presentation_percentile",
        ]

    def test_presentation_class_strong(self, tmp_path):
        """Percentile 0.3 (HLA-A*02:01) → strong."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01"],
        )
        df = pd.read_csv(out, sep="\t")
        assert (df["presentation_class"] == "strong").all()

    def test_presentation_class_weak(self, tmp_path):
        """Percentile 1.5 (HLA-B*07:02) → weak."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-B*07:02"],
        )
        df = pd.read_csv(out, sep="\t")
        assert (df["presentation_class"] == "weak").all()

    def test_presentation_class_non(self, tmp_path):
        """Unknown allele → percentile 5.0 → non."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-C*07:02"],
        )
        df = pd.read_csv(out, sep="\t")
        assert (df["presentation_class"] == "non").all()

    def test_best_allele_column_present(self, tmp_path):
        """best_allele from predictor is retained as the best_allele column."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01"],
        )
        df = pd.read_csv(out, sep="\t")
        assert "best_allele" in df.columns
        assert (df["best_allele"] == "HLA-A*02:01").all()

    def test_one_row_per_peptide(self, tmp_path):
        """Genotype prediction returns one row per peptide, not one per peptide×allele."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01", "HLA-B*07:02"],
        )
        df = pd.read_csv(out, sep="\t")
        assert len(df) == len(self.PEPTIDES)

    def test_per_allele_columns_present(self, tmp_path):
        """Output has one score+percentile column pair per allele."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01", "HLA-B*07:02"],
        )
        df = pd.read_csv(out, sep="\t")
        for allele in ("HLA-A*02:01", "HLA-B*07:02"):
            assert f"{allele}_presentation_score" in df.columns
            assert f"{allele}_presentation_percentile" in df.columns

    def test_genotype_presentation_score_computed(self, tmp_path):
        """genotype_presentation_score = 1 − ∏(1 − wᵢ·pᵢ); fake score=0.9, w=1.0 → GPS=0.9."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-A*02:01"],
        )
        df = pd.read_csv(out, sep="\t")
        # fake presentation_score = 0.9, w(A) = 1.0 → 1 - (1 - 0.9) = 0.9
        n = len(self.PEPTIDES)
        assert df["genotype_presentation_score"].to_list() == pytest.approx([0.9] * n, abs=1e-5)
        assert (df["n_strong_alleles"] == 1).all()   # percentile 0.3 <= 0.5%
        assert df["best_presentation_percentile"].to_list() == pytest.approx([0.3] * n, abs=1e-5)

    def test_hla_c_weight_applied(self, tmp_path):
        """HLA-C alleles use hla_c_weight (default 0.5), not 1.0."""
        out = tmp_path / "predictions.tsv"
        run_prediction(
            peptides_tsv=self._make_peptides_tsv(tmp_path),
            output_tsv=out,
            alleles=["HLA-C*07:02"],
            hla_c_weight=0.5,
        )
        df = pd.read_csv(out, sep="\t")
        # fake presentation_score = 0.9, w(C) = 0.5 → 1 - (1 - 0.5*0.9) = 1 - 0.55 = 0.45
        n = len(self.PEPTIDES)
        assert df["genotype_presentation_score"].to_list() == pytest.approx([0.45] * n, abs=1e-5)

    def test_load_predictor_populates_cache(self, monkeypatch):
        """_load_predictor must populate the module-level predictor cache."""
        import run_mhcflurry
        monkeypatch.setattr(run_mhcflurry, "_worker_predictor", None)
        run_mhcflurry._load_predictor()
        assert run_mhcflurry._worker_predictor is not None


# ---------------------------------------------------------------------------
# Percentile classification
# ---------------------------------------------------------------------------

class TestPercentileClassification:
    """Boundary behaviour: strong <= 0.5%, weak <= 2.0%, non otherwise."""

    @pytest.mark.parametrize("percentile,expected", [
        (0.1, "strong"),
        (0.4, "strong"),
        (0.5, "strong"),   # boundary: percentile <= 0.5 → strong
        (0.6, "weak"),
        (1.9, "weak"),
        (2.0, "weak"),     # boundary: percentile <= 2.0 → weak
        (2.1, "non"),
        (50.0, "non"),
    ])
    def test_classify_by_percentile_boundaries(self, percentile, expected):
        assert classify_by_percentile(percentile) == expected
