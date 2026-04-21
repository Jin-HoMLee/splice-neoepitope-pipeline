"""Tests for aggregate_hla_alleles.py — OptiType TSV parsing and aggregation."""

import csv
from pathlib import Path

import pytest

from aggregate_hla_alleles import (
    aggregate,
    load_optitype_result,
    load_serology,
    normalise_allele,
    _validate_vs_serology,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _write_optitype_tsv(path: Path, a1="A*02:01", a2="A*02:01",
                         b1="B*07:02", b2="B*07:02",
                         c1="C*07:02", c2="C*07:02",
                         reads=500) -> Path:
    """Write a minimal OptiType result TSV."""
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["", "A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"])
        writer.writerow([0, a1, a2, b1, b2, c1, c2, reads, 0.987654])
    return path


_SEROLOGY_FIELDS = [
    "serology_A1", "serology_A2",
    "serology_B1", "serology_B2",
    "serology_C1", "serology_C2",
]


def _write_samples_tsv(path: Path, rows: list[dict]) -> Path:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["patient_id", "sample_id", "sample_type", "fastq1", "fastq2"]
                       + _SEROLOGY_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


FALLBACK = {"A": "HLA-A*02:01", "B": "HLA-B*07:02", "C": "HLA-C*07:02"}


# ---------------------------------------------------------------------------
# normalise_allele
# ---------------------------------------------------------------------------

class TestNormaliseAllele:
    @pytest.mark.parametrize("raw,expected", [
        ("A*02:01",         "HLA-A*02:01"),
        ("HLA-A*02:01",     "HLA-A*02:01"),
        ("A*02:01:01:01",   "HLA-A*02:01"),   # truncate to 4-digit
        ("B*07:02",         "HLA-B*07:02"),
        ("C*07:02",         "HLA-C*07:02"),
        ("",                ""),
    ])
    def test_normalise(self, raw, expected):
        assert normalise_allele(raw) == expected


# ---------------------------------------------------------------------------
# load_optitype_result
# ---------------------------------------------------------------------------

class TestLoadOptitypeResult:
    def test_valid_tsv(self, tmp_path):
        tsv = _write_optitype_tsv(
            tmp_path / "sample_result.tsv",
            a1="A*02:01", a2="A*24:02",
            b1="B*07:02", b2="B*35:01",
            c1="C*07:01", c2="C*07:02",
            reads=1234,
        )
        result = load_optitype_result(str(tsv))
        assert result["A"] == (["HLA-A*02:01", "HLA-A*24:02"], 1234)
        assert result["B"] == (["HLA-B*07:02", "HLA-B*35:01"], 1234)
        assert result["C"] == (["HLA-C*07:01", "HLA-C*07:02"], 1234)

    def test_missing_file_returns_empty(self, tmp_path):
        assert load_optitype_result(str(tmp_path / "nonexistent.tsv")) == {}

    def test_header_only_returns_empty(self, tmp_path):
        tsv = tmp_path / "sample_result.tsv"
        tsv.write_text("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n")
        assert load_optitype_result(str(tsv)) == {}

    def test_homozygous_alleles(self, tmp_path):
        tsv = _write_optitype_tsv(
            tmp_path / "sample_result.tsv",
            a1="A*02:01", a2="A*02:01", reads=200,
        )
        result = load_optitype_result(str(tsv))
        assert result["A"] == (["HLA-A*02:01", "HLA-A*02:01"], 200)


# ---------------------------------------------------------------------------
# aggregate — priority order and serology
# ---------------------------------------------------------------------------

class TestAggregate:
    def test_tumor_first_policy(self, tmp_path):
        """Tumor sample call is preferred over normal."""
        normal_tsv = _write_optitype_tsv(
            tmp_path / "normal_result.tsv",
            a1="A*02:01", a2="A*02:01", reads=9999,
        )
        tumor_tsv = _write_optitype_tsv(
            tmp_path / "tumor_result.tsv",
            a1="A*03:01", a2="A*03:01", reads=100,
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Solid Tissue Normal", "fastq1": "x", "fastq2": ""},
            {"patient_id": "p1", "sample_id": "tumor",  "sample_type": "Primary Tumor",       "fastq1": "x", "fastq2": ""},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(normal_tsv), str(tumor_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        rows = list(csv.DictReader(out_alleles.open(), delimiter="\t"))
        a_row = next(r for r in rows if r["locus"] == "A")
        assert a_row["allele1"] == "HLA-A*03:01"

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "tumor"

    def test_serology_over_normal(self, tmp_path):
        """Serology is preferred over normal OptiType when tumor has no confident call."""
        normal_tsv = _write_optitype_tsv(
            tmp_path / "normal_result.tsv",
            a1="A*02:01", a2="A*02:01", reads=200,
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Solid Tissue Normal",
             "fastq1": "x", "fastq2": "", "serology_A1": "HLA-A*01:01", "serology_A2": "HLA-A*01:01"},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(normal_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        rows = list(csv.DictReader(out_alleles.open(), delimiter="\t"))
        a_row = next(r for r in rows if r["locus"] == "A")
        assert a_row["allele1"] == "HLA-A*01:01"

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "serology"

    def test_serology_not_used_when_tumor_available(self, tmp_path):
        """Tumor OptiType wins over serology even when serology is provided."""
        tumor_tsv = _write_optitype_tsv(
            tmp_path / "tumor_result.tsv",
            a1="A*03:01", a2="A*03:01", reads=200,
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "tumor", "sample_type": "Primary Tumor",
             "fastq1": "x", "fastq2": "", "serology_A1": "HLA-A*01:01", "serology_A2": "HLA-A*01:01"},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(tumor_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "tumor"
        assert a_qc["allele1"] == "HLA-A*03:01"

    def test_null_allele_excluded_from_prediction(self, tmp_path):
        """Null alleles in serology are not written to alleles.tsv."""
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Blood Derived Normal",
             "fastq1": "x", "fastq2": "",
             "serology_A1": "HLA-A*01:01", "serology_A2": "HLA-A*01:11N"},
        ])
        # Empty OptiType result so serology is used
        empty_tsv = tmp_path / "normal_result.tsv"
        empty_tsv.write_text("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n")

        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(empty_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        rows = list(csv.DictReader(out_alleles.open(), delimiter="\t"))
        a_row = next(r for r in rows if r["locus"] == "A")
        # Null allele excluded; both prediction slots use the expressed allele
        assert a_row["allele1"] == "HLA-A*01:01"
        assert a_row["allele2"] == "HLA-A*01:01"

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "serology"
        # Null allele still visible in QC
        assert a_qc["serology_allele2"] == "HLA-A*01:11N"

    def test_fallback_when_no_reads(self, tmp_path):
        """Fallback alleles used when all samples have empty result TSVs."""
        for name in ("normal_result.tsv", "tumor_result.tsv"):
            (tmp_path / name).write_text(
                "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
            )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Solid Tissue Normal", "fastq1": "x", "fastq2": ""},
            {"patient_id": "p1", "sample_id": "tumor",  "sample_type": "Primary Tumor",       "fastq1": "x", "fastq2": ""},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[
                str(tmp_path / "normal_result.tsv"),
                str(tmp_path / "tumor_result.tsv"),
            ],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        for row in qc_rows:
            assert row["source"] == "fallback"
            assert row["reads"] == "0"

    def test_low_reads_filtered_uses_fallback(self, tmp_path):
        """Samples below min_reads threshold are skipped → fallback used."""
        normal_tsv = _write_optitype_tsv(
            tmp_path / "normal_result.tsv",
            a1="A*03:01", a2="A*03:01", reads=5,   # below min_reads=30
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Solid Tissue Normal", "fastq1": "x", "fastq2": ""},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(normal_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "fallback"
        assert a_qc["allele1"] == "HLA-A*02:01"

    def test_tumor_primary_source(self, tmp_path):
        """Tumor call used when it is the only sample."""
        tumor_tsv = _write_optitype_tsv(
            tmp_path / "tumor_result.tsv",
            a1="A*11:01", a2="A*11:01", reads=200,
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "tumor", "sample_type": "Primary Tumor", "fastq1": "x", "fastq2": ""},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(tumor_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "tumor"
        assert a_qc["allele1"] == "HLA-A*11:01"

    def test_discrepancy_flagged(self, tmp_path):
        """Normal/tumor disagreement is recorded in hla_qc.tsv."""
        normal_tsv = _write_optitype_tsv(
            tmp_path / "normal_result.tsv",
            a1="A*02:01", a2="A*02:01", reads=200,
        )
        tumor_tsv = _write_optitype_tsv(
            tmp_path / "tumor_result.tsv",
            a1="A*03:01", a2="A*03:01", reads=200,
        )
        samples_tsv = _write_samples_tsv(tmp_path / "samples.tsv", [
            {"patient_id": "p1", "sample_id": "normal", "sample_type": "Solid Tissue Normal", "fastq1": "x", "fastq2": ""},
            {"patient_id": "p1", "sample_id": "tumor",  "sample_type": "Primary Tumor",       "fastq1": "x", "fastq2": ""},
        ])
        out_alleles = tmp_path / "alleles.tsv"
        out_qc      = tmp_path / "hla_qc.tsv"

        aggregate(
            result_tsvs=[str(normal_tsv), str(tumor_tsv)],
            samples_tsv=str(samples_tsv),
            min_reads=30,
            fallback=FALLBACK,
            out_alleles=str(out_alleles),
            out_qc=str(out_qc),
        )

        qc_rows = list(csv.DictReader(out_qc.open(), delimiter="\t"))
        a_qc = next(r for r in qc_rows if r["locus"] == "A")
        assert a_qc["source"] == "tumor"
        assert a_qc["discrepancy"] != ""
        assert "A*02:01" in a_qc["discrepancy"]
        assert "A*03:01" in a_qc["discrepancy"]


# ---------------------------------------------------------------------------
# _validate_vs_serology
# ---------------------------------------------------------------------------

class TestValidateVsSerology:
    def test_exact_match(self):
        result = _validate_vs_serology(
            ["HLA-A*01:01", "HLA-A*01:01"],
            ("HLA-A*01:01", "HLA-A*01:01"),
        )
        assert result == "match"

    def test_match_with_null_allele_excluded(self):
        result = _validate_vs_serology(
            ["HLA-A*01:01", "HLA-A*01:01"],
            ("HLA-A*01:01", "HLA-A*01:11N"),
        )
        assert result == "match (null allele excluded)"

    def test_mismatch(self):
        result = _validate_vs_serology(
            ["HLA-A*02:01", "HLA-A*02:01"],
            ("HLA-A*01:01", "HLA-A*01:01"),
        )
        assert result.startswith("mismatch")

    def test_no_optitype_call_returns_empty(self):
        assert _validate_vs_serology(None, ("HLA-A*01:01", "HLA-A*01:01")) == ""
        assert _validate_vs_serology([], ("HLA-A*01:01", "HLA-A*01:01")) == ""
