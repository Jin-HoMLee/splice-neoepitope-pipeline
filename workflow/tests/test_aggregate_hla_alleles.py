"""test_aggregate_hla_alleles.py — Unit tests for aggregate_hla_alleles.py.

Tests the normal-first policy, fallback logic, and QC discrepancy detection
without requiring Snakemake or arcasHLA to be installed.
"""

import csv
import json
import sys
from pathlib import Path

import pytest

# Add the scripts directory to the path so we can import the module directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from aggregate_hla_alleles import aggregate, normalise_allele, read_locus_reads


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_data(tmp_path):
    """Create a minimal matched tumor/normal sample set with mock arcasHLA outputs."""

    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "patient_id\tsample_id\tsample_type\tfastq1\n"
        "patient_001\ttumor_01\tPrimary Tumor\tdata/tumor.fq\n"
        "patient_001\tnormal_01\tSolid Tissue Normal\tdata/normal.fq\n"
    )

    def write_geno(sample_dir, alleles: dict):
        """Write a genotype.json for a sample."""
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / f"{sample_dir.name}.genotype.json").write_text(
            json.dumps(alleles)
        )

    def write_log(sample_dir, reads: dict):
        """Write a genotype.log with per-locus read counts."""
        lines = []
        for locus, n in reads.items():
            lines.append(f"[genotype] HLA-{locus}: {n} reads")
        (sample_dir / f"{sample_dir.name}.genotype.log").write_text("\n".join(lines))

    normal_dir = tmp_path / "normal_01"
    tumor_dir  = tmp_path / "tumor_01"

    write_geno(normal_dir, {
        "HLA-A": ["A*02:01:01", "A*24:02:01"],
        "HLA-B": ["B*07:02:01", "B*57:01:01"],
        "HLA-C": ["C*03:04:01", "C*06:02:01"],
    })
    write_log(normal_dir, {"A": 500, "B": 450, "C": 400})

    write_geno(tumor_dir, {
        "HLA-A": ["A*02:01:01", "A*24:02:01"],
        "HLA-B": ["B*07:02:01", "B*57:01:01"],
        "HLA-C": ["C*03:04:01", "C*06:02:01"],
    })
    write_log(tumor_dir, {"A": 300, "B": 280, "C": 250})

    jsons    = [str(normal_dir / "normal_01.genotype.json"),
                str(tumor_dir  / "tumor_01.genotype.json")]
    genologs = [str(normal_dir / "normal_01.genotype.log"),
                str(tumor_dir  / "tumor_01.genotype.log")]

    return {
        "tmp_path":   tmp_path,
        "samples_tsv": str(samples_tsv),
        "jsons":      jsons,
        "genologs":   genologs,
        "normal_dir": normal_dir,
        "tumor_dir":  tumor_dir,
    }


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------------------
# normalise_allele
# ---------------------------------------------------------------------------

def test_normalise_allele_trims_to_four_digit():
    assert normalise_allele("A*02:01:01:01") == "HLA-A*02:01"

def test_normalise_allele_adds_prefix():
    assert normalise_allele("A*02:01") == "HLA-A*02:01"

def test_normalise_allele_already_correct():
    assert normalise_allele("HLA-A*02:01") == "HLA-A*02:01"


# ---------------------------------------------------------------------------
# read_locus_reads
# ---------------------------------------------------------------------------

def test_read_locus_reads_found(tmp_path):
    log_file = tmp_path / "sample.genotype.log"
    log_file.write_text("[genotype] HLA-A: 1234 reads\n[genotype] HLA-B: 567 reads\n")
    assert read_locus_reads(str(log_file), "A") == 1234
    assert read_locus_reads(str(log_file), "B") == 567

def test_read_locus_reads_missing_locus(tmp_path):
    log_file = tmp_path / "sample.genotype.log"
    log_file.write_text("[genotype] HLA-A: 100 reads\n")
    assert read_locus_reads(str(log_file), "C") == 0

def test_read_locus_reads_missing_file(tmp_path):
    assert read_locus_reads(str(tmp_path / "nonexistent.log"), "A") == 0


# ---------------------------------------------------------------------------
# aggregate — normal-first policy
# ---------------------------------------------------------------------------

def test_normal_alleles_preferred(tmp_data):
    """Normal sample alleles should be used when reads meet the threshold."""
    out_alleles = str(tmp_data["tmp_path"] / "alleles.tsv")
    out_qc      = str(tmp_data["tmp_path"] / "qc.tsv")

    aggregate(
        jsons=tmp_data["jsons"],
        genologs=tmp_data["genologs"],
        samples_tsv=tmp_data["samples_tsv"],
        loci=["A", "B", "C"],
        min_reads=30,
        fallback={},
        out_alleles=out_alleles,
        out_qc=out_qc,
    )

    rows = read_tsv(out_alleles)
    assert len(rows) == 3
    a_row = next(r for r in rows if r["locus"] == "A")
    assert a_row["allele1"] == "HLA-A*02:01"
    assert a_row["allele2"] == "HLA-A*24:02"

    qc_rows = read_tsv(out_qc)
    assert all(r["source"] == "normal" for r in qc_rows)


def test_fallback_to_tumor_when_no_normal(tmp_data):
    """When normal reads are below min_reads, fall back to tumor."""
    # Overwrite normal log with very low read counts
    (tmp_data["normal_dir"] / "normal_01.genotype.log").write_text(
        "[genotype] HLA-A: 5 reads\n"
        "[genotype] HLA-B: 5 reads\n"
        "[genotype] HLA-C: 5 reads\n"
    )

    out_alleles = str(tmp_data["tmp_path"] / "alleles.tsv")
    out_qc      = str(tmp_data["tmp_path"] / "qc.tsv")

    aggregate(
        jsons=tmp_data["jsons"],
        genologs=tmp_data["genologs"],
        samples_tsv=tmp_data["samples_tsv"],
        loci=["A"],
        min_reads=30,
        fallback={},
        out_alleles=out_alleles,
        out_qc=out_qc,
    )

    qc_rows = read_tsv(out_qc)
    assert qc_rows[0]["source"] == "tumor"


def test_fallback_allele_used_when_no_calls(tmp_data):
    """When all samples have empty JSON, use the config fallback allele."""
    (tmp_data["normal_dir"] / "normal_01.genotype.json").write_text("{}")
    (tmp_data["tumor_dir"]  / "tumor_01.genotype.json").write_text("{}")

    out_alleles = str(tmp_data["tmp_path"] / "alleles.tsv")
    out_qc      = str(tmp_data["tmp_path"] / "qc.tsv")

    aggregate(
        jsons=tmp_data["jsons"],
        genologs=tmp_data["genologs"],
        samples_tsv=tmp_data["samples_tsv"],
        loci=["A"],
        min_reads=30,
        fallback={"A": "HLA-A*02:01"},
        out_alleles=out_alleles,
        out_qc=out_qc,
    )

    rows = read_tsv(out_alleles)
    assert rows[0]["allele1"] == "HLA-A*02:01"
    qc_rows = read_tsv(out_qc)
    assert qc_rows[0]["source"] == "fallback"


def test_discrepancy_flagged_in_qc(tmp_data):
    """QC TSV should flag when normal and tumor alleles differ."""
    (tmp_data["tumor_dir"] / "tumor_01.genotype.json").write_text(
        json.dumps({"HLA-A": ["A*01:01:01", "A*03:01:01"]})
    )

    out_alleles = str(tmp_data["tmp_path"] / "alleles.tsv")
    out_qc      = str(tmp_data["tmp_path"] / "qc.tsv")

    aggregate(
        jsons=tmp_data["jsons"],
        genologs=tmp_data["genologs"],
        samples_tsv=tmp_data["samples_tsv"],
        loci=["A"],
        min_reads=30,
        fallback={},
        out_alleles=out_alleles,
        out_qc=out_qc,
    )

    qc_rows = read_tsv(out_qc)
    assert qc_rows[0]["discrepancy"] != ""


def test_homozygous_allele_duplicated(tmp_data):
    """A homozygous call (single allele) should fill both allele1 and allele2."""
    (tmp_data["normal_dir"] / "normal_01.genotype.json").write_text(
        json.dumps({"HLA-A": ["A*02:01:01"]})
    )

    out_alleles = str(tmp_data["tmp_path"] / "alleles.tsv")
    out_qc      = str(tmp_data["tmp_path"] / "qc.tsv")

    aggregate(
        jsons=tmp_data["jsons"],
        genologs=tmp_data["genologs"],
        samples_tsv=tmp_data["samples_tsv"],
        loci=["A"],
        min_reads=30,
        fallback={},
        out_alleles=out_alleles,
        out_qc=out_qc,
    )

    rows = read_tsv(out_alleles)
    assert rows[0]["allele1"] == rows[0]["allele2"] == "HLA-A*02:01"
