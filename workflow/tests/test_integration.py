"""test_integration.py — Invariant checks on pipeline output files.

Reads results produced by a previous `snakemake --configfile config/test_config.yaml`
run and asserts structural and biological invariants that must hold regardless
of tool versions, model updates, or reference annotation changes.

Run after a local pipeline run:
    pytest workflow/tests/test_integration.py -v
"""

import csv
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _first_patient_id() -> str:
    tsv = REPO_ROOT / "config" / "samples" / "patient_001_test.tsv"
    with tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            if not pid or pid.startswith("#"):
                continue
            return pid
    raise RuntimeError(f"No patient_id found in {tsv}")


PATIENT_ID = _first_patient_id()
RESULTS = REPO_ROOT / "results" / PATIENT_ID

VALID_JUNCTION_ORIGINS = {"tumor_exclusive", "normal_shared"}
VALID_BINDER_CLASSES = {"strong", "weak", "non"}
VALID_SAMPLE_TYPES = {"Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _results_exist():
    return (RESULTS / "reports" / "report.html").exists()


# Skip the entire module if no pipeline outputs exist yet
pytestmark = pytest.mark.skipif(
    not _results_exist(),
    reason="No pipeline outputs found — run snakemake --configfile config/test_config.yaml first",
)


class TestJunctions:
    def test_novel_junctions_file_exists(self):
        assert (RESULTS / "junctions" / "novel_junctions.tsv").exists()

    def test_has_junctions(self):
        rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
        assert len(rows) > 0

    def test_junction_origins_are_valid(self):
        rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
        origins = {r["junction_origin"] for r in rows}
        assert origins <= VALID_JUNCTION_ORIGINS

    def test_all_junctions_have_required_columns(self):
        rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
        required = {"junction_id", "chrom", "start", "end", "strand",
                    "mapped_reads", "sample_id", "sample_type", "junction_origin"}
        assert required <= set(rows[0].keys())

    def test_coordinates_are_numeric(self):
        rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
        for r in rows:
            assert int(r["start"]) >= 0
            assert int(r["end"]) > int(r["start"])

    def test_normal_shared_filtering_applied(self):
        rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
        origins = {r["junction_origin"] for r in rows}
        # With a matched normal sample, we expect at least some normal_shared junctions
        assert "normal_shared" in origins, (
            "No normal_shared junctions found — normal filtering may not be working"
        )


class TestContigs:
    def test_contigs_file_exists(self):
        assert (RESULTS / "contigs" / "contigs.fa").exists()

    def test_has_contigs(self):
        fa = (RESULTS / "contigs" / "contigs.fa").read_text()
        assert fa.count(">") > 0

    def test_all_contigs_are_50nt(self):
        fa = (RESULTS / "contigs" / "contigs.fa").read_text()
        sequences = [
            line.strip() for line in fa.splitlines() if not line.startswith(">")
        ]
        for seq in sequences:
            assert len(seq) == 50, f"Contig length {len(seq)} != 50: {seq[:20]}..."

    def test_contigs_are_dna(self):
        fa = (RESULTS / "contigs" / "contigs.fa").read_text()
        sequences = [
            line.strip() for line in fa.splitlines() if not line.startswith(">")
        ]
        valid_bases = set("ACGTNacgtn")
        for seq in sequences:
            assert set(seq) <= valid_bases, f"Non-DNA characters in contig: {seq[:20]}..."


class TestPeptides:
    def test_peptides_file_exists(self):
        assert (RESULTS / "peptides" / "peptides.tsv").exists()

    def test_has_peptides(self):
        rows = _read_tsv(RESULTS / "peptides" / "peptides.tsv")
        assert len(rows) > 0

    def test_all_peptides_are_9mers(self):
        rows = _read_tsv(RESULTS / "peptides" / "peptides.tsv")
        for r in rows:
            assert len(r["peptide"]) == 9, f"Peptide length {len(r['peptide'])}: {r['peptide']}"

    def test_peptides_are_amino_acids(self):
        rows = _read_tsv(RESULTS / "peptides" / "peptides.tsv")
        for r in rows:
            assert set(r["peptide"]) <= AMINO_ACIDS, f"Non-AA characters: {r['peptide']}"

    def test_all_peptides_have_required_columns(self):
        rows = _read_tsv(RESULTS / "peptides" / "peptides.tsv")
        required = {"contig_key", "start_nt", "peptide"}
        assert required <= set(rows[0].keys())


class TestPredictions:
    def test_predictions_file_exists(self):
        assert (RESULTS / "predictions" / "mhc_affinity.tsv").exists()

    def test_has_predictions(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        assert len(rows) > 0

    def test_binder_classes_are_valid(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        classes = {r["binder_class"] for r in rows}
        assert classes <= VALID_BINDER_CLASSES

    def test_ic50_values_are_positive(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        for r in rows:
            assert float(r["ic50_nM"]) > 0, f"Non-positive IC50: {r['ic50_nM']}"

    def test_all_predictions_have_required_columns(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        required = {"peptide", "allele", "ic50_nM", "binder_class"}
        assert required <= set(rows[0].keys())

    def test_strong_binders_below_ic50_threshold(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        for r in rows:
            if r["binder_class"] == "strong":
                assert float(r["ic50_nM"]) <= 50, (
                    f"Strong binder above 50 nM threshold: {r['peptide']} {r['ic50_nM']}"
                )

    def test_non_binders_above_weak_threshold(self):
        rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
        for r in rows:
            if r["binder_class"] == "non":
                assert float(r["ic50_nM"]) > 500, (
                    f"Non-binder below 500 nM threshold: {r['peptide']} {r['ic50_nM']}"
                )


class TestReport:
    def test_report_html_exists(self):
        assert (RESULTS / "reports" / "report.html").exists()

    def test_report_is_not_empty(self):
        content = (RESULTS / "reports" / "report.html").read_text()
        assert len(content) > 100

    def test_report_is_valid_html(self):
        content = (RESULTS / "reports" / "report.html").read_text()
        assert "<html" in content.lower()
        assert "</html>" in content.lower()
