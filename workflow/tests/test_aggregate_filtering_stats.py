"""Tests for aggregate_filtering_stats — Issue #215 cross-step funnel."""

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aggregate_filtering_stats import aggregate, UNIFIED_COLUMNS


def _write_per_sample(path: Path, rows: list[dict]) -> None:
    """Write a per-sample stats TSV (junction_filter format)."""
    pd.DataFrame(rows, columns=["sample_id", "sample_type", "category", "count"]).to_csv(
        path, sep="\t", index=False,
    )


def _write_per_patient(path: Path, rows: list[tuple[str, float]]) -> None:
    """Write a per-patient stats TSV (2-column category, count)."""
    pd.DataFrame(rows, columns=["category", "count"]).to_csv(path, sep="\t", index=False)


@pytest.fixture
def stats_files(tmp_path):
    junction_filter = tmp_path / "junction_filter_stats.tsv"
    _write_per_sample(junction_filter, [
        {"sample_id": "tumor", "sample_type": "Primary Tumor", "category": "junctions_raw", "count": 100},
        {"sample_id": "tumor", "sample_type": "Primary Tumor", "category": "tumor_exclusive", "count": 30},
        {"sample_id": "tumor", "sample_type": "Primary Tumor", "category": "mean_reads", "count": 12.5},
    ])
    contig_assemble = tmp_path / "contig_assemble_stats.tsv"
    _write_per_patient(contig_assemble, [
        ("contigs_written", 25),
        ("skipped_softclip", 3),
        ("skipped_length", 2),
    ])
    translate = tmp_path / "translate_stats.tsv"
    _write_per_patient(translate, [("peptides_total", 1500)])
    proteome = tmp_path / "proteome_stats.tsv"
    _write_per_patient(proteome, [("novel", 1200), ("excluded", 300)])
    mhc = tmp_path / "mhc_stats.tsv"
    _write_per_patient(mhc, [("strong_presenters", 18), ("weak_presenters", 42)])
    return {
        "junction_filter": junction_filter,
        "contig_assemble": contig_assemble,
        "translate": translate,
        "proteome": proteome,
        "mhc": mhc,
    }


class TestAggregate:
    def test_unified_schema(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            proteome_tsv=stats_files["proteome"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == UNIFIED_COLUMNS
        assert (df["patient_id"] == "patient_001").all()

    def test_per_sample_step_preserves_sample_columns(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            proteome_tsv=stats_files["proteome"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t")
        junction_rows = df[df["step"] == "junction-filter"]
        assert (junction_rows["sample_id"] == "tumor").all()
        assert (junction_rows["sample_type"] == "Primary Tumor").all()

    def test_per_patient_steps_have_blank_sample_columns(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            proteome_tsv=stats_files["proteome"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t").fillna("")
        for step in ["contig-assemble", "peptide-translate", "proteome-filter", "mhc-affinity"]:
            step_rows = df[df["step"] == step]
            assert (step_rows["sample_id"] == "").all(), f"{step} should have empty sample_id"
            assert (step_rows["sample_type"] == "").all(), f"{step} should have empty sample_type"

    def test_includes_mhc_presenter_categories(self, tmp_path, stats_files):
        """Project vocabulary: 'presenter' not 'binder'."""
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            proteome_tsv=stats_files["proteome"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t")
        mhc_rows = df[df["step"] == "mhc-affinity"]
        assert set(mhc_rows["category"]) == {"strong_presenters", "weak_presenters"}

    def test_all_steps_present(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            proteome_tsv=stats_files["proteome"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t")
        assert set(df["step"]) == {
            "junction-filter", "contig-assemble", "peptide-translate",
            "proteome-filter", "mhc-affinity",
        }


class TestAggregateProteomeOptional:
    """Issue #215 follow-up: when ``proteome_filter.enabled: false`` in the
    pipeline config, the proteome step does not run and no stats file is
    produced. ``aggregate()`` must still succeed, simply omitting the
    proteome-filter rows from the audit trail."""

    def test_omits_proteome_when_argument_is_none(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
            proteome_tsv=None,
        )
        df = pd.read_csv(out, sep="\t")
        assert "proteome-filter" not in set(df["step"])
        assert set(df["step"]) == {
            "junction-filter", "contig-assemble", "peptide-translate", "mhc-affinity",
        }

    def test_omits_proteome_when_argument_is_omitted(self, tmp_path, stats_files):
        out = tmp_path / "filtering_stats.tsv"
        # proteome_tsv defaults to None — should behave identically to explicit None
        aggregate(
            patient_id="patient_001",
            junction_filter_tsv=stats_files["junction_filter"],
            contig_assemble_tsv=stats_files["contig_assemble"],
            translate_tsv=stats_files["translate"],
            mhc_tsv=stats_files["mhc"],
            output_tsv=out,
        )
        df = pd.read_csv(out, sep="\t")
        assert "proteome-filter" not in set(df["step"])


class TestAggregateColumnValidation:
    """Per-sample / per-patient stats files with missing required columns
    must raise a clear ``ValueError`` instead of propagating an opaque
    ``KeyError`` from the unified-schema reindex downstream."""

    def test_per_sample_missing_columns_raises(self, tmp_path, stats_files):
        bad = tmp_path / "bad_junction_filter_stats.tsv"
        # Missing sample_id + sample_type
        pd.DataFrame([{"category": "foo", "count": 1}]).to_csv(bad, sep="\t", index=False)
        with pytest.raises(ValueError, match="per-sample stats file missing"):
            aggregate(
                patient_id="patient_001",
                junction_filter_tsv=bad,
                contig_assemble_tsv=stats_files["contig_assemble"],
                translate_tsv=stats_files["translate"],
                proteome_tsv=stats_files["proteome"],
                mhc_tsv=stats_files["mhc"],
                output_tsv=tmp_path / "out.tsv",
            )

    def test_per_patient_missing_columns_raises(self, tmp_path, stats_files):
        bad = tmp_path / "bad_contig_stats.tsv"
        # Missing required category + count
        pd.DataFrame([{"foo": "x", "bar": 2}]).to_csv(bad, sep="\t", index=False)
        with pytest.raises(ValueError, match="per-patient stats file missing"):
            aggregate(
                patient_id="patient_001",
                junction_filter_tsv=stats_files["junction_filter"],
                contig_assemble_tsv=bad,
                translate_tsv=stats_files["translate"],
                proteome_tsv=stats_files["proteome"],
                mhc_tsv=stats_files["mhc"],
                output_tsv=tmp_path / "out.tsv",
            )
