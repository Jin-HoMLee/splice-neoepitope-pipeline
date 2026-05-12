"""Tests for aggregate_cohort — schema validation + concat correctness."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from aggregate_cohort import (
    REQUIRED_COLUMNS,
    aggregate,
    load_report,
    main,
    resolve_inputs,
)


def _write_report(path: Path, patient_id: str, rows: list[tuple[str, str, str, str]]) -> None:
    """rows: list of (stage, metric, value, notes) tuples."""
    df = pd.DataFrame(
        [{"patient_id": patient_id, "stage": s, "metric": m, "value": v, "notes": n} for s, m, v, n in rows],
        columns=REQUIRED_COLUMNS,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def test_aggregate_two_patients(tmp_path: Path):
    p1 = tmp_path / "p1.tsv"
    p2 = tmp_path / "p2.tsv"
    _write_report(p1, "patient_001", [("junction_filtering", "input", "100", "")])
    _write_report(p2, "patient_002", [("mhc_prediction", "strong_count", "12", "")])

    cohort = aggregate([p1, p2])
    assert len(cohort) == 2
    assert sorted(cohort["patient_id"].unique()) == ["patient_001", "patient_002"]
    assert list(cohort.columns) == list(REQUIRED_COLUMNS)


def test_aggregate_sorts_by_patient_stage_metric(tmp_path: Path):
    p1 = tmp_path / "p1.tsv"
    p2 = tmp_path / "p2.tsv"
    _write_report(p2, "patient_002", [
        ("mhc_prediction", "strong_count", "12", ""),
        ("junction_filtering", "input", "200", ""),
    ])
    _write_report(p1, "patient_001", [
        ("mhc_prediction", "strong_count", "5", ""),
        ("junction_filtering", "input", "100", ""),
    ])

    cohort = aggregate([p2, p1])
    assert cohort["patient_id"].tolist() == [
        "patient_001", "patient_001", "patient_002", "patient_002",
    ]
    assert cohort["stage"].tolist()[:2] == ["junction_filtering", "mhc_prediction"]


def test_load_report_rejects_missing_patient_id_column(tmp_path: Path):
    bad = tmp_path / "bad.tsv"
    pd.DataFrame({"stage": ["s"], "metric": ["m"], "value": ["1"], "notes": [""]}).to_csv(
        bad, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="missing required columns"):
        load_report(bad)


def test_load_report_raises_clear_error_for_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="report.tsv not found"):
        load_report(tmp_path / "does_not_exist.tsv")


def test_resolve_inputs_explicit_paths_passthrough():
    paths = resolve_inputs(["a.tsv", "b.tsv"], None, "results")
    assert paths == [Path("a.tsv"), Path("b.tsv")]


def test_resolve_inputs_patients_resolves_against_results_root():
    paths = resolve_inputs(None, ["patient_001", "patient_002"], "results")
    assert paths == [
        Path("results/patient_001/reports/report.tsv"),
        Path("results/patient_002/reports/report.tsv"),
    ]


def test_main_writes_cohort_tsv(tmp_path: Path):
    p1 = tmp_path / "p1.tsv"
    p2 = tmp_path / "p2.tsv"
    out = tmp_path / "out" / "cohort.tsv"
    _write_report(p1, "patient_001", [("junction_filtering", "input", "100", "")])
    _write_report(p2, "patient_002", [("mhc_prediction", "strong_count", "12", "")])

    rc = main(["--inputs", str(p1), str(p2), "--output", str(out)])
    assert rc == 0
    assert out.exists()
    cohort = pd.read_csv(out, sep="\t")
    assert sorted(cohort["patient_id"].unique()) == ["patient_001", "patient_002"]
