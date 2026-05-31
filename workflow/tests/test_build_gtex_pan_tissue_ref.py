"""Unit tests for workflow/scripts/build_gtex_pan_tissue_ref.py (Issue #211).

The builder is a standalone CLI (not a Snakemake `script:` entry), so it is
loaded by path via importlib. All tests run on synthetic Snaptron TSV strings
and a monkeypatched fetch — no network.
"""

import importlib.util
import io
import urllib.error
from collections import Counter
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_gtex_pan_tissue_ref",
    Path(__file__).parents[1] / "scripts" / "build_gtex_pan_tissue_ref.py",
)
gtex = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(gtex)

# A realistic Snaptron gtexv2 header. Column order is NOT relied upon — lookup is by name.
SNAPTRON_HEADER = (
    "DataSource:Type\tsnaptron_id\tchromosome\tstart\tend\tlength\tstrand\t"
    "annotated\tleft_motif\tright_motif\tleft_annotated\tright_annotated\t"
    "samples\tsamples_count\tcoverage_sum\tcoverage_avg\tcoverage_median\tsource_dataset_id"
)


def _synthetic_snaptron_tsv() -> list:
    # header + 3 junctions: samples_count 7, 1, 0 (the 0 must be dropped at min_samples>=1).
    return [
        SNAPTRON_HEADER,
        "GTEX:I\t1\tchr22\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1,2\t7\t9\t4.5\t4\t0",
        "GTEX:I\t2\tchr22\t301\t400\t99\t-\t0\tCT\tAC\t0\t0\t3\t1\t2\t2.0\t2\t0",
        "GTEX:I\t3\tchr22\t501\t600\t99\t+\t1\tGT\tAG\t1\t1\t\t0\t0\t0.0\t0\t0",
    ]


# ----- Task 2: build_col_index -----

def test_build_col_index_resolves_required_columns_by_name():
    idx = gtex.build_col_index(SNAPTRON_HEADER)
    assert idx["chromosome"] == 2
    assert idx["start"] == 3
    assert idx["end"] == 4
    assert idx["strand"] == 6
    assert idx["samples_count"] == 13


def test_build_col_index_rejects_missing_columns():
    with pytest.raises(ValueError, match="missing columns"):
        gtex.build_col_index("snaptron_id\tchromosome\tstart\tend")  # no strand / samples_count


# ----- Task 3: parse_snaptron_line -----

def test_parse_snaptron_line_applies_start_minus_one_transform():
    idx = gtex.build_col_index(SNAPTRON_HEADER)
    fields = (
        "GTEX:I\t42\tchr22\t16062315\t16063236\t921\t+\t"
        "1\tGT\tAG\t1\t1\t10,55,99\t7\t123\t17.6\t12\t0"
    ).split("\t")
    parsed = gtex.parse_snaptron_line(fields, idx)
    assert parsed == ("chr22", 16062314, 16063236, "+", 7)  # start-1, end passthrough


def test_parse_snaptron_line_rejects_malformed():
    idx = gtex.build_col_index(SNAPTRON_HEADER)
    assert gtex.parse_snaptron_line(["too", "few"], idx) is None
    bad = (
        "GTEX:I\t42\tchr22\tNOT_AN_INT\t16063236\t921\t+\t"
        "1\tGT\tAG\t1\t1\t10\t7\t123\t17.6\t12\t0"
    ).split("\t")
    assert gtex.parse_snaptron_line(bad, idx) is None


# ----- Task 4: accumulate_union -----

def test_accumulate_union_keeps_min_samples_and_transforms():
    keys, sweep = gtex.accumulate_union(_synthetic_snaptron_tsv(), min_samples=1)
    assert keys == {("chr22", 100, 200, "+"), ("chr22", 300, 400, "-")}
    assert sweep[1] == 2    # two junctions with samples_count >= 1
    assert sweep[5] == 1    # only the samples_count=7 junction
    assert sweep[10] == 0


def test_accumulate_union_min_samples_gate():
    keys, _ = gtex.accumulate_union(_synthetic_snaptron_tsv(), min_samples=5)
    assert keys == {("chr22", 100, 200, "+")}  # only samples_count=7 survives


def test_accumulate_union_restrict_chrom():
    lines = _synthetic_snaptron_tsv() + [
        "GTEX:I\t4\tchr1\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1\t9\t9\t9.0\t9\t0"
    ]
    keys, _ = gtex.accumulate_union(lines, min_samples=1, restrict_chrom="chr22")
    assert all(k[0] == "chr22" for k in keys)


def test_accumulate_union_empty_input():
    assert gtex.accumulate_union([], min_samples=1) == (set(), Counter())


# ----- Task 5: fetch_snaptron_region (monkeypatched urlopen — no real network) -----

def test_fetch_snaptron_region_yields_decoded_lines(monkeypatch):
    payload = (SNAPTRON_HEADER + "\n"
               "GTEX:I\t1\tchr22\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1\t7\t9\t4.5\t4\t0\n")

    class _FakeResp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            self.close()

    def _fake_urlopen(url, timeout=None):
        assert "regions=chr22:1-50818468" in url
        return _FakeResp(payload.encode("utf-8"))

    monkeypatch.setattr(gtex.urllib.request, "urlopen", _fake_urlopen)
    lines = list(gtex.fetch_snaptron_region("chr22:1-50818468"))
    assert lines[0] == SNAPTRON_HEADER
    assert lines[1].startswith("GTEX:I\t1\tchr22")


def test_fetch_snaptron_region_retries_then_raises(monkeypatch):
    calls = {"n": 0}

    def _always_fail(url, timeout=None):
        calls["n"] += 1
        raise urllib.error.URLError("boom")

    monkeypatch.setattr(gtex.urllib.request, "urlopen", _always_fail)
    monkeypatch.setattr(gtex.time, "sleep", lambda s: None)  # no real backoff wait
    with pytest.raises(urllib.error.URLError):
        list(gtex.fetch_snaptron_region("chr22:1-50818468", retries=3))
    assert calls["n"] == 3


# ----- Task 6: writers -----

def test_write_bed6_sorted_with_strand(tmp_path):
    keys = {
        ("chr22", 300, 400, "-"),
        ("chr22", 100, 200, "+"),
    }
    out = tmp_path / "panel.bed"
    gtex.write_bed6(keys, str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "chr22\t100\t200\tchr22:100-200:+\t0\t+"
    assert lines[1] == "chr22\t300\t400\tchr22:300-400:-\t0\t-"


def test_write_qc_sidecar_reports_sweep(tmp_path):
    sweep = Counter({1: 880769, 2: 500000, 5: 120000, 10: 40000, 20: 9000})
    out = tmp_path / "panel.qc.tsv"
    gtex.write_qc_sidecar(sweep, n_junctions=880769, path=str(out))
    text = out.read_text()
    assert "n_junctions\t880769" in text
    assert "min_samples_count\tn_junctions" in text
    assert "1\t880769" in text
    assert "20\t9000" in text


# ----- Task 7: build() orchestrator + restrict-chrom -----

def test_build_end_to_end_injected_fetch(tmp_path):
    region_lines = {
        "chr22:1-50818468": _synthetic_snaptron_tsv(),
        "chr21:1-46709983": [
            SNAPTRON_HEADER,
            "GTEX:I\t9\tchr21\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1\t12\t9\t4.5\t4\t0",
        ],
    }
    bed = tmp_path / "out.bed"
    qc = tmp_path / "out.qc.tsv"
    n = gtex.build(
        regions=["chr22:1-50818468", "chr21:1-46709983"],
        bed_path=str(bed), qc_path=str(qc), min_samples=1,
        line_source=lambda region: region_lines[region],
    )
    lines = bed.read_text().splitlines()
    assert lines[0] == "chr21\t100\t200\tchr21:100-200:+\t0\t+"
    assert "chr22\t100\t200\tchr22:100-200:+\t0\t+" in lines
    assert n == 3  # 2 from chr22 (sc 7,1) + 1 from chr21


def test_build_restrict_chrom_fixture(tmp_path):
    region_lines = {"chr22:1-50818468": _synthetic_snaptron_tsv()}
    bed = tmp_path / "out.chr22.bed"
    qc = tmp_path / "out.chr22.qc.tsv"
    gtex.build(
        regions=["chr22:1-50818468"], bed_path=str(bed), qc_path=str(qc),
        min_samples=1, restrict_chrom="chr22",
        line_source=lambda region: region_lines[region],
    )
    assert all(line.startswith("chr22\t") for line in bed.read_text().splitlines())
