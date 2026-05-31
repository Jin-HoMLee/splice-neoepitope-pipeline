# GTEx Pan-Tissue Junction Reference Set — Implementation Plan

> **✅ Scientist-approved + lab-owner governance-confirmed (2026-05-31).** Source is the
> **Snaptron `gtexv2`** endpoint (recount3 / GTEx v8 / hg38 / **19,214 samples** / ~33M junctions,
> novel + per-sample). The originally-scoped source (GTEx V10 portal `junctions.gct.gz`) was
> **rejected** — it is annotation-only (~99.7% annotated, chr1 + chr22 verified) and this gate
> only ever sees *novel* junctions, so it would have filtered ≈ 0 candidates (a silent no-op).
> Recon + rejection evidence: [`docs/gtex_pan_tissue_build.md`](../../gtex_pan_tissue_build.md);
> [recon comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587336618);
> [Scientist sign-off](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587458018);
> [governance confirm](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587902847).
> **Task 1 (recon) is complete.** Tasks 2–10 below are the live, Snaptron-targeted plan.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible, GCS-staged pan-tissue **novel-junction blacklist BED** from the Snaptron `gtexv2` compilation, consumed downstream by [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) as a population-level normal reference.

**Architecture:** A standalone one-shot script `workflow/scripts/build_gtex_pan_tissue_ref.py` (NOT part of the per-patient Snakemake DAG) queries the Snaptron `gtexv2` HTTP endpoint **per chromosome**, streams the headered TSV response line-by-line, keeps every junction with `samples_count >= min_samples` (default 1 = pan-tissue), and emits a sorted **BED6** blacklist plus a QC sidecar (a `samples_count` sensitivity sweep). The pure parse/accumulate functions are unit-tested on tiny synthetic TSV fixtures; the fetch layer is isolated behind an injectable line source so the end-to-end test needs no network. A chr22-restricted output is produced as a committed local-pipeline-test fixture and must reproduce #225's cached `chr22_gtex_panel.parquet`.

**Tech Stack:** Python 3.13 (test venv `workflow/tests/.venv`), **stdlib only** (`urllib`, `argparse`, `logging`, `collections`) — mirroring [build_reference_junctions.py](../../../workflow/scripts/build_reference_junctions.py). Network calls are isolated in one fetch function (ported from #225's `fetch_snaptron_chr22`); everything else is pure and tested. `gsutil`/`gcloud` for GCS staging. Pytest for unit tests.

---

## Key facts established during planning (verified, not assumed)

- **Source = Snaptron `gtexv2`** (`https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=<region>`). recount3 reprocessing of **GTEx v8**, hg38, **19,214 samples**, ~33M junctions, with per-junction `samples_count`. NOT the V10 portal (19,788 samples / annotation-only) and NOT the `/gtex/` recount2 endpoint (9,662 samples). Provenance citations: recount3 (Wilks et al., *Genome Biology* 2021) + Snaptron (Wilks et al., *Bioinformatics* 2018).
- **Response is a headered TSV with named columns.** Snaptron `gtexv2` returns a header row; the columns this build needs — `chromosome`, `start`, `end`, `strand`, `samples_count` — are addressed **by name** (not fixed position). Other columns (`snaptron_id`, `annotated`, `samples`, coverage stats, …) are ignored. Each junction is **one row** with its own `samples_count` across all samples — so per-region there is no within-region dedup to do.
- **Coordinate transform is the #1 silent-failure risk — and it is already pinned.** Snaptron `start`/`end` are 1-based-inclusive intron donor/acceptor → BED 0-based half-open via **`bed_start = start − 1`, `bed_end = end`**. This is the *same* transform #225's `snaptron_to_key_set()` uses and it was validated **259/259** (every matched-normal ∩ GENCODE ground-truth intron present in the GTEx panel under this normalisation, 2026-05-21). Empirically banked — do not re-derive.
- **Strand IS carried** (unlike the rejected GCT source). Snaptron emits a `strand` column, so the blacklist key is the full 4-tuple `(chrom, start, end, strand)` — a **true drop-in** for `filter_junctions.py`'s reference reader ([filter_junctions.py:56-71](../../../workflow/scripts/filter_junctions.py#L56-L71) keys on cols 1/2/3/6 = `(chrom,start,end,strand)`). **This removes the strand-agnostic special-case** the GCT plan needed — see the updated #212 cross-note at the end.
- **Output format is BED6**, byte-identical in shape to [build_reference_junctions.py:131-150](../../../workflow/scripts/build_reference_junctions.py#L131-L150): `chrom start end name score strand` with `name = chrom:start-end:strand`, score `0`. UCSC `chr*` naming (Snaptron `gtexv2` ships `chr`-prefixed, matching the pipeline's GENCODE primary assembly — CLAUDE.md).
- **Scale.** chr22 alone ≈ **880,769** junctions at `samples_count >= 1` (= #225's cached `chr22_gtex_panel.parquet`). Genome-wide at `min_samples=1` is ~33M rows → the BED is hundreds of MB; keep it out of git (GCS + `data_manifest.yaml`), commit only the chr22 fixture. The heavy genome-wide build runs on a highmem VM, not the M1.

## File Structure

- **Create** `workflow/scripts/build_gtex_pan_tissue_ref.py` — the one-shot builder. Pure functions (`build_col_index`, `parse_snaptron_line`, `accumulate_union`, `write_bed6`, `write_qc_sidecar`) + one network function (`fetch_snaptron_region`) + a `build()` orchestrator with an injectable line source + a `--restrict-chrom` / `--region` filter for the chr22 fixture + a CLI entry point. No Snakemake `script:` entry (not in the DAG), so a plain `if __name__ == "__main__"` CLI only.
- **Create** `workflow/tests/test_build_gtex_pan_tissue_ref.py` — unit tests on synthetic Snaptron TSV strings + a monkeypatched fetch.
- **Modify** `config/config.yaml` — add the `gtex_filter:` block (consumed by #212; declared here so the reference path + provenance live in one place).
- **Modify** `docs/gtex_pan_tissue_build.md` — append a `## Build command` runbook section (exact genome-wide invocation, GCS staging, refresh procedure). The recon findings + provenance + transform are already in this doc.
- **Produce (artifacts, chr22 fixture committed; genome-wide staged to GCS)** `gtex_gtexv2_pan_tissue_junctions.bed`, `gtex_gtexv2_pan_tissue_junctions.qc.tsv`, `gtex_gtexv2_pan_tissue_junctions.chr22.bed` → `gs://splice-neoepitope-project/resources/gtex/gtexv2/` + a committed `data_manifest.yaml`. The chr22 fixture (small) is committed under `resources/test/`.

---

### Task 1: Data reconnaissance — COMPLETE ✅

**Files:** `docs/gtex_pan_tissue_build.md` (findings).

Done (commits `de930f2`, `046e5a5`). Established: the V10 portal `.gct` is annotation-only (rejected); the correct source is Snaptron `gtexv2`; the coordinate transform (`start−1, end`) is validated 259/259 via #225; provenance pinned to 19,214 / v8 / recount3. No production code in this task. Evidence in [`docs/gtex_pan_tissue_build.md`](../../gtex_pan_tissue_build.md).

- [x] Source identified + rejected source documented
- [x] Coordinate transform pinned + empirically validated (259/259, #225)
- [x] Provenance + citations recorded

---

### Task 2: Map Snaptron header → column index

**Files:**
- Create: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
# workflow/tests/test_build_gtex_pan_tissue_ref.py
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_gtex_pan_tissue_ref",
    Path(__file__).parents[1] / "scripts" / "build_gtex_pan_tissue_ref.py",
)
gtex = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(gtex)

# A realistic Snaptron gtexv2 header (column order is NOT relied upon; lookup is by name).
SNAPTRON_HEADER = (
    "DataSource:Type\tsnaptron_id\tchromosome\tstart\tend\tlength\tstrand\t"
    "annotated\tleft_motif\tright_motif\tleft_annotated\tright_annotated\t"
    "samples\tsamples_count\tcoverage_sum\tcoverage_avg\tcoverage_median\tsource_dataset_id"
)


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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k build_col_index -v`
Expected: FAIL — `build_gtex_pan_tissue_ref.py` does not exist / `build_col_index` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""build_gtex_pan_tissue_ref.py — Build a pan-tissue novel-junction blacklist BED
from the Snaptron gtexv2 endpoint (recount3 / GTEx v8 / hg38 / 19,214 samples).

One-shot reference builder (NOT part of the per-patient Snakemake DAG). See
docs/gtex_pan_tissue_build.md for provenance + the coordinate transform, pinned
against real data and validated 259/259 via Issue #225's helpers (Issue #211).

Output BED is BED6 (chrom start end name score strand), drop-in with
references/reference_junctions.bed and filter_junctions.py's 4-tuple reference
reader. Snaptron carries strand, so (chrom,start,end,strand) is used end-to-end.
"""

import argparse
import logging
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SNAPTRON_GTEXV2_URL = "https://snaptron.cs.jhu.edu/gtexv2/snaptron"

# Columns this build addresses by name (Snaptron gtexv2 TSV carries a header row).
_REQUIRED_COLS = ("chromosome", "start", "end", "strand", "samples_count")

# samples_count thresholds reported in the QC sensitivity sweep.
QC_SWEEP_THRESHOLDS = (1, 2, 5, 10, 20)

# hg38 (UCSC) primary chromosome sizes — the genome-wide region set (one query each).
# chr22 = 50,818,468 matches Issue #225's SNAPTRON_CHR22_REGION exactly.
HG38_CHROM_SIZES: Dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def build_col_index(header_line: str) -> Dict[str, int]:
    """Map Snaptron column names -> positional indices.

    Lookup is by name so a future Snaptron column reorder does not silently
    shift the parse. Raises ValueError if any required column is absent.
    """
    cols = header_line.rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(cols)}
    missing = [c for c in _REQUIRED_COLS if c not in idx]
    if missing:
        raise ValueError(
            f"Snaptron header missing columns {missing}; got {cols}"
        )
    return idx
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k build_col_index -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): Snaptron header->column-index map (issue #211)"
```

---

### Task 3: Parse a Snaptron row → BED coords (the pinned transform)

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
def test_parse_snaptron_line_applies_start_minus_one_transform():
    idx = gtex.build_col_index(SNAPTRON_HEADER)
    # A row whose 1-based intron is 16,062,315..16,063,236 (+ strand), seen in 7 samples.
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k parse_snaptron_line -v`
Expected: FAIL — `parse_snaptron_line` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def parse_snaptron_line(
    fields: List[str], col_idx: Dict[str, int]
) -> Optional[Tuple[str, int, int, str, int]]:
    """Parse one Snaptron data row into (chrom, bed_start, bed_end, strand, samples_count).

    Snaptron start/end are 1-based inclusive intron donor/acceptor; convert to
    BED 0-based half-open via bed_start = start - 1, bed_end = end (the transform
    pinned + validated 259/259 in Issue #225's snaptron_to_key_set). Returns None
    on a short or non-integer row.
    """
    try:
        chrom = fields[col_idx["chromosome"]]
        start = int(fields[col_idx["start"]])
        end = int(fields[col_idx["end"]])
        strand = fields[col_idx["strand"]]
        samples_count = int(fields[col_idx["samples_count"]])
    except (IndexError, ValueError):
        return None
    return chrom, start - 1, end, strand, samples_count
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k parse_snaptron_line -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): parse_snaptron_line with start-1 transform (issue #211)"
```

---

### Task 4: Accumulate the per-region union + QC sweep

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
def _synthetic_snaptron_tsv() -> list:
    # header + 3 junctions: samples_count 7, 1, 0 (the 0 must be dropped at min_samples>=1).
    return [
        SNAPTRON_HEADER,
        "GTEX:I\t1\tchr22\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1,2\t7\t9\t4.5\t4\t0",
        "GTEX:I\t2\tchr22\t301\t400\t99\t-\t0\tCT\tAC\t0\t0\t3\t1\t2\t2.0\t2\t0",
        "GTEX:I\t3\tchr22\t501\t600\t99\t+\t1\tGT\tAG\t1\t1\t\t0\t0\t0.0\t0\t0",
    ]


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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k accumulate_union -v`
Expected: FAIL — `accumulate_union` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def accumulate_union(
    lines: Iterable[str],
    min_samples: int,
    restrict_chrom: Optional[str] = None,
) -> Tuple[Set[Tuple[str, int, int, str]], Counter]:
    """Consume Snaptron TSV lines (first line = header) into a junction-key set.

    Returns (keys, sweep):
      - keys: {(chrom, bed_start, bed_end, strand)} for samples_count >= min_samples
        (and chrom == restrict_chrom, if given).
      - sweep: Counter mapping each QC_SWEEP_THRESHOLDS value t -> number of
        (restrict-filtered) junctions with samples_count >= t. Drives the QC sidecar.

    Each Snaptron junction is one row, so within a single region there is no
    dedup; the set still guards against any Snaptron-side duplicate.
    """
    it = iter(lines)
    try:
        header = next(it)
    except StopIteration:
        return set(), Counter()
    col_idx = build_col_index(header)

    keys: Set[Tuple[str, int, int, str]] = set()
    sweep: Counter = Counter()
    n_rows = 0
    for line in it:
        if not line:
            continue
        parsed = parse_snaptron_line(line.split("\t"), col_idx)
        if parsed is None:
            continue
        chrom, bstart, bend, strand, sc = parsed
        if restrict_chrom is not None and chrom != restrict_chrom:
            continue
        n_rows += 1
        for t in QC_SWEEP_THRESHOLDS:
            if sc >= t:
                sweep[t] += 1
        if sc >= min_samples:
            keys.add((chrom, bstart, bend, strand))
    log.info(
        "Accumulated %d junctions (min_samples=%d) from %d parsed rows",
        len(keys), min_samples, n_rows,
    )
    return keys, sweep
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k accumulate_union -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): accumulate per-region union + QC sweep (issue #211)"
```

---

### Task 5: Fetch a Snaptron region (network, with retry) — ported from #225

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test (monkeypatch urlopen — no real network)**

```python
import io
import urllib.error


def test_fetch_snaptron_region_yields_decoded_lines(monkeypatch):
    payload = (SNAPTRON_HEADER + "\n"
               "GTEX:I\t1\tchr22\t101\t200\t99\t+\t1\tGT\tAG\t1\t1\t1\t7\t9\t4.5\t4\t0\n")

    class _FakeResp(io.BytesIO):
        def __enter__(self): return self
        def __exit__(self, *a): self.close()

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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k fetch_snaptron_region -v`
Expected: FAIL — `fetch_snaptron_region` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def fetch_snaptron_region(
    region: str,
    endpoint: str = SNAPTRON_GTEXV2_URL,
    timeout_s: int = 180,
    retries: int = 2,
) -> Iterator[str]:
    """Stream a Snaptron region query response as decoded, newline-stripped lines.

    Ported from Issue #225's fetch_snaptron_chr22: connect with a bounded retry
    (transient URLError -> 5s backoff), then yield lines lazily so a whole
    chromosome's TSV is never fully held in memory. The connection is retried
    only at open time; a mid-stream failure surfaces to the caller.
    """
    url = f"{endpoint}?regions={region}"
    attempt = 0
    resp = None
    while True:
        attempt += 1
        try:
            log.info("Snaptron query (attempt %d/%d): %s", attempt, retries, url)
            resp = urllib.request.urlopen(url, timeout=timeout_s)
            break
        except urllib.error.URLError as e:
            if attempt >= retries:
                raise
            log.warning("Snaptron query failed (%s); retrying in 5s", e)
            time.sleep(5)
    with resp:
        for raw in resp:
            yield raw.decode("utf-8").rstrip("\n")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k fetch_snaptron_region -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): fetch_snaptron_region with retry (issue #211)"
```

---

### Task 6: Write BED6 + QC sidecar

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k "write_bed6 or write_qc" -v`
Expected: FAIL — writers undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def write_bed6(keys: Set[Tuple[str, int, int, str]], path: str) -> None:
    """Write the pan-tissue blacklist as sorted BED6 (chrom,start,end,name,0,strand).

    name = 'chrom:start-end:strand' and score 0 — identical in shape to
    build_reference_junctions.py, so filter_junctions.py's 4-tuple reference
    reader consumes it unchanged.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(keys, key=lambda k: (k[0], k[1], k[2], k[3]))
    with out.open("w") as fh:
        for chrom, start, end, strand in ordered:
            name = f"{chrom}:{start}-{end}:{strand}"
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    log.info("Wrote %d junctions to %s", len(ordered), out)


def write_qc_sidecar(sweep: Counter, n_junctions: int, path: str) -> None:
    """Write a QC TSV: total union size + the samples_count sensitivity sweep."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"n_junctions\t{n_junctions}\n")
        fh.write("min_samples_count\tn_junctions\n")
        for t in QC_SWEEP_THRESHOLDS:
            fh.write(f"{t}\t{sweep.get(t, 0)}\n")
    log.info("Wrote QC sidecar to %s", out)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k "write_bed6 or write_qc" -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): BED6 + QC sweep writers (issue #211)"
```

---

### Task 7: `build()` orchestrator + CLI (`--region` / `--restrict-chrom`, genome-wide default)

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test (end-to-end, fetch injected — no network)**

```python
def test_build_end_to_end_injected_fetch(tmp_path):
    # Inject a fetcher so build() needs no network. Two regions, one row each.
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
    # chr21 sorts before chr22; both kept junctions present; strand carried.
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k build_ -v`
Expected: FAIL — `build` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def build(
    regions: List[str],
    bed_path: str,
    qc_path: str,
    min_samples: int = 1,
    restrict_chrom: Optional[str] = None,
    endpoint: str = SNAPTRON_GTEXV2_URL,
    line_source=None,
) -> int:
    """Query each region, union the kept junctions, write BED6 + QC. Returns union size.

    line_source(region) -> iterable[str] is injectable for tests; default is the
    live Snaptron fetch. Per-region QC sweeps are summed (regions are disjoint
    per-chromosome queries, so no double counting).
    """
    if line_source is None:
        line_source = lambda region: fetch_snaptron_region(region, endpoint=endpoint)

    union: Set[Tuple[str, int, int, str]] = set()
    sweep: Counter = Counter()
    for region in regions:
        keys, region_sweep = accumulate_union(
            line_source(region), min_samples=min_samples, restrict_chrom=restrict_chrom
        )
        union |= keys
        sweep += region_sweep
        log.info("Region %s: +%d junctions (running total %d)", region, len(keys), len(union))

    write_bed6(union, bed_path)
    write_qc_sidecar(sweep, n_junctions=len(union), path=qc_path)
    return len(union)


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a Snaptron gtexv2 pan-tissue novel-junction blacklist BED."
    )
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-qc", required=True)
    parser.add_argument("--min-samples", type=int, default=1,
                        help="Keep junctions seen in >= this many samples (default 1).")
    parser.add_argument("--endpoint", default=SNAPTRON_GTEXV2_URL)
    parser.add_argument(
        "--region", action="append", default=None,
        help="Snaptron region (e.g. chr22:1-50818468). Repeatable. "
             "Default: all hg38 primary chromosomes (genome-wide).")
    parser.add_argument("--restrict-chrom", default=None,
                        help="Emit only this chromosome (chr22 fixture build).")
    args = parser.parse_args()

    if args.region:
        regions = args.region
    else:
        regions = [f"{c}:1-{size}" for c, size in HG38_CHROM_SIZES.items()]

    n = build(
        regions=regions, bed_path=args.output_bed, qc_path=args.output_qc,
        min_samples=args.min_samples, restrict_chrom=args.restrict_chrom,
        endpoint=args.endpoint,
    )
    log.info("Done: %d junctions in the pan-tissue union", n)


if __name__ == "__main__":
    _cli_main()
```

- [ ] **Step 4: Run the full test file**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -v`
Expected: PASS (all tests across Tasks 2–7)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): build() orchestrator + genome-wide CLI (issue #211)"
```

---

### Task 8: Config schema + runbook

**Files:**
- Modify: `config/config.yaml`
- Modify: `docs/gtex_pan_tissue_build.md`

- [ ] **Step 1: Add the `gtex_filter` block to config**

Add near the existing `reference:` neighbourhood in `config/config.yaml` (consumed by #212; `reference_bed` is the Task-9 GCS staging target):
```yaml
gtex_filter:
  enabled: false                 # opt-in here in #211; #212 flips the default at integration
  source: snaptron_gtexv2        # recount3 / GTEx v8 / hg38 / 19,214 samples
  endpoint: "https://snaptron.cs.jhu.edu/gtexv2/snaptron"
  min_samples: 1                 # most aggressive — precision over recall (vaccine safety, #126)
  reference_bed: "gs://splice-neoepitope-project/resources/gtex/gtexv2/gtex_gtexv2_pan_tissue_junctions.bed"
  reference_bed_chr22: "resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed"
```

- [ ] **Step 2: Verify config still parses**

Run: `conda activate snakemake && python -c "import yaml; yaml.safe_load(open('config/config.yaml')); print('config OK')"`
Expected: `config OK`

- [ ] **Step 3: Write the build runbook**

Append to `docs/gtex_pan_tissue_build.md` a `## Build command` section: the exact genome-wide invocation, the chr22-fixture invocation, the GCS staging + `data_manifest.yaml` commands (Task 9), and an explicit note that this is a manual one-shot (not in the per-patient DAG) + how to refresh on a new Snaptron/recount release.

- [ ] **Step 4: Commit**

```bash
git add config/config.yaml docs/gtex_pan_tissue_build.md
git commit -m "feat(gtex): gtex_filter config block + build runbook (issue #211)"
```

---

### Task 9: Production genome-wide build + GCS staging + chr22 fixture (heavy; on a VM)

**Files:**
- Produce + stage artifacts (genome-wide BED staged to GCS, chr22 fixture committed, `data_manifest.yaml` committed)

This task runs the real genome-wide Snaptron build. It is heavy (24 chromosome queries, ~33M junctions, hundreds of MB) and memory-bounded by the union set — run on the highmem VM, not the M1. **Needs explicit go-ahead** (network + VM time).

- [ ] **Step 1: chr22 fixture FIRST — it gates correctness**

```bash
conda activate snakemake
python workflow/scripts/build_gtex_pan_tissue_ref.py \
  --region chr22:1-50818468 --restrict-chrom chr22 \
  --output-bed resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed \
  --output-qc /tmp/gtex_gtexv2_pan_tissue_junctions.chr22.qc.tsv \
  --min-samples 1
wc -l resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed   # expect ~880,769
```

- [ ] **Step 2: Reproduce #225's cached panel (consistency check)**

The chr22 BED must contain the same junction keys as #225's `chr22_gtex_panel.parquet` (~880,769 @ `samples_count >= 1`). Compare key sets:
```bash
research/.venv/bin/python - <<'PY'
import pandas as pd
panel = pd.read_parquet("research/experiments/issue_225_normal_junction_filter_strength/outputs/chr22_gtex_panel.parquet")
panel = panel[(panel.samples_count >= 1) & (panel.chromosome == "chr22")]
ref = {(c, s-1, e, st) for c, s, e, st in zip(panel.chromosome, panel.start, panel.end, panel.strand)}
bed = set()
for line in open("resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed"):
    c, s, e, _n, _sc, st = line.rstrip("\n").split("\t")
    bed.add((c, int(s), int(e), st))
print(f"parquet keys: {len(ref):,}  bed keys: {len(bed):,}  identical: {ref == bed}")
print(f"only-in-parquet: {len(ref - bed):,}  only-in-bed: {len(bed - ref):,}")
PY
```
Expected: `identical: True` (or a near-zero symmetric diff explained by a Snaptron refresh). **A large diff means the transform/parse regressed — STOP and diff a few keys.**

- [ ] **Step 3: Genome-wide build**

```bash
python workflow/scripts/build_gtex_pan_tissue_ref.py \
  --output-bed /tmp/gtex_gtexv2_pan_tissue_junctions.bed \
  --output-qc /tmp/gtex_gtexv2_pan_tissue_junctions.qc.tsv \
  --min-samples 1
head /tmp/gtex_gtexv2_pan_tissue_junctions.bed
cat  /tmp/gtex_gtexv2_pan_tissue_junctions.qc.tsv   # sanity-check the sweep monotonically decreases
```

- [ ] **Step 4: Coordinate cross-check against the GENCODE reference (the #148/#370 silent-empty guard)**

```bash
conda activate snakemake
python workflow/scripts/build_reference_junctions.py \
  --gtf references/gencode.v47.annotation.gtf.gz --output /tmp/ref_full.bed
comm -12 \
  <(grep -P '^chr22\t' /tmp/gtex_gtexv2_pan_tissue_junctions.bed | cut -f1-3 | sort -u) \
  <(grep -P '^chr22\t' /tmp/ref_full.bed | cut -f1-3 | sort -u) | wc -l
```
Expected: a **substantial non-zero** chr22 overlap (most annotated junctions are expressed somewhere in GTEx). **Zero ⇒ coordinate transform wrong — STOP.**

- [ ] **Step 5: Stage to GCS + commit the chr22 fixture + `data_manifest.yaml`**

```bash
gsutil cp /tmp/gtex_gtexv2_pan_tissue_junctions.bed    gs://splice-neoepitope-project/resources/gtex/gtexv2/
gsutil cp /tmp/gtex_gtexv2_pan_tissue_junctions.qc.tsv gs://splice-neoepitope-project/resources/gtex/gtexv2/
# data_manifest.yaml: artifact GCS paths + sha256 + the exact build command (CLAUDE.md >100 MB rule)
git add resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed \
        research/experiments/issue_211_gtex_pan_tissue_ref/data_manifest.yaml
git commit -m "feat(gtex): stage gtexv2 pan-tissue BED to GCS + chr22 fixture + manifest (issue #211)"
```

---

### Task 10: Re-run Issue #225 notebook §2(c) — NO-GO verdict guard (AC)

**Files:**
- Run: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` §2(c)

- [ ] **Step 1: Re-run §2(c) against the production panel from Task 9**

Open the notebook in `research/.venv`; point the GTEx-panel cell at the staged production BED (or the chr22 fixture for a fast check) and re-run §2(c). Because the production panel is the genome-wide extension of #225's exact chr22 source, this is now a **clean consistency check**, not a new experiment.

- [ ] **Step 2: Confirm the verdict holds**

Expected: the AlphaGenome-as-3rd-filter **NO-GO** verdict still holds (the production panel is at least as inclusive as the chr22 proxy at F1-max τ). Record before/after numbers.

- [ ] **Step 3: Lab notebook + commit (post-review, per the lab-notebook timing rule)**

Capture the §2(c) re-run conclusion in `research/lab_notebook/developer.md` AFTER review feedback is incorporated, before merge.

---

## Acceptance-criteria coverage (self-review against Issue #211)

- [x] Snaptron `gtexv2` provenance pinned in config (recount3 / v8 / hg38 / 19,214 samples + citations) → Task 8 + doc
- [x] Genome-wide region-query union, clean exit + retry on endpoint error → Tasks 5+7 (per-chromosome `--region` default; `fetch_snaptron_region` retry)
- [x] `min_samples` exposed as config param, default 1 → Task 7 CLI + Task 8 config
- [x] BED6, sorted/dedup'd, UCSC `chr*`, strand-carrying, transform `start-1,end` → Tasks 3+6
- [x] QC sidecar reports the `samples_count` sweep (≥ {1,2,5,10,20}); count-only → Tasks 4+6
- [x] BED + QC staged to GCS; `data_manifest.yaml` committed → Task 9 Step 5
- [x] chr22 fixture produced AND reproduces #225's `chr22_gtex_panel.parquet` (~880,769) → Task 9 Steps 1–2
- [x] Pytest covers Snaptron parse + union construction on synthetic input → Tasks 2–7 tests
- [x] METHODS caveat (recount3/STAR vs HISAT2/STAR → conservative under-filtering) → doc (already in `gtex_pan_tissue_build.md`)
- [x] Re-run #225 §2(c) against production panel; NO-GO holds → Task 10

## Cross-issue note for [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) (integration)

**Updated for Snaptron (strand IS carried — supersedes the GCT plan's strand-agnostic note).** The blacklist BED carries real strand, so #212 can intersect tumour junctions on the full 4-tuple `(chrom, start, end, strand)` — i.e. **reuse `filter_junctions.py`'s existing `_load_reference_junctions` reader verbatim**; no strand-agnostic special-case is needed. Tag filtered junctions `origin = gtex_pantissue_shared`. (The one remaining caveat: recount3/STAR vs our HISAT2/STAR may disagree by ±1 bp on a minority of junctions → conservative under-filtering, documented as a METHODS caveat — not a coordinate bug.) The **immune-privileged-tissue / cancer-testis-antigen exemption** (whether to keep testis-only junctions) is a #212 design decision — captured on [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212#issuecomment-4587476577); it needs per-tissue provenance, which is exactly why #211 ships count-only.

## Out of scope (this plan)
- `filter_junctions.py` integration / stacking with matched-normal → #212
- patient_002 / patient_001 validation runs → #212
- Immune-privileged-tissue (CTA) exemption → #212
- AlphaGenome predicted-normal axis → #203
