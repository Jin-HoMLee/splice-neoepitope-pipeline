# GTEx Pan-Tissue Junction Reference Set — Implementation Plan

> **✅ Scientist-APPROVED redirect (2026-05-31); PM governance confirm pending (non-blocking).**
> Task 1 (recon) is **complete**: the originally-scoped source (**GTEx V10 portal `junctions.gct.gz`**)
> is annotation-only → a silent no-op; the verified replacement is the **Snaptron `gtexv2`** endpoint
> ([Scientist sign-off](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587458018)). The only open item is a PM read (in-scope-without-re-commitment +
> milestone feasibility) — governance, not technical. See the **"Redirect addendum"** immediately
> below for the revised architecture + the Scientist's confirmed parameters. **Do not execute
> Tasks 2–10 as written** — they target the rejected GCT source; they need a per-task re-write onto
> Snaptron (design fully specified in the addendum). Recon evidence + redirect rationale:
> [`docs/gtex_pan_tissue_build.md`](../../gtex_pan_tissue_build.md).

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

---

## Redirect addendum (2026-05-31) — Snaptron `gtexv2`, supersedes the GCT-based tasks below

**Why:** the pipeline's normal-junction gate only sees **novel/unannotated** junctions (annotated
ones are discarded one step earlier). The portal `.gct` is ~99.7% annotated (chr1 + chr22 verified)
→ it would filter ≈ 0 candidates. The gate needs a **raw, novel-containing** population reference.
Full rationale: [`docs/gtex_pan_tissue_build.md`](../../gtex_pan_tissue_build.md) + [#211 recon comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587336618).

**Revised architecture:** `workflow/scripts/build_gtex_pan_tissue_ref.py` queries the **Snaptron
`gtexv2`** HTTP endpoint (`https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=<region>`; recount3 /
GTEx v8 / hg38 / ~19,788 samples / ~33M junctions, novel + per-sample) region-by-region across the
genome, keeps junctions with `samples_count >= min_samples` (default 1 = pan-tissue), and writes the
same sorted **BED6** blacklist. Reuse #225's helpers: `fetch_snaptron_chr22(region)` +
`snaptron_to_key_set(df, min_samples)` from `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` §2(c).

**What carries over from the recon / old plan (unchanged):**
- Coordinate transform **`bed_start = start − 1`, `bed_end = end`, strand passthrough** — Snaptron uses the *same* 1-based-inclusive-intron convention (validated 259/259 in #225). Empirically banked.
- **BED6** output (drop-in with [filter_junctions.py:56-71](../../../workflow/scripts/filter_junctions.py#L56-L71)); strand IS carried.
- The `--restrict-chrom` chr22-fixture path, the BED6 + QC writers (Task 6), the CLI shape (Task 7), the #225 §2(c) re-run AC (Task 10).

**What changes vs the GCT-based tasks below:**
- **Task 2 (`parse_sample_attributes`) — likely DROPPED.** With `min_samples=1` ("any sample"), no sample→tissue map is needed. The "seen in N tissues" QC would require mapping Snaptron sample IDs → GTEx tissue via recount3 metadata — **open question for Sci/PM** (keep the QC, or ship count-only?).
- **Task 3 (`parse_junction_name`) — simplified.** Snaptron returns parsed `chromosome/start/end/strand` columns (18-col TSV: col 8 `annotated`, col 13 `samples`, col 14 `samples_count`) — no `chr:start-end:strand` string parsing.
- **Task 4 (GCT streaming) — REPLACED** by paginated Snaptron region queries (per-chromosome or windowed; handle response size + retries, per #225's `fetch` helper).
- **Task 5 (accumulate) — simplified** to `samples_count >= min_samples` per record; optionally retain the `annotated` flag for QC.
- **Task 8 (config)** — keys become `gtex_filter.source: snaptron_gtexv2`, `snaptron_endpoint`, `min_samples` (not `min_read_count`/GCS `.gct` path).
- **Task 9 (build)** — genome-wide Snaptron query (chr22 alone ≈ 880k junctions at `samples_count≥1`; the genome-wide BED is tens of millions of rows → size for GCS, keep out of git, `data_manifest.yaml`).

**Resolved by [Scientist sign-off](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211#issuecomment-4587458018) (2026-05-31):**
- **(a) Snaptron `gtexv2` — APPROVED** over controlled-access raw V10 (novel-junction containment is the gate's whole purpose; v10 recency would forfeit it).
- **(b) `min_samples = 1` — CONFIRMED**, two conditions: expose as config param `gtex_filter.min_samples` (not hard-coded); QC sidecar must report a **`samples_count` sensitivity sweep** (union size at `>= {1, 2, 5, 10, 20}`).
- **(c) QC — COUNT-ONLY** for this slice; per-tissue provenance deferred to [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) (where the cancer-testis-antigen / immune-privileged-tissue exemption is decided).

**Additional Scientist conditions folded into the rewritten [#211 ACs](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211):**
- Provenance pinned to **19,214 samples / GTEx v8 / recount3** (NOT 19,788 portal / 9,662 recount2) + recount3 (Wilks 2021) & Snaptron (Wilks 2018) citations.
- chr22 build slice must **reproduce #225's cached `chr22_gtex_panel.parquet`** (~880,769 @ `samples_count >= 1`) — consistency check.
- METHODS caveat: recount3/STAR vs our HISAT2/STAR aligner difference → conservative under-filtering on 1-bp coord disagreements.

**Still pending (non-blocking): PM read** — (i) confirm in-scope without re-commitment (source refinement, still M); (ii) genome-wide-build feasibility vs the `i2-S3` milestone (due 2026-06-03).

---

**Goal:** Build a reproducible, GCS-staged pan-tissue splice-junction blacklist BED from GTEx V10, consumed downstream by [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) as a population-level normal reference.

**Architecture:** A standalone one-shot script `workflow/scripts/build_gtex_pan_tissue_ref.py` (NOT part of the per-patient Snakemake DAG) streams the single GTEx V10 junction-counts GCT row-by-row, joins each sample column to its tissue via the `SampleAttributesDS` annotation, and emits a sorted **BED6** blacklist of every junction with `≥ min_read_count` reads in any tissue, plus a QC sidecar. The pure parse/accumulate functions are unit-tested on tiny synthetic GCT/annotation fixtures; a chr22-restricted output is produced for local pipeline tests.

**Tech Stack:** Python 3.13 (test venv `workflow/tests/.venv`), stdlib only (`gzip`, `csv`, `argparse`, `logging`) — mirroring [build_reference_junctions.py](workflow/scripts/build_reference_junctions.py). `gsutil`/`gcloud` for GCS staging. Pytest for unit tests.

---

## Key facts established during planning (verified, not assumed)

- **GTEx V10 distributes ONE junction matrix, not 54 per-tissue files.** A single `*_junctions.gct.gz` (rows = junctions, cols = ~20k samples) + a `*_Annotations_SampleAttributesDS.txt` mapping `SAMPID → SMTSD` (tissue). The issue's "download ~54 per-tissue files" framing is superseded: stream the one GCT, group columns by tissue via the annotation.
- **Output format is BED6**, drop-in with the existing reference: [build_reference_junctions.py:131-150](workflow/scripts/build_reference_junctions.py#L131-L150) writes `chrom start end name score strand` (`name = chrom:start-end:strand`, score `0`), and [filter_junctions.py:56-71](workflow/scripts/filter_junctions.py#L56-L71) reads reference BEDs keying on cols 1/2/3/6 = `(chrom, start, end, strand)`.
- **Coordinate convention is the #1 silent-failure risk.** reference_junctions.bed uses **0-based half-open** intron coords (`start` = upstream exon end 0-based, `end` = downstream exon start 0-based — see [build_reference_junctions.py:62-128](workflow/scripts/build_reference_junctions.py#L62-L128)). GTEx junction-ID coords must be normalized to the *same* space or the downstream intersect in #212 silently matches nothing (the Issue #148 / #370 empty-result bug class). Task 1 pins GTEx's convention against a known GENCODE junction before any parser code is trusted.
- **GTEx junction IDs carry no strand.** Historically `chr_start_end`. The reference reader expects a strand in col 6. Decision (this plan): emit strand `.` and treat the GTEx blacklist as **strand-agnostic** (positional) — which is also the *more conservative / safer* choice for a vaccine off-tumour-toxicity filter. This is a **cross-issue note for #212** (see end of plan): #212 must intersect the GTEx source on `(chrom,start,end)` only, NOT reuse `_load_reference_junctions`'s 4-tuple key verbatim.

## File Structure

- **Create** `workflow/scripts/build_gtex_pan_tissue_ref.py` — the one-shot builder. Pure functions (`parse_sample_attributes`, `parse_junction_name`, `iter_gct_rows`, `accumulate_pan_tissue_union`, `write_bed6`, `write_qc_sidecar`) + a `--restrict-chrom` filter for the chr22 fixture + CLI entry point. No Snakemake `script:` entry (not in the DAG), so a plain `if __name__ == "__main__"` CLI only.
- **Create** `workflow/tests/test_build_gtex_pan_tissue_ref.py` — unit tests on synthetic GCT + annotation strings.
- **Modify** `config/config.yaml` — add the `gtex_filter:` block (consumed by #212; declared here so the reference path is documented in one place).
- **Create** `docs/gtex_pan_tissue_build.md` — short runbook: exact V10 URLs, build command, GCS staging command, refresh procedure (the build is manual/one-shot, so the runbook IS the reproducibility record).
- **Produce (artifacts, not committed if >100 MB)** `gtex_v10_pan_tissue_junctions.bed`, `gtex_v10_pan_tissue_junctions.qc.tsv`, `gtex_v10_pan_tissue_junctions.chr22.bed` → staged to `gs://splice-neoepitope-project/resources/gtex/v10/`. The chr22 fixture (small) is committed under `resources/test/`.

---

### Task 1: Data reconnaissance — pin real V10 filenames + junction-ID semantics

**Files:**
- Create: `docs/gtex_pan_tissue_build.md` (findings section)

This task writes NO production code. It de-risks the parser by recording ground truth from the real files. Do it before any parsing code so the unit-test fixtures mirror reality.

- [ ] **Step 1: Locate the V10 open-access bucket paths**

Run (the GTEx open-access bulk files live on a public GCS bucket; list to find exact names):
```bash
gsutil ls -r 'gs://adult-gtex/bulk-gex/v10/rna-seq/**' 2>/dev/null | grep -iE 'junction|SampleAttributes' || \
  echo "bucket path differs — fall back to the portal Downloads page: https://gtexportal.org/home/downloads/adult-gtex"
```
Expected: an exact `...v10..._junctions.gct.gz` URL and a `GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt` URL. If `gsutil ls` fails (auth/path), open the portal Downloads → "Bulk tissue expression" and copy the two HTTPS links. Record both exact URLs in `docs/gtex_pan_tissue_build.md` under a `## Source files (V10)` heading.

- [ ] **Step 2: Inspect the GCT header + first data row WITHOUT downloading the whole file**

Run (stream just the head of the gzip):
```bash
curl -s '<JUNCTIONS_GCT_URL>' | zcat 2>/dev/null | head -4 | cut -c1-300
```
Expected shape:
```
#1.2
<nrows>	<ncols>
Name	Description	<SAMPID_1>	<SAMPID_2>	...
<junction_id>	<gene_id>	<count>	<count>	...
```
Record in the findings doc: (a) the literal junction-ID format in the `Name` column (e.g. `chr1_14829_14969` vs `chr1:14829-14969`), (b) whether a strand token is present, (c) the column-header sample-ID format.

- [ ] **Step 3: Pin the coordinate convention against a known GENCODE junction**

Pick any junction that exists in the committed `references/reference_junctions.bed` (build it first if absent: `conda activate snakemake && python workflow/scripts/build_reference_junctions.py --gtf references/gencode.v47.annotation.gtf.gz --output /tmp/ref.bed`, then `head /tmp/ref.bed`). Take its `chrom/start/end`. Then grep the GTEx GCT `Name` column for the same locus:
```bash
curl -s '<JUNCTIONS_GCT_URL>' | zcat | cut -f1 | grep -E '<chrom>[_:]<approx_start>' | head
```
Determine whether GTEx coords equal the BED 0-based-half-open `start/end`, are 1-based, or are off-by-one/anchor-shifted. **Record the exact transform** (e.g. "GTEx start is 1-based intron start → BED start = gtex_start − 1; GTEx end = BED end") in the findings doc. This transform becomes the body of `parse_junction_name` in Task 3.

- [ ] **Step 4: Inspect the sample-attributes columns**

Run:
```bash
curl -s '<SAMPLE_ATTRIBUTES_URL>' | head -2
```
Confirm the column names: `SAMPID` (sample id, matches GCT column headers) and `SMTSD` (detailed tissue). Also note `SMAFRZE` (analysis-freeze flag; value `RNASEQ` marks the bulk RNA-seq samples to keep). Record the exact header in the findings doc.

- [ ] **Step 5: Commit the findings doc**

```bash
git add docs/gtex_pan_tissue_build.md
git commit -m "docs(gtex): pin V10 junction GCT + sample-attributes format (issue #211 recon)"
```

---

### Task 2: Parse sample attributes → sample-id → tissue map

**Files:**
- Create: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
# workflow/tests/test_build_gtex_pan_tissue_ref.py
import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "build_gtex_pan_tissue_ref",
    Path(__file__).parents[1] / "scripts" / "build_gtex_pan_tissue_ref.py",
)
gtex = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(gtex)


def test_parse_sample_attributes_maps_sampid_to_tissue(tmp_path):
    p = tmp_path / "attrs.txt"
    p.write_text(
        "SAMPID\tSMTSD\tSMAFRZE\n"
        "S-1\tWhole Blood\tRNASEQ\n"
        "S-2\tLung\tRNASEQ\n"
        "S-3\tLung\tEXCLUDE\n"          # non-RNASEQ freeze → dropped
    )
    mapping = gtex.parse_sample_attributes(str(p))
    assert mapping == {"S-1": "Whole Blood", "S-2": "Lung"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py::test_parse_sample_attributes_maps_sampid_to_tissue -v`
Expected: FAIL — `build_gtex_pan_tissue_ref.py` does not exist / `parse_sample_attributes` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""build_gtex_pan_tissue_ref.py — Build a pan-tissue splice-junction blacklist
BED from the GTEx V10 junction read-counts GCT.

One-shot reference builder (NOT part of the per-patient Snakemake DAG). See
docs/gtex_pan_tissue_build.md for source URLs + the coordinate transform pinned
against real V10 data (Issue #211).

Output BED is BED6 (chrom start end name score strand), drop-in with
references/reference_junctions.bed. GTEx junctions carry no strand → strand
column is '.', and the blacklist is positional (strand-agnostic) — the safer
choice for a vaccine off-tumour-toxicity filter. See Issue #212 cross-note.
"""

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def _open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_sample_attributes(path: str) -> Dict[str, str]:
    """Map SAMPID -> SMTSD (detailed tissue) for RNA-seq freeze samples only.

    Rows whose SMAFRZE column is not 'RNASEQ' are skipped.
    """
    mapping: Dict[str, str] = {}
    with _open_maybe_gzip(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("SMAFRZE") != "RNASEQ":
                continue
            sampid = row["SAMPID"]
            mapping[sampid] = row["SMTSD"]
    log.info("Parsed %d RNA-seq samples across %d tissues",
             len(mapping), len(set(mapping.values())))
    return mapping
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py::test_parse_sample_attributes_maps_sampid_to_tissue -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): parse_sample_attributes (issue #211)"
```

---

### Task 3: Parse a GTEx junction ID → BED coords

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

> **NOTE:** the `start`/`end` arithmetic below assumes the most common GTEx convention (junction Name = `chrom_start_end`, 1-based intron start). **If Task 1 Step 3 pinned a different transform, edit the two arithmetic lines to match the recorded transform** — the test values must reflect the real convention so the keys land in the same 0-based-half-open space as reference_junctions.bed.

- [ ] **Step 1: Write the failing test**

```python
def test_parse_junction_name_underscore_to_bed_coords():
    # GTEx Name 'chr1_14830_14969' (1-based intron start) -> BED 0-based start
    chrom, start, end = gtex.parse_junction_name("chr1_14830_14969")
    assert chrom == "chr1"
    assert start == 14829   # 1-based 14830 -> 0-based 14829
    assert end == 14969


def test_parse_junction_name_rejects_malformed():
    assert gtex.parse_junction_name("not_a_junction") is None
    assert gtex.parse_junction_name("chr1_abc_def") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k parse_junction_name -v`
Expected: FAIL — `parse_junction_name` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def parse_junction_name(name: str) -> "Tuple[str, int, int] | None":
    """Parse a GTEx junction ID into 0-based half-open BED coords.

    GTEx V10 junction Name format: 'chrom_start_end' (1-based intron start).
    Converted to BED 0-based half-open to match references/reference_junctions.bed.
    Returns (chrom, start, end) or None on malformed input.

    NOTE: GENCODE primary assembly uses UCSC 'chr' naming (CLAUDE.md). GTEx V10
    also ships 'chr'-prefixed, so no chr/MT normalization is needed; if a future
    release drops the prefix, normalize here.
    """
    parts = name.split("_")
    if len(parts) != 3:
        return None
    chrom, raw_start, raw_end = parts
    try:
        start = int(raw_start) - 1   # 1-based -> 0-based half-open start
        end = int(raw_end)
    except ValueError:
        return None
    return chrom, start, end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k parse_junction_name -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): parse_junction_name -> BED coords (issue #211)"
```

---

### Task 4: Stream GCT rows (header sample list + per-junction counts)

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
def _synthetic_gct() -> str:
    # GCT: line1 version, line2 dims, line3 header, then data rows.
    return (
        "#1.2\n"
        "2\t3\n"
        "Name\tDescription\tS-1\tS-2\tS-3\n"
        "chr1_101_200\tGENE_A\t0\t5\t0\n"
        "chr1_301_400\tGENE_B\t0\t0\t0\n"
    )


def test_iter_gct_rows_yields_samples_then_rows(tmp_path):
    p = tmp_path / "j.gct"
    p.write_text(_synthetic_gct())
    samples, rows = gtex.iter_gct_rows(str(p))
    rows = list(rows)
    assert samples == ["S-1", "S-2", "S-3"]
    assert rows[0] == ("chr1_101_200", [0, 5, 0])
    assert rows[1] == ("chr1_301_400", [0, 0, 0])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k iter_gct_rows -v`
Expected: FAIL — `iter_gct_rows` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def iter_gct_rows(path: str) -> "Tuple[List[str], Iterator[Tuple[str, List[int]]]]":
    """Open a (optionally gzipped) GCT and return (sample_ids, row_iterator).

    GCT layout: line 1 '#1.2', line 2 '<nrows>\\t<ncols>', line 3 column header
    'Name\\tDescription\\t<sample>...'. Data rows stream lazily so the full
    junctions x ~20k-sample matrix is never held in memory.
    """
    fh = _open_maybe_gzip(path)
    fh.readline()                      # '#1.2'
    fh.readline()                      # dims
    header = fh.readline().rstrip("\n").split("\t")
    samples = header[2:]               # drop Name, Description

    def _rows() -> Iterator[Tuple[str, List[int]]]:
        try:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                name = fields[0]
                counts = [int(float(x)) for x in fields[2:]]
                yield name, counts
        finally:
            fh.close()

    return samples, _rows()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k iter_gct_rows -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): stream GCT rows (issue #211)"
```

---

### Task 5: Accumulate the pan-tissue union (junction → set of tissues)

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
def test_accumulate_pan_tissue_union(tmp_path):
    p = tmp_path / "j.gct"
    p.write_text(_synthetic_gct())     # S-2 has 5 reads on chr1_101_200
    samples, rows = gtex.iter_gct_rows(str(p))
    sample_tissue = {"S-1": "Blood", "S-2": "Lung", "S-3": "Lung"}
    union = gtex.accumulate_pan_tissue_union(samples, rows, sample_tissue, min_read_count=1)
    # chr1_101_200 present (>=1 read in Lung); chr1_301_400 all-zero -> absent.
    assert ("chr1", 100, 200) in union          # 1-based 101 -> 0-based 100
    assert union[("chr1", 100, 200)] == {"Lung"}
    assert ("chr1", 300, 400) not in union


def test_accumulate_respects_min_read_count(tmp_path):
    p = tmp_path / "j.gct"
    p.write_text(_synthetic_gct())
    samples, rows = gtex.iter_gct_rows(str(p))
    sample_tissue = {"S-1": "Blood", "S-2": "Lung", "S-3": "Lung"}
    union = gtex.accumulate_pan_tissue_union(samples, rows, sample_tissue, min_read_count=6)
    assert union == {}                  # max count is 5 < 6
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k accumulate -v`
Expected: FAIL — `accumulate_pan_tissue_union` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def accumulate_pan_tissue_union(
    samples: List[str],
    rows: "Iterator[Tuple[str, List[int]]]",
    sample_tissue: Dict[str, str],
    min_read_count: int,
) -> Dict[Tuple[str, int, int], Set[str]]:
    """Build {(chrom, start, end): {tissues}} for every junction with
    >= min_read_count reads in at least one sample of that tissue.

    Samples absent from sample_tissue (non-RNASEQ freeze) are ignored. The
    tissue set per junction drives the QC sidecar ('seen in N tissues').
    """
    # Precompute column index -> tissue (skip columns we don't track).
    col_tissue: List["str | None"] = [sample_tissue.get(s) for s in samples]

    union: Dict[Tuple[str, int, int], Set[str]] = {}
    n_seen = 0
    n_kept = 0
    for name, counts in rows:
        n_seen += 1
        tissues: Set[str] = set()
        for idx, c in enumerate(counts):
            if c >= min_read_count:
                t = col_tissue[idx]
                if t is not None:
                    tissues.add(t)
        if not tissues:
            continue
        parsed = parse_junction_name(name)
        if parsed is None:
            continue
        union[parsed] = tissues
        n_kept += 1
    log.info("Pan-tissue union: %d/%d junctions kept (min_read_count=%d)",
             n_kept, n_seen, min_read_count)
    return union
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k accumulate -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): accumulate pan-tissue union (issue #211)"
```

---

### Task 6: Write BED6 + QC sidecar

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test**

```python
def test_write_bed6_sorted_strand_agnostic(tmp_path):
    union = {
        ("chr1", 300, 400): {"Lung"},
        ("chr1", 100, 200): {"Lung", "Blood"},
    }
    out = tmp_path / "panel.bed"
    gtex.write_bed6(union, str(out))
    lines = out.read_text().splitlines()
    # sorted by chrom,start,end; strand '.'; name chrom:start-end:.
    assert lines[0] == "chr1\t100\t200\tchr1:100-200:.\t0\t."
    assert lines[1] == "chr1\t300\t400\tchr1:300-400:.\t0\t."


def test_write_qc_sidecar_counts(tmp_path):
    union = {
        ("chr1", 100, 200): {"Lung", "Blood"},
        ("chr1", 300, 400): {"Lung"},
    }
    out = tmp_path / "panel.qc.tsv"
    gtex.write_qc_sidecar(union, str(out))
    text = out.read_text()
    assert "n_junctions\t2" in text
    # distribution: 1 junction in 2 tissues, 1 junction in 1 tissue
    assert "n_tissues\tn_junctions" in text
    assert "1\t1" in text
    assert "2\t1" in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k "write_bed6 or write_qc" -v`
Expected: FAIL — writers undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def write_bed6(union: Dict[Tuple[str, int, int], Set[str]], path: str) -> None:
    """Write the pan-tissue blacklist as sorted BED6 with strand '.'.

    Strand is '.' because GTEx junction IDs carry no strand; the blacklist is
    positional (strand-agnostic). See Issue #212 cross-note — downstream
    intersection must key on (chrom,start,end), not the 4-tuple incl. strand.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted(union.keys(), key=lambda k: (k[0], k[1], k[2]))
    with out.open("w") as fh:
        for chrom, start, end in keys:
            name = f"{chrom}:{start}-{end}:."
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\n")
    log.info("Wrote %d junctions to %s", len(keys), out)


def write_qc_sidecar(union: Dict[Tuple[str, int, int], Set[str]], path: str) -> None:
    """Write a QC TSV: total junction count + a 'seen in N tissues' histogram."""
    from collections import Counter
    hist = Counter(len(tissues) for tissues in union.values())
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"n_junctions\t{len(union)}\n")
        fh.write("n_tissues\tn_junctions\n")
        for n_tissues in sorted(hist):
            fh.write(f"{n_tissues}\t{hist[n_tissues]}\n")
    log.info("Wrote QC sidecar to %s", out)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k "write_bed6 or write_qc" -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): BED6 + QC sidecar writers (issue #211)"
```

---

### Task 7: CLI entry point + `--restrict-chrom` (chr22 fixture support)

**Files:**
- Modify: `workflow/scripts/build_gtex_pan_tissue_ref.py`
- Test: `workflow/tests/test_build_gtex_pan_tissue_ref.py`

- [ ] **Step 1: Write the failing test (end-to-end on synthetic files, incl. chrom restriction)**

```python
def test_build_end_to_end_with_chrom_restriction(tmp_path):
    gct = tmp_path / "j.gct"
    gct.write_text(
        "#1.2\n"
        "2\t3\n"
        "Name\tDescription\tS-1\tS-2\tS-3\n"
        "chr22_101_200\tG_A\t0\t5\t0\n"
        "chr1_101_200\tG_B\t9\t0\t0\n"     # different chrom -> excluded by restrict
    )
    attrs = tmp_path / "attrs.txt"
    attrs.write_text(
        "SAMPID\tSMTSD\tSMAFRZE\n"
        "S-1\tBlood\tRNASEQ\n"
        "S-2\tLung\tRNASEQ\n"
        "S-3\tLung\tRNASEQ\n"
    )
    bed = tmp_path / "out.bed"
    qc = tmp_path / "out.qc.tsv"
    gtex.build(
        gct_path=str(gct), attributes_path=str(attrs),
        bed_path=str(bed), qc_path=str(qc),
        min_read_count=1, restrict_chrom="chr22",
    )
    lines = bed.read_text().splitlines()
    assert lines == ["chr22\t100\t200\tchr22:100-200:.\t0\t."]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -k end_to_end -v`
Expected: FAIL — `build` undefined.

- [ ] **Step 3: Write minimal implementation**

```python
def build(
    gct_path: str,
    attributes_path: str,
    bed_path: str,
    qc_path: str,
    min_read_count: int = 1,
    restrict_chrom: "str | None" = None,
) -> None:
    """Full build: attributes -> GCT stream -> union -> BED6 + QC."""
    sample_tissue = parse_sample_attributes(attributes_path)
    samples, rows = iter_gct_rows(gct_path)
    union = accumulate_pan_tissue_union(samples, rows, sample_tissue, min_read_count)
    if restrict_chrom is not None:
        union = {k: v for k, v in union.items() if k[0] == restrict_chrom}
        log.info("Restricted to %s: %d junctions", restrict_chrom, len(union))
    write_bed6(union, bed_path)
    write_qc_sidecar(union, qc_path)


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Build GTEx V10 pan-tissue junction blacklist BED."
    )
    parser.add_argument("--junctions-gct", required=True, help="GTEx junctions .gct[.gz]")
    parser.add_argument("--sample-attributes", required=True, help="SampleAttributesDS .txt")
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-qc", required=True)
    parser.add_argument("--min-read-count", type=int, default=1)
    parser.add_argument("--restrict-chrom", default=None,
                        help="Emit only this chromosome (for the chr22 test fixture)")
    args = parser.parse_args()
    build(
        gct_path=args.junctions_gct,
        attributes_path=args.sample_attributes,
        bed_path=args.output_bed,
        qc_path=args.output_qc,
        min_read_count=args.min_read_count,
        restrict_chrom=args.restrict_chrom,
    )


if __name__ == "__main__":
    _cli_main()
```

- [ ] **Step 4: Run the full test file**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_build_gtex_pan_tissue_ref.py -v`
Expected: PASS (all tests across Tasks 2–7)

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/build_gtex_pan_tissue_ref.py workflow/tests/test_build_gtex_pan_tissue_ref.py
git commit -m "feat(gtex): build() + CLI with --restrict-chrom (issue #211)"
```

---

### Task 8: Config schema + runbook

**Files:**
- Modify: `config/config.yaml`
- Modify: `docs/gtex_pan_tissue_build.md`

- [ ] **Step 1: Add the `gtex_filter` block to config**

Add under the existing `reference:` neighbourhood in `config/config.yaml` (consumed by #212; the `reference_bed` path is the GCS staging target from Task 9):
```yaml
gtex_filter:
  enabled: false          # opt-in here in #211; #212 flips the default to true at integration
  release: v10
  min_read_count: 1       # most aggressive — precision over recall (vaccine safety, per #126)
  reference_bed: "gs://splice-neoepitope-project/resources/gtex/v10/gtex_v10_pan_tissue_junctions.bed"
  reference_bed_chr22: "resources/test/gtex_v10_pan_tissue_junctions.chr22.bed"
```

- [ ] **Step 2: Verify config still parses**

Run: `conda activate snakemake && python -c "import yaml; yaml.safe_load(open('config/config.yaml')); print('config OK')"`
Expected: `config OK`

- [ ] **Step 3: Write the build runbook**

Append to `docs/gtex_pan_tissue_build.md` a `## Build command` section with the exact production invocation (using the URLs pinned in Task 1) and the GCS staging commands from Task 9. State explicitly that this is a manual one-shot — not in the per-patient DAG — and how to refresh on a new GTEx release.

- [ ] **Step 4: Commit**

```bash
git add config/config.yaml docs/gtex_pan_tissue_build.md
git commit -m "feat(gtex): config block + build runbook (issue #211)"
```

---

### Task 9: Production build + GCS staging + chr22 fixture

**Files:**
- Produce + stage artifacts (not committed if >100 MB; chr22 fixture committed)

This task runs the real build against full V10 data. It is heavy (multi-GB stream) — run on a small VM or a machine with bandwidth, not necessarily the M1.

- [ ] **Step 1: Run the full pan-tissue build**

```bash
conda activate snakemake
# URLs from Task 1 findings doc:
curl -s '<JUNCTIONS_GCT_URL>' -o /tmp/gtex_v10_junctions.gct.gz
curl -s '<SAMPLE_ATTRIBUTES_URL>' -o /tmp/gtex_v10_sample_attrs.txt
python workflow/scripts/build_gtex_pan_tissue_ref.py \
  --junctions-gct /tmp/gtex_v10_junctions.gct.gz \
  --sample-attributes /tmp/gtex_v10_sample_attrs.txt \
  --output-bed /tmp/gtex_v10_pan_tissue_junctions.bed \
  --output-qc /tmp/gtex_v10_pan_tissue_junctions.qc.tsv \
  --min-read-count 1
```
Expected: logs reporting samples/tissues parsed and `N/M junctions kept`. Sanity-check `head /tmp/gtex_v10_pan_tissue_junctions.bed` shows `chr*` BED6 rows; `cat /tmp/...qc.tsv` shows a plausible total + N-tissues histogram.

- [ ] **Step 2: Build the chr22 fixture**

```bash
python workflow/scripts/build_gtex_pan_tissue_ref.py \
  --junctions-gct /tmp/gtex_v10_junctions.gct.gz \
  --sample-attributes /tmp/gtex_v10_sample_attrs.txt \
  --output-bed resources/test/gtex_v10_pan_tissue_junctions.chr22.bed \
  --output-qc /tmp/gtex_v10_pan_tissue_junctions.chr22.qc.tsv \
  --min-read-count 1 --restrict-chrom chr22
```
Expected: a small BED restricted to chr22. Confirm size is commit-friendly (`wc -l` + `du -h`).

- [ ] **Step 3: Coordinate-space cross-check (the #148/#370 silent-empty guard)**

Confirm GTEx coords land in the same space as `reference_junctions.bed` by checking overlap is non-trivial on chr22:
```bash
conda activate snakemake
python workflow/scripts/build_reference_junctions.py \
  --gtf <(zcat references/gencode.v47.annotation.gtf.gz | awk '$1=="chr22"') \
  --output /tmp/ref_chr22.bed 2>/dev/null || \
  python workflow/scripts/build_reference_junctions.py --gtf references/gencode.v47.annotation.gtf.gz --output /tmp/ref_full.bed
# Expect a substantial intersection on (chrom,start,end) — NOT zero:
cut -f1-3 resources/test/gtex_v10_pan_tissue_junctions.chr22.bed | sort -u > /tmp/gtex_keys.txt
grep -P '^chr22\t' /tmp/ref_*.bed | cut -f1-3 | sort -u > /tmp/ref_keys.txt
echo "shared chr22 junctions: $(comm -12 /tmp/gtex_keys.txt /tmp/ref_keys.txt | wc -l)"
```
Expected: a **non-zero, substantial** shared count (most annotated chr22 junctions are expressed somewhere in GTEx). **Zero means the coordinate transform in `parse_junction_name` is wrong — STOP and re-pin Task 1 Step 3.**

- [ ] **Step 4: Stage to GCS + commit the chr22 fixture**

```bash
gsutil cp /tmp/gtex_v10_pan_tissue_junctions.bed     gs://splice-neoepitope-project/resources/gtex/v10/
gsutil cp /tmp/gtex_v10_pan_tissue_junctions.qc.tsv  gs://splice-neoepitope-project/resources/gtex/v10/
git add resources/test/gtex_v10_pan_tissue_junctions.chr22.bed
git commit -m "feat(gtex): stage V10 pan-tissue BED to GCS + commit chr22 fixture (issue #211)"
```

---

### Task 10: Re-run Issue #225 notebook §2(c) — NO-GO verdict guard (AC)

**Files:**
- Run (no code change expected): `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` §2(c)

- [ ] **Step 1: Re-run §2(c) against the production V10 panel from Task 9**

Open the notebook in the research venv (`research/.venv`), point the GTEx-panel construction cell at the staged `gtex_v10_pan_tissue_junctions.bed` (or its chr22 fixture for a fast check), and re-run §2(c).

- [ ] **Step 2: Confirm the verdict holds**

Expected: the AlphaGenome-as-3rd-filter **NO-GO** verdict still holds (the production panel is at least as inclusive as the Snaptron chr22 proxy at F1-max τ). Record the before/after numbers.

- [ ] **Step 3: Document in lab notebook + commit (post-review, per lab-notebook timing rule)**

Capture the §2(c) re-run conclusion in `research/lab_notebook/developer.md` AFTER review feedback is incorporated (lab-notebook-entry-after-review rule), before merge.

---

## Acceptance-criteria coverage (self-review against Issue #211)

- [x] GTEx V10 release version pinned in config → Task 8 (`gtex_filter.release: v10`)
- [x] All tissues' junction read-counts ingested reproducibly, clean exit on missing/changed file → Tasks 4+7 (single GCT stream; CLI errors on missing path) + Task 1 (pinned URLs in runbook)
- [x] Pan-tissue union with `min_read_count=1` default → Task 5 + Task 7 CLI default
- [x] BED uses GRCh38 `chr*` naming consistent with pipeline → Task 3 (`parse_junction_name` note) + Task 9 Step 3 cross-check
- [x] BED + QC TSV staged to GCS at documented path → Task 9 Step 4
- [x] chr22 test fixture produced → Task 7 (`--restrict-chrom`) + Task 9 Step 2
- [x] Pytest covers parse + union construction on synthetic input → Tasks 2–7 tests
- [x] Re-run #225 notebook §2(c), confirm NO-GO holds → Task 10

## Cross-issue note for [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) (integration)

The GTEx blacklist BED uses strand `.` (positional/strand-agnostic). When #212 wires it into `filter_junctions.py`, it MUST intersect tumour junctions on `(chrom, start, end)` only — do NOT reuse `_load_reference_junctions`'s `(chrom,start,end,strand)` 4-tuple key for the GTEx source (a `+`/`-` tumour junction would never match a `.` blacklist row). Add a strand-agnostic loader/intersection path for `gtex_pantissue_shared`. Tag filtered junctions `origin = gtex_pantissue_shared`.

## Out of scope (this plan)
- `filter_junctions.py` integration / stacking with matched-normal → #212
- patient_002 / patient_001 validation runs → #212
- AlphaGenome predicted-normal axis → #203
