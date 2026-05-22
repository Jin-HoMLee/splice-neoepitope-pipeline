# Issue #204 VDJdb Panel + stitchr Setup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204): a new Snakemake rule `fetch_vdjdb_panel` that produces a per-patient reference TCR panel (top-10 paired α/β TCRs per HLA allele) sourced from VDJdb 2026-05-16 and reconstructed via stitchr.

**Architecture:** New conda env `workflow/envs/vdjdb.yaml` (stitchr + IMGTgeneDL) for local-runnable, non-GPU panel construction. Two new download rules in `download.smk` (VDJdb release with SHA256 verification + IMGT germlines via `stitchrdl`) feed a new rule `fetch_vdjdb_panel` in NEW `workflow/rules/tcr_panel.smk`, which runs `workflow/scripts/fetch_vdjdb_panel.py` per patient. Outputs land in the forward-compatible `results/{p}/tcr_panel/vdjdb/panel.tsv` + `panel_qc.tsv` layout (room for future TRUST4/ProTCR sources from [Issue #24](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/24)).

**Tech Stack:** Snakemake 8, Python 3.11, pandas, stitchr 1.x, IMGTgeneDL 0.x, pytest. Conda env-driven; no Docker/GPU.

**Reference:** Design spec at [docs/superpowers/specs/2026-05-20-vdjdb-panel-design.md](../specs/2026-05-20-vdjdb-panel-design.md).

**Branch already set up:** `204-feattcrdock-fetch_vdjdb_panel-snakemake-rule-+-stitchrvdjdb-setup` (checked out), spec committed at `b978897`. Issue #204 Status: In Progress on the project board.

---

## File Structure

**Files to create:**
- `workflow/envs/vdjdb.yaml` — conda env: stitchr, IMGTgeneDL, pandas, requests
- `workflow/rules/tcr_panel.smk` — new rule module; one rule `fetch_vdjdb_panel`
- `workflow/scripts/fetch_vdjdb_panel.py` — implementation; lazy-imports pandas/stitchr
- `workflow/tests/test_fetch_vdjdb_panel.py` — unit tests
- `workflow/tests/fixtures/vdjdb_mini.tsv` — synthetic ~10-row VDJdb fixture (35 columns)

**Files to modify:**
- `Snakefile` — add `include: "workflow/rules/tcr_panel.smk"` line
- `workflow/rules/download.smk` — add `download_vdjdb_release` + `download_imgt_germlines` rules
- `config/config.yaml` — extend `tcrdock:` section with `vdjdb_release`, `vdjdb_sha256`, `vdjdb_min_score`, `vdjdb_panel_size`
- `research/glossary.md` — add new `## K` section with KIR entry

**Files NOT touched in this PR (per spec out-of-scope):**
- `workflow/scripts/generate_report.py` — sub-issue [#206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206)
- `workflow/rules/structure.smk` — sub-issue [#205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205) wires selection in
- `results/{p}/predictions/` flatten — [Issue #435](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/435)

---

## Task 1: Add KIR glossary entry

**Files:**
- Modify: `research/glossary.md` (insert new `## K` section between `## I` and `## L`)

- [ ] **Step 1: Insert KIR entry**

Open `research/glossary.md`. Find the section divider between `## I` (around line 89) and `## L` (around line 97). Insert this block immediately before the `## L` heading:

```markdown
## K

**KIR** — Killer-cell Immunoglobulin-like Receptor. Family of NK-cell receptors whose dominant ligands are HLA-C alleles; HLA-C–KIR engagement inhibits NK cytotoxicity under "missing-self" surveillance. Relevant here because HLA-C is studied primarily through this NK biology rather than T-cell presentation, which explains why TCR repertoire databases (VDJdb, McPAS-TCR) have thin paired α/β coverage for C-locus alleles. See [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) empirical coverage check. *Domain: bio.*

```

- [ ] **Step 2: Verify section ordering**

Run:
```bash
grep -n '^## ' research/glossary.md
```
Expected: `## A → ## C → ## D → ## E → ## G → ## H → ## I → ## K → ## L → ## M → ## N → ## O → ## P → ## R → ## S → ## T → ## U → ## W` (K now between I and L).

- [ ] **Step 3: Commit**

```bash
git add research/glossary.md
git commit -m "docs(glossary): add KIR entry — motivates HLA-C thin coverage in VDJdb (Issue #204)"
```

---

## Task 2: Add VDJdb config keys

**Files:**
- Modify: `config/config.yaml` (extend existing `tcrdock:` block)

- [ ] **Step 1: Read current tcrdock: block to find insertion point**

Run:
```bash
grep -n '^tcrdock:' config/config.yaml
```
Expected: one line, around line ~120-150 of `config/config.yaml`.

- [ ] **Step 2: Add four new keys under the existing `tcrdock:` block**

The current block contains `enabled: false`. Add these keys immediately after `enabled: false`, preserving 2-space indentation and the existing comment block above the section:

```yaml
tcrdock:
  enabled: false
  # VDJdb panel fetch (Issue #204) — runs whenever hla.enabled is true,
  # independent of the tcrdock.enabled GPU gate.
  vdjdb_release: "2026-05-16"
  vdjdb_sha256: "0dce79ec55c109000da10b7bc72300e352ffb7df92d5de30682de20bab35a366"
  vdjdb_min_score: 2
  vdjdb_panel_size: 10
```

- [ ] **Step 3: Verify yaml is parseable**

Run:
```bash
python -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['tcrdock'])"
```
Expected: a dict with keys `enabled, vdjdb_release, vdjdb_sha256, vdjdb_min_score, vdjdb_panel_size`. No exceptions.

- [ ] **Step 4: Commit**

```bash
git add config/config.yaml
git commit -m "config(tcrdock): add VDJdb panel keys (release pin, SHA256, score threshold, panel size) (Issue #204)"
```

---

## Task 3: Create vdjdb.yaml conda env

**Files:**
- Create: `workflow/envs/vdjdb.yaml`

- [ ] **Step 1: Create the env file**

Write `workflow/envs/vdjdb.yaml` with this content:

```yaml
name: vdjdb
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - pandas
  - requests
  - pip
  - pip:
    - stitchr>=1.1.0
    - IMGTgeneDL>=0.7.0
```

- [ ] **Step 2: Verify the env solves (dry-run)**

Run:
```bash
conda env create -n vdjdb-probe -f workflow/envs/vdjdb.yaml --dry-run 2>&1 | tail -20
```
Expected: solver completes without conflicts; final line indicates "Dry run" success. If solver fails on pip packages, that's fine for dry-run — `pip:` resolves at install time, not solve time.

If solver fails on the conda side, capture the error and surface to user before continuing.

- [ ] **Step 3: Commit**

```bash
git add workflow/envs/vdjdb.yaml
git commit -m "env(vdjdb): add conda env for stitchr + IMGTgeneDL (Issue #204)"
```

---

## Task 4: Create synthetic VDJdb fixture

**Files:**
- Create: `workflow/tests/fixtures/vdjdb_mini.tsv`

- [ ] **Step 1: Confirm fixtures dir doesn't exist yet and create it**

Run:
```bash
mkdir -p workflow/tests/fixtures
```

- [ ] **Step 2: Write the fixture**

Write `workflow/tests/fixtures/vdjdb_mini.tsv` with 10 data rows + header. The schema matches the real `vdjdb_full.txt` (35 columns from VDJdb 2026-05-16, paired-α/β format). For brevity, only the columns we filter on (`species`, `mhc.a`, `mhc.class`, `vdjdb.score`) and the columns we emit (`cdr3.alpha`, `v.alpha`, `j.alpha`, `cdr3.beta`, `v.beta`, `j.beta`, `meta.subject.id`) carry meaningful values. Other columns are blank or sentinel.

Tabs are required between fields. Use this exact content:

```tsv
cdr3.alpha	v.alpha	j.alpha	cdr3.beta	v.beta	d.beta	j.beta	species	mhc.a	mhc.b	mhc.class	antigen.epitope	antigen.gene	antigen.species	reference.id	method.identification	method.frequency	method.singlecell	method.sequencing	method.verification	meta.study.id	meta.cell.subset	meta.subject.cohort	meta.subject.id	meta.replica.id	meta.clone.id	meta.epitope.id	meta.tissue	meta.donor.MHC	meta.donor.MHC.method	meta.structure.id	cdr3fix.alpha	cdr3fix.beta	vdjdb.score	TCR_hash
CAVNFGGGKLI	TRAV12-2	TRAJ21	CASSLAGGRPEQYF	TRBV6-5	TRBD1	TRBJ2-7	HomoSapiens	HLA-A*02:01	B2M	MHCI	NLVPMVATV	pp65	CMV	PMID:11111	ID-1	freq-1	sc-1	seq-1	verif-1	S1	CD8	cohort-A	donor-001	r1	c1	e1	blood	A02:01	pcr	1aa1	fixA	fixB	3	hash-001
CAVNFGGGKLI	TRAV12-2	TRAJ21	CASSLAGGRPEQYF	TRBV6-5	TRBD1	TRBJ2-7	HomoSapiens	HLA-A*02:01:110	B2M	MHCI	NLVPMVATV	pp65	CMV	PMID:11112	ID-2	freq-2	sc-2	seq-2	verif-2	S2	CD8	cohort-A	donor-002	r2	c2	e2	blood	A02:01:110	pcr	1bb1	fixA	fixB	3	hash-002
CASSDR	TRAV1-1	TRAJ34	CASSDR	TRBV1	TRBD1	TRBJ2-3	HomoSapiens	HLA-A*02:01	B2M	MHCI	XXXX	geneX	speciesY	PMID:22221	ID-3	freq-3	sc-3	seq-3	verif-3	S3	CD8	cohort-B	donor-003	r3	c3	e3	blood	A02:01	pcr	2aa	fixA	fixB	2	hash-003
CASSDR	TRAV1-1	TRAJ34	CASSDR	TRBV1	TRBD1	TRBJ2-3	HomoSapiens	HLA-A*02:01	B2M	MHCI	YYYY	geneY	speciesZ	PMID:22222	ID-4	freq-4	sc-4	seq-4	verif-4	S4	CD8	cohort-B	donor-004	r4	c4	e4	blood	A02:01	pcr	2bb	fixA	fixB	0	hash-004
CASSDR	TRAV1-1	TRAJ34	CASSDR	TRBV1	TRBD1	TRBJ2-3	HomoSapiens	HLA-A*02:02	B2M	MHCI	ZZZZ	geneZ	speciesW	PMID:22223	ID-5	freq-5	sc-5	seq-5	verif-5	S5	CD8	cohort-C	donor-005	r5	c5	e5	blood	A02:02	pcr	2cc	fixA	fixB	3	hash-005
CASSDR	TRAV2	TRAJ35	CASSDR	TRBV2	TRBD1	TRBJ2-4	HomoSapiens	HLA-B*08:01	B2M	MHCI	FLY	geneM	speciesN	PMID:33331	ID-6	freq-6	sc-6	seq-6	verif-6	S6	CD8	cohort-D	donor-006	r6	c6	e6	blood	B08:01	pcr	3aa	fixA	fixB	2	hash-006
CASSDR	TRAV3	TRAJ36	CASSDR	TRBV3	TRBD1	TRBJ2-5	MusMusculus	H2-Kb	B2M	MHCI	SIINFEKL	OVA	chicken	PMID:44441	ID-7	freq-7	sc-7	seq-7	verif-7	S7	CD8	cohort-E	donor-007	r7	c7	e7	blood	H2-Kb	pcr	4aa	fixA	fixB	3	hash-007
CASSDR	TRAV4	TRAJ37	CASSDR	TRBV4	TRBD1	TRBJ2-6	HomoSapiens	HLA-DRB1*01:01	HLA-DRA*01:01	MHCII	PEPTIDE2	geneP	speciesQ	PMID:55551	ID-8	freq-8	sc-8	seq-8	verif-8	S8	CD4	cohort-F	donor-008	r8	c8	e8	blood	DRB1*01:01	pcr	5aa	fixA	fixB	3	hash-008
CASSDR	TRAV5	TRAJ38	CASSDR	TRBV5	TRBD1	TRBJ2-1	HomoSapiens	HLA-A*02:01	B2M	MHCI	NLVPMVATV	pp65	CMV	PMID:66661	ID-9	freq-9	sc-9	seq-9	verif-9	S9	CD8	cohort-G	donor-009	r9	c9	e9	blood	A02:01	pcr	6aa	fixA	fixB	3	hash-009
CASSDR	TRAV6	TRAJ39	CASSDR	TRBV6	TRBD1	TRBJ2-2	HomoSapiens	HLA-A*02:01	B2M	MHCI	OTHER	geneO	speciesR	PMID:77771	ID-10	freq-10	sc-10	seq-10	verif-10	S10	CD8	cohort-H	donor-010	r10	c10	e10	blood	A02:01	pcr	7aa	fixA	fixB	3	hash-010
```

Row inventory (used by tests later):
- Rows 1, 9, 10: `HomoSapiens`, `MHCI`, `HLA-A*02:01`, score 3 (passes filters; 3 entries on A*02:01)
- Row 2: `HomoSapiens`, `MHCI`, `HLA-A*02:01:110` (6-digit, normalizes to A*02:01), score 3 → also matches A*02:01 (4 total)
- Rows 3: `HomoSapiens`, `MHCI`, `HLA-A*02:01`, score 2 (passes) (5 total)
- Row 4: `HomoSapiens`, `MHCI`, `HLA-A*02:01`, score 0 (FAILS score filter)
- Row 5: `HomoSapiens`, `MHCI`, `HLA-A*02:02`, score 3 (passes; different allele)
- Row 6: `HomoSapiens`, `MHCI`, `HLA-B*08:01`, score 2 (passes; B locus)
- Row 7: `MusMusculus`, `MHCI` (FAILS species filter)
- Row 8: `HomoSapiens`, `MHCII` (FAILS class filter)

Expected `HLA-A*02:01` exact-match-after-4-digit-normalize with score≥2: 5 rows. Useful for `test_top_n_deterministic_tiebreak`, `test_exact_match_per_allele`, `test_filter_*`.

- [ ] **Step 3: Verify the fixture loads cleanly**

Run:
```bash
python -c "import pandas as pd; df = pd.read_csv('workflow/tests/fixtures/vdjdb_mini.tsv', sep='\t'); print(df.shape); print(df.columns.tolist())"
```
Expected: `(10, 35)`, full column list matching the VDJdb 2026-05-16 schema.

- [ ] **Step 4: Commit**

```bash
git add workflow/tests/fixtures/vdjdb_mini.tsv
git commit -m "test(fixture): add synthetic vdjdb_mini.tsv for fetch_vdjdb_panel tests (Issue #204)"
```

---

## Task 5: TDD — normalize_allele_to_4digit

**Files:**
- Create: `workflow/scripts/fetch_vdjdb_panel.py` (module skeleton + first function)
- Create: `workflow/tests/test_fetch_vdjdb_panel.py` (first test)

- [ ] **Step 1: Write the failing test**

Create `workflow/tests/test_fetch_vdjdb_panel.py` with:

```python
"""Tests for fetch_vdjdb_panel.py — VDJdb panel construction."""

import pytest

from fetch_vdjdb_panel import normalize_allele_to_4digit


class TestNormalizeAlleleTo4Digit:
    @pytest.mark.parametrize("inp,expected", [
        ("HLA-A*02:01", "HLA-A*02:01"),       # already 4-digit
        ("HLA-A*02:01:110", "HLA-A*02:01"),   # 6-digit truncates
        ("HLA-B*15:63", "HLA-B*15:63"),       # rare allele, already 4-digit
        ("HLA-C*07:01:01:03", "HLA-C*07:01"), # 8-digit truncates
    ])
    def test_4digit_or_longer_normalizes(self, inp, expected):
        assert normalize_allele_to_4digit(inp) == expected

    @pytest.mark.parametrize("inp", [
        "HLA-A*02",      # 2-digit only — excluded
        "HLA-A",         # no allele subtype at all
        "",              # empty
    ])
    def test_2digit_or_invalid_returns_none(self, inp):
        assert normalize_allele_to_4digit(inp) is None
```

- [ ] **Step 2: Run the test to confirm it fails**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v
```
Expected: collection failure — `ImportError: cannot import name 'normalize_allele_to_4digit' from 'fetch_vdjdb_panel'` (or module not found).

- [ ] **Step 3: Create the script skeleton + function**

Create `workflow/scripts/fetch_vdjdb_panel.py`:

```python
"""Fetch a per-patient VDJdb TCR panel (Issue #204).

Reads a pinned VDJdb release, filters to HomoSapiens + MHCI + paired α/β +
vdjdb.score >= threshold, exact-matches by 4-digit HLA allele to a patient's
alleles.tsv, picks the top-N TCRs per allele (ranked by score, deterministic
tiebreak by donor ID), reconstructs full α/β chains via stitchr, and writes:

    results/{patient_id}/tcr_panel/vdjdb/panel.tsv    — TCRs with full sequences
    results/{patient_id}/tcr_panel/vdjdb/panel_qc.tsv — per-allele coverage status

Lazy-imports pandas + stitchr at first use (per PR #428 pattern) so pytest
collection stays fast.
"""

from __future__ import annotations

import logging
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Allele normalization
# ---------------------------------------------------------------------------

def normalize_allele_to_4digit(allele: str) -> Optional[str]:
    """Truncate an HLA allele string to its 4-digit form (`HLA-X*GG:PP`).

    Anything shorter than 4-digit returns None — we require exact 4-digit
    matching, and shorter forms (e.g. `HLA-A*02`) are excluded.

    Examples:
        HLA-A*02:01      -> HLA-A*02:01
        HLA-A*02:01:110  -> HLA-A*02:01
        HLA-A*02         -> None
        ""               -> None
    """
    if not allele or ":" not in allele:
        return None
    parts = allele.split(":")
    if len(parts) < 2:
        return None
    return f"{parts[0]}:{parts[1]}"
```

- [ ] **Step 4: Run the test to confirm it passes**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestNormalizeAlleleTo4Digit -v
```
Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add normalize_allele_to_4digit + script skeleton (Issue #204)

Test + impl committed together per TDD red-green cycle."
```

---

## Task 6: TDD — load_and_filter_vdjdb

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py` (add `load_and_filter_vdjdb`)
- Modify: `workflow/tests/test_fetch_vdjdb_panel.py` (add tests)

- [ ] **Step 1: Add failing tests**

Append to `workflow/tests/test_fetch_vdjdb_panel.py`:

```python
from pathlib import Path

from fetch_vdjdb_panel import load_and_filter_vdjdb


FIXTURE_PATH = Path(__file__).parent / "fixtures" / "vdjdb_mini.tsv"


class TestLoadAndFilterVdjdb:
    def test_filters_to_homosapiens_mhci_minscore(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        # Expected pass rows: 1, 2, 3, 5, 6, 9, 10 — i.e. 7 rows
        # (row 4 fails score=0; row 7 fails species=Mouse; row 8 fails MHCII)
        assert len(df) == 7
        assert (df["species"] == "HomoSapiens").all()
        assert (df["mhc.class"] == "MHCI").all()
        assert (df["vdjdb.score"] >= 2).all()

    def test_min_score_threshold_is_inclusive(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=3)
        # Drops the score=2 rows (rows 3 and 6 in fixture)
        assert (df["vdjdb.score"] >= 3).all()
        assert len(df) == 5

    def test_normalizes_mhc_a_inplace(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        # Row 2's HLA-A*02:01:110 must be normalized to HLA-A*02:01
        a0201_rows = df[df["mhc.a_4digit"] == "HLA-A*02:01"]
        assert len(a0201_rows) == 5  # rows 1, 2, 3, 9, 10

    def test_drops_rows_with_unnormalizable_allele(self, tmp_path):
        # Synthetic fixture with a 2-digit-only allele — should be dropped
        tsv = tmp_path / "vdjdb_partial.tsv"
        header = open(FIXTURE_PATH).readline()
        # 35 cols total: 7 empty (1-7) + 4 set (8-11: species, mhc.a, mhc.b, mhc.class)
        # + 22 empty (12-33) + vdjdb.score (34) + TCR_hash (35).
        bad_row = "\t".join([""] * 7 + ["HomoSapiens", "HLA-A*02", "B2M", "MHCI"] + [""] * 22 + ["3", ""]) + "\n"
        tsv.write_text(header + bad_row)
        df = load_and_filter_vdjdb(tsv, min_score=2)
        assert len(df) == 0
```

- [ ] **Step 2: Run the new tests to confirm they fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestLoadAndFilterVdjdb -v
```
Expected: `ImportError` on `load_and_filter_vdjdb`.

- [ ] **Step 3: Add the implementation**

Append to `workflow/scripts/fetch_vdjdb_panel.py`:

```python
# ---------------------------------------------------------------------------
# VDJdb load + filter
# ---------------------------------------------------------------------------

def load_and_filter_vdjdb(vdjdb_full_tsv, min_score: int):
    """Load vdjdb_full.txt and apply HomoSapiens + MHCI + score filters.

    Adds a `mhc.a_4digit` column with 4-digit-normalized allele.
    Drops rows where the allele cannot be normalized to 4-digit.
    Returns a pandas DataFrame.
    """
    import pandas as pd  # lazy import

    df = pd.read_csv(vdjdb_full_tsv, sep="\t", dtype=str)
    df["vdjdb.score"] = pd.to_numeric(df["vdjdb.score"], errors="coerce")
    df = df[
        (df["species"] == "HomoSapiens")
        & (df["mhc.class"] == "MHCI")
        & (df["vdjdb.score"] >= min_score)
    ].copy()
    df["mhc.a_4digit"] = df["mhc.a"].apply(normalize_allele_to_4digit)
    df = df[df["mhc.a_4digit"].notna()].copy()
    log.info(
        "Loaded VDJdb (%s): %d rows after filters (HomoSapiens + MHCI + score>=%d)",
        vdjdb_full_tsv, len(df), min_score,
    )
    return df
```

- [ ] **Step 4: Run all tests so far**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v
```
Expected: all green (the 7 normalize tests + 4 new load_and_filter tests).

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add load_and_filter_vdjdb (HomoSapiens + MHCI + score) (Issue #204)"
```

---

## Task 7: TDD — select_top_n_for_allele (panel ranking with deterministic tiebreak)

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py`
- Modify: `workflow/tests/test_fetch_vdjdb_panel.py`

- [ ] **Step 1: Add failing tests**

Append to `workflow/tests/test_fetch_vdjdb_panel.py`:

```python
from fetch_vdjdb_panel import select_top_n_for_allele


class TestSelectTopNForAllele:
    def test_exact_match_filter(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=10)
        # rows 1, 2, 3, 9, 10 — 5 total (includes the 6-digit row 2 after normalization)
        assert len(result) == 5
        assert (result["mhc.a_4digit"] == "HLA-A*02:01").all()

    def test_does_not_match_different_allele(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:02", n=10)
        assert len(result) == 1  # only row 5 in fixture

    def test_returns_top_n_when_more_available(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=2)
        assert len(result) == 2

    def test_returns_fewer_when_fewer_available(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-B*08:01", n=10)
        assert len(result) == 1  # only row 6 matches B*08:01

    def test_returns_empty_when_no_matches(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*31:01", n=10)
        assert len(result) == 0

    def test_sorted_by_score_then_donor_id(self):
        df = load_and_filter_vdjdb(FIXTURE_PATH, min_score=2)
        result = select_top_n_for_allele(df, allele="HLA-A*02:01", n=10)
        # Of the 5 matching rows: 4 with score=3 (donors 001, 002, 009, 010), 1 with score=2 (donor 003).
        # Score 3 rows come first; among them, donor IDs sorted ascending lexicographically.
        scores = result["vdjdb.score"].tolist()
        assert scores == [3, 3, 3, 3, 2]
        donor_order = result["meta.subject.id"].tolist()
        # Top 4 (score=3) sorted ascending: 001, 002, 009, 010. Last is score=2 (donor 003).
        assert donor_order == ["donor-001", "donor-002", "donor-009", "donor-010", "donor-003"]
```

- [ ] **Step 2: Run the new tests to confirm they fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestSelectTopNForAllele -v
```
Expected: `ImportError` on `select_top_n_for_allele`.

- [ ] **Step 3: Add the implementation**

Append to `workflow/scripts/fetch_vdjdb_panel.py`:

```python
# ---------------------------------------------------------------------------
# Per-allele top-N selection
# ---------------------------------------------------------------------------

def select_top_n_for_allele(df, allele: str, n: int):
    """Filter `df` to exact 4-digit matches for `allele`, sort by
    (vdjdb.score DESC, meta.subject.id ASC) for deterministic tiebreak,
    return top `n` rows.
    """
    matched = df[df["mhc.a_4digit"] == allele].copy()
    matched = matched.sort_values(
        ["vdjdb.score", "meta.subject.id"],
        ascending=[False, True],
        kind="mergesort",  # stable
    )
    return matched.head(n).reset_index(drop=True)
```

- [ ] **Step 4: Run all tests so far**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v
```
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add select_top_n_for_allele with deterministic tiebreak (Issue #204)"
```

---

## Task 8: TDD — classify_panel_status (ok / low_coverage / empty)

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py`
- Modify: `workflow/tests/test_fetch_vdjdb_panel.py`

- [ ] **Step 1: Add failing tests**

Append:

```python
from fetch_vdjdb_panel import classify_panel_status


class TestClassifyPanelStatus:
    @pytest.mark.parametrize("n_in_panel,target,expected", [
        (10, 10, "ok"),
        (15, 10, "ok"),   # clipped to top-N before reaching here; defensive case
        (9, 10, "low_coverage"),
        (1, 10, "low_coverage"),
        (0, 10, "empty"),
    ])
    def test_classification(self, n_in_panel, target, expected):
        assert classify_panel_status(n_in_panel, target_size=target) == expected
```

- [ ] **Step 2: Run the new tests to confirm they fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestClassifyPanelStatus -v
```
Expected: `ImportError` on `classify_panel_status`.

- [ ] **Step 3: Add the implementation**

Append:

```python
# ---------------------------------------------------------------------------
# Panel status classification
# ---------------------------------------------------------------------------

def classify_panel_status(n_in_panel: int, target_size: int) -> str:
    """Return 'ok' if panel reached target_size, 'low_coverage' if 1..target_size-1,
    'empty' if 0.
    """
    if n_in_panel <= 0:
        return "empty"
    if n_in_panel < target_size:
        return "low_coverage"
    return "ok"
```

- [ ] **Step 4: Run all tests so far**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v
```
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add classify_panel_status (ok/low_coverage/empty) (Issue #204)"
```

---

## Task 9: TDD — stitch_chain (stitchr smoke test with real IMGT)

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py`
- Modify: `workflow/tests/test_fetch_vdjdb_panel.py`

**Note:** This task requires `stitchr` + IMGT germline data installed locally. Set up before running tests:

```bash
workflow/tests/.venv/bin/pip install stitchr IMGTgeneDL
workflow/tests/.venv/bin/stitchrdl --species HUMAN  # one-time IMGT download
```

If `workflow/tests/.venv` doesn't have stitchr (the conda env in Task 3 has it, but the test venv may not), install with the commands above. The test is marked `@pytest.mark.network` so CI can skip it if IMGT is unavailable, and a `skipif shutil.which("stitchr") is None` gate ensures local runs without stitchr silently skip.

- [ ] **Step 1: Add failing test**

Append:

```python
import shutil

from fetch_vdjdb_panel import stitch_chain


@pytest.mark.network
@pytest.mark.skipif(shutil.which("stitchr") is None, reason="stitchr CLI not on PATH")
class TestStitchChain:
    """Smoke test against the real stitchr CLI + cached IMGT data.

    Skipped automatically when the `stitchr` CLI is not on PATH (e.g. CI runners
    without the vdjdb conda env activated). The `network` marker is informational
    — actual skip is driven by the `skipif` so CI works regardless of `-m` flag.
    """
    def test_dmf5_alpha_chain(self):
        # DMF5 (the existing single-fallback TCR in config/gpu_config.yaml)
        alpha = stitch_chain(
            v_gene="TRAV12-2",
            j_gene="TRAJ21",
            cdr3="CAVNFGGGKLI",
            chain="A",
        )
        # Stitched chain should be a non-trivial protein sequence containing the CDR3
        assert alpha is not None
        assert "CAVNFGGGKLI" in alpha
        # TRAV12-2 framework begins with a known leader; check for FW1 motif "QSV"
        # (allowing some flexibility — sanity check, not equality)
        assert len(alpha) > 100  # full Vα is ~100-110 aa
```

- [ ] **Step 2: Run the test to confirm it fails**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestStitchChain -v -m network
```
Expected: `ImportError` on `stitch_chain`.

- [ ] **Step 3: Add the implementation**

Append:

```python
# ---------------------------------------------------------------------------
# stitchr wrapper
# ---------------------------------------------------------------------------

def stitch_chain(v_gene: str, j_gene: str, cdr3: str, chain: str) -> Optional[str]:
    """Reconstruct a full TCR chain (V + CDR3 + J + framework) via stitchr.

    Returns the AA sequence string, or None if stitchr fails for any reason
    (caller is expected to log and skip).

    `chain` is 'A' for alpha (TRA) or 'B' for beta (TRB).
    """
    try:
        from Stitchr import stitchrfunctions, stitchr  # lazy import
    except ImportError:
        log.error("stitchr not installed in this environment")
        return None

    try:
        # stitchr's Python API: stitchr.stitch(specific_args, chain_data, ...)
        # Practical interface: use stitchrfunctions.stitch_a_chain or similar.
        # If the Python API is unstable, fall back to subprocess: `stitchr -v TRAV12-2 -j TRAJ21 -cdr3 CAVNFGGGKLI -species HUMAN`.
        import subprocess
        import shlex
        chain_arg = "TRA" if chain == "A" else "TRB"
        cmd = [
            "stitchr",
            "-v", v_gene,
            "-j", j_gene,
            "-cdr3", cdr3,
            "-species", "HUMAN",
            "-m", "AA",  # amino acid output
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if proc.returncode != 0:
            log.error("stitchr failed (chain=%s, V=%s, J=%s, CDR3=%s): %s",
                      chain, v_gene, j_gene, cdr3, proc.stderr.strip())
            return None
        # stitchr stdout is a FASTA-like blob; extract the sequence
        seq_lines = [ln.strip() for ln in proc.stdout.splitlines() if ln and not ln.startswith(">")]
        return "".join(seq_lines) if seq_lines else None
    except Exception as exc:
        log.error("stitchr crashed (chain=%s, V=%s, J=%s, CDR3=%s): %s",
                  chain, v_gene, j_gene, cdr3, exc)
        return None
```

- [ ] **Step 4: Configure pytest to register the `network` marker**

If `workflow/tests/pytest.ini` or `pyproject.toml` doesn't already register the `network` marker, add to `pyproject.toml` (or create `workflow/tests/pytest.ini`):

```ini
[pytest]
markers =
    network: tests requiring external network or IMGT data
```

- [ ] **Step 5: Run the stitchr test**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestStitchChain -v -m network
```
Expected: green if stitchr + IMGT data installed; skipped otherwise.

If the stitchr CLI invocation doesn't work as written (different flag names, output format), debug interactively first with `stitchr --help`, adjust the `cmd` list, and re-test.

- [ ] **Step 6: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py pyproject.toml
git commit -m "feat(vdjdb): add stitch_chain (stitchr subprocess wrapper) (Issue #204)"
```

---

## Task 10: TDD — build_panel (orchestrate + write panel.tsv + panel_qc.tsv)

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py`
- Modify: `workflow/tests/test_fetch_vdjdb_panel.py`

This task wires the helpers together into a single function that writes the two output TSVs from an alleles list + VDJdb path. The test uses a mock stitchr (monkeypatch) since the real one is exercised in Task 9.

- [ ] **Step 1: Add failing tests**

Append:

```python
from fetch_vdjdb_panel import build_panel


class TestBuildPanel:
    def test_writes_panel_and_qc_with_correct_schema(self, tmp_path, monkeypatch):
        # Monkeypatch stitch_chain so we don't depend on real stitchr in this test
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}_{cdr3}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        build_panel(
            vdjdb_full_tsv=FIXTURE_PATH,
            alleles=["HLA-A*02:01", "HLA-A*31:01", "HLA-B*08:01"],
            output_panel=panel_tsv,
            output_qc=qc_tsv,
            min_score=2,
            panel_size=10,
        )

        import pandas as pd
        panel = pd.read_csv(panel_tsv, sep="\t")
        qc = pd.read_csv(qc_tsv, sep="\t")

        # Panel schema
        assert list(panel.columns) == [
            "allele", "va_gene", "ja_gene", "cdr3a",
            "vb_gene", "jb_gene", "cdr3b",
            "alpha_seq", "beta_seq",
            "vdjdb_score", "vdjdb_donor_id",
        ]
        # QC schema
        assert list(qc.columns) == [
            "allele", "n_exact_matches", "n_in_panel", "panel_status"
        ]

    def test_qc_status_classification(self, tmp_path, monkeypatch):
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        build_panel(
            vdjdb_full_tsv=FIXTURE_PATH,
            alleles=["HLA-A*02:01", "HLA-A*31:01", "HLA-B*08:01"],
            output_panel=panel_tsv,
            output_qc=qc_tsv,
            min_score=2,
            panel_size=10,
        )

        import pandas as pd
        qc = pd.read_csv(qc_tsv, sep="\t").set_index("allele")
        # HLA-A*02:01: 5 matches, panel_size=10 → low_coverage
        assert qc.loc["HLA-A*02:01", "n_in_panel"] == 5
        assert qc.loc["HLA-A*02:01", "panel_status"] == "low_coverage"
        # HLA-A*31:01: 0 matches in fixture → empty
        assert qc.loc["HLA-A*31:01", "n_in_panel"] == 0
        assert qc.loc["HLA-A*31:01", "panel_status"] == "empty"
        # HLA-B*08:01: 1 match → low_coverage
        assert qc.loc["HLA-B*08:01", "n_in_panel"] == 1
        assert qc.loc["HLA-B*08:01", "panel_status"] == "low_coverage"

    def test_empty_allele_emits_warning(self, tmp_path, monkeypatch, caplog):
        import fetch_vdjdb_panel as mod
        monkeypatch.setattr(mod, "stitch_chain",
                            lambda v_gene, j_gene, cdr3, chain: f"MOCK_{chain}")

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"

        with caplog.at_level(logging.WARNING):
            build_panel(
                vdjdb_full_tsv=FIXTURE_PATH,
                alleles=["HLA-A*31:01"],
                output_panel=panel_tsv,
                output_qc=qc_tsv,
                min_score=2,
                panel_size=10,
            )

        warning_msgs = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("HLA-A*31:01" in m and "empty" in m.lower() for m in warning_msgs)

    def test_stitch_failure_skips_row_continues(self, tmp_path, monkeypatch, caplog):
        # stitch_chain returns None for the second call only — simulates one-row stitch failure
        import fetch_vdjdb_panel as mod
        call_count = {"n": 0}
        def flaky_stitch(v_gene, j_gene, cdr3, chain):
            call_count["n"] += 1
            return None if call_count["n"] == 2 else f"MOCK_{chain}_{cdr3}"
        monkeypatch.setattr(mod, "stitch_chain", flaky_stitch)

        panel_tsv = tmp_path / "panel.tsv"
        qc_tsv = tmp_path / "panel_qc.tsv"
        with caplog.at_level(logging.ERROR):
            build_panel(
                vdjdb_full_tsv=FIXTURE_PATH,
                alleles=["HLA-A*02:01"],
                output_panel=panel_tsv,
                output_qc=qc_tsv,
                min_score=2,
                panel_size=3,  # need 3, will fall back when one fails
            )
        import pandas as pd
        panel = pd.read_csv(panel_tsv, sep="\t")
        # Should still get 3 rows because we keep iterating after a stitch failure
        assert len(panel) == 3
        # ERROR was logged for the failed row
        assert any("stitchr" in r.message.lower() or "stitch" in r.message.lower()
                   for r in caplog.records if r.levelno == logging.ERROR)
```

- [ ] **Step 2: Run new tests to confirm failure**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestBuildPanel -v
```
Expected: `ImportError` on `build_panel`.

- [ ] **Step 3: Add the implementation**

Append to `workflow/scripts/fetch_vdjdb_panel.py`:

```python
# ---------------------------------------------------------------------------
# Build panel — orchestrate + write outputs
# ---------------------------------------------------------------------------

PANEL_COLUMNS = [
    "allele", "va_gene", "ja_gene", "cdr3a",
    "vb_gene", "jb_gene", "cdr3b",
    "alpha_seq", "beta_seq",
    "vdjdb_score", "vdjdb_donor_id",
]

QC_COLUMNS = ["allele", "n_exact_matches", "n_in_panel", "panel_status"]


def build_panel(
    vdjdb_full_tsv,
    alleles: list,
    output_panel,
    output_qc,
    min_score: int,
    panel_size: int,
) -> None:
    """Build a per-patient VDJdb panel and write panel.tsv + panel_qc.tsv.

    For each allele in `alleles`:
      - exact-match filter (4-digit normalized)
      - sort by score DESC + donor_id ASC
      - iterate top-down: call stitchr per row, accumulate successes,
        skip+log on stitch failures, stop at panel_size or exhaustion
      - record n_exact_matches + n_in_panel + panel_status in QC
    """
    import pandas as pd  # lazy import
    from pathlib import Path

    output_panel = Path(output_panel)
    output_qc = Path(output_qc)
    output_panel.parent.mkdir(parents=True, exist_ok=True)
    output_qc.parent.mkdir(parents=True, exist_ok=True)

    df = load_and_filter_vdjdb(vdjdb_full_tsv, min_score=min_score)

    panel_rows = []
    qc_rows = []
    for allele in alleles:
        candidates = select_top_n_for_allele(df, allele=allele, n=len(df))
        n_exact = len(candidates)
        n_in_panel = 0
        for _, row in candidates.iterrows():
            if n_in_panel >= panel_size:
                break
            alpha = stitch_chain(
                v_gene=row["v.alpha"], j_gene=row["j.alpha"],
                cdr3=row["cdr3.alpha"], chain="A",
            )
            beta = stitch_chain(
                v_gene=row["v.beta"], j_gene=row["j.beta"],
                cdr3=row["cdr3.beta"], chain="B",
            )
            if alpha is None or beta is None:
                log.error("Skipping VDJdb row for allele %s (donor %s) — stitch failed",
                          allele, row["meta.subject.id"])
                continue
            panel_rows.append({
                "allele": allele,
                "va_gene": row["v.alpha"], "ja_gene": row["j.alpha"], "cdr3a": row["cdr3.alpha"],
                "vb_gene": row["v.beta"],  "jb_gene": row["j.beta"],  "cdr3b": row["cdr3.beta"],
                "alpha_seq": alpha, "beta_seq": beta,
                "vdjdb_score": int(row["vdjdb.score"]),
                "vdjdb_donor_id": row["meta.subject.id"],
            })
            n_in_panel += 1

        status = classify_panel_status(n_in_panel, target_size=panel_size)
        qc_rows.append({
            "allele": allele,
            "n_exact_matches": n_exact,
            "n_in_panel": n_in_panel,
            "panel_status": status,
        })
        if status == "empty":
            log.warning("Allele %s: zero VDJdb entries (panel_status=empty)", allele)
        elif status == "low_coverage":
            log.warning("Allele %s: only %d entries available (panel_status=low_coverage)",
                        allele, n_in_panel)

    pd.DataFrame(panel_rows, columns=PANEL_COLUMNS).to_csv(output_panel, sep="\t", index=False)
    pd.DataFrame(qc_rows, columns=QC_COLUMNS).to_csv(output_qc, sep="\t", index=False)

    n_ok = sum(1 for r in qc_rows if r["panel_status"] == "ok")
    n_low = sum(1 for r in qc_rows if r["panel_status"] == "low_coverage")
    n_empty = sum(1 for r in qc_rows if r["panel_status"] == "empty")
    log.info(
        "VDJdb panel built — %d alleles: %d ok, %d low_coverage, %d empty. Wrote %s and %s.",
        len(alleles), n_ok, n_low, n_empty, output_panel, output_qc,
    )
```

- [ ] **Step 4: Run all tests**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v -m "not network"
```
Expected: all green (network-marked stitchr test skipped).

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add build_panel orchestrator + panel.tsv/panel_qc.tsv writers (Issue #204)"
```

---

## Task 11: Add Snakemake + CLI entry point

**Files:**
- Modify: `workflow/scripts/fetch_vdjdb_panel.py`

- [ ] **Step 1: Append CLI + Snakemake entry**

Append to `workflow/scripts/fetch_vdjdb_panel.py`:

```python
# ---------------------------------------------------------------------------
# Allele loading from alleles.tsv (HLA typing output)
# ---------------------------------------------------------------------------

def load_alleles_tsv(alleles_tsv) -> list:
    """Load unique 4-digit alleles from alleles.tsv produced by aggregate_hla_alleles.

    alleles.tsv schema (existing): rows per locus (A/B/C) with allele1, allele2 columns.
    Returns a deduplicated list of allele strings, all 4-digit.
    """
    import pandas as pd  # lazy import
    df = pd.read_csv(alleles_tsv, sep="\t", dtype=str)
    alleles: set = set()
    for col in ("allele1", "allele2"):
        if col not in df.columns:
            continue
        for a in df[col].dropna().unique():
            normalized = normalize_allele_to_4digit(a)
            if normalized:
                alleles.add(normalized)
    return sorted(alleles)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    """Entry point when called via Snakemake `script:` directive."""
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s — %(message)s")
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    alleles = load_alleles_tsv(snakemake.input.alleles_tsv)  # type: ignore[name-defined]  # noqa: F821
    build_panel(
        vdjdb_full_tsv=snakemake.input.vdjdb_tsv,  # type: ignore[name-defined]  # noqa: F821
        alleles=alleles,
        output_panel=snakemake.output.panel,  # type: ignore[name-defined]  # noqa: F821
        output_qc=snakemake.output.qc,  # type: ignore[name-defined]  # noqa: F821
        min_score=snakemake.params.min_score,  # type: ignore[name-defined]  # noqa: F821
        panel_size=snakemake.params.panel_size,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    """Entry point when called from the command line (for development/debug)."""
    import argparse

    parser = argparse.ArgumentParser(description="Build a VDJdb TCR panel for one patient.")
    parser.add_argument("--vdjdb-tsv", required=True, help="Path to vdjdb_full.txt")
    parser.add_argument("--alleles-tsv", required=True, help="Path to alleles.tsv")
    parser.add_argument("--output-panel", required=True, help="Output panel.tsv path")
    parser.add_argument("--output-qc", required=True, help="Output panel_qc.tsv path")
    parser.add_argument("--min-score", type=int, default=2)
    parser.add_argument("--panel-size", type=int, default=10)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s — %(message)s")
    alleles = load_alleles_tsv(args.alleles_tsv)
    build_panel(
        vdjdb_full_tsv=args.vdjdb_tsv, alleles=alleles,
        output_panel=args.output_panel, output_qc=args.output_qc,
        min_score=args.min_score, panel_size=args.panel_size,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
    except NameError:
        _cli_main()
    else:
        _snakemake_main()
```

- [ ] **Step 2: Add a test for load_alleles_tsv**

Append to `workflow/tests/test_fetch_vdjdb_panel.py`:

```python
from fetch_vdjdb_panel import load_alleles_tsv


class TestLoadAllelesTsv:
    def test_extracts_unique_4digit_alleles(self, tmp_path):
        tsv = tmp_path / "alleles.tsv"
        tsv.write_text(
            "locus\tallele1\tallele2\n"
            "A\tHLA-A*02:01\tHLA-A*31:01\n"
            "B\tHLA-B*08:01\tHLA-B*08:01\n"  # homozygous → dedupe
            "C\tHLA-C*07:01\tHLA-C*03:03\n"
        )
        alleles = load_alleles_tsv(tsv)
        assert sorted(alleles) == [
            "HLA-A*02:01", "HLA-A*31:01",
            "HLA-B*08:01",
            "HLA-C*03:03", "HLA-C*07:01",
        ]
```

- [ ] **Step 3: Run all tests**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py -v -m "not network"
```
Expected: all green.

- [ ] **Step 4: Smoke-test the CLI**

Run a fixture-based CLI smoke test (no real stitchr — just verify args plumbing). Create a stub alleles file:

```bash
echo -e "locus\tallele1\tallele2\nA\tHLA-A*02:01\tHLA-A*31:01" > /tmp/test_alleles.tsv
python workflow/scripts/fetch_vdjdb_panel.py \
  --vdjdb-tsv workflow/tests/fixtures/vdjdb_mini.tsv \
  --alleles-tsv /tmp/test_alleles.tsv \
  --output-panel /tmp/panel.tsv \
  --output-qc /tmp/panel_qc.tsv \
  --min-score 2 --panel-size 10 2>&1 | tail -10
```
Expected: log lines including the summary `VDJdb panel built — 2 alleles:`. Files exist at `/tmp/panel.tsv` (may be empty because real stitchr is needed) and `/tmp/panel_qc.tsv` (will have 2 rows: A*02:01 = empty if stitch fails or n=5 if works, A*31:01 = empty).

Stitch failures here are expected if stitchr isn't installed in `.venv` — that doesn't break this smoke test, which just confirms argparse + entry plumbing works.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/fetch_vdjdb_panel.py workflow/tests/test_fetch_vdjdb_panel.py
git commit -m "feat(vdjdb): add Snakemake + CLI entry points + load_alleles_tsv (Issue #204)"
```

---

## Task 12: Add tcr_panel.smk rule module + Snakefile include

**Files:**
- Create: `workflow/rules/tcr_panel.smk`
- Modify: `Snakefile`

- [ ] **Step 1: Create the rule module**

Write `workflow/rules/tcr_panel.smk`:

```python
# =============================================================================
# Rule module: TCR panel construction (Issue #204)
# =============================================================================
#
# Per-patient: filters VDJdb for HLA-matched paired α/β TCRs, reconstructs
# full chain sequences via stitchr, and writes a reference TCR panel.
#
# Outputs
# -------
#   results/{patient_id}/tcr_panel/vdjdb/panel.tsv      — top-10 per allele
#   results/{patient_id}/tcr_panel/vdjdb/panel_qc.tsv   — per-allele coverage
#
# Gating
# ------
# Gated on config[hla][enabled]. When HLA typing is disabled, alleles.tsv
# isn't produced and this rule is excluded from the DAG. Independent of
# config[tcrdock][enabled] — the panel is a reference set, not a GPU output.

_HLA_ENABLED = config.get("hla", {}).get("enabled", False)


if _HLA_ENABLED:

    rule fetch_vdjdb_panel:
        """Build the per-patient VDJdb TCR panel."""
        input:
            vdjdb_tsv = lambda wc: f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
            vdjdb_sentinel = lambda wc: f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
            imgt_sentinel = "resources/imgt_germlines/.download.done",
            alleles_tsv = os.path.join(_RES, "{patient_id}", "hla_typing", "alleles.tsv"),
        output:
            panel = os.path.join(_RES, "{patient_id}", "tcr_panel", "vdjdb", "panel.tsv"),
            qc = os.path.join(_RES, "{patient_id}", "tcr_panel", "vdjdb", "panel_qc.tsv"),
        log:
            os.path.join(_LOGS, "{patient_id}", "tcr_panel", "fetch_vdjdb_panel.log"),
        params:
            min_score = config["tcrdock"]["vdjdb_min_score"],
            panel_size = config["tcrdock"]["vdjdb_panel_size"],
        conda:
            "../envs/vdjdb.yaml"
        script:
            "../scripts/fetch_vdjdb_panel.py"
```

- [ ] **Step 2: Add include line to Snakefile**

Open `Snakefile`. Find the existing include block (currently 10 lines from `common.smk` → `structure.smk`). Add a new line for `tcr_panel.smk` in a sensible position — between `mhc_affinity.smk` and `analysis.smk` keeps it close to other prediction-adjacent rules:

```diff
 include: "workflow/rules/common.smk"
 include: "workflow/rules/download.smk"
 include: "workflow/rules/alignment.smk"
 include: "workflow/rules/hla_typing.smk"
 include: "workflow/rules/filter_junctions.smk"
 include: "workflow/rules/assemble_contigs.smk"
 include: "workflow/rules/translate_peptides.smk"
 include: "workflow/rules/proteome_filter.smk"
 include: "workflow/rules/mhc_affinity.smk"
+include: "workflow/rules/tcr_panel.smk"
 include: "workflow/rules/analysis.smk"
 include: "workflow/rules/structure.smk"
```

- [ ] **Step 3: Snakemake parse check**

The download rules don't exist yet (Task 13 + 14 add them), so a full dry-run will fail. Run a parse-only check:

```bash
python -c "
import snakemake
# Snakemake 8 parses .smk on workflow construction; just check syntax.
print('parsing Snakefile...')
" 2>&1
python -c "
from snakemake.api import SnakemakeApi
print('OK if no SyntaxError below')
" 2>&1
```

A safer parse check is just to ensure the rule module syntax is valid Python+Snakemake DSL:

```bash
python -c "
# Faux-parse: just ensure the file is syntactically importable as text Python
src = open('workflow/rules/tcr_panel.smk').read()
compile(src, 'tcr_panel.smk', 'exec')
print('OK')
"
```
Expected: `OK` — no SyntaxError. Full DAG resolution comes after Tasks 13 + 14.

- [ ] **Step 4: Commit**

```bash
git add workflow/rules/tcr_panel.smk Snakefile
git commit -m "feat(vdjdb): add tcr_panel.smk module with fetch_vdjdb_panel rule + Snakefile include (Issue #204)

Multi-file commit: the rule module + its Snakefile include are atomically
related — neither is useful alone."
```

---

## Task 13: Add download_vdjdb_release rule

**Files:**
- Modify: `workflow/rules/download.smk`

- [ ] **Step 1: Read download.smk to find the right insertion point**

Run:
```bash
grep -n '^rule\|^# ====' workflow/rules/download.smk | head -30
```

- [ ] **Step 2: Add the rule**

Append to `workflow/rules/download.smk` (at the end of the file):

```python
# =============================================================================
# VDJdb release download — Issue #204
# =============================================================================

rule download_vdjdb_release:
    """Download a pinned VDJdb release, verify SHA256, extract vdjdb_full.txt.

    Sentinel-gated for idempotency. Re-runs are no-ops once the sentinel exists.
    """
    output:
        vdjdb_tsv = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
        sentinel = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
    log:
        f"logs/download/vdjdb_{config['tcrdock']['vdjdb_release']}.log",
    params:
        release = config["tcrdock"]["vdjdb_release"],
        sha256 = config["tcrdock"]["vdjdb_sha256"],
    conda:
        "../envs/vdjdb.yaml"
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.sentinel})
        mkdir -p "$DIR"
        ZIP="$DIR/vdjdb-{params.release}.zip"
        URL="https://github.com/antigenomics/vdjdb-db/releases/download/{params.release}/vdjdb-{params.release}.zip"

        echo "Downloading VDJdb {params.release} from $URL" >> {log} 2>&1
        curl -fsSL "$URL" -o "$ZIP" >> {log} 2>&1

        echo "Verifying SHA256..." >> {log} 2>&1
        ACTUAL=$(shasum -a 256 "$ZIP" | awk '{{print $1}}')
        if [ "$ACTUAL" != "{params.sha256}" ]; then
            echo "SHA256 mismatch! expected={params.sha256} actual=$ACTUAL" >> {log} 2>&1
            exit 1
        fi
        echo "SHA256 OK" >> {log} 2>&1

        echo "Extracting..." >> {log} 2>&1
        unzip -o "$ZIP" -d "$DIR/extracted" >> {log} 2>&1
        cp "$DIR/extracted/vdjdb-{params.release}/vdjdb_full.txt" {output.vdjdb_tsv}

        touch {output.sentinel}
        echo "Done." >> {log} 2>&1
        """
```

- [ ] **Step 3: Snakemake parse check**

```bash
python -c "
src = open('workflow/rules/download.smk').read()
compile(src, 'download.smk', 'exec')
print('OK')
"
```
Expected: `OK`.

- [ ] **Step 4: Commit**

```bash
git add workflow/rules/download.smk
git commit -m "feat(vdjdb): add download_vdjdb_release rule (SHA256-verified) (Issue #204)"
```

---

## Task 14: Add download_imgt_germlines rule

**Files:**
- Modify: `workflow/rules/download.smk`

**Note before writing:** Confirm where `stitchrdl` writes IMGT data by default. Run (with the vdjdb conda env active, since stitchr lives there per Task 3):

```bash
conda activate vdjdb && stitchrdl --help 2>&1 | head -30
```
If `stitchrdl` has a `--output-dir` (or similar) flag, use it to redirect output into `resources/imgt_germlines/`. If it writes to a fixed cache location (e.g. `~/.local/share/stitchr/`), the rule can `cp -r` or `ln -s` after the fact. Adjust the rule accordingly.

- [ ] **Step 1: Add the rule**

Append to `workflow/rules/download.smk`:

```python
# =============================================================================
# IMGT germline download via stitchrdl — Issue #204
# =============================================================================

rule download_imgt_germlines:
    """Download IMGT germline reference data via stitchrdl.

    Sentinel-gated for idempotency. stitchrdl pulls latest IMGT release; the
    rule pins by sentinel only (IMGT itself doesn't expose tagged releases —
    see Known limitations in the design spec).
    """
    output:
        sentinel = "resources/imgt_germlines/.download.done",
    log:
        "logs/download/imgt_germlines.log",
    conda:
        "../envs/vdjdb.yaml"
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.sentinel})
        mkdir -p "$DIR"
        echo "Downloading IMGT germline data via stitchrdl..." >> {log} 2>&1
        # stitchrdl writes to its default cache; we copy into resources/ for reproducibility
        stitchrdl --species HUMAN >> {log} 2>&1
        # Discover stitchr's cache dir from python (preferred over hard-coding)
        STITCHR_CACHE=$(python -c "import IMGTgeneDL; from pathlib import Path; print(Path(IMGTgeneDL.__file__).parent / 'data' / 'HUMAN')")
        if [ -d "$STITCHR_CACHE" ]; then
            cp -r "$STITCHR_CACHE" "$DIR/HUMAN"
            echo "Copied IMGT HUMAN data to $DIR/HUMAN" >> {log} 2>&1
        fi
        touch {output.sentinel}
        echo "Done." >> {log} 2>&1
        """
```

**Caveat:** the exact cache-discovery line is best-effort. If `IMGTgeneDL` doesn't expose a discoverable path, you may need to query `stitchrdl --help` for an `--output-dir` flag, or adjust the path logic. Surface to user if the cache location is unclear.

- [ ] **Step 2: Snakemake parse check**

```bash
python -c "
src = open('workflow/rules/download.smk').read()
compile(src, 'download.smk', 'exec')
print('OK')
"
```
Expected: `OK`.

- [ ] **Step 3: Commit**

```bash
git add workflow/rules/download.smk
git commit -m "feat(vdjdb): add download_imgt_germlines rule (stitchrdl wrapper) (Issue #204)"
```

---

## Task 15: Full Snakemake dry-run validation

**Files:** none modified — verification only.

- [ ] **Step 1: Run dry-run on the production config**

Run (note the `--` terminator per CLAUDE.md Snakemake 8 gotcha):

```bash
conda activate snakemake
snakemake -n --use-conda --configfile config/config.yaml -- results/test/patient_001/tcr_panel/vdjdb/panel.tsv 2>&1 | tail -40
```

Expected: Snakemake resolves the DAG including:
- `download_vdjdb_release`
- `download_imgt_germlines`
- `aggregate_hla_alleles` (existing)
- `fetch_vdjdb_panel`

No `MissingInputException`, no `AmbiguousRuleException`.

- [ ] **Step 2: If dry-run fails, diagnose and fix**

Common failures and fixes:
- `MissingInputException` on alleles.tsv → ensure `config['hla']['enabled']` is `true` in `config.yaml` (or use `test_config.yaml` if it is)
- `KeyError` on `config['tcrdock']['vdjdb_release']` → Task 2 didn't add keys correctly; re-read `config/config.yaml`
- Snakefile include order error → `tcr_panel.smk` must come after `hla_typing.smk` (in Task 12 we placed it after `mhc_affinity.smk`, which is fine)

- [ ] **Step 3: Verify both target outputs are resolvable**

Run a second dry-run targeting both outputs simultaneously:

```bash
snakemake -n --use-conda --configfile config/config.yaml -- \
  results/test/patient_001/tcr_panel/vdjdb/panel.tsv \
  results/test/patient_001/tcr_panel/vdjdb/panel_qc.tsv 2>&1 | tail -20
```

Expected: same DAG; both outputs produced by a single invocation of `fetch_vdjdb_panel`.

- [ ] **Step 4: No commit; this is validation only.**

If the dry-run revealed issues that required code edits, those edits should have gone into the appropriate prior tasks' files and been committed there. This task is verify-only.

---

## Task 16: Full pytest sweep + manual smoke test on chr22 fixture

**Files:** none modified — verification only.

- [ ] **Step 1: Full pytest sweep (excluding network-marked stitchr)**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/ -v --tb=short -m "not network" 2>&1 | tail -30
```
Expected: all green. Ensure no regressions in pre-existing tests (e.g., `test_run_tcrdock.py`, `test_run_mhcflurry.py`).

- [ ] **Step 2: Network-marked stitchr smoke test**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_fetch_vdjdb_panel.py::TestStitchChain -v -m network 2>&1 | tail -10
```
Expected: green if stitchr + IMGT data are installed locally. If skipped, that's acceptable for CI but verify locally before opening PR.

- [ ] **Step 3: Render the DAG (per `feedback_dag_visualization` memory)**

The Snakemake rule changes require regenerating the DAG diagram:

```bash
bash scripts/visualize_dag.sh 2>&1 | tail -10
```
Expected: `dag.svg` regenerated (the script's default output filename — see [scripts/visualize_dag.sh:43](../../../scripts/visualize_dag.sh#L43)); new `download_vdjdb_release`, `download_imgt_germlines`, `fetch_vdjdb_panel` nodes appear.

- [ ] **Step 4: Commit the regenerated DAG**

```bash
git add dag.svg
git commit -m "chore(dag): regenerate after Issue #204 rule additions"
```

- [ ] **Step 5: Surface the manual local run command to the user**

Per `feedback_local_pipeline_run` memory, the user always runs snakemake locally themselves. Present the command:

```bash
# Manual real-run command for local M1 verification (do NOT execute):
conda activate snakemake
snakemake --cores 4 --use-conda --configfile config/config.yaml -- \
  results/test/patient_001/tcr_panel/vdjdb/panel.tsv \
  results/test/patient_001/tcr_panel/vdjdb/panel_qc.tsv
```

Expected when user runs: VDJdb 30 MB download + stitchrdl IMGT download (one-time), then per-patient panel + QC files produced.

---

## Task 17: PR open — title, body, board status

**Files:** none modified — GitHub orchestration only.

- [ ] **Step 1: Push branch**

```bash
git push -u origin 204-feattcrdock-fetch_vdjdb_panel-snakemake-rule-+-stitchrvdjdb-setup
```

- [ ] **Step 2: Open the PR**

Use the canonical body format. Per repository convention (Created by: Developer at top of body):

```bash
gh pr create \
  --title "feat(vdjdb): fetch_vdjdb_panel Snakemake rule + stitchr/VDJdb setup (Issue #204)" \
  --body "$(cat <<'EOF'
**Created by:** Developer

Closes [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) (closes).
Parent epic: [Issue #86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86) (epic).
Spinoff: [Issue #435](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/435) (predictions/ flatten, deferred).

## Summary

Adds a per-patient VDJdb TCR panel rule. Pipeline produces `results/{p}/tcr_panel/vdjdb/panel.tsv` (top-10 paired α/β TCRs per HLA allele, with full chain sequences from stitchr) + `panel_qc.tsv` (per-allele coverage status: ok/low_coverage/empty).

Foundation for sub-issues [Issue #205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205) (selection) and [Issue #206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206) (report).

## Empirical pre-step

Performed before any code per the issue's pre-implementation step. Results in [Issue #204 comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204#issuecomment-4502381646):
- patient_001: 3/6 alleles empty in VDJdb (motivates deferred supertype-fallback follow-up)
- patient_002: 0 empty, mostly low-coverage on HLA-C

## Test plan

- [ ] Pytest sweep green (incl. new `test_fetch_vdjdb_panel.py`)
- [ ] Snakemake dry-run green on `config/config.yaml`
- [ ] DAG diagram regenerated
- [ ] Stitchr smoke test on DMF5 alpha chain green (real IMGT data)
- [ ] Local end-to-end run on chr22 test config: `panel.tsv` + `panel_qc.tsv` produced

## Out of scope

- Report wiring → [Issue #206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206)
- TCR selection from panel → [Issue #205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205)
- Supertype fallback for empty alleles → follow-up issue motivated by the pre-step data
- Flattening `results/{p}/predictions/` wrapper → [Issue #435](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/435)

## Cross-references

- Design spec: [docs/superpowers/specs/2026-05-20-vdjdb-panel-design.md](docs/superpowers/specs/2026-05-20-vdjdb-panel-design.md)
- Implementation plan: [docs/superpowers/plans/2026-05-20-vdjdb-panel-plan.md](docs/superpowers/plans/2026-05-20-vdjdb-panel-plan.md)
EOF
)"
```

- [ ] **Step 3: Set Status: Ready for review on the board**

```bash
gh api graphql -f query='mutation { updateProjectV2ItemFieldValue(input: {projectId: "PVT_kwHOB17eGc4BSomP", itemId: "PVTI_lAHOB17eGc4BSomPzgrdSfU", fieldId: "PVTSSF_lAHOB17eGc4BSomPzhAHFf8", value: {singleSelectOptionId: "8bf9192f"}}) { projectV2Item { id } } }'
```

Per memory: "When opening a PR, immediately set Status to 'Ready for review' (`8bf9192f`) — never 'In review' (`df73e18b`), that's the reviewer's transition."

- [ ] **Step 4: After review feedback incorporated, write lab notebook entry**

Per the always-in-effect memory: "Lab notebook entry comes AFTER review, before merge — not before commit." This step happens AFTER review comments have been addressed and the PR is rebased to current main. Use the existing `research/lab_notebook/developer.md` format.

(Do NOT do this step pre-review; that violates the timing rule.)

- [ ] **Step 5: Merge via the closure-ritual gate script**

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER> --squash --delete-branch
```

Per CLAUDE.md: `audit_and_merge.sh` enforces that all `- [ ]` boxes on the PR Test plan + linked Issue Acceptance criteria are ticked before merging. Do NOT use bare `gh pr merge`.

---

## Self-review checklist (performed during plan-writing)

- ✅ **Spec coverage:** All design sections mapped — empirical pre-step (done), env (Task 3), download rules (13, 14), `fetch_vdjdb_panel` rule (12), implementation (5–11), config keys (2), glossary KIR (1), output layout (12), error handling (built into 9, 10), tests (5–10), acceptance criteria (mirrored in PR body Task 17).
- ✅ **Placeholder scan:** No TBDs, TODOs, vague "add error handling" steps. Every code step has actual code; every command has expected output.
- ✅ **Type consistency:** `normalize_allele_to_4digit` returns `Optional[str]`, `load_and_filter_vdjdb` returns DataFrame with `mhc.a_4digit` column, `select_top_n_for_allele` consumes that column, `build_panel` writes the 11-col panel + 4-col QC schema consistently across tests and implementation. `panel_status` strings match between `classify_panel_status` and the schema assertions.
- ✅ **One known caveat surfaced inline:** `stitchrdl` cache discovery in Task 14 is best-effort; deferred to runtime adjustment if needed (better than a hard placeholder).
