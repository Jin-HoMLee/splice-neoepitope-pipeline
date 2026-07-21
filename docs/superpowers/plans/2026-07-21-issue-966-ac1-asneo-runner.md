# #966 AC-1 ASNEO End-to-End Runner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the caller-benchmark harness invoke ASNEO end-to-end on a chr22 `SJ.out.tab` and collect its output into the common schema, satisfying #966 AC-1.

**Architecture:** Add a `RUNNERS` registry beside the existing `ADAPTERS` registry; a runner shells out to the option-B-patched ASNEO caller (reusing the proven #965 smoke sequence) and an adapter parses its peptide-only output. A new `record_level` schema discriminator + a legal-combination validator represent ASNEO's peptide-level records honestly.

**Tech Stack:** Python 3.13 (harness venv), pytest, `conda run -n asneo` for the caller, dataclasses, plain TSV.

## Global Constraints

- Harness lives in `research/experiments/issue_679_caller_benchmark/harness/`; run pytest from that directory (conftest puts the dir on `sys.path`, so imports are `from common_schema import ...` and `from adapters.X import ...`).
- Test interpreter: `workflow/tests/.venv/bin/python -m pytest` is the repo default, but this harness has no external deps; `python -m pytest` from the harness dir with any 3.13 works. Use `conda activate snakemake` for ad-hoc runs.
- No em-dash / en-dash in any authored text (plain hyphen only).
- ASNEO caller runs only in the `asneo` conda env (exists at `~/miniforge3/envs/asneo`); scratch (ASNEO clone, hg19 chr22 FASTA) stays OUT of the repo tree (default `${TMPDIR}`), so no >100 MB artifact is committed and AC-4 stays vacuous.
- `record_level` is a required field with values `"junction"` | `"peptide"`, orthogonal to `event_type` (splice-event kind).

---

### Task 1: Schema - `record_level` discriminator + legal-combination validator

**Files:**
- Modify: `research/experiments/issue_679_caller_benchmark/harness/common_schema.py`
- Test: `research/experiments/issue_679_caller_benchmark/harness/tests/test_common_schema.py`

**Interfaces:**
- Produces: `CommonRecord` gains a required `record_level: str` field; `junction_id` and `strand` become `Optional[str]` (still required to pass, may be `None`). New `validate(record: CommonRecord) -> None`.

- [ ] **Step 1: Write the failing test** (append to `test_common_schema.py`)

```python
import pytest
from common_schema import CommonRecord, validate


def _peptide_record(**overrides):
    base = dict(caller="asneo", junction_id=None, genome_build="hg19",
                strand=None, peptide="LEQGTHPKFQ", record_level="peptide")
    base.update(overrides)
    return CommonRecord(**base)


def test_validate_peptide_level_allows_null_junction():
    validate(_peptide_record())  # must not raise


def test_validate_junction_level_requires_junction_id():
    # Matched pair: flip only record_level; opposite outcomes.
    bad = _peptide_record(record_level="junction")   # junction_id/strand still None
    with pytest.raises(ValueError):
        validate(bad)
    good = _peptide_record(record_level="junction",
                           junction_id="chr22:100-200:+", strand="+")
    validate(good)  # must not raise


def test_validate_peptide_level_rejects_a_junction_id():
    with pytest.raises(ValueError):
        validate(_peptide_record(junction_id="chr22:100-200:+", strand="+"))


def test_validate_unknown_record_level_raises():
    with pytest.raises(ValueError):
        validate(_peptide_record(record_level="isoform"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_common_schema.py -v`
Expected: FAIL (`ImportError: cannot import name 'validate'` / `TypeError` on `record_level`).

- [ ] **Step 3: Write minimal implementation** (edit `common_schema.py`: insert `record_level`, retype `junction_id`/`strand`, add `validate`)

Change the dataclass field block to:

```python
    caller: str
    junction_id: Optional[str]
    genome_build: str
    strand: Optional[str]
    peptide: str
    record_level: str  # "junction" | "peptide"; orthogonal to event_type
    frame_shift: Optional[bool] = None
    event_type: Optional[str] = None
    transcript_id: Optional[str] = None
    presentation_class: Optional[str] = None
    presentation_score: Optional[float] = None
    presentation_percentile: Optional[float] = None
    provenance: dict = field(default_factory=dict)
```

Add at module level:

```python
def validate(record: "CommonRecord") -> None:
    """Enforce the legal record_level / field combinations.

    A ``junction``-level record must carry junction coordinates; a
    ``peptide``-level record must not (its junction linkage is unavailable by
    construction - see #1258). Raise loudly on any illegal combination so an
    adapter bug cannot smuggle a mislabeled record downstream.
    """
    if record.record_level == "junction":
        if record.junction_id is None or record.strand is None:
            raise ValueError(
                f"junction-level record for {record.caller!r} needs junction_id "
                f"and strand; got junction_id={record.junction_id!r}, strand={record.strand!r}"
            )
    elif record.record_level == "peptide":
        if record.junction_id is not None or record.strand is not None:
            raise ValueError(
                f"peptide-level record for {record.caller!r} must have null junction_id "
                f"and strand; got junction_id={record.junction_id!r}, strand={record.strand!r}"
            )
    else:
        raise ValueError(
            f"unknown record_level {record.record_level!r} for {record.caller!r}; "
            f"expected 'junction' or 'peptide'"
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_common_schema.py -v`
Expected: PASS (all validate tests green).

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_679_caller_benchmark/harness/common_schema.py \
        research/experiments/issue_679_caller_benchmark/harness/tests/test_common_schema.py
git commit -m "feat(benchmark): add record_level discriminator + legal-combination validator (#966)"
```

---

### Task 2: splice2neo adapter declares `record_level="junction"`

**Files:**
- Modify: `research/experiments/issue_679_caller_benchmark/harness/adapters/splice2neo.py:47-62`
- Test: `research/experiments/issue_679_caller_benchmark/harness/tests/test_splice2neo_adapter.py`

**Interfaces:**
- Consumes: `CommonRecord` (Task 1). Produces: splice2neo records now carry `record_level="junction"`.

- [ ] **Step 1: Write the failing test** (add to `test_splice2neo_adapter.py`)

```python
def test_splice2neo_records_are_junction_level(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(
        "junc_id\ttx_id\tframe_shift\tcts_seq_len\tpeptide_context\n"
        "chr2:152389996-152392205:-\tENST00000409198\tTRUE\t400\tINRHFKYATQLMNEIC\n"
    )
    (r,) = parse_splice2neo(tsv, genome_build="hg19")
    assert r.record_level == "junction"
    assert r.junction_id is not None and r.strand is not None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_splice2neo_adapter.py::test_splice2neo_records_are_junction_level -v`
Expected: FAIL (`TypeError: __init__() missing 1 required positional argument: 'record_level'`).

- [ ] **Step 3: Write minimal implementation** (in `parse_splice2neo`, add `record_level="junction"` to the `CommonRecord(...)` call, right after `peptide=peptide,`)

```python
                    peptide=peptide,
                    record_level="junction",
                    frame_shift=_to_bool(row.get("frame_shift", "")),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_splice2neo_adapter.py -v`
Expected: PASS (all splice2neo tests green).

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_679_caller_benchmark/harness/adapters/splice2neo.py \
        research/experiments/issue_679_caller_benchmark/harness/tests/test_splice2neo_adapter.py
git commit -m "feat(benchmark): splice2neo adapter emits record_level=junction (#966)"
```

---

### Task 3: ASNEO adapter - `parse_asneo`

**Files:**
- Create: `research/experiments/issue_679_caller_benchmark/harness/adapters/asneo.py`
- Test: `research/experiments/issue_679_caller_benchmark/harness/tests/test_asneo_adapter.py`

**Interfaces:**
- Consumes: `CommonRecord` (Task 1). Produces: `parse_asneo(path, genome_build="hg19", provenance=None) -> list[CommonRecord]`, all `record_level="peptide"`.

- [ ] **Step 1: Write the failing test** (`test_asneo_adapter.py`)

```python
from adapters.asneo import parse_asneo


def test_parses_bare_peptides_as_peptide_level(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text(
        "# ASNEO chr22 smoke - first N candidate peptides\n"
        "LEQGTHPKFQ\n"
        "\n"
        "ILQPKPVD\n"
    )
    records = parse_asneo(f, genome_build="hg19",
                          provenance={"sj_tab_digest": "abc123", "thresholds": "reads2 psi0.05"})
    assert [r.peptide for r in records] == ["LEQGTHPKFQ", "ILQPKPVD"]
    r0 = records[0]
    assert r0.caller == "asneo"
    assert r0.record_level == "peptide"
    assert r0.junction_id is None and r0.strand is None
    assert r0.genome_build == "hg19"
    assert r0.event_type == "junction"
    assert r0.provenance["source"] == "asneo"
    assert r0.provenance["sj_tab_digest"] == "abc123"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_asneo_adapter.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'adapters.asneo'`).

- [ ] **Step 3: Write minimal implementation** (`adapters/asneo.py`)

```python
"""ASNEO adapter (#966).

ASNEO (bm2-lab, Apache-2.0) is an end-to-end splice-neoantigen caller. Its
open-only stop-point output (``putative_peptide.txt``, ``ASNEO.py:246-253``)
is a normal-subtracted SET of k-mer peptide strings - one bare peptide per
line. That set operation discards the junction -> peptide linkage, so these
records are peptide-level (``record_level="peptide"``, null junction fields).
Recovering the linkage is tracked in #1258 (a patch to ASNEO), NOT #1100.
"""

from typing import Optional

from common_schema import CommonRecord

CALLER = "asneo"


def parse_asneo(path, genome_build: str = "hg19",
                provenance: Optional[dict] = None) -> list[CommonRecord]:
    """Parse an ASNEO ``putative_peptide.txt`` (one bare peptide per line).

    Blank lines and ``#`` comment lines (our smoke annotations) are skipped.
    Each peptide becomes a peptide-level record with null junction fields.
    """
    base_prov = {"source": CALLER}
    if provenance:
        base_prov.update(provenance)
    records: list[CommonRecord] = []
    with open(path) as fh:
        for line in fh:
            peptide = line.strip()
            if not peptide or peptide.startswith("#"):
                continue
            records.append(
                CommonRecord(
                    caller=CALLER,
                    junction_id=None,
                    genome_build=genome_build,
                    strand=None,
                    peptide=peptide,
                    record_level="peptide",
                    event_type="junction",
                    provenance=dict(base_prov),
                )
            )
    return records
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_asneo_adapter.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_679_caller_benchmark/harness/adapters/asneo.py \
        research/experiments/issue_679_caller_benchmark/harness/tests/test_asneo_adapter.py
git commit -m "feat(benchmark): ASNEO peptide-only adapter (#966)"
```

---

### Task 4: ASNEO runner - subprocess wrapper + pure validation

**Files:**
- Create: `research/experiments/issue_679_caller_benchmark/harness/runners/__init__.py` (empty)
- Create: `research/experiments/issue_679_caller_benchmark/harness/runners/asneo.py`
- Create: `research/experiments/issue_679_caller_benchmark/harness/runners/run_asneo.sh`
- Test: `research/experiments/issue_679_caller_benchmark/harness/tests/test_asneo_runner.py`

**Interfaces:**
- Produces: `run_asneo(sj_tab: str, workdir: str) -> str` (path to `putative_peptide.txt`); pure helpers `peptide_output_path(workdir) -> str` and `validate_output(path) -> None`.

- [ ] **Step 1: Write the failing test** (`test_asneo_runner.py` - covers only the pure parts; the live subprocess run is Task 6)

```python
import pytest
from runners.asneo import peptide_output_path, validate_output


def test_peptide_output_path_is_under_workdir(tmp_path):
    assert peptide_output_path(str(tmp_path)).endswith("putative_peptide.txt")


def test_validate_output_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_output(str(tmp_path / "nope.txt"))


def test_validate_output_zero_peptides_raises(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text("# only a comment, no peptides\n\n")
    with pytest.raises(ValueError):
        validate_output(str(f))


def test_validate_output_with_peptides_passes(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text("# header\nLEQGTHPKFQ\n")
    validate_output(str(f))  # must not raise
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_asneo_runner.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'runners'`).

- [ ] **Step 3: Write minimal implementation**

`runners/__init__.py`: empty file.

`runners/asneo.py`:

```python
"""ASNEO runner (#966, AC-1).

Invokes the open-only (option-B patched) ASNEO caller on a chr22-scale
SJ.out.tab in the ``asneo`` conda env, returning the path to the produced
putative_peptide.txt. The clone / patch / genome / run sequence lives in the
proven ``run_asneo.sh`` (adapted from the #965 smoke); this module shells out
to it and validates the result. Command orchestration is bash; the output
contract is checked here in pure, unit-tested Python.
"""

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).with_name("run_asneo.sh")


def peptide_output_path(workdir: str) -> str:
    """Where run_asneo.sh writes the caller's candidate peptides (pure)."""
    return str(Path(workdir) / "out" / "putative_peptide.txt")


def validate_output(path: str) -> None:
    """Fail loud if ASNEO produced no candidate peptides (pure)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"ASNEO produced no output at {path}")
    peptides = [ln.strip() for ln in p.read_text().splitlines()
                if ln.strip() and not ln.strip().startswith("#")]
    if not peptides:
        raise ValueError(f"ASNEO output {path} has 0 candidate peptides")


def run_asneo(sj_tab: str, workdir: str) -> str:
    """Run patched ASNEO on ``sj_tab``; return the putative_peptide.txt path."""
    subprocess.run(["bash", str(_SCRIPT), sj_tab, workdir], check=True)
    out = peptide_output_path(workdir)
    validate_output(out)
    return out
```

`runners/run_asneo.sh` (adapt the #965 `asneo_smoke.sh`: take the input `SJ.out.tab` as `$1`, workdir as `$2`, and write `putative_peptide.txt` to `$2/out/`):

```bash
#!/usr/bin/env bash
# ASNEO runner (#966 AC-1): run the option-B-patched ASNEO caller on a given
# chr22-scale SJ.out.tab, no STAR, against the hg19 chr22 FASTA. Reuses the
# proven #965 smoke sequence but parametrized on an input junction table.
# Usage: bash run_asneo.sh <sj_out_tab> <workdir>
set -euo pipefail
SJ_TAB="$1"
WORK="${2:?workdir required}"
REPO="$(git rev-parse --show-toplevel)"
PATCH="$REPO/research/experiments/issue_566_asneo_crosscheck/apply_optionB_patch.py"
mkdir -p "$WORK" "$WORK/out"

[ -d "$WORK/ASNEO" ] || git clone --depth 1 https://github.com/bm2-lab/ASNEO.git "$WORK/ASNEO"
if [ ! -f "$WORK/chr22.fa" ]; then
  curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz" -o "$WORK/chr22.fa.gz"
  gunzip -f "$WORK/chr22.fa.gz"
fi
grep -q "option B" "$WORK/ASNEO/ASNEO.py" || conda run -n asneo python "$PATCH" "$WORK/ASNEO/ASNEO.py"

cd "$WORK/ASNEO"
conda run -n asneo python ASNEO.py -j "$SJ_TAB" -a HLA-A02:01 -g "$WORK/chr22.fa" \
  -o "$WORK/out" -l 8,9,10,11 --reads 2 --psi 0.05
# The option-B patch copies putative_peptide.txt into -o (outdir) before tmp/ cleanup.
test -f "$WORK/out/putative_peptide.txt"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_asneo_runner.py -v`
Expected: PASS (pure helpers; no subprocess invoked in tests).

- [ ] **Step 5: Commit**

```bash
chmod +x research/experiments/issue_679_caller_benchmark/harness/runners/run_asneo.sh
git add research/experiments/issue_679_caller_benchmark/harness/runners/
git add research/experiments/issue_679_caller_benchmark/harness/tests/test_asneo_runner.py
git commit -m "feat(benchmark): ASNEO runner (subprocess wrapper + validated output contract) (#966)"
```

---

### Task 5: `collect.py` - `RUNNERS` registry, `--run` mode, `validate()` on merge

**Files:**
- Modify: `research/experiments/issue_679_caller_benchmark/harness/collect.py`
- Test: `research/experiments/issue_679_caller_benchmark/harness/tests/test_collect.py`

**Interfaces:**
- Consumes: `parse_asneo` (Task 3), `run_asneo` (Task 4), `validate` (Task 1). Produces: `RUNNERS` dict; `collect(specs)` calls `validate` on every record; CLI gains `--run CALLER:SJ_TAB`.

- [ ] **Step 1: Write the failing test** (add to `test_collect.py`)

```python
def test_collect_validates_records(monkeypatch):
    # An adapter that emits an illegal record (peptide-level with a junction_id)
    # must be caught at collect time, not silently written.
    import collect as collect_mod
    from common_schema import CommonRecord

    def bad_adapter(path):
        return [CommonRecord(caller="x", junction_id="chr1:1-2:+", genome_build="hg19",
                             strand="+", peptide="MEIC", record_level="peptide")]
    monkeypatch.setitem(collect_mod.ADAPTERS, "bad", bad_adapter)
    import pytest
    with pytest.raises(ValueError):
        collect_mod.collect([("bad", "ignored")])


def test_runners_registry_has_asneo():
    import collect as collect_mod
    assert "asneo" in collect_mod.RUNNERS
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/test_collect.py -v`
Expected: FAIL (`AttributeError: module 'collect' has no attribute 'RUNNERS'`, and the bad record is not yet rejected).

- [ ] **Step 3: Write minimal implementation** (edit `collect.py`)

Add imports near the top:

```python
from adapters.asneo import parse_asneo
from common_schema import CommonRecord, validate
from runners.asneo import run_asneo
```

Add the adapter + runner registries:

```python
ADAPTERS = {
    "splice2neo": parse_splice2neo,
    "asneo": parse_asneo,
}

# The runner registry: caller -> a function (sj_tab, workdir) -> output_path.
# One line per caller alongside its adapter is the whole extension point (AC-5).
RUNNERS = {
    "asneo": run_asneo,
}
```

Make `collect` validate every record:

```python
def collect(specs):
    """Dispatch each ``(caller, path)`` spec to its adapter, validate, merge."""
    records = []
    for caller, path in specs:
        for rec in get_adapter(caller)(path):
            validate(rec)
            records.append(rec)
    return records
```

Add a `--run` argument and dispatch (in `main`, before `collect(...)`):

```python
    parser.add_argument(
        "--run", dest="runs", action="append", default=[], type=_parse_input_spec,
        metavar="CALLER:SJ_TAB",
        help="run a caller on an SJ.out.tab, then ingest its output (repeatable)",
    )
```

and in `main`, resolve runs into input specs before collecting:

```python
    specs = list(args.inputs or [])
    for caller, sj_tab in args.runs:
        try:
            runner = RUNNERS[caller]
        except KeyError:
            known = ", ".join(sorted(RUNNERS))
            raise ValueError(f"no runner for {caller!r}; registered runners: {known}")
        import tempfile
        out_path = runner(sj_tab, tempfile.mkdtemp(prefix=f"{caller}_run_"))
        specs.append((caller, out_path))
    records = collect(specs)
```

Change `--input` to `required=False` (since `--run` can supply specs) and add a guard that at least one of `--input`/`--run` is given:

```python
    parser.add_argument("--input", dest="inputs", action="append", default=[],
                        type=_parse_input_spec, metavar="CALLER:PATH",
                        help="a caller name and its output path (repeatable)")
    ...
    if not args.inputs and not args.runs:
        parser.error("provide at least one --input or --run")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_679_caller_benchmark/harness && python -m pytest tests/ -v`
Expected: PASS (full suite green, including the two new collect tests).

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_679_caller_benchmark/harness/collect.py \
        research/experiments/issue_679_caller_benchmark/harness/tests/test_collect.py
git commit -m "feat(benchmark): RUNNERS registry + --run mode + validate-on-merge (#966)"
```

---

### Task 6: Live chr22 ASNEO run (ticks AC-1) + fixture + README

**Files:**
- Create: `research/experiments/issue_679_caller_benchmark/inputs/chr22_SRR2660032.SJ.out.tab` (committed fixture)
- Create: `research/experiments/issue_679_caller_benchmark/harness/outputs/asneo_chr22_unified.tsv` (the AC-1 deliverable)
- Modify: `research/experiments/issue_679_caller_benchmark/README.md`

**Interfaces:** consumes everything above; produces the recorded end-to-end artifact.

- [ ] **Step 1: Generate + commit the chr22 input fixture**

```bash
WORK="${TMPDIR:-/tmp}/asneo_fixture"; mkdir -p "$WORK"
[ -d "$WORK/ASNEO" ] || git clone --depth 1 https://github.com/bm2-lab/ASNEO.git "$WORK/ASNEO"
mkdir -p research/experiments/issue_679_caller_benchmark/inputs
awk '$1=="chr22"' "$WORK/ASNEO/test/SRR2660032.SJ.out.tab" \
  > research/experiments/issue_679_caller_benchmark/inputs/chr22_SRR2660032.SJ.out.tab
wc -l research/experiments/issue_679_caller_benchmark/inputs/chr22_SRR2660032.SJ.out.tab   # expect ~6194
```

- [ ] **Step 2: Run the harness end-to-end (the AC-1 smoke)**

Run:
```bash
cd research/experiments/issue_679_caller_benchmark/harness
conda activate snakemake  # for the python that drives collect.py; the runner uses `conda run -n asneo` internally
python collect.py --run asneo:../inputs/chr22_SRR2660032.SJ.out.tab --out outputs/asneo_chr22_unified.tsv
```
Expected: prints `collected N records from 1 input(s) -> outputs/asneo_chr22_unified.tsv` with N in the hundreds (the relaxed-threshold chr22 run emits ~800 candidate peptides).

- [ ] **Step 3: Sanity-check the output**

Run:
```bash
head -3 outputs/asneo_chr22_unified.tsv
awk -F'\t' 'NR>1{print $NF}' outputs/asneo_chr22_unified.tsv | head   # provenance JSON present
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)print i,$i}' outputs/asneo_chr22_unified.tsv | grep record_level
```
Expected: header includes `record_level`; every ASNEO row has `record_level=peptide`, empty `junction_id`/`strand`, and a provenance JSON cell.

- [ ] **Step 4: Update the README** (add a `## Running a caller end-to-end (AC-1)` section)

```markdown
## Running a caller end-to-end (AC-1)

ASNEO (open-only, option-B patched) runs locally on a chr22 `SJ.out.tab`:

    cd harness
    python collect.py --run asneo:../inputs/chr22_SRR2660032.SJ.out.tab \
      --out outputs/asneo_chr22_unified.tsv

The runner clones + patches ASNEO and fetches the hg19 chr22 FASTA into an
out-of-tree scratch dir (no >100 MB artifact committed), runs the caller, and
ingests its peptide-only output. ASNEO records are `record_level=peptide`
(null junction fields); per-peptide junction linkage is tracked in #1258.

Adding a caller = one `RUNNERS` + one `ADAPTERS` entry (AC-5).
```

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_679_caller_benchmark/inputs/chr22_SRR2660032.SJ.out.tab \
        research/experiments/issue_679_caller_benchmark/harness/outputs/asneo_chr22_unified.tsv \
        research/experiments/issue_679_caller_benchmark/README.md
git commit -m "feat(benchmark): ASNEO end-to-end chr22 run wired into the harness (closes AC-1, #966)"
```

---

## Post-implementation (not a task, do after Task 6)

- Regenerate `outputs/unified_smoke.tsv` if it is still referenced, since the column order changed (added `record_level`): re-run its documented `--input` invocation, or note it as superseded by `asneo_chr22_unified.tsv`.
- Resolve PR #1251's project-map-atlas conflict (merge main, take main's atlas, regenerate) before pushing - same procedure used for #1250.
- Tick #966 AC-1 (and re-verify AC-2/3/4/5) on the Issue body, then run `scripts/audit_and_merge.sh 1251` (do NOT merge without Jin-Ho's go).
- Lab-notebook entry for the AC-1 completion before merge.

## Self-Review

- **Spec coverage:** AC-1 -> Tasks 4+6; schema/record_level -> Task 1; splice2neo consistency -> Task 2; ASNEO adapter -> Task 3; RUNNERS/--run/validate -> Task 5; README/AC-3 -> Task 6; AC-4 vacuous (out-of-tree scratch, Global Constraints); AC-5 -> Tasks 3/4/5 registries. `record_level` discriminator + validator, deferral to #1258-not-#1100, per-caller conda env, live-integration smoke all present.
- **Placeholder scan:** no TBD/TODO; all code shown. The one bash file (`run_asneo.sh`) is given in full and marked as an adaptation of the verified #965 `asneo_smoke.sh`.
- **Type consistency:** `parse_asneo`, `run_asneo(sj_tab, workdir)`, `peptide_output_path`, `validate_output`, `validate`, `record_level`, `RUNNERS`/`ADAPTERS` used consistently across tasks.
