# ORF-stretch FASTA Emitter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Emit per-junction 3-frame stop-split ORF stretches crossing the breakpoint as a FASTA of mini-proteins, consumable by a nonspecific MS search, to unblock the [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) re-search.

**Architecture:** Reuse `assemble_contigs.py` at 30-nt flanks (decoupled from the 27-nt MHC path) to cut 60-nt contigs; a new `orf_fasta_from_contigs.py` translates all 3 frames, splits at stop codons, keeps the single breakpoint-crossing stretch, and writes it as FASTA. A new optional Snakemake rule wires it; the production MHCflurry k-mer path is untouched.

**Tech Stack:** Python 3.10+, Biopython (`Bio.Seq`), Snakemake 8, pytest, Sage (MS search, manual smoke).

## Global Constraints

- Header format, exact: `>{junction_id}|{frame}|{chrom}:{start}-{end}:{strand}` (drop the contig header's trailing `|{sample_type}`).
- Flank width 30 nt each side; breakpoint boundary at nt 30 (0-based); `min_peptide_len` default 8.
- Crossing condition, exact: keep a stop-free stretch spanning codons `[c0, c1)` iff `frame + 3*c0 < 30 < frame + 3*c1`.
- Drop any crossing stretch containing `X` or shorter than `min_peptide_len`.
- No `from __future__ import annotations` in any `script:`-invoked file (Snakemake wrapper prepends a preamble; `__future__` then raises `SyntaxError`). PEP 604 unions (`str | None`) are fine on 3.10+.
- Pytest: `workflow/tests/.venv/bin/python -m pytest`.
- No em-dash or en-dash anywhere; plain hyphen only.
- Snakemake dry-run: single `--configfile A B`, and `--` before any positional target.
- `translation.peptide_lengths` stays `[8,9,10]`; do not touch the MHCflurry path.

## File Structure

- Create `workflow/scripts/orf_fasta_from_contigs.py` - the emitter (core function + FASTA I/O + CLI/Snakemake entry).
- Create `workflow/rules/ms_search_db.smk` - a wide (30-nt) contig-assembly rule reusing `assemble_contigs.py`, plus the emit rule.
- Modify `config/config.yaml` - add the `ms_search_db:` block.
- Modify `workflow/Snakefile` - `include:` the new rule file.
- Create `workflow/tests/test_orf_fasta_from_contigs.py` - unit tests + structural/multi-chr check.
- Create `workflow/tests/data/orf_fasta/` - multi-chromosome contig fixture + tiny canonical proteome.
- Create `research/experiments/issue_1204_orf_fasta/` - Sage smoke fixtures, config, run script, README.

---

### Task 1: Core ORF-stretch extraction

**Files:**
- Create: `workflow/scripts/orf_fasta_from_contigs.py`
- Test: `workflow/tests/test_orf_fasta_from_contigs.py`

**Interfaces:**
- Produces: `crossing_orf_stretch(contig_seq: str, frame: int, breakpoint_nt: int = 30) -> str | None` - the amino-acid string of the single stop-free stretch that crosses the breakpoint, or `None` if none crosses in that frame.

- [ ] **Step 1: Write the failing tests**

```python
# workflow/tests/test_orf_fasta_from_contigs.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from orf_fasta_from_contigs import crossing_orf_stretch


def _codon_seq(aas, frame=0, upstream_codons=10):
    # Build a 60-nt contig whose frame-`frame` translation is `aas`.
    # One unambiguous codon per residue; breakpoint at nt 30.
    table = {"M": "ATG", "K": "AAA", "L": "CTG", "F": "TTT", "*": "TAA", "P": "CCC"}
    body = "".join(table[a] for a in aas)
    pad = "G" * frame
    seq = (pad + body)
    return (seq + "C" * 60)[:60]


def test_clean_crossing_stretch_frame0():
    # 20 codons, no stop: one stretch spanning nt [0,60), crosses 30.
    contig = "".join(["AAA"] * 20)  # all Lys, no stop
    assert crossing_orf_stretch(contig, 0) == "K" * 20


def test_stop_at_breakpoint_yields_none():
    # Stop codon occupies codon index 10 -> nt [30,33). Upstream run ends at
    # nt 30 (entirely upstream), downstream run starts at nt 33: neither crosses.
    contig = "".join(["AAA"] * 10 + ["TAA"] + ["AAA"] * 9)
    assert crossing_orf_stretch(contig, 0) is None


def test_upstream_only_stretch_dropped():
    # Stop at codon 5 (nt 15) then run to end: the surviving run nt [18,60)
    # starts downstream? No - it starts at 18 < 30 < 60, so it DOES cross.
    # Construct an upstream-only run: stop at codon 12 (nt 36), first run
    # nt [0,36) crosses; make first run end BEFORE 30 with a stop at codon 8.
    contig = "".join(["AAA"] * 8 + ["TAA"] + ["AAA"] * 11)
    # first run nt [0,24) entirely upstream (24 < 30) -> dropped;
    # second run nt [27,60) crosses (27 < 30 < 60) -> returned.
    assert crossing_orf_stretch(contig, 0) == "K" * 11


def test_frame_offset_shifts_breakpoint_math():
    # frame 1: codon i covers nt [1+3i, 1+3i+3); 19 codons over nt [1,58).
    contig = "G" + "".join(["AAA"] * 19)
    assert crossing_orf_stretch(contig, 1) == "K" * 19
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py -v`
Expected: FAIL with `ModuleNotFoundError` / `ImportError: cannot import name 'crossing_orf_stretch'`.

- [ ] **Step 3: Write the minimal implementation**

```python
#!/usr/bin/env python3
"""orf_fasta_from_contigs.py - emit breakpoint-crossing ORF stretches as FASTA.

For each junction contig (30 nt upstream + 30 nt downstream of the breakpoint),
translate all three reading frames, split each translation at stop codons, and
keep the single stop-free ORF stretch that crosses the breakpoint (nt 30). Each
kept stretch is written as a FASTA mini-protein for a nonspecific MS search.
"""

import argparse
import logging
from pathlib import Path

from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_CODON = 3


def crossing_orf_stretch(contig_seq, frame, breakpoint_nt=30):
    """Return the amino-acid string of the stop-free stretch crossing the
    breakpoint in this frame, or None.

    Codon i (0-based) covers nt [frame + 3i, frame + 3i + 3). A stretch
    spanning codons [c0, c1) covers nt [frame + 3*c0, frame + 3*c1) and
    crosses iff frame + 3*c0 < breakpoint_nt < frame + 3*c1.
    """
    usable = contig_seq[frame:]
    usable = usable[: (len(usable) // _CODON) * _CODON]
    aa = str(Seq(usable).translate())
    c0 = 0
    for i, residue in enumerate(list(aa) + ["*"]):  # sentinel flushes last run
        if residue == "*":
            if i > c0:
                nt_start = frame + _CODON * c0
                nt_end = frame + _CODON * i
                if nt_start < breakpoint_nt < nt_end:
                    return aa[c0:i]
            c0 = i + 1
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py -v`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/orf_fasta_from_contigs.py workflow/tests/test_orf_fasta_from_contigs.py
git commit -m "feat(translation): breakpoint-crossing ORF-stretch extraction core"
```

---

### Task 2: FASTA emitter (I/O, header reformat, filters, CLI)

**Files:**
- Modify: `workflow/scripts/orf_fasta_from_contigs.py`
- Test: `workflow/tests/test_orf_fasta_from_contigs.py`

**Interfaces:**
- Consumes: `crossing_orf_stretch` (Task 1).
- Produces: `emit_orf_fasta(contigs_fasta, output_fasta, flank_nt: int = 30, min_peptide_len: int = 8, frames=(0, 1, 2)) -> None` - reads a contig FASTA, writes the ORF-stretch FASTA; `_orf_header(contig_header: str, frame: int) -> str` - reformats `{junc}|{coords}|{sample}` to `{junc}|{frame}|{coords}`.

- [ ] **Step 1: Write the failing tests**

```python
# append to workflow/tests/test_orf_fasta_from_contigs.py
from orf_fasta_from_contigs import _orf_header, emit_orf_fasta


def test_orf_header_reformat():
    h = "JUNC1|chr7:100-200:+|tumor"
    assert _orf_header(h, 2) == "JUNC1|2|chr7:100-200:+"


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def test_emit_drops_x_and_short_and_writes_headers(tmp_path):
    contigs = tmp_path / "contigs.fa"
    out = tmp_path / "orf.fa"
    # J1: clean 20-Lys crossing stretch in frame 0.
    clean = "".join(["AAA"] * 20)
    # J2: contains an ambiguous base (N) inside a crossing run -> 'X' -> dropped.
    ambig = "".join(["AAA"] * 9) + "ANA" + "".join(["AAA"] * 10)
    _write_fasta(contigs, [
        ("J1|chr1:10-20:+|tumor", clean),
        ("J2|chr2:30-40:-|tumor", ambig),
    ])
    emit_orf_fasta(contigs, out, flank_nt=30, min_peptide_len=8)
    text = out.read_text()
    assert ">J1|0|chr1:10-20:+" in text
    assert "K" * 20 in text
    assert "J2|" not in text  # X-containing stretch dropped
    # headers unique
    headers = [ln for ln in text.splitlines() if ln.startswith(">")]
    assert len(headers) == len(set(headers))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py -v`
Expected: FAIL with `ImportError: cannot import name '_orf_header'`.

- [ ] **Step 3: Write the minimal implementation**

```python
# add to workflow/scripts/orf_fasta_from_contigs.py

def _parse_fasta(fasta_path):
    records = []
    header = None
    parts = []
    with Path(fasta_path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def _orf_header(contig_header, frame):
    # contig header: {junction_id}|{chrom}:{start}-{end}:{strand}|{sample_type}
    fields = contig_header.split("|")
    junction_id, coords = fields[0], fields[1]
    return f"{junction_id}|{frame}|{coords}"


def emit_orf_fasta(contigs_fasta, output_fasta, flank_nt=30, min_peptide_len=8,
                   frames=(0, 1, 2)):
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    records = _parse_fasta(contigs_fasta)

    n_written = n_drop_x = n_drop_short = n_none = 0
    with output_fasta.open("w") as out:
        for header, seq in records:
            seq = seq.upper()
            if len(seq) != 2 * flank_nt:
                log.warning("Skipping %s: length %d != %d", header, len(seq),
                            2 * flank_nt)
                continue
            for frame in frames:
                stretch = crossing_orf_stretch(seq, frame, breakpoint_nt=flank_nt)
                if stretch is None:
                    n_none += 1
                    continue
                if "X" in stretch:
                    n_drop_x += 1
                    continue
                if len(stretch) < min_peptide_len:
                    n_drop_short += 1
                    continue
                out.write(f">{_orf_header(header, frame)}\n{stretch}\n")
                n_written += 1

    log.info(
        "ORF FASTA: %d written, %d dropped (X), %d dropped (<%d aa), "
        "%d frames with no crossing stretch",
        n_written, n_drop_x, n_drop_short, min_peptide_len, n_none,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py -v`
Expected: PASS (6 tests).

- [ ] **Step 5: Add the CLI / Snakemake entry point**

```python
# add to workflow/scripts/orf_fasta_from_contigs.py

def _snakemake_main():
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))
    emit_orf_fasta(
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_fasta=snakemake.output.orf_fasta,  # type: ignore[name-defined]  # noqa: F821
        flank_nt=snakemake.params.flank_nt,  # type: ignore[name-defined]  # noqa: F821
        min_peptide_len=snakemake.params.min_peptide_len,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main():
    p = argparse.ArgumentParser(description="Emit breakpoint-crossing ORF-stretch FASTA.")
    p.add_argument("--contigs-fasta", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--flank-nt", type=int, default=30)
    p.add_argument("--min-peptide-len", type=int, default=8)
    args = p.parse_args()
    emit_orf_fasta(args.contigs_fasta, args.output, args.flank_nt, args.min_peptide_len)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
```

- [ ] **Step 6: Verify CLI runs on a fixture**

Run: `workflow/tests/.venv/bin/python workflow/scripts/orf_fasta_from_contigs.py --help`
Expected: argparse usage printed, exit 0.

- [ ] **Step 7: Commit**

```bash
git add workflow/scripts/orf_fasta_from_contigs.py workflow/tests/test_orf_fasta_from_contigs.py
git commit -m "feat(translation): ORF-FASTA emitter I/O, header reformat, X/min-len filters, CLI"
```

---

### Task 3: Snakemake wiring (config + rules, decoupled 30-nt path)

**Files:**
- Modify: `config/config.yaml`
- Create: `workflow/rules/ms_search_db.smk`
- Modify: `workflow/Snakefile`

**Interfaces:**
- Consumes: `assemble_contigs.py` (existing, reused at 30-nt params), `orf_fasta_from_contigs.py` (Tasks 1-2).
- Produces: an optional target `results/{sample_dir}/ms_search_db/junction_orf.fasta`.

- [ ] **Step 1: Add the config block**

```yaml
# append under an appropriate section in config/config.yaml
# -----------------------------------------------------------------
# MS search database (optional; junction-spanning ORF FASTA for #1176)
# -----------------------------------------------------------------
ms_search_db:
  flank_nt: 30            # 3 * (11 - 1); all class-I 8-11mers span the breakpoint
  min_peptide_len: 8      # drop crossing stretches shorter than this
```

- [ ] **Step 2: Write the rule file**

Inspect an existing rule for the exact wildcard/path pattern first:
Run: `sed -n '1,60p' workflow/rules/assemble_contigs.smk`

Then create `workflow/rules/ms_search_db.smk`, mirroring the `assemble_contigs` rule's input/output/wildcard conventions (adjust the paths to match what `assemble_contigs.smk` actually uses):

```python
# workflow/rules/ms_search_db.smk
# Decoupled 30-nt contig assembly (reuses assemble_contigs.py) feeding the
# ORF-stretch FASTA emitter. Optional target; not part of `rule all`.

rule assemble_contigs_wide:
    input:
        novel_junctions=<same as assemble_contigs rule>,
        genome_fasta=<same as assemble_contigs rule>,
    output:
        contigs_fasta="results/{sample_dir}/ms_search_db/contigs_wide.fa",
    params:
        upstream_nt=config["ms_search_db"]["flank_nt"],
        downstream_nt=config["ms_search_db"]["flank_nt"],
    log:
        "logs/{sample_dir}/ms_search_db/assemble_contigs_wide.log",
    conda:
        "../envs/<same env as assemble_contigs>.yaml"
    script:
        "../scripts/assemble_contigs.py"


rule emit_orf_fasta:
    input:
        contigs_fasta="results/{sample_dir}/ms_search_db/contigs_wide.fa",
    output:
        orf_fasta="results/{sample_dir}/ms_search_db/junction_orf.fasta",
    params:
        flank_nt=config["ms_search_db"]["flank_nt"],
        min_peptide_len=config["ms_search_db"]["min_peptide_len"],
    log:
        "logs/{sample_dir}/ms_search_db/emit_orf_fasta.log",
    conda:
        "../envs/<same env as translate_peptides>.yaml"
    script:
        "../scripts/orf_fasta_from_contigs.py"
```

- [ ] **Step 3: Include the rule file**

Add to `workflow/Snakefile` next to the other `include:` lines:

```python
include: "rules/ms_search_db.smk"
```

- [ ] **Step 4: Dry-run to verify the DAG parses and the target resolves**

Run (note the single `--configfile` and the `--` before the target, per the Snakemake 8 gotchas):
```bash
conda activate snakemake
snakemake -n --configfile config/test_config.yaml -- results/test/patient_001/ms_search_db/junction_orf.fasta
```
Expected: a dry-run plan showing `assemble_contigs_wide` then `emit_orf_fasta`, no `NameError`/`FileNotFoundError`. Adjust the `results/...` wildcard path to match the test config's sample-dir convention.

- [ ] **Step 5: Render the DAG (memory rule for rule changes)**

Run: `bash scripts/visualize_dag.sh` (after the dry-run passes). Confirm the two new rules appear.

- [ ] **Step 6: Commit**

```bash
git add config/config.yaml workflow/rules/ms_search_db.smk workflow/Snakefile
git commit -m "feat(translation): optional ms_search_db rule - 30-nt contigs + ORF-FASTA emit"
```

---

### Task 4: Structural validation + multi-chromosome capability check (AC-1, AC-3 structural)

**Files:**
- Create: `workflow/tests/data/orf_fasta/contigs_multichr.fa`
- Create: `workflow/tests/data/orf_fasta/canonical_tiny.fasta`
- Modify: `workflow/tests/test_orf_fasta_from_contigs.py`

**Interfaces:**
- Consumes: `emit_orf_fasta` (Task 2).

- [ ] **Step 1: Create the multi-chromosome contig fixture**

Junctions on three different chromosomes prove the path is not chr22-scoped. Each contig is 60 nt with a clean crossing ORF in at least one frame:

```
# workflow/tests/data/orf_fasta/contigs_multichr.fa
>JCHR1|chr1:1000-2000:+|tumor
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>JCHR7|chr7:5000-6000:-|tumor
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
>JCHRX|chrX:9000-9500:+|tumor
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

```
# workflow/tests/data/orf_fasta/canonical_tiny.fasta
>sp|TEST1|canonical decoy target
MKLPQRSTVWYACDEFGHIKLMNPQRST
```

- [ ] **Step 2: Write the failing structural test**

```python
# append to workflow/tests/test_orf_fasta_from_contigs.py
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def test_multichr_capability_and_structural_validity(tmp_path):
    data = Path(__file__).parent / "data" / "orf_fasta"
    out = tmp_path / "orf.fa"
    emit_orf_fasta(data / "contigs_multichr.fa", out, flank_nt=30, min_peptide_len=8)
    records = _parse_fasta(out)  # reuse the emitter's parser

    # AC-1: entries from all three chromosomes (not chr22-scoped).
    chroms = {h.split("|")[2].split(":")[0] for h, _ in records}
    assert {"chr1", "chr7", "chrX"} <= chroms

    # Structural: unique headers, valid alphabet, no stop/X.
    headers = [h for h, _ in records]
    assert len(headers) == len(set(headers))
    for _, seq in records:
        assert set(seq) <= _VALID_AA

    # Concatenation with the canonical proteome yields no duplicate IDs.
    canon = _parse_fasta(data / "canonical_tiny.fasta")
    all_ids = [h.split()[0] for h, _ in records + canon]
    assert len(all_ids) == len(set(all_ids))
```

- [ ] **Step 3: Run to verify it fails, then passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py::test_multichr_capability_and_structural_validity -v`
Expected first: FAIL if fixtures/paths are off; fix fixture content until PASS. (The emitter already exists, so this validates fixtures + structural invariants, not new code.)

- [ ] **Step 4: Run the full test module**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_orf_fasta_from_contigs.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Commit**

```bash
git add workflow/tests/data/orf_fasta/ workflow/tests/test_orf_fasta_from_contigs.py
git commit -m "test(translation): multi-chr capability + structural validity for ORF FASTA"
```

---

### Task 5: Sage smoke harness (AC-3 live)

**Files:**
- Create: `research/experiments/issue_1204_orf_fasta/README.md`
- Create: `research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh`
- Create: `research/experiments/issue_1204_orf_fasta/sage_config.json`
- Create: `research/experiments/issue_1204_orf_fasta/synthetic.mgf` (or `.mzML` if MGF unsupported by the pinned Sage)

**Interfaces:**
- Consumes: an emitted `junction_orf.fasta` (Task 3 output) + `workflow/tests/data/orf_fasta/canonical_tiny.fasta` (Task 4).

- [ ] **Step 1: Write a minimal nonspecific Sage config**

```json
{
  "database": {
    "enzyme": { "cleave_at": "" },
    "fragment_min_mz": 100.0,
    "fragment_max_mz": 2000.0,
    "peptide_min_len": 8,
    "peptide_max_len": 11,
    "generate_decoys": true,
    "fasta": "combined.fasta"
  },
  "precursor_tol": { "ppm": [-20.0, 20.0] },
  "fragment_tol": { "ppm": [-20.0, 20.0] },
  "output_directory": "sage_out"
}
```

- [ ] **Step 2: Write a synthetic spectra file**

Author a minimal one-spectrum `synthetic.mgf` whose precursor/fragments correspond to a peptide present in the emitted FASTA (a substring of a `JCHR*` stretch), so the search has something to match. Keep it to a single `BEGIN IONS ... END IONS` block. (If the pinned Sage build rejects MGF, convert once to `synthetic.mzML` via `pyteomics`/`psims` and commit that instead; the loadability proof does not depend on the format.)

- [ ] **Step 3: Write the run script (bash-3.2-safe)**

```bash
#!/bin/bash
# research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh
# One-time local loadability smoke: does Sage load the junction-partition FASTA
# concatenated with the canonical proteome and complete a nonspecific search?
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ORF_FASTA="${1:?path to emitted junction_orf.fasta}"
CANON="$HERE/../../../workflow/tests/data/orf_fasta/canonical_tiny.fasta"
cat "$CANON" "$ORF_FASTA" > "$HERE/combined.fasta"
# sage must be on PATH (brew install sage-proteomics, or the github.com/lazear/sage release binary)
sage "$HERE/sage_config.json" "$HERE/synthetic.mgf"
echo "SAGE SMOKE OK"
```

- [ ] **Step 4: Write the README (procedure + what "pass" means)**

```markdown
# Issue #1204 - Sage loadability smoke

One-time pre-merge verification that the emitted ORF-stretch FASTA is consumable
by a nonspecific MS search (AC-3). NOT run in CI; the CI-side check is the
structural test in `workflow/tests/test_orf_fasta_from_contigs.py`.

## Run

    conda activate snakemake
    snakemake --configfile config/test_config.yaml -- results/test/patient_001/ms_search_db/junction_orf.fasta
    bash research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh \
        results/test/patient_001/ms_search_db/junction_orf.fasta

## Pass criteria

Sage exits 0, writes `sage_out/`, and logs the loaded target+decoy count without
a FASTA parse error. Recovering the synthetic peptide is a bonus, not required;
the deliverable is loadability, not a hit (real hits are #1176 AC-6).
```

- [ ] **Step 5: Make the script executable and commit**

```bash
chmod +x research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh
git add research/experiments/issue_1204_orf_fasta/
git commit -m "test(translation): Sage loadability smoke harness for the ORF FASTA (#1204 AC-3)"
```

- [ ] **Step 6: Run the smoke once (manual, pre-merge; record in the PR Test plan)**

Install Sage locally (arm64 macOS binary), run the README procedure, and paste the `SAGE SMOKE OK` result into the PR Test plan. This is the live-integration-smoke evidence.

---

### Task 6: Integration run + PR prep

**Files:** none (verification + PR).

- [ ] **Step 1: Run the chr22 integration for the new rule (memory rule for new Snakemake rules)**

```bash
conda activate snakemake
snakemake --cores 4 --use-conda --configfile config/test_config.yaml -- results/test/patient_001/ms_search_db/junction_orf.fasta
```
Expected: `assemble_contigs_wide` + `emit_orf_fasta` run green; `junction_orf.fasta` is non-empty and its headers match `{junc}|{frame}|{coords}`.

- [ ] **Step 2: Run all five CI pytest dirs locally (memory: pytest does not recurse)**

```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests tools/ci tools/news tools/project_map scripts/tests -q
```
Expected: all pass.

- [ ] **Step 3: Open the PR**

Follow the PR-open 4-step checklist (project add, `**Created by:** Developer`, Status -> Ready for review, mirror on review request). Test plan lists: unit tests, structural test, DAG render, chr22 integration output, and the pasted Sage smoke result. Body carries the closure-ritual AC ticks for [#1204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1204).

---

## Self-Review

**Spec coverage:**
- Goal / why-new-emitter -> Tasks 1-2. Architecture (reuse stage 1, new stage 2) -> Tasks 2-3. ORF-stretch algorithm -> Task 1. Header format -> Task 2. Fork 1-A decoupled 30-nt -> Task 3 (`assemble_contigs_wide` distinct output, `peptide_lengths` untouched). Fork 2-both -> Task 4 (structural CI) + Task 5 (Sage smoke). Genome-wide capability (AC-1) -> Task 4 multi-chr fixture. Testing -> Tasks 1,2,4,5. Out-of-scope respected (no MHC-path change; #1176/#1255 untouched). All covered.
- AC-1 (genome-wide FASTA from tumor_exclusive) -> Task 3 reuses `assemble_contigs.py`'s tumor_exclusive filter + Task 4 multi-chr capability. AC-2 (30-nt flanks, traceable headers) -> Task 2 header + Task 3 flank_nt=30. AC-3 (consumable by no-enzyme MS search) -> Tasks 4+5.

**Placeholder scan:** the rule file (Task 3 Step 2) has intentional `<same as ...>` markers for the input paths + conda env, because the exact wildcard/path convention must be read from `assemble_contigs.smk` at implementation time (Step 2 begins with that `sed`); this is a read-then-fill instruction, not an unresolved placeholder. All code steps otherwise contain complete content.

**Type consistency:** `crossing_orf_stretch(contig_seq, frame, breakpoint_nt=30)` and `emit_orf_fasta(contigs_fasta, output_fasta, flank_nt=30, min_peptide_len=8, frames=(0,1,2))` and `_orf_header(contig_header, frame)` are used consistently across Tasks 1, 2, 4. The Snakemake `params.flank_nt` feeds `emit_orf_fasta(flank_nt=...)` which is passed as `breakpoint_nt` to the core - consistent (symmetric flank means breakpoint == flank_nt).
