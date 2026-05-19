# Issue #17 — Switch production aligner to STAR (minimum-viable migration)

**Date:** 2026-05-19
**Issue:** [Issue #17](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17) (STAR aligner)
**Author:** Developer

## Goal

Switch the production aligner default from HISAT2 to STAR, with the minimum-viable set of STAR parameter changes that make the switch correctness-equivalent to the Nature paper's command. Defer paper-fidelity sensitivity tuning to a follow-up Issue with per-flag benchmarks.

## Motivation

STAR consistently outperforms HISAT2 for novel/unannotated splice junction detection — the critical step for tumor-neoepitope discovery. The current `star_align` rule exists and is used in PoC runs, but the production default is `hisat2` and the existing STAR command silently drops multi-mappers and runs single-pass — both deviations from the Nature paper's reference recipe that materially affect novel-junction recall.

## Scope decision — middle ground, not full paper-fidelity

The [2026-04-21 implementation comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17#issuecomment) on the Issue lists ~14 parameter changes vs the paper's command. We deliberately limit THIS PR to the changes that are **correctness issues** rather than tunables:

| Change | Reason for inclusion / exclusion |
|---|---|
| `aligner: "hisat2"` → `"star"` in `config/config.yaml` | The headline change; ships the switch |
| Add `--twopassMode Basic` | THE flag that gives STAR its novel-junction edge — without it, switching aligners is mostly cosmetic |
| Add `--limitSjdbInsertNsj 2000000` | Defensive raise from default 1M; 2-pass + multi-mappers can hit the cap and cause a fatal abort mid-run |
| Remove `--outSJfilterReads Unique` | Silently drops multi-mapper-derived junctions vs paper |
| Remove `--outSJfilterCountUniqueMin 1 1 1 1` | Silently throttles output; paper uses defaults |
| Remove `--outSJfilterCountTotalMin 1 1 1 1` | Same |
| (deferred) Sensitivity flags: `--outFilterMatchNminOverLread 0.33`, `--alignIntronMax 500000`, `--alignSJDBoverhangMin 1`, `--outFilterMultimapNmax 20`, `--alignMatesGapMax 1000000` | Tunables — each should be justified with a benchmark, not copy-pasted |
| (deferred) Remove `--sjdbGTFfile` from `star_index` | Stylistic; GTF-in-index doesn't hurt novel-junction discovery in 2-pass mode |
| (deferred) BAM emission (`--outSAMtype BAM SortedByCoordinate` + SAM-formatting flags) | Downstream consumes `SJ.out.tab`, not BAM; current `--outSAMtype None` saves disk + memory |
| (deferred) Pre-built STAR index download | No widely-standard GRCh38 STAR index URL with matching settings; build-from-genome works |

The 2026-04-21 comment's point 3 ("switch from awk parsing SJ.out.tab to regtools on BAM") is **already obsolete** — [PR #402](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402) wired `workflow/scripts/star_sj_to_junctions.py` to parse `SJ.out.tab` directly with strand rescue.

## Changes — four files

### 1. `config/config.yaml`
Single-line change:
```yaml
alignment:
  aligner: "star"   # was "hisat2"
```
Comment in the file already documents both options; no doc edit needed.

### 2. `config/test_config.yaml`
`test_config.yaml` inherits `aligner` from `config.yaml` via Snakemake's recursive merge. Without an explicit override, the M1 chr22 local dev environment would inherit `star` after change 1 and attempt to build a 32 GB STAR index on 8 GB RAM. Add an explicit pin to the `alignment:` block:

```yaml
alignment:
  aligner: "hisat2"   # explicit: M1 chr22 dev can't fit STAR's 32 GB index
  hisat2_index_dir: "resources/test/hisat2_index"
  ...
```

### 3. `workflow/rules/alignment.smk` — `star_align` rule

In the shell block at [workflow/rules/alignment.smk:338-349](../../../workflow/rules/alignment.smk#L338-L349), modify the STAR command:

```diff
 STAR \
     --runMode alignReads \
     --runThreadN {threads} \
     --genomeDir {input.index_dir} \
     --readFilesIn $FASTQ_FILES \
     $READCMD \
     --outFileNamePrefix {params.output_prefix} \
     --outSAMtype None \
-    --outSJfilterReads Unique \
-    --outSJfilterCountUniqueMin 1 1 1 1 \
-    --outSJfilterCountTotalMin 1 1 1 1 \
+    --twopassMode Basic \
+    --limitSjdbInsertNsj 2000000 \
     2>&1 | tee {log}
```

### 4. `workflow/tests/test_alignment_star_command.py` (new file)

Programmatic snapshot test of the rendered `star_align` shell command. Approximate shape:

```python
import subprocess
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

@pytest.fixture
def star_dry_run_output(tmp_path):
    # Write a minimal override config that pins aligner=star and points
    # at a stub samples.tsv with one patient_id + one sample_id.
    # Run: snakemake -n --printshellcmds --configfile config/test_config.yaml <override> -- <star_align target>
    # Capture stdout, return as a single string for substring assertions below.
    ...

def test_star_align_uses_twopass_mode(star_dry_run_output):
    assert "--twopassMode Basic" in star_dry_run_output

def test_star_align_raises_sjdb_insert_limit(star_dry_run_output):
    assert "--limitSjdbInsertNsj 2000000" in star_dry_run_output

def test_star_align_does_not_filter_to_unique(star_dry_run_output):
    assert "--outSJfilterReads Unique" not in star_dry_run_output

def test_star_align_does_not_throttle_min_counts(star_dry_run_output):
    assert "--outSJfilterCountUniqueMin" not in star_dry_run_output
    assert "--outSJfilterCountTotalMin" not in star_dry_run_output
```

This is the first shell-level rule snapshot test in the codebase — existing tests cover pure Python helpers (`test_strandness.py`, `test_bed12_to_junctions.py`, etc.). The pattern is worth establishing: `star_sj_to_junctions.py` from [PR #402](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402) will eventually want similar rule-level coverage.

## Verification

### In-PR (gates merge)
- `pytest` (existing CI check) — new `test_alignment_star_command.py` passes
- `snakemake -n` with both `aligner: hisat2` and `aligner: star` (existing CI check) — parse-time validation of the new flags

chr22 cannot exercise STAR end-to-end: its 32 GB index won't fit on the M1 8 GB local dev environment. End-to-end STAR verification happens only on cloud production VMs.

### Post-merge (via [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378))
[Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) is already queued as a full patient_002 production re-run to verify the [PR #402](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402) strand fix. Adding the new STAR params on top means **one** cloud run verifies both changes simultaneously — saves duplicate cloud cost vs running a separate single-sample smoke test as part of THIS PR.

This PR's wrap-up adds three explicit ACs to [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378):

1. `SJ.out.tab` contains rows with `multi_reads > 0 ∧ unique_reads == 0` — confirms multi-mapper inclusion is wired through end-to-end
2. 2-pass log line present in `star_align.log` (e.g. `STAR: 2nd pass:`)
3. No `limitSjdbInsertNsj exceeded` fatal in stderr — confirms the 2M raise was sufficient

Adding these is a procedural follow-on, not a scope-growth of [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378).

## Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Bad STAR flag syntax causes startup abort | low | `snakemake -n` parse + new pytest catch at PR time |
| 2-pass mode wall time roughly 2× single-pass | certain | Accepted — this IS the value of 2-pass; cloud runs are not time-critical |
| [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) finds an unexpected regression after merge | low | Recoverable: flip `aligner` back to `hisat2` in config OR remove `--twopassMode Basic` alone to isolate which change regressed |
| Novel-junction count explodes and exceeds `limitSjdbInsertNsj 2000000` | very low | Paper authors bumped to 2M as defensive insurance; would surface as a fatal log line, easily caught and raised further if needed |

## Follow-up Issue to file (after THIS PR merges)

Title: *"STAR sensitivity tuning — benchmark Nature paper command flags individually"*
Priority: P2
Size: M
Scope: Per-flag benchmark on patient_002 `SJ.out.tab` (junction count, % annotated, runtime) for each deferred flag — adopt the flag if its delta justifies its complexity. Avoids cargo-culting the Nature paper recipe wholesale.

Flags to evaluate:
- `--outFilterMatchNminOverLread 0.33` (vs default 0.66)
- `--outFilterScoreMinOverLread 0.33` (vs default 0.66)
- `--outFilterMultimapNmax 20` (vs default 10)
- `--alignIntronMax 500000` (vs STAR default)
- `--alignMatesGapMax 1000000` (vs STAR default)
- `--alignSJDBoverhangMin 1` (vs default 3)
- `--sjdbGTFfile` removal from `star_index` (vs keep)
- BAM emission via `--outSAMtype BAM SortedByCoordinate` (vs current None)

## Out of scope

- Any change to the HISAT2 alignment rule (params, behavior)
- Pre-built STAR index download / sourcing
- Pipeline-wide BAM consumers (no current downstream rule consumes the alignment BAM beyond the HISAT2 path itself)
