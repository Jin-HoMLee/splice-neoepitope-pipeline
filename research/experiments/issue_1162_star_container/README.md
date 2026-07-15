# Issue #1162 - containerized linux-64 STAR path, chr22 end-to-end

Experiment record for the STAR-via-container local path (ADR-0003) and its AC4
STAR-vs-HISAT2 characterization. Fixture: chr22 test pair (SRR9143066 tumor /
SRR9143065 normal), 500K single-end reads each.

## The fix (AC1/AC2): vfork ENOMEM, not the sjdb cap

The 2026-07-14 run died `Failed vforking readFilesCommand / 12: Cannot allocate
memory` at the **start of 2nd-pass mapping**, after junction insertion had already
*succeeded*. Every subprocess spawn (even STAR's internal `ls`) failed - a
fork-accounting failure, not the `--limitSjdbInsertNsj 2000000` array (the log
shows only ~2.1 Mbp / ~10K junctions were actually inserted). The colima `vz` VM
has **no swap**, so a 4 GB VM refuses the fork's memory reservation under Linux'
default heuristic overcommit, even though a `vfork` shares the address space and
immediately execs tiny `zcat`.

Fix: `vm.overcommit_memory=1` in the VM, wired into `scripts/run_local_linux64.sh`
(idempotent, every run; host-RAM-neutral - no `--memory` bump that would starve
the 8 GB host). Controlled result: same command / inputs / 4 GB VM as the failing
run, only overcommit changed → STAR cleared 2nd-pass on both samples and the DAG
completed `3 of 3 steps (100%)`.

Counts (knob-off = default, count-all):

| | raw junctions | tumor_exclusive |
|---|---|---|
| STAR   | tumor 3379 / normal 3105 | 319 |
| HISAT2 | tumor 1872 | 122 |

## AC4 - STAR-vs-HISAT2 difference, attributed to the aligner not the semantic

Two comparisons; artifacts under `star_off/`, `star_on/`, `hisat2_off/`.

### 1. Knob matched-pair (STAR off vs on) - the knob IS wired

Regenerated from the same `SJ.out.tab` with/without `--unique-only`:

- tumor: knob-off 3379 → knob-on 1581; the **1798** dropped are exactly the
  multimapper-only junctions (`supported only by multimappers: 1798 (dropped)`),
  knob-on ⊂ knob-off (on-only = 0).
- of the 1581 shared, **254** carry higher counts under knob-off (multimapper
  reads added; e.g. `chr22:42926812:43177170:-` 1612 vs 296, +1316), 1327 are
  pure-unique and unchanged.

Confirms AC4's predicted values in both membership and per-junction counts. The
two runs are **not** identical, so the knob is wired (not a no-op).

### 2. Aligner difference (STAR-off vs HISAT2-off, identical semantic)

Both knob-off, so read-support semantic is held fixed and any difference is the
aligner.

- **No coordinate artifact**: 499 exact raw matches, a single (-2,0) near-miss -
  no systematic off-by-one. The divergence is real aligner biology, not a
  BED12/SJ.out.tab coordinate-convention bug.
- **Confidence-dependent overlap** (raw tumor, shared as % of HISAT2):
  reads≥1 26% · reads≥2 41% · reads≥3 60% · reads≥5 66%. Agreement is strong on
  high-confidence junctions; disagreement lives in the single-read noise floor,
  where calls are inherently aligner-specific.
- STAR finds ~1.8x more raw junctions (two-pass novel-junction sensitivity).
- **Classified tumor_exclusive overlap is only 4/122**, far below the raw 499:
  108 of 122 HISAT2-TE junctions are reads=2-level calls absent from STAR's set
  entirely. The TE set is dominated by low-count junctions, so the aligners' TE
  sets barely overlap - a downstream amplification of low-count instability, not
  a semantic difference.

### Caveat

The TE comparison sits downstream of the pipeline's mean-reads depth gate
(Issue #1161, Scientist-owned) and the low-count instability interacts with it.
Which low-count junctions to trust is the Scientist's call (#1122/#1161), out of
scope for #1162, which only had to characterize the difference and attribute it.

## Reproduce

```bash
# STAR path (container):
bash scripts/run_local_linux64.sh -- --cores 4 --use-conda --rerun-triggers mtime \
  --configfile config/test_config.yaml config/test_star_config.yaml \
  -- results/patient_001_test/junctions/novel_junctions.tsv

# knob matched-pair from an existing SJ.out.tab:
python3 workflow/scripts/star_sj_to_junctions.py \
  --input <sample>_SJ.out.tab --output on.tsv --unique-only

# comparisons:
python3 research/experiments/issue_1162_star_container/compare_junctions.py --help
```
