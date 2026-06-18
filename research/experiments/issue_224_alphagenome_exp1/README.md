# Issue #224 — AlphaGenome predicted-normal validation (Exp 1, patient_001)

**Issue:** [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)
**Parent:** [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) — rethink normal filtering with population panel + AlphaGenome fallback (Experiment 1)
**Status:** complete — GREEN-with-caveats verdict (chr22 PoC, 2026-05-18; metrics corrected post-[PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398) review)
**Ship PR:** [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398)

## Goal

Test whether AlphaGenome predicts known tissue-expressed splicing from the GRCh38 reference alone (no patient variants). Run AlphaGenome on chr22 for gastric/stomach tracks and compare its predicted-junction set against the gold-standard matched-normal RNA-seq for patient_001. Outputs the Exp 1 row of #203's decision-rule table.

## Inputs

- **Matched-normal junctions:** `results/patient_001_test/alignment/SRR9143065_test/junctions.tsv` (produced by the chr22 test pipeline; not tracked in git)
- **GENCODE chr22 annotation:** `resources/test/chr22.gtf.gz` (built by `scripts/prepare_test_data.sh`)
- **AlphaGenome API:** `ALPHAGENOME_API_KEY` in `.env` (env `workflow/envs/alphagenome.yaml`, shipped via [PR #386](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/386))

## Outputs

- `outputs/chr22_stomach_predicted_junctions.parquet` — dedup'd AlphaGenome chr22 stomach-track predictions (2.6 M junctions, ~17 MB). **Not tracked in git** (API-derived; regenerate via the notebook §4 sweep cell — delete the parquet to force a fresh sweep). Cross-experiment dependency for #225 and #393.

## Cross-experiment deps

- **Consumed by [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)** (filter-strength comparison) and **[Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393)** (slide-deck PoC), both reading `outputs/chr22_stomach_predicted_junctions.parquet` via explicit path reference per the cross-experiment-sharing convention in CLAUDE.md (single producer here; multiple read-only consumers).

## Headline result (chr22 PoC)

Framing 1 (annotated-chr22-introns universe; 7,731 universe / 259 positives / 7,472 negatives): **AP = 0.214** (6.3× random baseline), **best F1 = 0.300** at τ = 3.16, **recall = 0.405** (the honest headline — depth confounder doesn't bias recall). Verdict: AG is a viable *secondary* evidence stream, not a standalone matched-normal replacement at this scale. See notebook §6 for the full caveat list (chr22-only, 500K-read depth, single stomach track, annotated-only ground truth, n=259 bootstrap CI [0.258, 0.333]).

## Caveats

- **chr22 PoC scope** — full-genome generalisation unvalidated; expect different runtime + AG cost.
- **Experiment 2 (germline-aware AG) deferred** to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS acquisition.
- **Exp 3 next step** — comparative filter strength is [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225).

## Conda env

`splice-neoepitope-alphagenome` (defined in `workflow/envs/alphagenome.yaml`).
