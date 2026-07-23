# Normal-filter decision for the Courcelles cohort (Issue #1176)

**Date:** 2026-07-23
**Author:** Scientist
**Issue:** [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176), the non-canonical public-spectra re-search (best-bet 3, epic [#681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681))
**Status:** decided, with a stated cost and a named re-open trigger

## The problem

Our junction-origin classification filters tumor junctions against a **matched normal**.
The Courcelles CRC cohort has none: all 26 RNA-seq runs in `PRJNA1372708` / `GSE312236` carry the sample title `Flash-frozen biopsy`, and there is no normal arm in the deposit (verified against the ENA run table, 2026-07-23).

With no normal, `classify_junctions` labels every surviving unannotated tumor junction `tumor_exclusive` and emits a warning, leaving the **GTEx pan-tissue population filter as the sole specificity control**.
That matters beyond taxonomy: an inflated `tumor_exclusive` set inflates the search database, and database inflation is the documented mechanism that degrades group-specific FDR, which is this Issue's central methodological risk.

The authors solved the same problem differently, with **eight mTEC controls** plus relaxed k-mer occurrence filtering.

## Decision

**Proceed with the GTEx pan-tissue filter as the sole normal control, and report the substitution explicitly as a limitation.**
Do not attempt to synthesise a pseudo-normal from the cohort, and do not treat `tumor_exclusive` on this cohort as equivalent to `tumor_exclusive` on a matched-pair patient.

## Evidence

A matched-pair control on data where we **do** have a matched normal (the chr22 gastric tumor/normal pair), running the pipeline's own `classify_junctions` twice with exactly one variable flipped.
Regenerate with `normal_vs_gtex_control.py` in this directory; raw numbers in `outputs/normal_vs_gtex.json`.

| category | A: matched normal + GTEx | B: GTEx only (the Courcelles case) |
|---|---|---|
| `junctions_raw` | 1872 | 1872 |
| `mean_reads_filtered` | 1638 | 1638 |
| `annotated_discarded` | 79 | 79 |
| `normal_shared` | 4 | 0 |
| `gtex_pantissue_shared` | 29 | 32 |
| **`tumor_exclusive`** | **122** | **123** |

Both funnels reconcile against `junctions_raw` (asserted in the script, not eyeballed).

**The result: of the 4 junctions the matched normal removes, GTEx independently blacklists 3.**
Exactly one is lost, `chr22:35501382:35584228:+` at 5 reads, and the candidate set inflates by **+1 junction, 0.8%**.

This check could have failed and would have been informative either way.
Had GTEx recovered none of the four, the population filter would have been shown to be blind to precisely the class a matched normal exists to catch, and the honest response would have been to refuse the cohort.

## Why the loss is small, and what it is

A population reference is structurally blind to **patient-private** non-tumor junctions: germline variants creating private splice sites, or individual-specific expression that no GTEx donor happens to share.
A matched normal sees exactly those.
So the residual is not random noise, it is a specific and predictable class, and `chr22:35501382:35584228:+` is an instance of it.

Two structural points worth carrying forward:

1. **GTEx is applied *after* normal subtraction in our code.**
   The `gtex_pantissue_shared` counts in our ordinary matched-pair runs are therefore only GTEx's *marginal* contribution, and they systematically **understate** what it does when operating alone.
   Reading those historical counts as an estimate of standalone GTEx performance would have been wrong; that is why this measurement flips the arm rather than reusing an existing run.
2. **On this data GTEx does roughly seven times more filtering work than the matched normal.**
   Of the 234 junctions surviving the depth gate, the matched normal removes 4 (1.7%) and GTEx removes 29 (12.4%).
   Note this is a within-run-A partition, so the 29 is GTEx's *marginal* count with the normal also active; standalone in run B it is 32, per point 1 above.
   The matched normal is the finer instrument, not the heavier one.

## Honest limitations

These are the reasons this decision is "proceed with a stated cost" rather than "proceed, the substitution is free".

- **n = 4.**
  The entire matched-normal contribution on this data is four junctions.
  A 3-of-4 recovery rate is directionally reassuring and is **not** a calibrated estimate; the binomial interval on 3/4 spans most of the unit interval.
  Do not quote "GTEx recovers 75% of matched-normal hits" as a number.
- **The depth gate dominates everything.**
  The per-file mean-reads filter removes 1638 of 1872 junctions (87.5%) *before* any normal or GTEx logic runs.
  At this depth, specificity is overwhelmingly set by the depth gate, not by the normal-versus-GTEx choice.
  This interacts directly with [#1161](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1161) (the floating mean-depth gate), and the measurement here describes behaviour on the surviving 12.5%.
- **Wrong tissue, wrong depth.**
  The control pair is gastric cancer at 500K reads on chr22; Courcelles is colorectal at full depth, genome-wide.
  GTEx pan-tissue includes colon so the reference remains applicable, but the private-junction fraction may scale differently with depth, and deeper data admits more low-support junctions to the classification step.
- **One patient.**
  Patient-private junction burden is a per-individual property; a single pair cannot bound its variance.

## What would re-open this decision

- A measured `tumor_exclusive` count on Courcelles that is implausibly large relative to the chr22 scaling, which would suggest the private-junction fraction grows with depth rather than staying near 1%.
- Any recovered hit whose only support is a junction that GTEx would not have blacklisted and that has no independent RNA-seq witness in a second patient.
- Discovery of a matched normal arm for this cohort, or a suitable colon-specific normal panel, either of which is strictly better than the pan-tissue substitute.

## Consequence for the write-up

The `tier=presentation-prevalence` rows this Issue produces must carry the substitution in their provenance, not just in prose.
A reader comparing our recovery rate against Courcelles' own must be able to see that our specificity control was a population reference while theirs was eight mTEC controls, because that difference plausibly moves the candidate count in our favour and would otherwise read as sensitivity.

## Reproducibility note

The GTEx chr22 blacklist is gitignored and regenerable.
Two traps cost time on 2026-07-23 and are recorded so the next person does not re-pay them:

- **`--restrict-chrom` does not restrict the fetch.**
  It filters rows after download; the query set defaults to every hg38 primary chromosome.
  Use **`--region chr22`** to scope the tabix query.
  The genome-wide default spent 17 minutes on chr1 alone before a broken pipe.
- **The committed manifest's `gcs_path` is dead.**
  `resources/test/gtex_gtexv2_data_manifest.yaml` records `gcs_staged: true` against `gs://splice-neoepitope-project/...`, a bucket deleted in the June GCP decommission.
  The artifact is fully recoverable by rebuild: the regenerated BED matched the manifest's recorded sha256, byte count, and junction count exactly (880769 junctions, 47561526 bytes, `d727ce5c...`), so this is a stale pointer rather than data loss.
