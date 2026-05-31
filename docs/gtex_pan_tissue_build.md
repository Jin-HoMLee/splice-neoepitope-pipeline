# GTEx Pan-Tissue Junction Reference — Recon + Build Runbook

> **⚠️ STATUS (2026-05-31): source redirected.** This document began as a runbook for the GTEx
> **V10 portal `junctions.gct.gz`**, but Task-1 recon proved that file is **annotation-only**
> (~99.7% annotated on chr1 + chr22) and therefore **cannot** serve as a novel-junction normal
> filter — see "BLOCKING FINDING" below. The verified source is the **Snaptron `gtexv2`
> endpoint** (see "Redirect"). The portal-`.gct` sections are retained as recon evidence, not as
> the production recipe. Parent: [Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126); this slice: [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211).

Reproducibility record + source-format reconnaissance for the GTEx pan-tissue junction blacklist
(parent [Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126)). The build is a **manual one-shot** (NOT part of the per-patient Snakemake DAG); this
document IS the reproducibility record.

## Source files (V10) — verified against the public GCS bucket 2026-05-31

| Role | Path |
|------|------|
| Junction read-counts GCT | `gs://adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz` |
| Sample attributes (sample→tissue) | `gs://adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt` |

HTTPS equivalents: prefix `https://storage.googleapis.com/` to the bucket path (the
bucket is public; no auth needed — verified `gsutil ls` + `curl` both work).

**v11 annotations already exist** on the bucket (`annotations/v11/`), but bulk-gex
RNA-seq junctions are published at **v10** — V10 is the pinned release per #126/#211.

## GCT format (verified by streaming the header + rows)

```
#1.2
392955	19788                         <- 392,955 junctions × 19,788 samples
Name	Description	GTEX-1117F-0005-SM-HL9SH	...   <- col1 Name, col2 gene, then sample columns
chr1:12058-12178:+	ENSG00000223972.5	0	0	...
```

- **Compressed size: 3.85 GiB.** Stream row-by-row — never load the full matrix.
- **Junction `Name` format: `chrom:start-end:strand`** (e.g. `chr1:12228-12612:+`).
  - **Strand IS present** (`:+` / `:-`) — earlier planning assumed it was absent; that was wrong. The GTEx BED can carry real strand, so it is drop-in with the existing `(chrom,start,end,strand)` reference reader in [filter_junctions.py:56-71](../workflow/scripts/filter_junctions.py#L56-L71). The "strand-agnostic / cross-issue note for #212" caveat is **retracted**.
  - This is the **same string format** as `build_reference_junctions.py`'s `name` column (`{chrom}:{start}-{end}:{strand}`).
  - Chromosomes are **UCSC `chr`-prefixed** (matches the pipeline's GENCODE primary assembly convention — no chr/MT normalization needed).
- **Coordinate convention:** GTEx junction coords are **1-based inclusive intron boundaries** (donor first intron base .. acceptor last intron base; STARv2.7.10a-derived). See the verified transform below.

## Coordinate transform (GTEx Name → BED 0-based half-open) — EMPIRICALLY CONFIRMED

```
bed_start = gtex_start - 1      # 1-based donor base -> 0-based half-open start
bed_end   = gtex_end            # 1-based 'last intron base' == 0-based exclusive end (numerically equal)
strand    = gtex_strand         # passthrough (verified: 0 disagreements)
```

**Verification (2026-05-31, chr1):** 5 transform hypotheses were tested by intersecting
the full chr1 GTEx junction set (37,548) against chr1 GENCODE v47 reference junctions
(49,780, built via `build_reference_junctions.py`). Result is a sharp cliff — only the
transform above matches:

| hypothesis | ref junctions matched | recall |
|---|---|---|
| **`start-1, end`** | **37,452** | **75.2%** |
| `start, end` | 2 | 0.0% |
| `start-1, end-1` | 1 | 0.0% |
| `start, end+1` | 1 | 0.0% |
| `start-1, end+1` | 1 | 0.0% |

Strand check under the winning transform: **37,452 agree, 0 disagree** — strand parsing
of the `:+`/`:-` token is confirmed. Worked example (DDX11L1 first intron): GTEx
`chr1:12228-12612:+` → BED `(chr1, 12227, 12612, +)`.

## ⚠️ BLOCKING FINDING — the portal `.gct` is the WRONG source (annotation-only)

The same chr1 cross-check revealed that **`GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz`
contains essentially no novel junctions** and therefore **cannot serve as the pan-tissue
normal filter** this issue (#211 / parent #126) needs:

Corroborated on two independent chromosomes (both vs the **full** GENCODE v47 annotation):

| chromosome | GTEx junctions | annotated | novel |
|---|---|---|---|
| chr1 | 37,548 | 99.74% | **0.26%** |
| chr22 | 8,347 | 99.69% | **0.31%** |

A junction set that is ~99.7% annotated (and, on chr1, *smaller* than the annotation:
37,548 < 49,780) cannot contain meaningful novel content — it is a *subset* of annotated
(collapsed model + a detection threshold). The ~0.3% residual is likely v39↔v47 version
drift, not true novel discovery.

> **Caveat noted in passing:** the chr22 test fixture `resources/test/chr22.gtf.gz` is a
> *reduced* annotation (7,731 junctions vs 10,970 in full-v47 chr22). Comparing GTEx against
> the fixture spuriously inflated "novel" to 24.6%; always compare against the full annotation.

**Why this breaks the filter:** the normal-junction gate only ever sees *novel /
unannotated* junctions — annotated ones are discarded one step earlier
([filter_junctions.py:357](../workflow/scripts/filter_junctions.py#L357)). The matched-normal
gate works on novels because it is built from **raw** alignment junctions
([filter_junctions.py:279](../workflow/scripts/filter_junctions.py#L279),
[:316](../workflow/scripts/filter_junctions.py#L316)). The GTEx gate is the population-scale
analogue of the matched-normal gate, so it **must also be raw / novel-containing**. The
portal `.gct` is the analogue of `reference_junctions.bed` (annotation) instead — a category
error. Intersecting an annotation-only set against all-novel candidates removes ≈ 0 → a silent
no-op filter.

**Redirect → Snaptron `gtexv2` endpoint** (verified 2026-05-31 by a 5-agent recon workflow,
incl. two independent adversarial live probes — both `refuted=false`, high confidence):

- **Source:** `https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=<region>` — **recount3
  reprocessing of GTEx v8, hg38, 19,214 samples / 972 donors, ~33M junctions.** Returns novel +
  annotated junctions with per-sample coverage. Sample count verified live:
  `gtexv2/samples?all=1` → 19,214 rows. Cite **recount3 (Wilks et al., *Genome Biology* 2021)**
  + **Snaptron (Wilks et al., *Bioinformatics* 2018)**.
  - ⚠️ **Do not confuse three different numbers:** `19,214` = Snaptron **`gtexv2`** (recount3/v8 — the
    source we use); `9,662` = the older **`/gtex/`** endpoint (recount2/v6 — `research/scripts/issue_299/`,
    verified `gtex/samples?all=1` → 9,662); `19,788` = the **GTEx V10 portal `.gct`** sample columns
    (the rejected source, unrelated to Snaptron). An earlier draft of this doc conflated these.
- **Live-probe proof (mechanism):** the two adversarial verifiers probed the `/gtex/` (recount2)
  endpoint and found chr1:1,000,000–1,300,000 → 30,134 novel (`annotated=0`) vs 259 annotated,
  every novel record carrying valid `sample_id:count` pairs (reproduced on chr12: 17,957 novel,
  read-counts >1). This confirms the **mechanism** (novel + per-sample) that `/gtexv2/` shares; the
  production build uses `/gtexv2/` (19,214 samples), not `/gtex/`. The portal `.gct` cannot do this.
- **Response schema (18-col TSV):** col 3/4/5/7 = chrom/start/end/strand, **col 8 = `annotated`**
  (0=novel, 1=annotated), **col 13 = `samples`** (RLE `id:count` list), **col 14 = `samples_count`**.
- **Coordinate convention is identical** to the portal: Snaptron `start` is 1-based inclusive
  intron donor → `bed_start = start - 1`, `bed_end = end`. So the transform proven above
  transfers directly (already validated 259/259 against chr22 ground truth in #225's notebook).
- **Reusable code:** `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`
  §2(c) — `fetch_snaptron_chr22(region)` + `snaptron_to_key_set(df, min_samples=1)`. #225's panel
  was the chr22 proxy; #211 generalizes it genome-wide (chr22 alone was 880,769 junctions at
  `samples_count >= 1`, so the genome-wide blacklist is tens of millions of rows — size it for GCS).

**Raw GTEx V10 `SJ.out.tab` ruled out for open access:** exhaustive `gsutil ls -r gs://adult-gtex/**`
returns **zero** `SJ.out.tab` files; per-sample V10 junctions exist only behind dbGaP
phs000424.v10.p1 controlled access (DAR/DUOS + re-aligning ~19,788 samples). Track as a separate
controlled-access subproject only if V10-recency novel junctions ever become mission-critical.
The V10 portal `.gct` may still serve as a *complementary annotated-junction* presence reference,
but NOT as the novel-junction normal filter this issue needs.

### Methods caveats for the Snaptron source (Scientist sign-off, 2026-05-31)

- **Aligner difference → conservative under-filtering.** recount3/Snaptron junctions are
  **STAR**-derived (Monorail pipeline); our pipeline calls junctions with **HISAT2** (local) /
  **STAR** (prod). Any 1-bp coordinate disagreement between aligners means some true normal
  junctions are *missed* by the GTEx filter — i.e. it **under**-filters, never over-filters. That
  is conservative for vaccine safety (we don't wrongly drop a candidate), but it should get a
  sentence in METHODS. The `start-1, end` transform was validated against GENCODE-annotated
  junctions; novel-junction coordinate agreement across aligners is the residual risk.
- **`samples_count >= 1` is intentionally aggressive.** At the `min_samples=1` default the union is
  dominated by single-sample recount artifacts (~880k on chr22 alone). This is the precision-over-
  recall / vaccine-safety choice (parent [Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126)), but the QC sidecar must report a
  `samples_count` sensitivity sweep (union size at `>= {1, 2, 5, 10, …}`) so the threshold can be
  revisited if it costs real candidate yield. `min_samples` ships as a config param, not hard-coded.
- **Immune-privileged tissues (deferred to filter integration, [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)):** a pan-tissue
  blacklist at `min_samples >= 1` excludes junctions seen *only* in testis — but testis is
  immune-privileged and testis-restricted reactivation defines the **cancer-testis antigen** class,
  a validated vaccine-target family. Whether to exempt immune-privileged tissues is a sub-issue-B
  design decision (needs per-tissue provenance, which this construction slice ships count-only).

## Sample attributes format (verified)

Tab-separated; header columns include `SAMPID`, `SMTS` (tissue area), `SMTSD` (detailed
tissue), `SMAFRZE` (analysis-freeze flag). 48,231 sample rows total.

`SMAFRZE` value distribution:

| SMAFRZE | count |
|---------|-------|
| RNASEQ | 19,788 |
| SMLRNA | 16,761 |
| EXCLUDE | 4,750 |
| (blank) | 4,258 |
| WES | 979 |
| WGS | 944 |
| DEEPWGS | 301 |
| OMNI | 450 |

- **Keep `SMAFRZE == 'RNASEQ'` only** (19,788 samples) — these are the bulk RNA-seq
  samples whose `SAMPID` matches the GCT column headers.
- **68 distinct RNASEQ tissues** (`SMTSD`) — not "~54" as the issue estimated.

## counts-by-tissue/ — NOT a junction source

`gs://adult-gtex/bulk-gex/v10/rna-seq/counts-by-tissue/` contains **gene-level** files
only (`gene_reads_v10_<tissue>.gct.gz`). There is no per-tissue junction breakdown, so
the single `..._junctions.gct.gz` is the sole junction source. The issue's "download ~54
per-tissue junction files" framing is superseded: stream the one GCT, group sample
columns by tissue via `SampleAttributesDS`.

## Build command

<!-- Filled in at Task 8/9 with the exact production invocation + GCS staging commands. -->

## Refresh procedure (new GTEx release)

1. Re-run the `Source files` `gsutil ls` to find the new release's junction GCT + SampleAttributesDS paths.
2. Re-verify the GCT `Name` format + coordinate convention (this doc's two sections) — STAR version bumps can change conventions.
3. Re-run the build command with the new URLs + `--release vNN`; re-stage to `gs://splice-neoepitope-project/resources/gtex/vNN/`.
