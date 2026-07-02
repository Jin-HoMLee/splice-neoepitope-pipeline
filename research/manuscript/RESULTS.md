# Results

*Living document — intended to be refined into a publication results section.*

---

## Patient_001 — Gastric Cancer

### Dataset

Matched tumor/normal RNA-seq pair from a gastric cancer patient. SRR9143066 (Primary
Tumor, gastric cancer surgical section) and SRR9143065 (Solid Tissue Normal, adjacent
stomach tissue). Illumina HiSeq 3000, single-end. Obtained from the European Nucleotide
Archive (ENA). Processed with **STAR** (2-pass), the production default aligner, on the
2026-06-23 cohort run. The numbers below supersede an earlier HISAT2-path draft and reflect
two corrections: the switch to STAR, and the fix of the [#370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370)
anchor-outer coordinate bug that had inflated the earlier HISAT2 junction counts by orders of
magnitude (the pre-#370 draft reported 27,348 tumor-exclusive junctions; see the corrected
funnel below). A GTEx pan-tissue population-normal filter stage
([#211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211)/[#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212))
is now applied after the matched-normal subtraction.

### HLA Typing

OptiType typed six HLA Class I alleles from the RNA-seq data. The valid run used
**tumor-source** calls (per `report.tsv` `hla_typing` rows):

| Locus | Allele 1 | Allele 2 | Source | Notes |
|-------|----------|----------|--------|-------|
| HLA-A | A\*31:01 | A\*26:01 | tumor | normal call concordant |
| HLA-B | B\*18:01 | B\*15:63 | tumor | normal call discrepant: B\*15:01 / B\*18:02 |
| HLA-C | C\*07:01 | C\*03:03 | tumor | normal call concordant |

HLA-B shows a known noise pattern between tumor and normal calls (second-field shifts
B\*15:01 ↔ B\*15:63 and B\*18:02 ↔ B\*18:01), consistent with allelic typing variability
in tumor RNA-seq rather than true LOH. The MHCflurry predictions in this run used the
tumor calls listed above. No ground-truth HLA alleles are available for this patient;
typed alleles are not independently validated (Red Cross serology validation is only
available for patient_002).

### Junction Filtering

| Stage | Count | Passing |
|-------|-------|---------|
| Raw junctions extracted (tumor) | 131,321 | 100.0% |
| Removed - insufficient read support | 99,918 | - |
| Passing read-support filter | 31,403 | 23.9% |
| Annotated (GENCODE v47, discarded) | 31,262 | - |
| Unannotated | 141 | 0.11% |
| Normal-shared (matched normal, filtered) | 94 | - |
| GTEx pan-tissue shared (filtered) | 39 | - |
| **Tumor-exclusive candidates** | **8** | **0.006%** |

After the read-support filter, 31,262 of the 31,403 surviving junctions were annotated in
GENCODE v47 and discarded, leaving 141 unannotated junctions. Of those, 94 were present in
the matched Solid Tissue Normal (germline or tissue-specific splicing) and a further 39 were
seen in the GTEx pan-tissue population normal, leaving **8 tumor-exclusive candidates**. This
supersedes the pre-#370 draft's 27,348 - a ~3,400x inflation that was an artifact of the
anchor-outer coordinate bug ([#370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370):
shifted coordinates failed exact GENCODE matching, so annotated junctions silently leaked
through the discard filter), not biology.

Read support across the 8 `tumor_exclusive` junctions: min 20, median 28, mean 61, max 173.
These are shallower than patient_002's tumor-exclusive set (min 164, median 373), but that
comparison is confounded by a different tumor and a materially weaker matched normal for
patient_002 (see the patient_002 matched-normal caveat below).

Counts are read from `junctions/junction_filter_stats.tsv` on Cloudflare R2 (the 2026-06-23
STAR run).

### Peptide Translation

Junction-spanning peptides were extracted from contigs assembled at each
tumor-exclusive junction, across all three reading frames, with the
complete-codon junction-spanning filter applied. `peptide_lengths: [8, 9, 10]`
(PR #99), flank size derived automatically as `3 x (max_length - 1) = 27 nt` per side.
Across the 8 tumor-exclusive junctions, 418 peptides were translated; 23 were excluded
(incomplete codons / stop-interrupted frames), leaving **395 novel peptides**:

| Length | Count |
|--------|-------|
| 8-mer | 113 |
| 9-mer | 133 |
| 10-mer | 149 |
| **Total novel** | **395** |
| **Unique sequences** | **395** |

All 395 novel peptides are unique sequences. Counts are computed in
`research/notebooks/patient_001_results.ipynb` §3 from `peptides_novel.tsv` (2026-06-23
STAR run). MHCflurry is run on the unique sequences and joined back to each source row
([`run_mhcflurry.py:475`](workflow/scripts/run_mhcflurry.py#L475)) to retain the
`contig_key`/`start_nt` mapping.

### MHC Presentation Predictions

MHCflurry 2.x `Class1PresentationPredictor` was run for all six expressed HLA alleles across
the 395 novel peptides. Presentation class is defined by `presentation_percentile`:
strong <= 0.5%, weak <= 2.0%.

| Presentation class | Count | % of total |
|---|---|---|
| Strong (percentile <= 0.5%) | 17 | 4.3% |
| Weak (0.5% < percentile <= 2.0%) | 22 | 5.6% |
| Non (percentile > 2.0%) | 356 | 90.1% |
| **Total predictions** | **395** | - |

Strong presenters by best-allele attribution (computed in notebook §4.2):

| Allele | Strong-presenting peptides |
|--------|---------------------------|
| HLA-C\*03:03 | 9 |
| HLA-C\*07:01 | 3 |
| HLA-A\*31:01 | 3 |
| HLA-B\*18:01 | 1 |
| HLA-B\*15:63 | 1 |

HLA-C alleles combined account for **12 of 17 (71%)** strong presenters, recapitulating the
HLA-C dominance reported for patient_002 despite different HLA-C alleles (patient_001:
C\*03:03 / C\*07:01; patient_002: C\*01:02 / C\*07:01); this signal survived the #370
correction. At the corrected (much smaller) candidate count the per-allele counts are single-
to low-double-digit, so these proportions are directional rather than precise.

### Genotype Presentation Score (GPS)

Each peptide was scored with the Genotype Presentation Score:

$$\text{GPS} = 1 - \prod_{i} (1 - w_i \cdot p_i)$$

where $p_i$ is per-allele `presentation_score` and $w_i$ is the locus weight
(HLA-A/B = 1.0, HLA-C = 0.5, reflecting ~50% lower surface density). GPS estimates
the probability that at least one allele in the patient's genotype presents the peptide.

GPS distribution across all 395 predictions: mean 0.072, median 0.020.
2 peptides (0.51%) had GPS > 0.9, both with `n_strong_alleles >= 1`. The GPS inflation
edge case (GPS > 0.9 with `n_strong_alleles = 0`) is therefore **0** in the corrected
patient_001 run, down from 25 in the pre-#370 draft. This matters for ranking: the inflated
pre-#370 pool of 27,348 mostly-spurious junctions was far more likely to *contain* a
near-ceiling-GPS peptide, which is how the pre-#370 draft's top candidate reached GPS 0.9999;
the corrected top candidate sits at a realistic GPS 0.942 (below).

### Top Neoepitope Candidate

**SQVTRGLAM / HLA-B\*15:63**

| Metric | Value |
|--------|-------|
| IC50 | 64.7 nM |
| Presentation percentile | 0.257% |
| GPS | 0.942 |
| n_strong_alleles | 3 of 6 |
| Presentation class | strong |

SQVTRGLAM ranks first by GPS and is presented as strong by 3 of 6 patient alleles. Unlike the
pre-#370 draft's top pick (SQIPRTHSY / HLA-C\*07:01, GPS 0.9999, percentile 0.005%), the
corrected top candidate has non-ceiling scores - a realistic GPS 0.942 and a 0.26% presentation
percentile - consistent with removal of the inflated-pool artifact rather than a peptide-specific
outlier. The second-ranked candidate is RTVLQSLWFR / HLA-A\*31:01 (GPS 0.903).

The rank change is a direct consequence of the #370 correction plus the STAR/GTEx pipeline: the
pre-#370 draft's SQIPRTHSY was drawn from the spurious inflated junction pool. One residual
quality flag remains in the corrected set - ranks 3, 5, 6, 7, and 9 (the FFNVGPVLLR / FFNVGPVL
family) derive from a single low-complexity poly-T junction (chr19:39227510) whose assembled
contig is a long poly-T/poly-A run, a known alignment-artifact risk; these are flagged for
scrutiny or complexity-filtering before any downstream use.

### Structural Validation (TCRdock)

TCRdock was run on SQVTRGLAM / HLA-B\*15:63, the top GPS-ranked candidate. The predicted
TCR-peptide-MHC ternary complex was generated (model pLDDT 92.4, PAE 6.44) and is available in
`results/patient_001/tcrdock/`, rendered interactively using Mol\* 4.x in the pipeline HTML
report with chain labels A=MHC, B=peptide, C=TCR-α, D=TCR-β.

**Caveat:** no patient-matched or epitope-matched TCR was available, so the docking used a
fallback TCR (DMF5; `panel_status = dmf5_fallback`, no VDJdb match). The confidence metrics
therefore describe the predicted geometry of a generic TCR against this peptide-MHC, not a
validated cognate pairing - read as a foldability/geometry check, not evidence of
patient-specific recognition.

---

## Patient_002 — Osteosarcoma

### Dataset

BG003082 T0 tumor sample (Boston Gene, Nov 2022), paired-end RNA-seq (~10 GB). Processed with
**STAR** (the production default aligner) on the 2026-06-23 cohort run. The matched normal is a
**CD3+ T-cell PBMC scRNA-seq** sample (Hudson Lab, Jan 2025; sample-sheet issue #277), which
replaced an earlier WES blood normal. Critically, **no tissue-matched normal** (osteoblast,
mesenchymal, or adjacent bone) is available for this patient - a limitation that governs how its
junction-level results can be interpreted (see below).

### HLA Typing

OptiType typed six HLA Class I alleles from the tumor RNA-seq data:

| Locus | Allele 1 | Allele 2 | Source | Discrepancy |
|-------|----------|----------|--------|-------------|
| HLA-A | A\*01:01 | A\*01:01 | tumor | concordant |
| HLA-B | B\*08:01 | B\*27:05 | tumor | concordant |
| HLA-C | C\*07:01 | C\*01:02 | tumor | concordant |

**HLA typing validated against Red Cross serology** (A\*01:01/A\*01:11N,
B\*08:01/B\*27:05, C\*01:02/C\*07:01). All six alleles were an exact match.
A\*01:11N is a null allele; OptiType correctly identified A\*01:01 as the expressed
allele and called it homozygous. This is the first direct validation of OptiType
accuracy in this pipeline.

### Junction filtering and the matched-normal limitation

Unlike patient_001 (a tissue-matched tumor / adjacent-normal pair), patient_002 has **no
tissue-matched normal**. Its only normal available for junction subtraction is the CD3+ T-cell
PBMC transcriptome, and because the pipeline subtracts junctions seen in *any* normal-typed
sample ([#940](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/940)), this
blood/immune-lineage sample drove the subtraction in the 2026-06-23 STAR run: of the tumor's
442 unannotated junctions, 154 were labeled normal-shared, 246 were removed by the GTEx
pan-tissue population filter, leaving 42 nominal tumor-exclusive junctions.

That set cannot be read as tissue-specificity-controlled. A CD3+ T-cell transcriptome is a poor
junction-subtraction normal for a bone/osteosarcoma tumor: T cells express a narrow,
lineage-specific splicing repertoire, so junctions that are normal for bone or mesenchymal
tissue but simply absent from T cells are never subtracted and survive as spurious
"tumor-exclusive" calls. The GTEx pan-tissue filter partially compensates - it removes 246
junctions here versus 39 for patient_001 - but bulk pan-tissue coverage does not fully capture
bone- or osteosarcoma-microenvironment splicing. (The tumor-exclusive junctions also run much
deeper here, median 373 mapped reads vs 28 for patient_001, reflecting a different tumor and
library rather than stronger tumor-specificity.)

We therefore treat patient_002 as a **demonstration of a methodological boundary condition**, not
a second validated result. Its 42 tumor-exclusive junctions are **not comparable** to
patient_001's 8, and the downstream figures (peptide translation, MHC presentation, GPS, and the
top-GPS candidate VPQVRVTVL / HLA-B\*08:01) are regenerated for the record in
`research/notebooks/patient_002_results.ipynb` but are deliberately **omitted from the
comparative results here**. The generalizable lesson: for tumors lacking a tissue-matched normal,
junction-level tumor-specificity cannot be established from a blood normal alone, and
population-level filters are load-bearing but insufficient. A tissue-appropriate normal - or an
explicit population-only mode with the blood normal excluded from subtraction
([#940](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/940)) - is the
prerequisite for a defensible patient_002 candidate set. The HLA typing above is unaffected and
remains the pipeline's first serology-validated OptiType result.

*Downstream figures for patient_002 (peptide translation, MHC presentation, GPS, top
candidate, and TCRdock on VPQVRVTVL / HLA-B\*08:01) are regenerated in
`research/notebooks/patient_002_results.ipynb` for the record but are intentionally not
tabulated here - see the matched-normal limitation above.*
