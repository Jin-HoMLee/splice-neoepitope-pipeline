# Results

*Living document — intended to be refined into a publication results section.*

---

## Patient_001 — Gastric Cancer

### Dataset

Matched tumor/normal RNA-seq pair from a gastric cancer patient. SRR9143066 (Primary
Tumor, gastric cancer surgical section) and SRR9143065 (Solid Tissue Normal, adjacent
stomach tissue). Illumina HiSeq 3000, single-end. Obtained from the European Nucleotide
Archive (ENA). Processed with HISAT2 alignment on GCP. Numbers below are from the valid
post-#148 run (HISAT2 chr-naming bugfix).

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

| Stage | Count | % of total extracted |
|-------|-------|----------------------|
| Total junctions extracted (tumor) | 146,647 | 100.0% |
| Annotated (GENCODE v47, discarded) | 116,618 | 79.5% |
| Unannotated | 30,029 | 20.5% |
| Normal-shared (filtered out) | 2,681 | 1.8% |
| **Tumor-exclusive candidates** | **27,348** | **18.6%** |

Of all extracted junctions, 79.5% were annotated in GENCODE v47 and discarded. A
further 1.8% were unannotated but present in the matched normal (germline or
tissue-specific splicing), leaving 18.6% of extracted junctions as tumor-exclusive
neoepitope candidates. Read support across `tumor_exclusive` junctions: min 16,
median 26, mean 49, max 82,098 — substantially shallower than patient_002 (min 174,
median 425, mean 953).

The unannotated, normal-shared, and tumor-exclusive counts are read from `report.tsv`.
The total-extracted and annotated-discarded rows are derived from the line count of
`alignment/SRR9143066/junctions.tsv` (post-#148 run); these will be added directly to
`report.tsv` once [Issue #104](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104) lands.

### Peptide Translation

Junction-spanning peptides were extracted from contigs assembled at each
tumor-exclusive junction, across all three reading frames, with the
complete-codon junction-spanning filter applied. `peptide_lengths: [8, 9, 10]`
(PR #99), flank size derived automatically as `3 × (max_length − 1) = 27 nt` per side.

| Length | Count |
|--------|-------|
| 8-mer | 382,792 |
| 9-mer | 430,178 |
| 10-mer | 473,522 |
| **Total** | **1,286,492** |
| **Unique sequences** | **1,260,074** |

97.9% of peptides are unique sequences. Peptide totals are computed in
`research/notebooks/patient_001_results.ipynb` §3 from `peptides_novel.tsv`;
the total matches `report.tsv` `mhc_prediction.total_predictions` (one prediction
per unique peptide).

### MHC Presentation Predictions

MHCflurry 2.x `Class1PresentationPredictor` was run for all six expressed HLA alleles.
Presentation class is defined by `presentation_percentile`: strong ≤ 0.5%, weak ≤ 2.0%.

| Presentation class | Count | % of total |
|---|---|---|
| Strong (percentile ≤ 0.5%) | 44,916 | 3.5% |
| Weak (0.5% < percentile ≤ 2.0%) | 125,775 | 9.8% |
| Non (percentile > 2.0%) | 1,115,801 | 86.7% |
| **Total predictions** | **1,286,492** | — |

Strong presenters per allele (best-allele attribution, computed in notebook §4.2):

| Allele | Strong-presenting peptides | % of total strong |
|--------|---------------------------|-------------------|
| HLA-C\*03:03 | 15,123 | 33.7% |
| HLA-A\*31:01 | 12,685 | 28.2% |
| HLA-C\*07:01 | 10,547 | 23.5% |
| HLA-B\*18:01 | 4,068 | 9.1% |
| HLA-A\*26:01 | 1,344 | 3.0% |
| HLA-B\*15:63 | 1,149 | 2.6% |

HLA-C alleles combined account for 57.2% of strong presenters — partial recapitulation
of the HLA-C dominance seen in patient_002 (~69%) despite different HLA-C alleles
(patient_001: C\*03:03 / C\*07:01; patient_002: C\*01:02 / C\*07:01). HLA-A\*31:01 is
a substantial contributor (28.2%), unlike patient_002's HLA-A\*01:01 which was nearly
silent. Median per-allele percentile among genotype-strong presenters is 1.25 for
HLA-C\*03:03 and 10.66 for HLA-B\*18:01, identifying B\*18:01 as the patient's
"passenger" allele.

### Genotype Presentation Score (GPS)

Each peptide was scored with the Genotype Presentation Score:

$$\text{GPS} = 1 - \prod_{i} (1 - w_i \cdot p_i)$$

where $p_i$ is per-allele `presentation_score` and $w_i$ is the locus weight
(HLA-A/B = 1.0, HLA-C = 0.5, reflecting ~50% lower surface density). GPS estimates
the probability that at least one allele in the patient's genotype presents the peptide.

GPS distribution across all 1,286,492 predictions: mean 0.089, median 0.023.
Only 9,276 peptides (0.72%) had GPS > 0.9 — slightly tighter than patient_002 (1.1%),
confirming the metric is discriminating across patients. The GPS inflation edge case
(GPS > 0.9 with `n_strong_alleles = 0`) is markedly smaller in patient_001:
**25 candidates** vs 174 in patient_002 (~7× drop), likely reflecting fewer
borderline-allele-breadth cases when individual allele strengths separate more cleanly.
The proposed `n_strong_alleles ≥ 1` quality gate would remove these 25 from the
strong-presenter pool.

### Top Neoepitope Candidate

**SQIPRTHSY / HLA-C\*07:01**

| Metric | Value |
|--------|-------|
| IC50 | 33.6 nM |
| Presentation percentile | 0.0052% |
| GPS | 0.9999 |
| n_strong_alleles | 5 of 6 |
| Presentation class | strong |

SQIPRTHSY ranks first by GPS and is presented as strong by 5 of 6 patient alleles —
broader breadth than patient_002's FADLRPLLL (4 of 5 = 80% vs SQIPRTHSY 5 of 6 ≈ 83%).
Its presentation percentile places it in the top 0.005% of all HLA-C\*07:01-presented
peptides, indicating high intrinsic affinity rather than a borderline call. Notably,
both patients' top candidates land at GPS ≈ 0.9999 / IC50 ≈ 33–34 nM despite different
splice contexts and different best-allele HLA-C alleles — suggestive of a structural
ceiling around the GPS-best slot rather than a peptide-specific outlier.

HLA-C\*07:01 is shared between patient_001 and patient_002 (patient_001: C\*07:01 /
C\*03:03; patient_002: C\*01:02 / C\*07:01), but the top candidates rest on different
HLA-C alleles per patient (C\*07:01 here, C\*01:02 in patient_002). Patient_001
benefits from a matched RNA-seq normal applied at the junction-filtering step — a
stronger tumor-exclusivity claim than patient_002's WES-only normal.

### Structural Validation (TCRdock)

TCRdock was run on SQIPRTHSY / HLA-C\*07:01, the top GPS-ranked candidate, on a GCP
Spot GPU VM (NVIDIA P100, 16 GB VRAM). The predicted TCR–peptide–MHC ternary complex
structure was successfully generated (`pdb_available: true` in `report.tsv`) and is
available in `results/patient_001/predictions/tcrdock/`, rendered interactively using
Mol\* 4.x in the pipeline HTML report with chain labels A=MHC, B=peptide, C=TCR-α,
D=TCR-β.

*Structural interpretation (pLDDT, CDR3 contacts) to be added after Developer
review of the PDB output, in parallel with patient_002.*

---

## Patient_002 — Osteosarcoma

### Dataset

BG003082 T0 tumor sample (Boston Gene, Nov 2022), paired-end RNA-seq (~10 GB).
No matched RNA-seq normal was available. BG003082_N0_WES (whole-exome sequencing, DNA)
was used as the normal input for HLA typing (successful) but contributes no junctions
to normal subtraction — WES contains no RNA splice junctions by design. Processed with
HISAT2 alignment on GCP (n2-highmem-8).

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

### Junction Filtering

| Stage | Count | % of total extracted |
|-------|-------|----------------------|
| Total junctions extracted (tumor) | 364,168 | 100.0% |
| Annotated (GENCODE v47, discarded) | 305,254 | 83.8% |
| Unannotated | 58,914 | 16.2% |
| Normal-shared (WES normal) | 0 | 0.0% |
| **Tumor-exclusive candidates** | **58,914** | **16.2%** |

BG003082_N0_WES was used as the normal input and yielded valid HLA typing results.
However, WES data contains no RNA splice junctions, so `normal_shared = 0` is
expected — junction-level normal subtraction was not effective. Without a matched
RNA-seq normal, a fraction of the 58,914 tumor-exclusive junctions may represent
tissue-specific or germline splicing rather than tumor-specific events.

### Peptide Translation

Junction-spanning peptides were extracted from contigs assembled at each
tumor-exclusive junction, across all three reading frames, with the
complete-codon junction-spanning filter applied.

| Length | Count |
|--------|-------|
| 8-mer | 703,106 |
| 9-mer | 781,159 |
| 10-mer | 846,422 |
| **Total** | **2,330,687** |
| **Unique sequences** | **2,313,700** |

99.3% of peptides are unique sequences. A previous run (older pipeline version)
reported 781,424 9-mers; the 265-peptide discrepancy is untracked — likely due to
the HISAT2 chr-naming fix (#148) altering junction calls and/or the difference in
normal input (WES proxy vs. no normal). Motivates the run registry (Issue TBD).

### MHC Presentation Predictions

MHCflurry 2.x `Class1PresentationPredictor` was run for all five expressed HLA alleles
(HLA-A homozygous A\*01:01 counted once). Presentation class is defined by
`presentation_percentile`: strong ≤ 0.5%, weak ≤ 2.0%.

| Presentation class | Count | % of total |
|---|---|---|
| Strong (percentile ≤ 0.5%) | 67,935 | 2.9% |
| Weak (0.5% < percentile ≤ 2.0%) | 222,823 | 9.6% |
| Non (percentile > 2.0%) | 2,039,929 | 87.5% |
| **Total predictions** | **2,330,687** | — |

Strong presenters per allele (best-allele attribution):

| Allele | Strong-presenting peptides | % of total strong |
|--------|---------------------------|-------------------|
| HLA-C\*01:02 | 26,236 | 38.6% |
| HLA-C\*07:01 | 20,443 | 30.1% |
| HLA-B\*08:01 | 11,477 | 16.9% |
| HLA-B\*27:05 | 5,193 | 7.6% |
| HLA-A\*01:01 | 4,586 | 6.8% |

HLA-C alleles dominated the strong-presenter set (~69% combined). HLA-A\*01:01 was
nearly silent (median presentation percentile 8.5% among strong presenters).

### Genotype Presentation Score (GPS)

Each peptide was scored with the Genotype Presentation Score:

$$\text{GPS} = 1 - \prod_{i} (1 - w_i \cdot p_i)$$

where $p_i$ is per-allele `presentation_score` and $w_i$ is the locus weight
(HLA-A/B = 1.0, HLA-C = 0.5, reflecting ~50% lower surface density). GPS estimates
the probability that at least one allele in the patient's genotype presents the peptide.

GPS distribution across all 2,330,687 predictions: mean 0.101, median 0.026.
Only 24,961 peptides (1.1%) had GPS > 0.9, confirming the metric is
discriminating. An inflation edge case was identified: 174 candidates (0.7% of
GPS > 0.9) had `n_strong_alleles = 0` with `best_presentation_percentile` ~0.5–0.55%,
just above the strong threshold. The current quality gate does not catch these;
an additional `n_strong_alleles ≥ 1` filter is under consideration.

### Top Neoepitope Candidate

**FADLRPLLL / HLA-C\*01:02**

| Metric | Value |
|--------|-------|
| IC50 | 33.2 nM |
| Presentation percentile | 0.0045% |
| GPS | 0.9999 |
| n_strong_alleles | 4 of 5 |
| Presentation class | strong |

FADLRPLLL ranks first by GPS and is presented as strong by 4 of 5 patient alleles.
Its presentation percentile places it in the top 0.005% of all HLA-C\*01:02-presented
peptides. HLA-C\*01:02 was also the dominant allele in a pre-#148 (invalidated) run,
suggesting it is a consistently strong presenter for splice-junction-derived 9-mers
in this patient.

Without a matched RNA-seq normal, the originating junction cannot be confirmed as
absent from healthy osteoblasts or mesenchymal stem cells.

### Structural Validation (TCRdock)

TCRdock was run on FADLRPLLL / HLA-C\*01:02, the top GPS-ranked candidate.
The predicted TCR–peptide–MHC ternary complex structure was successfully generated
and is available in `results/patient_002/predictions/tcrdock/`.

*Structural interpretation (pLDDT, CDR3 contacts) to be added after Developer
review of the PDB output.*
