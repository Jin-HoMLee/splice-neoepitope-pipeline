# Results

*Living document — intended to be refined into a publication results section.*

---

## Patient_001 — Gastric Cancer

### Dataset

Matched tumor/normal RNA-seq pair from a gastric cancer patient. SRR9143066 (Primary
Tumor, gastric cancer surgical section) and SRR9143065 (Solid Tissue Normal, adjacent
stomach tissue). Illumina HiSeq 3000, single-end. Obtained from the European Nucleotide
Archive (ENA). Processed with HISAT2 alignment on GCP (n1-standard-8).

### HLA Typing

OptiType typed six HLA Class I alleles from the RNA-seq data (normal-first policy):

| Locus | Allele 1 | Allele 2 | Source | Discrepancy |
|-------|----------|----------|--------|-------------|
| HLA-A | A\*31:01 | A\*26:01 | normal | concordant |
| HLA-B | B\*15:01 | B\*18:02 | normal | normal=[B\*15:01, B\*18:02] / tumor=[B\*15:63, B\*18:01] |
| HLA-C | C\*07:01 | C\*03:03 | normal | concordant |

HLA-B showed a discrepancy between normal and tumor calls, consistent with allelic
variation in tumor typing rather than true LOH — the normal calls are used per
normal-first policy. No ground-truth HLA alleles are available for this patient;
typed alleles cannot be independently validated. Patient_002 (osteosarcoma) will
provide the first direct validation opportunity.

### Junction Filtering

| Stage | Count | % of total extracted |
|-------|-------|----------------------|
| Total junctions extracted (tumor) | 146,644 | 100.0% |
| Annotated (GENCODE v47, discarded) | 116,615 | 79.5% |
| Unannotated | 30,029 | 20.5% |
| Normal-shared (filtered out) | 2,682 | 1.8% |
| **Tumor-exclusive candidates** | **27,347** | **18.6%** |

Of all extracted junctions, 79.5% were annotated in GENCODE v47 and discarded. A
further 1.8% were unannotated but present in the matched normal (germline or
tissue-specific splicing), leaving 18.6% of extracted junctions as tumor-exclusive
neoepitope candidates.

### Peptide Translation

Junction-spanning peptides were extracted from contigs assembled at each tumor-exclusive
junction, across all three reading frames, with the complete-codon junction-spanning
filter applied.

| Metric | Run 1 (9-mers only) | Run 2 (8/9/10-mers) |
|--------|---------------------|---------------------|
| Contig length | 50 nt | 54 nt |
| Total peptides | 433,129 | 1,292,977 |
| Unique peptide sequences | 424,133 | 1,286,492 |

Run 2 uses `peptide_lengths: [8, 9, 10]` (PR #99), with flank size derived automatically
as `3 × (max_length − 1) = 27 nt` per side.

### MHC Binding Predictions

Two runs are compared: Run 1 used `Class1AffinityPredictor` with per-allele IC50
thresholds; Run 2 used `Class1PresentationPredictor` with a genotype-level call
returning one best-allele prediction per peptide, classified by `presentation_percentile`.

| Metric | Run 1 (AffinityPredictor, IC50) | Run 2 (PresentationPredictor, percentile) |
|---|---|---|
| Total predictions | 2,598,774 (6 alleles × peptides) | 1,286,492 (1 per unique peptide) |
| Strong | 14,990 (IC50 ≤ 50 nM, 0.58%) | 44,916 (percentile ≤ 0.5%, 3.49%) |
| Weak | 99,444 (IC50 ≤ 500 nM, 3.83%) | 125,775 (percentile ≤ 2%, 9.78%) |
| Non | 2,484,340 (95.6%) | 1,115,801 (86.7%) |

The 6× reduction in total predictions reflects the switch from per-allele rows to one
row per peptide. The increase in strong presenters (14,990 → 44,916) reflects the
broader, biologically better-calibrated coverage of the percentile-based threshold
versus the IC50 cutoff.

### Top Neoepitope Candidate

| | Run 1 (AffinityPredictor) | Run 2 (PresentationPredictor) |
|---|---|---|
| Peptide | EVAEYNASF | FAFPFAQTL |
| Allele | HLA-A\*26:01 | HLA-C\*03:03 |
| IC50 (nM) | 16.5 | 19.98 |
| Presentation percentile | — | 0.0007% |
| Class | strong binder | strong presenter |

The change in top candidate reflects ranking by `presentation_percentile` (composite
affinity + processing) rather than IC50 alone. FAFPFAQTL ranks highest by presentation
likelihood across the patient's HLA-C\*03:03 allele; EVAEYNASF ranked highest by raw
binding affinity on HLA-A\*26:01.

### Structural Validation (TCRdock)

TCRdock was run on the top candidate from each run on a GCP Spot GPU VM (NVIDIA P100,
16 GB VRAM). Run 1 modelled EVAEYNASF / HLA-A\*26:01; Run 2 modelled FAFPFAQTL /
HLA-C\*03:03. Both predicted TCR–peptide–MHC ternary complex structures were
successfully generated and rendered interactively using Mol\* 4.x in the pipeline HTML
report, with chain labels A=MHC, B=peptide, C=TCR-α, D=TCR-β.

---

## Patient_002 — Osteosarcoma

### Dataset

BG003082 T0 tumor sample (Boston Gene, Nov 2022), paired-end RNA-seq (~10 GB).
No matched RNA-seq normal was available. BG003082_N0_WES (whole-exome sequencing, DNA)
was used as the normal input; WES-derived spliced HISAT2 alignments are largely
alignment artifacts — see *Discussions* for implications. Processed with HISAT2
alignment on GCP (n2-highmem-8).

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
| Total junctions extracted (tumor) | 347,046 | 100.0% |
| Annotated (GENCODE v47, discarded) | 291,131 | 83.9% |
| Unannotated | 55,915 | 16.1% |
| Normal-shared (WES proxy, discarded) | 3 | ~0.0% |
| **Tumor-exclusive candidates** | **55,912** | **16.1%** |

The WES normal produced only 3 overlapping junctions (all likely alignment artifacts).
Without RNA-seq-quality normal subtraction, a fraction of the 55,912 tumor-exclusive
junctions may represent patient-specific but non-tumor splicing.

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
