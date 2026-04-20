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

Junction-spanning 9-mers were extracted from 50 nt contigs assembled at each
tumor-exclusive junction, across all three reading frames. The junction-spanning filter
(complete upstream and downstream codon requirement) was applied.

| Metric | Count |
|--------|-------|
| Total 9-mers | 433,129 |
| Unique 9-mer sequences | 424,133 |

### MHC Binding Predictions

MHCflurry 2.x was run for all six OptiType-called HLA alleles.

| Binder class | Count | % of total |
|---|---|---|
| Strong (IC50 ≤ 50 nM) | 14,990 | 0.58% |
| Weak (IC50 ≤ 500 nM) | 99,444 | 3.83% |
| Non-binder (IC50 > 500 nM) | 2,484,340 | 95.6% |
| **Total predictions** | **2,598,774** | — |

### Top Neoepitope Candidate

**EVAEYNASF / HLA-A\*26:01, IC50 = 16.5 nM (strong binder)**

This peptide is the highest-affinity predicted binder across all alleles. An IC50 of
16.5 nM is well below the 50 nM strong-binder threshold and represents a compelling
MHC class I binding prediction.

### Structural Validation (TCRdock)

TCRdock was run on the top candidate (EVAEYNASF / HLA-A\*26:01) on a GCP Spot GPU VM
(NVIDIA P100, 16 GB VRAM). The predicted TCR–peptide–MHC ternary complex structure was
successfully generated. The structure was rendered interactively using Mol\* 4.x in the
pipeline HTML report, with chain labels A=MHC, B=peptide, C=TCR-α, D=TCR-β.

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

Junction-spanning 9-mers were extracted from 50 nt contigs across all three reading
frames with the complete-codon junction-spanning filter applied.

| Metric | Count |
|--------|-------|
| Total 9-mers | 781,424 |
| Unique 9-mer sequences | 775,440 |

### MHC Binding Predictions

MHCflurry 2.x was run for all six OptiType-called HLA alleles.

| Binder class | Count | % of total |
|---|---|---|
| Strong (IC50 ≤ 50 nM) | 12,430 | 0.32% |
| Weak (IC50 ≤ 500 nM) | 211,418 | 5.41% |
| Non-binder (IC50 > 500 nM) | 3,683,272 | 94.27% |
| **Total predictions** | **3,907,120** | — |

Best strong binder per allele:

| Allele | Peptide | IC50 |
|--------|---------|------|
| HLA-A\*01:01 | TTDPVQALY | 23.9 nM |
| HLA-B\*08:01 | HAYTKIHSL | 29.5 nM |
| HLA-B\*27:05 | GRFSKVHTF | 45.0 nM |
| HLA-C\*01:02 | AAPPHPLSL | 17.9 nM |
| HLA-C\*07:01 | YRIDRTLSL | 22.7 nM |

### Top Neoepitope Candidate

**TTDPVQALY / HLA-A\*01:01, IC50 = 23.9 nM (strong binder)**

Top predicted binder for the expressed HLA-A allele. HLA-C\*01:02 produced the
highest-affinity prediction overall (AAPPHPLSL, 17.9 nM), consistent with
HLA-C\*01:02 peptide-binding preferences.

### Structural Validation (TCRdock)

TCRdock was run on this sample. Due to a sentinel bug (since fixed in issue #65),
the structural prediction used a fallback allele (HLA-A\*02:01) rather than the
patient's actual HLA-A\*01:01, and the submitted peptide (FMSGFLYFV) is a non-binder
for all six patient alleles (IC50 > 2,200 nM across all alleles). The predicted
structure (pLDDT 92.25) is retained as a technical demonstration of the GPU
infrastructure but should not be interpreted as a biologically meaningful prediction
for this patient. A corrected TCRdock run with TTDPVQALY / HLA-A\*01:01 is planned.
