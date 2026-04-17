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

## Patient_002 — Osteosarcoma (planned)

Patient IPISRC044 from the osteosarcoma public dataset (https://osteosarc.com/).
T0 Boston Gene tumor RNA-seq (Nov 2022), paired-end. Pipeline run planned; results to
be added here after completion.

**Planned validation:** OptiType HLA typing results will be compared against the known
Red Cross serology results (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01) to
evaluate the accuracy of RNA-seq-based HLA typing.
