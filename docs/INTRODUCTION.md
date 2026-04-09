# Introduction

Draft introduction for the splice neoepitope prediction pipeline.
Intended as a living document to be refined into a publication introduction section.

---

## Cancer Neoepitopes

The immune system can recognize and eliminate tumor cells through cytotoxic T lymphocytes
(CTLs) that detect tumor-specific peptides presented on MHC class I molecules at the cell
surface. These peptides — neoepitopes — arise from somatic alterations in tumor cells that
produce protein sequences absent from normal tissue and therefore not subject to central
immune tolerance.

Neoepitope-based immunotherapy has emerged as a promising approach in cancer treatment,
underpinning the development of personalized cancer vaccines and adoptive T-cell therapies.
Most neoepitope discovery pipelines focus on single-nucleotide variants (SNVs) and small
insertions/deletions (indels) as the source of altered peptide sequences. However, a
broader class of tumor-specific alterations — aberrant splicing — has received comparatively
less attention despite its potential to generate highly immunogenic neoepitopes.

---

## Aberrant Splicing in Cancer

Alternative splicing is a fundamental mechanism of gene regulation, and its dysregulation
is a hallmark of cancer. Tumors frequently exhibit splice junctions that are absent from
matched normal tissue, arising from mutations in splicing factor genes, alterations at
splice donor/acceptor sites, or broader transcriptional dysregulation. These tumor-specific
junctions can produce novel open reading frames encoding peptide sequences that are
entirely absent from the normal proteome.

Because splice-junction-derived neoepitopes arise from intronic or exon-boundary sequences
that are never translated in normal tissue, they represent a class of truly tumor-specific
antigens with strong potential for immune recognition.

---

## The Original 2015 Pipeline

This project is a modernised reimplementation of a neoepitope prediction pipeline first
developed in 2015 at Seoul National University (Jin-Ho Lee). The original pipeline
identified tumor-specific splice junctions from RNA-seq data and predicted MHC-binding
neoepitopes from the resulting novel peptide sequences. The core biological logic —
filtering junctions against matched normal tissue and retaining only those absent from
the normal transcriptome — remains unchanged.

The reimplementation updates the toolchain to current standards:

| Component | Original (2015) | Current |
|---|---|---|
| Aligner | TopHat | HISAT2 / STAR |
| Junction extraction | TopHat output | regtools |
| MHC binding prediction | NetMHCPan | MHCflurry 2.x |
| Workflow management | custom scripts | Snakemake |
| Reference genome | GRCh37/hg19 | GRCh38/hg38 |
| Gene annotation | Ensembl | GENCODE v47 |

---

## Study Design

This pipeline is applied to matched tumor/normal RNA-seq pairs. The matched normal sample
serves as a patient-specific baseline: any splice junction present in the normal is
considered a germline or tissue-specific splicing event and is excluded from neoepitope
prediction. Only junctions that are (i) absent from the GENCODE annotation and (ii) absent
from the matched normal sample are carried forward as tumor-specific candidates.

The current study uses a matched gastric cancer tumor/normal pair (SRR9143066 / SRR9143065,
gastric cancer surgical section and adjacent stomach tissue respectively; Illumina HiSeq
3000, single-end). Data were obtained from the European Nucleotide Archive (ENA).
