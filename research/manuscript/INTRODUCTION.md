# Introduction

*Living document — intended to be refined into a publication introduction section.*

---

## Cancer Neoepitopes

The immune system can recognize and eliminate tumor cells through cytotoxic T lymphocytes
(CTLs) that detect tumor-specific peptides presented on MHC class I molecules at the cell
surface. These peptides — neoepitopes — arise from somatic alterations in tumor cells that
produce protein sequences absent from normal tissue and therefore not subject to central
immune tolerance.

Neoepitope-based immunotherapy has emerged as a promising approach in cancer treatment,
underpinning the development of personalized cancer vaccines and adoptive T-cell therapies.
Phase 1 trials of personalized mRNA neoantigen vaccines have established clinical
feasibility — most notably autogene cevumeran in pancreatic ductal adenocarcinoma, where
8 of 16 patients developed vaccine-induced T-cell responses (Rojas et al., *Nature* 2023).
Most neoepitope discovery pipelines focus on single-nucleotide variants (SNVs) and small
insertions/deletions (indels) as the source of altered peptide sequences. However, a
broader class of tumor-specific alterations — aberrant splicing — has received comparatively
less attention despite its potential to generate highly immunogenic neoepitopes.

Recent neoantigen prediction pipelines have advanced rapidly. Contemporary deep learning
approaches such as CNNeoPP (Cai et al., *Frontiers in Immunology* 2026) integrate large
language model-derived sequence representations with multi-modal feature fusion to push
SNV-driven prediction beyond the limits of earlier tools. These advances remain confined
to two axes: peptide sequence as the antigen source, and binding-prediction confidence as
the selection criterion. The pipeline described here differs along both axes — splice
junctions as a complementary antigen source not surfaced by SNV-driven sequence
prediction, and an explicit TCR–pMHC docking stage as a downstream structural confidence
layer.

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

This premise is increasingly supported by independent cohorts. A transcript-targeted
antigen-mapping study across 587 glioma patients built per-patient repertoires of
tumor-enriched isoform antigens and their HLA-I-presented peptides, and showed that the
periostin isoform POSTN-203 — associated with poor patient survival — carries multiple
HLA-restricted splice-junction epitopes, one of which (an HLA-A11-restricted peptide)
elicited antigen-specific T-cell responses against POSTN-203-expressing glioma cells
(Xiong et al., *Genes Immun* 2025). That patient-specific, HLA-presented splice-junction
neoepitopes are recoverable from tumor RNA-seq at cohort scale is precisely the premise on
which the pipeline described here is built.

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
| Aligner | TopHat | STAR / HISAT2 (low-memory) |
| Junction extraction | TopHat output | STAR SJ.out.tab / regtools |
| MHC binding prediction | NetMHCPan | MHCflurry 2.x |
| Structural validation | — | TCRdock (Bradley, *eLife* 2023) |
| Workflow management | custom scripts | Snakemake |
| Reference genome | GRCh37/hg19 | GRCh38/hg38 |
| Gene annotation | Ensembl | GENCODE v47 |

Structural validation is a new stage absent from the 2015 pipeline. MHC binding prediction
identifies peptides with the thermodynamic potential to occupy the MHC groove, but surface
presentation is necessary rather than sufficient for T-cell activation: a predicted binder
must also form a productive complex with a T-cell receptor (TCR) to elicit a cytotoxic
response. TCRdock models the three-dimensional TCR-pMHC complex using AlphaFold2 fine-tuned
on crystal structures, providing a docking confidence score as a complementary filter
downstream of MHCflurry.

---

## HLA Genotype, Surface Expression, and Allele Breadth

Each individual carries up to six HLA class I alleles — two from each of the HLA-A, HLA-B,
and HLA-C loci — that are co-dominantly expressed on the surface of nucleated cells,
including tumor cells. Each allele independently samples the intracellular peptide pool and
presents a distinct but overlapping repertoire of peptides to CD8+ T cells. A neoepitope
prediction pipeline that reports only the single best-scoring allele per peptide therefore
discards biologically relevant information: a peptide that is presented at moderate affinity
by all six of a patient's alleles may be a stronger immunotherapy target than a peptide
presented at high affinity by only one.

### Genotype presentation score as a ranking criterion

For a neoepitope to trigger a cytotoxic T cell response, at least one of the patient's HLA
alleles must present it on the tumor cell surface. This is naturally framed as a
complementary probability:

```
P(presented by ≥ 1 allele) = 1 − ∏ᵢ (1 − pᵢ)
```

where `pᵢ` is the presentation probability for allele `i` (e.g. MHCflurry 2.x
`presentation_score`). This genotype presentation score (GPS) rewards candidates that are broadly
presented across the patient's genotype and rises steeply as additional alleles contribute,
while each marginal allele at low affinity adds diminishing probability.

Allele breadth is also clinically relevant beyond the initial ranking:

- **Robustness to HLA loss of heterozygosity (LOH).** Tumors frequently downregulate or
  delete individual HLA alleles as an immune evasion mechanism. A neoepitope presented by
  multiple alleles remains visible to the immune system even if one allele is silenced.
- **Broader T cell recruitment.** Each HLA allele constitutes an independent antigen
  presentation pathway that can prime a distinct T cell clone. A peptide presented by
  three alleles has three independent opportunities to engage a cognate T cell in the
  patient's repertoire.

The integer count of alleles for which a peptide qualifies as a strong presenter
(`n_strong_alleles`) provides a secondary signal: among peptides with equal
GPS, a peptide classified as strong by three alleles is a
biologically cleaner winner than one with one strong and five near-zero predictions.

A complementary but distinct biological force is **immunodominance** — the well-documented
hierarchy in which the strongest-binding epitope for the dominant allele tends to drive the
T cell response while weaker epitopes become subdominant or cryptic (Yewdell & Bennink,
*Annu Rev Immunol* 1999). In natural anti-tumor immunity, a peptide with one exceptionally
strong allele can elicit a more potent focused T cell response than a peptide with moderate
breadth, through intramolecular competition of peptides for the available MHC grooves of
that allele and through immunodomination — dominant T cell clones monopolizing
antigen-presenting cells and suppressing priming of clones restricted to weaker alleles
(Chen & McCluskey, *Adv Cancer Res* 2006). GPS alone does not capture this.

In the context of therapeutic cancer vaccination, however, immunodomination is largely
bypassed: each neoepitope is delivered as a discrete immunogen, allowing independent
priming of T cell clones without the APC-level competition that enforces immunodominance
in natural immunity. The degree of bypass depends on vaccine format — short peptides
competing directly with the endogenous peptidome for empty MHC grooves achieve a more
complete bypass than mRNA vaccines, which re-enter intracellular antigen processing inside
the APC — but the suppression of subdominant T cell clones through APC killing is
alleviated in either case. Combined with the threat of HLA loss of heterozygosity under
immune pressure and the limited number of peptide slots in a personalized vaccine
formulation (typically 10–20 candidates in current clinical trials; Sahin et al.,
*Nature* 2017; Ott et al., *Nature* 2017; Sahin et al., *Nature* 2026), GPS becomes the committed primary
ranking criterion for our application. `best_presentation_percentile` is retained not as a
competing ranking dimension but as a minimum quality gate — at least one allele must reach
strong or weak presenter threshold to ensure sufficient pMHC density for vaccine-primed T
cells to recognize the tumor at the site of disease.

### Differential surface expression across HLA loci

HLA-A and HLA-B are expressed at broadly similar surface densities. HLA-C is
constitutively expressed at approximately 10–50% of HLA-A/B levels across most cell types
(Zemmour & Parham 1992; Trolle et al. 2016). Two mechanisms contribute: HLA-C has lower
intrinsic affinity for β2-microglobulin, reducing pMHC complex stability, and HLA-C
allotypes serve as ligands for killer immunoglobulin-like receptors (KIRs) on NK cells,
with lower surface expression reflecting a balance between NK regulation and CD8+ T cell
priming.

For neoepitope ranking, this means that equivalent `presentation_score` values carry
different expected surface densities depending on the allele's locus. Concretely, an
HLA-C-restricted peptide predicted at the same score as an HLA-A-restricted peptide is
expected to generate fewer pMHC complexes per cell and thus a lower T cell activation
signal. Locus-level weights — HLA-A and HLA-B at 1.0, HLA-C at approximately 0.5 —
scale the effective contribution of each allele to the GPS. Within a locus,
expression differences between individual alleles are comparatively small and are not
captured by current predictors; intra-locus alleles are therefore weighted equally. The
HLA-C weight is treated as a configurable parameter, as published estimates of the
HLA-A/B:HLA-C expression ratio vary across cell types and between studies.

---

## Study Design

This pipeline is applied to matched tumor/normal RNA-seq pairs where available. The matched
normal sample serves as a patient-specific baseline: any splice junction present in the
normal is considered a germline or tissue-specific splicing event and is excluded from
neoepitope prediction. Only junctions that are (i) absent from the GENCODE annotation and
(ii) absent from the matched normal sample are carried forward as tumor-specific candidates.
When no matched normal is available, all unannotated junctions are treated as
tumor-exclusive candidates (see Discussion).

The pipeline is currently applied to two patients:

**Patient_001 — Gastric cancer (matched tumor/normal)**
SRR9143066 (Primary Tumor) / SRR9143065 (Solid Tissue Normal), gastric cancer surgical
section and adjacent stomach tissue. Illumina HiSeq 3000, single-end. Obtained from the
European Nucleotide Archive (ENA). Matched normal available; full junction-level filtering
applied. HLA alleles typed by OptiType from RNA-seq reads; no ground-truth alleles
available for validation.

**Patient_002 — Osteosarcoma (tumor only)**
Patient IPISRC044 (BostonGene sample BG003082) from the publicly available osteosarcoma
dataset (https://osteosarc.com/). T0 tumor RNA-seq (Nov 2022), paired-end Illumina. No
tissue-matched normal is available; the only normal is a CD3+ T-cell PBMC sample, so
patient_002 serves as a matched-normal limitation case rather than a second validated result
(see Results). Its Red Cross serology-typed HLA Class I alleles (A\*01:01/01:11N,
B\*08:01/27:05, C\*01:02/07:01) provide a direct validation opportunity for the OptiType HLA
typing step. Longitudinal tumor samples (T0, T1, T2) enable future analysis of neoepitope
evolution over the course of treatment.
