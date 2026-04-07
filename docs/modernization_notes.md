# Modernisation Notes

## Splice Neoepitope Pipeline — Changes from the 2015 Original

This document provides a detailed, technical account of every substantive
change made when reimplementing the pipeline originally described in Jin-Ho
Lee's 2015 thesis work (Seoul National University).  It is intended for
reviewers and future maintainers who need to understand *why* changes were
made and what behavioural differences to expect.

---

## 1. Reference Genome: hg19 → GRCh38/hg38

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Assembly | GRCh37/hg19 | GRCh38/hg38 |
| Annotation source | UCSC RefSeq hg19 | GENCODE v47 (GRCh38) |

### Rationale
GRCh38 is the current human reference genome released by the Genome Reference
Consortium.  It contains substantially more resolved sequences (especially in
centromeric and pericentromeric regions), updated chromosome names (e.g., ``chr1``
is still ``chr1`` but scaffold naming has changed), and improved annotation.

### Behavioural differences
- **Coordinate shift**: Almost all genomic positions are offset between hg19
  and GRCh38 (liftover is required to compare results between assemblies).
- **Novel chromosomes**: GRCh38 contains unplaced scaffolds (``chrUn_…``) that
  are absent from hg19.  These are included in the reference junction list but
  may rarely appear in TCGA RNA-Seq data.
- **Sample counts differ**: The original paper analysed BRCA (729 tumour, 87
  normal), LUAD (488 tumour, 53 normal), and LAML (154 tumour, 0 normal).  The
  GDC now hosts updated and additional samples; the exact counts will differ.

---

## 2. Data Access: TCGA HTTP Directory → GDC Data Portal API

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Source | TCGA Open-Access HTTP Directory | GDC Data Portal REST API |
| URL base | `https://tcga-data.nci.nih.gov/…` | `https://api.gdc.cancer.gov/` |
| File type | Level 3 junction quantification | Splice Junction Quantification (RNA-Seq) |
| Aligner | TopHat2 (hg19) | STAR (GRCh38, GDC harmonised) |

### Rationale
The legacy TCGA HTTP directory was retired in 2016 when all data was migrated
to the Genomic Data Commons (GDC).  The GDC provides a modern REST API with
structured metadata queries, enabling reproducible programmatic access.

### Behavioural differences
- **Re-aligned reads**: GDC RNA-Seq data was re-processed using STAR against
  GRCh38.  Junction quantification files use STAR's ``SJ.out.tab`` format
  (or its derivative), which has slightly different column structure from
  TopHat2's junction outputs.
- **Junction ID format**: STAR reports junctions as 1-based genomic intervals
  with strand information encoded differently from TopHat2.  The
  `filter_junctions.py` script handles both ``chr:start:end:strand`` and
  ``chr:start-end:strand`` formats.
- **File count**: The GDC may contain more (or fewer) files per cancer type
  as data continues to be curated.

---

## 3. Splice-Junction Reference: Manual hg19 List → GENCODE GTF Derived

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Reference | Manual list of known hg19 junctions | Derived programmatically from GENCODE v47 GTF |
| Script | Not available | `build_reference_junctions.py` |

### Rationale
GENCODE provides a comprehensive, regularly updated annotation for GRCh38
that covers all GENCODE-annotated transcripts.  Deriving the reference junction
list programmatically ensures reproducibility and allows easy updates to newer
GENCODE releases.

### Behavioural differences
- **Coverage**: GENCODE v47 is substantially more comprehensive than the
  hg19 reference list, covering more alternative transcripts.  This means
  *fewer* junctions will pass the novelty filter compared to the original,
  which may reduce the number of neoepitope candidates (more true annotation
  → less "novel" signal).
- **Reproducibility**: The reference junction list is now derived
  programmatically and versioned alongside the pipeline configuration.

---

## 4. Epitope Predictor: NetMHCPan 2.8 → MHCflurry 2.x

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Version | NetMHCPan 2.8 | MHCflurry 2.x |
| License | Academic registration required | Open source (Apache 2.0) |
| Installation | Manual download + registration | `pip install mhcflurry` |
| Neural network | Shallow ANN | Deep neural network with antigen processing |
| Pan-allele coverage | Limited | >16,000 alleles |
| Training data | ~170,000 peptides | >850,000 MS ligandomics peptides |

### Rationale
MHCflurry 2.x is an open-source MHC-I binding predictor that achieves
state-of-the-art performance comparable to NetMHCPan 4.1.  Unlike NetMHCPan,
it does not require academic registration or institutional email addresses,
making it accessible to all users.  MHCflurry can be installed via pip and
its models are downloaded automatically.

### Behavioural differences
- **Different IC50 values**: Absolute IC50 predictions will differ between
  MHCflurry and NetMHCPan.  Peptides classified as strong binders by one tool
  may not all be classified identically by the other.
- **Output format**: MHCflurry output includes affinity (IC50 in nM) and
  percentile rank; the `run_mhcflurry.py` script parses this format.
- **Model download**: MHCflurry requires a one-time model download
  (`mhcflurry-downloads fetch`) which caches ~1 GB of trained models.
- **No registration**: Unlike NetMHCPan, MHCflurry requires no registration
  or licence agreement.

### Reference
> O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
> of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
> *Cell Systems*, 11(1), 42-48.e7.

---

## 5. Biopython API: Bio.Alphabet (removed) → Modern Bio.Seq

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| API | `Bio.Alphabet.generic_dna` | `Bio.Seq.Seq` (no alphabet) |
| Biopython version | <1.78 | ≥1.78 |

### Rationale
`Bio.Alphabet` was deprecated in Biopython 1.78 (released 2020) and removed
in later versions.  The modern API uses `Seq` objects directly without
alphabet specification.  All translation in `translate_peptides.py` uses
`Seq(sequence).translate()` without any alphabet argument.

### Behavioural differences
- None for translation results.  The only difference is API compatibility
  with current Biopython.

---

## 6. Workflow Management: Manual Shell Scripts → Snakemake

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Orchestration | Manual/ad-hoc shell scripts | Snakemake 7+ |
| Reproducibility | Low | High (rule-based DAG, conda envs) |
| Parallelism | Manual | Automatic (``--cores N``) |

### Rationale
Snakemake provides explicit dependency tracking, automatic re-execution of
invalidated steps, and integration with conda environments.  This makes the
pipeline reproducible and easier to extend.

---

## 7. Environment Management: None → Conda

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Dependencies | System-installed | Conda environments (per rule) |
| Reproducibility | None | Pinned via `workflow/envs/*.yaml` |

### Rationale
Conda environments eliminate "works on my machine" issues by precisely
specifying all software versions.

---

## 8. Contig Construction: Unchanged

The contig assembly logic (26 nt upstream + 24 nt downstream = 50 nt) and
the three-reading-frame translation to 16-mer peptides are unchanged from the
original pipeline description.

---

## 9. Statistical Analysis: Equivalent

The one-tailed Fisher's exact test for tumour-vs-normal epitope enrichment is
equivalent to the original.  The implementation uses ``scipy.stats.fisher_exact``
with ``alternative="greater"``.

---

## Summary Change Table

| Component | Original (2015) | Modernised | Reason |
|-----------|----------------|------------|--------|
| Reference genome | hg19 | GRCh38/hg38 | Current standard |
| Reference annotation | UCSC RefSeq hg19 | GENCODE v47 GRCh38 | Comprehensive, reproducible |
| Data source | TCGA HTTP (retired) | GDC Data Portal API | TCGA HTTP was retired in 2016 |
| Aligner | TopHat2 | STAR (GDC harmonised) | GDC re-aligned all data |
| Epitope predictor | NetMHCPan 2.8 | MHCflurry 2.x | Open source, no registration, SOTA |
| Biopython API | Bio.Alphabet | Bio.Seq only | Bio.Alphabet removed in ≥1.78 |
| Workflow | Manual scripts | Snakemake | Reproducibility, parallelism |
| Environments | None | Conda | Reproducibility |

---

## Citation

> Jin-Ho Lee. "Identification of Cancer-Specific Neoepitopes Arising from
> Alternative Splicing Detected by RNA-Seq." Seoul National University, 2015.
