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
  are absent from hg19.  These are included in the reference junction list.

---

## 2. Splice-Junction Reference: Manual hg19 List → GENCODE GTF Derived

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

## 3. Aligner: TopHat2 → HISAT2 or STAR

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Aligner | TopHat2 (hg19) | HISAT2 or STAR (GRCh38), user-configurable |
| Input | Pre-downloaded junction files | FASTQ files aligned locally |

### Rationale
TopHat2 is unmaintained and not compatible with GRCh38 workflows.  HISAT2 is
its successor from the same group, uses substantially less memory (~8 GB vs
TopHat2's multi-pass approach), and produces comparable accuracy.  STAR is an
alternative for users with more RAM (~32 GB) who need maximum sensitivity.

### Behavioural differences
- **Junction ID format**: HISAT2 uses regtools to extract junctions (BED-like
  format); STAR uses ``SJ.out.tab``.  Both are normalised to
  ``chr:start:end:strand`` by the alignment rules before downstream processing.
- **Memory**: HISAT2 ~8 GB; STAR ~32 GB for full GRCh38 index.

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

## 5. HLA Typing: Hardcoded Allele → OptiType Per-Patient

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| HLA alleles | Single hardcoded allele | OptiType HLA typing per patient |
| Source | Not documented | FASTQs (same files used for alignment) |

### Rationale
Patient-specific HLA alleles are required for accurate neoepitope prediction.
The original pipeline used a single representative allele for all samples.
OptiType infers HLA-A/B/C genotype from RNA-Seq reads and is run automatically
when `hla.enabled: true` in the config.

### Behavioural differences
- **More alleles per patient**: MHCflurry is run against up to 6 alleles
  (HLA-A/B/C heterozygous pair), increasing prediction count proportionally.
- **Normal-first policy**: If a normal sample is present, its HLA calls are
  used (germline is more reliable than tumour, which may have LOH).

---

## 6. Biopython API: Bio.Alphabet (removed) → Modern Bio.Seq

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

## 7. Workflow Management: Manual Shell Scripts → Snakemake

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

## 8. Environment Management: None → Conda

| Aspect | Original (2015) | Modern |
|--------|----------------|--------|
| Dependencies | System-installed | Conda environments (per rule) |
| Reproducibility | None | Pinned via `workflow/envs/*.yaml` |

### Rationale
Conda environments eliminate "works on my machine" issues by precisely
specifying all software versions.

---

## 9. Contig Construction: Unchanged

The contig assembly logic (26 nt upstream + 24 nt downstream = 50 nt) and
the three-reading-frame translation to 9-mer peptides are unchanged from the
original pipeline description.

---

## Summary Change Table

| Component | Original (2015) | Modernised | Reason |
|-----------|----------------|------------|--------|
| Reference genome | hg19 | GRCh38/hg38 | Current standard |
| Reference annotation | UCSC RefSeq hg19 | GENCODE v47 GRCh38 | Comprehensive, reproducible |
| Input | Pre-downloaded junction files | FASTQ files | No controlled-access requirement |
| Aligner | TopHat2 | HISAT2 or STAR | TopHat2 unmaintained; successors more accurate |
| HLA alleles | Hardcoded | OptiType per patient | Patient-specific improves recall |
| Epitope predictor | NetMHCPan 2.8 | MHCflurry 2.x | Open source, no registration, SOTA |
| Biopython API | Bio.Alphabet | Bio.Seq only | Bio.Alphabet removed in ≥1.78 |
| Workflow | Manual scripts | Snakemake | Reproducibility, parallelism |
| Environments | None | Conda | Reproducibility |

---

## Citation

> Jin-Ho Lee. "Identification of Cancer-Specific Neoepitopes Arising from
> Alternative Splicing Detected by RNA-Seq." Seoul National University, 2015.
