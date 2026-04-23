# Methods

*Living document — intended to be refined into a publication methods section.*

---

## Overview

The pipeline identifies tumor-specific splice junctions from RNA-seq data and predicts
peptides derived from those junctions that may bind MHC class I molecules (neoepitopes).
The workflow is implemented in Snakemake and consists of seven major stages: alignment,
junction extraction, junction classification, HLA typing, contig assembly and translation,
MHC binding prediction, and structural validation.

---

## 1. RNA-seq Alignment

Reads are aligned to the GRCh38 (hg38) reference genome using HISAT2 (v2.x), a
splice-aware aligner. Alignments are output as coordinate-sorted BAM files using samtools.
Single-end data uses the `-U` flag; paired-end data uses `-1`/`-2`.

---

## 2. Splice Junction Extraction

Splice junctions are extracted from the aligned BAM files using regtools
(`junctions extract`), with the following parameters:

- Minimum anchor length: 8 nt (`-a 8`)
- Minimum intron length: 50 nt (`-m 50`)
- Maximum intron length: 500,000 nt (`-M 500000`)
- Strand specificity inferred from the XS tag (`-s XS`)

Output is a BED file of junction coordinates with read support counts.

---

## 3. Junction Classification

Junctions are classified into three categories:

```
all junctions
  └─ annotated          (present in GENCODE annotation)   → discard
  └─ unannotated        (absent from GENCODE annotation)
       ├─ normal_shared   (also found in matched normal) → excluded from prediction (retained in TSV)
       └─ tumor_exclusive (absent in matched normal)     → neoepitope prediction
```

**Annotated junctions** are filtered against the GENCODE v47 GTF reference junction set.

**normal_shared junctions** are identified by comparing tumor junctions against the
matched normal sample. A junction present in the normal at ≥ 2 reads is classified
as `normal_shared` (germline or tissue-specific splicing) and excluded from prediction.
These are retained in the output TSV for reporting.

**tumor_exclusive junctions** — unannotated and absent from the matched normal — are
carried forward to neoepitope prediction.

When no matched normal sample is available, all unannotated junctions are labeled
`tumor_exclusive` with a warning.

---

## 4. HLA Typing

Patient HLA-A, -B, and -C alleles are typed using OptiType (Szolek et al., 2014).
OptiType aligns RNA-seq reads to an HLA-specific reference (IMGT/HLA) using razers3
and solves an integer linear programme to call the most likely genotype.

A **normal-first policy** is applied: calls from the Solid Tissue Normal sample are
preferred over the Primary Tumor, as HLA alleles are germline and the tumor may carry
loss of heterozygosity (LOH) at the HLA locus. When the normal sample provides
insufficient reads (< 30 per locus), the tumor sample is used. If neither sample
yields a confident call, per-locus fallback alleles from `config.mhcflurry.fallback_alleles`
are substituted with a warning.

The aggregated per-patient alleles (`alleles.tsv`) and a QC file flagging the source
of each call and any normal/tumor discrepancies (`hla_qc.tsv`) are written to
`results/hla_typing/{patient_id}/`.

---

## 5. Contig Assembly and Translation

For each tumor-specific junction, a nucleotide contig is assembled by joining flanking
sequences extracted from the GRCh38 reference genome (via `bedtools getfasta`). The flank
size is `3 × (max_peptide_length − 1)` nucleotides on each side; with peptide lengths
8–10 this gives 27 nt symmetric flanks (54 nt contigs).

Contigs containing soft-clipped (lower-case) bases are excluded.

For each contig, junction-spanning peptides are extracted in all three reading frames
(offsets 0, 1, 2) and written to a TSV file (`contig_key`, `start_nt`, `peptide`).

### Canonical reading frame annotation

For each tumor-exclusive junction, the GENCODE v47 CDS annotation is used to derive the
**canonical reading frame** at the splice donor. Only `transcript_type "protein_coding"`
CDS records are considered; NMD and other non-canonical isoforms are excluded.

For a junction whose splice **donor** (5′ splice site) exactly matches the end of a
protein-coding CDS exon, the frame offset is computed from the GTF `frame` field:

```
phase_at_donor = (exon_length − gtf_frame) mod 3
frame_offset   = (−phase_at_donor)          mod 3
```

where `gtf_frame` encodes how many bases at the 5′ start of the CDS exon complete a codon
from the previous exon. If a donor is covered by multiple protein-coding transcripts that
disagree on the reading frame, the union of attested frame offsets is recorded. If the
donor does not match any protein-coding CDS exon end, the `reading_frame` field is left
empty.

This annotation is written to the `reading_frame` column of `novel_junctions.tsv` and is
available for downstream filtering and biological interpretation. It is **not** used to
restrict the set of frames translated, for the reasons discussed in the Discussion.

### Junction-spanning filter

Only 27 nt windows whose first codon is fully upstream and last codon is fully downstream
are retained. For a window starting at nucleotide `start`:

```
Keep if:   upstream_nt - (window_size - 1) × 3   ≤   start   ≤   upstream_nt - 3

With defaults (upstream_nt = 26, window_size = 9):
           2   ≤   start   ≤   23
```

This filter prevents false positives from 9-mers that fall entirely within one exon
and are therefore identical to normal-proteome sequences. See Discussion for the
chimeric codon rationale.

---

## 6. MHC Binding Prediction

Junction-spanning peptides are scored using MHCflurry 2.x `Class1PresentationPredictor`,
a composite model that integrates MHC class I binding affinity with antigen processing
predictions (proteasomal cleavage, TAP transport). Predictions are run for each
patient-specific HLA allele resolved by OptiType, producing one row per peptide × allele
combination. All alleles are run sequentially in the main process to avoid competing CUDA
contexts on GPU and OOM errors from loading multiple model copies on CPU.

Each prediction row contains four informative scores:

| Column | Description |
|---|---|
| `ic50_nM` | Binding affinity (nM); informational |
| `processing_score` | Predicted antigen processing efficiency (0–1) |
| `presentation_score` | Composite presentation probability (0–1); integrates affinity + processing |
| `presentation_percentile` | Percentile rank of presentation_score among random peptides for that allele |

Peptides are assigned a single classification label based on `presentation_percentile`
(lower = better, per allele):

- **`presentation_class`** — strong (≤ 0.5%), weak (≤ 2%), non (> 2%)

The 0.5% threshold is consistent with Jiang et al. (2024, *Communications Biology*)
and analogous to the conventional IC50 ≤ 50 nM cutoff.

Epitopes are ranked in the output report by `presentation_percentile` (ascending).

---

## 7. Structural Validation (TCRdock)

The top strong-binding candidate is submitted to TCRdock (Bradley et al.) for structural
prediction of the TCR–peptide–MHC ternary complex. TCRdock uses a modified AlphaFold v2
multimer backend adapted specifically for TCR:pMHC modelling. Predictions run on a GCP
Spot GPU VM (NVIDIA P100, 16 GB) inside a Docker container to isolate CUDA/cuDNN
dependencies from the host environment.

The predicted PDB structure has chain IDs reassigned post-prediction (A=MHC heavy chain,
B=peptide, C=TCR-α, D=TCR-β) for compatibility with molecular viewers. The structure and
docking geometry metrics are embedded in the HTML summary report for interactive
visualisation via Mol\* 4.x.

---

## Output

| File | Contents |
|------|----------|
| `results/hla_typing/{patient_id}/alleles.tsv` | Patient HLA-A/B/C alleles (normal-first) |
| `results/hla_typing/{patient_id}/hla_qc.tsv` | Per-locus source, read counts, discrepancies |
| `results/junctions/{patient_id}/novel_junctions.tsv` | All unannotated junctions with origin labels |
| `results/peptides/{patient_id}/peptides.tsv` | Junction-spanning 9-mers (contig_key, start_nt, peptide) |
| `results/predictions/{patient_id}/mhc_affinity.tsv` | All peptide × allele predictions: ic50_nM, processing_score, presentation_score, presentation_percentile, presentation_class |
| `results/predictions/{patient_id}/tcrdock/top_candidate.pdb` | Predicted TCR–peptide–MHC ternary complex (PDB) |
| `results/predictions/{patient_id}/tcrdock/docking_scores.tsv` | TCRdock geometry metrics |
| `results/reports/{patient_id}/report.html` | Summary HTML report with HLA QC and Mol\* viewer |

---

## Known Limitations and Future Work

- **Reference-based contig assembly:** flanking exon sequences are extracted from the GRCh38 reference genome rather than assembled from patient reads. Somatic SNVs or indels within the exon flanks are therefore not captured, and predicted peptides assume a wild-type exonic context. This is a standard simplification; the dominant neoepitope signal arises from the novel exon–exon junction sequence itself. In hypermutated tumors (e.g. MSI-high, POLE-mutant), read-based local assembly would improve accuracy.
- **No proteome filter:** peptides are not cross-referenced against the full human
  proteome beyond the junction-spanning filter. A BLAST or exact-match check against
  the reference proteome would catch any remaining false positives.
- **HISAT2 vs STAR:** STAR has been shown to detect more novel splice junctions in
  benchmarks. Future production runs will compare results between the two aligners
  (issue #17).
- **Pre-built genome index:** the HISAT2 `genome_tran` index (GRCh38 + GENCODE splice sites) is downloaded from the HISAT2 S3 mirror at pipeline start rather than built from scratch, saving ~60–90 min per VM. The chr22 test configuration still builds from the local FASTA.
- **HLA typing from RNA-seq:** OptiType is run on RNA-seq reads, which may have lower
  HLA coverage than WES/WGS. Low read depth (< 30 reads per locus) triggers fallback
  to configured default alleles with a warning.
- **No matched RNA-seq normal for patient_002:** without a normal sample, all unannotated
  junctions are treated as tumor-exclusive. Based on patient_001 statistics (~8.9%
  normal-shared rate), a small fraction of candidates will be false positives arising
  from patient-specific but non-tumor splicing.
