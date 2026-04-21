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

For each tumor-specific junction, a 50 nt nucleotide contig is assembled by joining:

- 26 nt immediately upstream of the junction (last 26 nt of the upstream exon)
- 24 nt immediately downstream of the junction (first 24 nt of the downstream exon)

Contigs containing soft-clipped (lower-case) bases are excluded.

For each contig, junction-spanning 9-mers are extracted directly in all three reading
frames (offsets 0, 1, 2) and written to a TSV file (`contig_key`, `start_nt`, `peptide`).

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

Junction-spanning 9-mers are scored for MHC class I binding affinity using MHCflurry
2.x (`Class1AffinityPredictor`). Predictions are run for all patient-specific HLA alleles
resolved by OptiType, producing one row per 9-mer × allele combination. Predictions for
multiple alleles are run in parallel using Python's `ProcessPoolExecutor`; each worker
process loads its own predictor instance with TF/BLAS thread counts constrained to 1
to prevent CPU oversubscription.

Peptides are classified as:

- **Strong binder:** IC50 ≤ 50 nM
- **Weak binder:** IC50 ≤ 500 nM
- **Non-binder:** IC50 > 500 nM

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
| `results/predictions/{patient_id}/mhc_affinity.tsv` | All 9-mer × allele IC50 predictions with binder class |
| `results/predictions/{patient_id}/tcrdock/top_candidate.pdb` | Predicted TCR–peptide–MHC ternary complex (PDB) |
| `results/predictions/{patient_id}/tcrdock/docking_scores.tsv` | TCRdock geometry metrics |
| `results/reports/{patient_id}/report.html` | Summary HTML report with HLA QC and Mol\* viewer |

---

## Known Limitations and Future Work

- **No proteome filter:** peptides are not cross-referenced against the full human
  proteome beyond the junction-spanning filter. A BLAST or exact-match check against
  the reference proteome would catch any remaining false positives.
- **HISAT2 vs STAR:** STAR has been shown to detect more novel splice junctions in
  benchmarks. Future production runs will compare results between the two aligners
  (issue #17).
- **Pre-built genome index:** the HISAT2 index is currently built from scratch on each
  new VM run (~60–90 min). Using a pre-built index (issue #16) would reduce startup time.
- **HLA typing from RNA-seq:** OptiType is run on RNA-seq reads, which may have lower
  HLA coverage than WES/WGS. Low read depth (< 30 reads per locus) triggers fallback
  to configured default alleles with a warning.
- **No matched RNA-seq normal for patient_002:** without a normal sample, all unannotated
  junctions are treated as tumor-exclusive. Based on patient_001 statistics (~8.9%
  normal-shared rate), a small fraction of candidates will be false positives arising
  from patient-specific but non-tumor splicing.
