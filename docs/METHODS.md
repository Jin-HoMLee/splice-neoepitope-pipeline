# Pipeline Methods

Draft methods section for the splice neoepitope prediction pipeline.
Intended as a living document to be refined into a publication methods section.

---

## Overview

The pipeline identifies tumor-specific splice junctions from RNA-seq data and predicts
peptides derived from those junctions that may bind MHC class I molecules (neoepitopes).
The workflow is implemented in Snakemake and consists of five major stages: alignment,
junction extraction, junction classification, contig assembly and translation, and
MHC binding prediction.

---

## 1. RNA-seq Alignment

Reads are aligned to the GRCh38 (hg38) reference genome using HISAT2 (v2.x), a
splice-aware aligner. Alignments are output as coordinate-sorted BAM files using
samtools. For single-end data, HISAT2 is run with the `-U` flag.

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
       ├─ patient_specific   (also found in matched normal) → discard
       └─ tumor_specific     (absent in matched normal)     → neoepitope prediction
```

**Annotated junctions** are filtered against the GENCODE v47 GTF reference junction set.

**Patient-specific junctions** are filtered by comparing tumor junctions against the
matched normal sample. A junction present in the normal at ≥ 2 reads is considered
patient-specific (germline or tissue-specific splicing) and discarded. These are
retained in the output TSV for reference but excluded from prediction.

**Tumor-specific junctions** — unannotated and absent from the matched normal — are
carried forward to neoepitope prediction.

When no matched normal sample is available, all unannotated junctions are labeled
`tumor_specific` with a warning.

---

## 4. Contig Assembly and Translation

For each tumor-specific junction, a 50 nt nucleotide contig is assembled by joining:

- 26 nt immediately upstream of the junction (last 26 nt of the upstream exon)
- 24 nt immediately downstream of the junction (first 24 nt of the downstream exon)

Contigs containing soft-clipped (lower-case) bases are excluded.

Each contig is translated in all three reading frames (offsets 0, 1, 2), yielding
up to three 16-mer peptide sequences per junction.

---

## 5. MHC Binding Prediction and Junction-Spanning Filter

### 5.1 The junction-spanning filter

A 16-mer peptide is sliced into 9-mers using a sliding window (8 windows per 16-mer).
**Not all 9-mers are novel sequences**: those that fall entirely within one exon are
identical to sequences present in the normal proteome and are therefore not valid
neoepitope candidates.

The diagram below illustrates the issue for reading frame 1 (offset 0):

```
50 nt contig
|<---------- 26 nt upstream ---------->|<------- 24 nt downstream ------->|
 nt: 0                                 26                                  49
     ─────────────────────────────────┬──────────────────────────────────
                                      ^
                               junction breakpoint

Codons (frame offset = 0, each codon spans 3 nt):
  AA:  0     1     2     3  ...  7     8    ...
  nt: 0-2   3-5   6-8  9-11 ... 21-23 24-26 ...
                                        ^^^
                                 codon spans junction (nt 24, 25 upstream; nt 26 downstream)

9-mer sliding window positions (start_nt = i × 3 for frame offset 0):

  i=0  start_nt= 0  covers nt  0–26
       AA 0–7: nt  0–23  entirely upstream (8 normal amino acids)
       AA 8:   nt 24–26  chimeric codon: 2 nt upstream + 1 nt downstream
                                         ✗ FALSE POSITIVE — see note below

  i=1  start_nt= 3  covers nt  3–29
       AA 0:   nt  3– 5  entirely upstream (complete upstream codon) ✓
       ...
       AA 8:   nt 27–29  entirely downstream (complete downstream codon) ✓

  i=2  start_nt= 6  covers nt  6–32  ✓
  ...
  i=7  start_nt=21  covers nt 21–47  ✓
```

**Why a chimeric codon is not enough (i=0 subtlety).**
At i=0, the 9th amino acid's codon straddles the junction (2 upstream nt + 1 downstream
nt).  In principle this codon could encode a novel amino acid not seen in the normal
protein.  In practice it frequently does not — the downstream exon may contribute a
nucleotide that still encodes the same amino acid by chance (codon degeneracy).  This is
exactly what happened with `YLADLYHFV`: the last codon still encoded valine (V), so all
9 amino acids matched the normal SH3BP1 sequence.  A chimeric codon gives no reliable
guarantee of novelty; only a **complete downstream codon** does.

### 5.2 Spanning condition

Only 9-mers that contain **at least one complete codon from each side** of the junction
are kept. For a 9-mer at 0-indexed position `i` with reading frame offset `f`:

```
start_nt = f + i × 3

Keep if:   upstream_nt - (window_size - 1) × 3   ≤   start_nt   ≤   upstream_nt - 3
           ───────────────────────────────────        ───────────────────────────────
           ensures last AA is fully downstream       ensures first AA is fully upstream
           (complete downstream codon)               (complete upstream codon)

With defaults (upstream_nt = 26, window_size = 9):
           2   ≤   start_nt   ≤   23
```

This filter was present in the original 2015 pipeline and is critical for avoiding
false positives. Without it, the top-ranked strong binder in the gastric cancer
production run (`YLADLYHFV`, IC50 = 9.4 nM) was found by BLAST to be residues 209–217
of the normal SH3BP1 protein — a peptide the immune system is already tolerant to.

### 5.3 MHCflurry prediction

Junction-spanning 9-mers are scored for MHC class I binding affinity using MHCflurry
2.x (`Class1AffinityPredictor`). Peptides are classified as:

- **Strong binder:** IC50 < 50 nM
- **Weak binder:** IC50 < 500 nM
- **Non-binder:** IC50 ≥ 500 nM

The default HLA allele is HLA-A\*02:01 (the most common allele in populations of
European descent). Future work will incorporate patient-specific HLA typing.

---

## 6. Output

The pipeline produces:

| File | Contents |
|------|----------|
| `results/junctions/<cancer_type>/novel_junctions.tsv` | All unannotated junctions with origin labels |
| `results/predictions/<cancer_type>/predictions.tsv` | All 9-mer predictions with IC50 and binder class |
| `results/reports/<cancer_type>/report.html` | Summary HTML report |

---

## Known Limitations and Future Work

- **Single HLA allele:** currently only HLA-A\*02:01 is scored. Patient-specific HLA
  typing (e.g. via OptiType or HLA-HD) would improve clinical relevance.
- **No proteome filter:** peptides are not cross-referenced against the full human
  proteome beyond the junction-spanning filter. A BLAST or exact-match check against
  the reference proteome would catch any remaining false positives.
- **HISAT2 vs STAR:** STAR has been shown to detect more novel splice junctions in
  benchmarks. Future production runs will compare results between the two aligners
  (see issue #17).
- **Pre-built genome index:** the HISAT2 index is currently built from scratch on each
  new VM run (~60–90 min). Using a pre-built index (issue #16) would reduce startup time.
