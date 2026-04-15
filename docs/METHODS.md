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
       ├─ normal_shared   (also found in matched normal) → discard
       └─ tumor_exclusive     (absent in matched normal)     → neoepitope prediction
```

**Annotated junctions** are filtered against the GENCODE v47 GTF reference junction set.

**Patient-specific junctions** are filtered by comparing tumor junctions against the
matched normal sample. A junction present in the normal at ≥ 2 reads is considered
patient-specific (germline or tissue-specific splicing) and discarded. These are
retained in the output TSV for reference but excluded from prediction.

**Tumor-specific junctions** — unannotated and absent from the matched normal — are
carried forward to neoepitope prediction.

When no matched normal sample is available, all unannotated junctions are labeled
`tumor_exclusive` with a warning.

---

## 4. Contig Assembly and Translation

For each tumor-specific junction, a 50 nt nucleotide contig is assembled by joining:

- 26 nt immediately upstream of the junction (last 26 nt of the upstream exon)
- 24 nt immediately downstream of the junction (first 24 nt of the downstream exon)

Contigs containing soft-clipped (lower-case) bases are excluded.

For each contig, junction-spanning 9-mers are extracted directly in all three reading
frames (offsets 0, 1, 2) and written to a TSV file (`contig_key`, `start_nt`, `peptide`).
This replaces the earlier approach of translating full 16-mers and filtering 9-mers
during prediction.

---

## 5. MHC Binding Prediction and Junction-Spanning Filter

### 5.1 The junction-spanning filter

**Not all possible 9-mers from a contig are novel sequences**: those that fall entirely
within one exon are identical to sequences present in the normal proteome and are
therefore not valid neoepitope candidates.

The contig below shows why. `X` = upstream nucleotide, `O` = downstream nucleotide:

```
Contig:  [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXO][OOO][OOO][OOO][OOO][OOO][OOO][OOO][OO]

Legend:
  start  — 0-based nucleotide start of the 27 nt window
  step   — nucleotide step between windows (= 3, one codon)
  X      — upstream nucleotide
  O      — downstream nucleotide

━━━ Frame 0 (start=0, step=3) ━━━━━━━━━━━━━━━━━━━━━━━━━━━

start= 0: [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXO]  ✗ last codon straddles junction
start= 3: [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXO][OOO]  ✓ first valid
start=21: [XXX][XXO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✓ last valid
start=24: [XXO][OOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✗ first codon straddles junction

━━━ Frame 1 (start=1, step=3) ━━━━━━━━━━━━━━━━━━━━━━━━━━━

start= 1: [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XXX][XOO]  ✗ last codon straddles junction
start= 4: [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XOO][OOO]  ✓ first valid
start=22: [XXX][XOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✓ last valid
start=25: [XOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✗ first codon straddles junction

━━━ Frame 2 (start=2, step=3) ━━━━━━━━━━━━━━━━━━━━━━━━━━━

start= 2: [XXX][XXX][XXX][XXX][XXX][XXX][XXX][XOO][OOO]  ✓ first position already valid
start=23: [XXX][OOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✓ last valid
start=26: [OOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO][OOO]  ✗ first codon fully downstream
```

**Why a chimeric codon is not enough.**
A codon that straddles the junction (e.g. 2 upstream nt + 1 downstream nt) could in
principle encode a novel amino acid — but in practice frequently does not due to codon
degeneracy. This was observed directly in the gastric cancer production run: the chimeric
last codon of `YLADLYHFV` still encoded valine (V), making the entire 9-mer identical to
SH3BP1 residues 209–217. A chimeric codon gives no reliable guarantee of novelty; only a
**complete downstream codon** does.

### 5.2 Spanning condition

Only 27 nt windows whose first codon is fully upstream and last codon is fully downstream
are retained. For a window starting at nucleotide `start`:

```
Keep if:   upstream_nt - (window_size - 1) × 3   ≤   start   ≤   upstream_nt - 3
           ───────────────────────────────────        ─────────────────────────────
           ensures last AA is fully downstream       ensures first AA is fully upstream
           (complete downstream codon)               (complete upstream codon)

With defaults (upstream_nt = 26, window_size = 9):
           2   ≤   start   ≤   23
```

This filter is critical for avoiding false positives. Without it, the top-ranked strong
binder in the gastric cancer production run (`YLADLYHFV`, IC50 = 9.4 nM) was found by
BLAST to be residues 209–217 of the normal SH3BP1 protein — a peptide the immune system
is already tolerant to.

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
