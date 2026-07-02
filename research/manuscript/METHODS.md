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

Reads are aligned to the GRCh38 (hg38) reference genome with a splice-aware aligner,
selected via the `config.alignment.aligner` switch. The production default is **STAR**
(Dobin et al., 2013), run in two-pass mode (`--twopassMode Basic`) so that novel junctions
found in the first pass are incorporated into the second-pass splice-junction database. The
genome index is built with the GENCODE v47 annotation as its splice-junction database
(`--sjdbGTFfile`, `--sjdbOverhang 100`), so annotated junctions are known to the aligner
while two-pass discovery recovers unannotated ones. The novel-junction output filters
(`--outSJfilter*`) are not overridden, so STAR's permissive defaults apply and discovered
junctions are not down-filtered, and the novel-junction insertion limit is raised
(`--limitSjdbInsertNsj 2000000`) to retain the large set of unannotated junctions of interest. STAR is run without writing alignments (`--outSAMtype None`); splice
junctions are taken directly from its `SJ.out.tab` output (§2).

For memory-constrained environments — for example local development on a macOS M1 with 8 GB
RAM, where STAR's ~32 GB genome index does not fit — the aligner can be switched to **HISAT2**
(Kim et al., 2019), which writes coordinate-sorted BAM files via samtools (single-end data
uses `-U`; paired-end `-1`/`-2`). Both aligners feed a common junction table, so all
downstream stages are aligner-agnostic.

The patient_001 and patient_002 results reported in this manuscript were generated with the
HISAT2 path, the pipeline's earlier default; STAR is the current production default for new runs.

---

## 2. Splice Junction Extraction

Splice junctions are extracted in an aligner-specific manner and normalized to a common
per-sample junction table; downstream stages do not depend on which aligner produced it.

On the **STAR** path, junctions are read directly from STAR's `SJ.out.tab` output. Each
junction carries intron donor/acceptor coordinates, an intron-motif code, and read support;
only uniquely-mapping reads are counted toward support.

On the **HISAT2** path, junctions are extracted from the coordinate-sorted BAM with regtools
(`junctions extract`; Cotto et al., 2023):

- Minimum anchor length: 8 nt (`-a 8`)
- Minimum intron length: 50 nt (`-m 50`)
- Maximum intron length: 500,000 nt (`-M 500000`)
- Strand specificity inferred from the XS tag (`-s XS`)

The resulting BED12 records are converted to the common junction table; the two coordinate
conventions are reconciled as described below. Read support is counted slightly differently by
the two paths — uniquely-mapping reads on the STAR path, total spliced reads on the regtools
path — so a read-count threshold is not strictly comparable across aligners.

### Junction coordinate handling

The two extraction paths reach the same intron coordinates by different routes, and each
requires a path-specific correction.

regtools reports its junctions in BED12, where the start/end columns are the **outer
boundaries of the spliced-read anchors**, not the intron donor and acceptor. The true intron
boundaries are recovered from the block geometry (the donor is the start plus the first block
size; the acceptor is the start plus the second block offset) rather than read from the
start/end columns directly — which would otherwise displace every junction by the anchor
length (typically 100–150 nt per side) and, during classification, silently fail to match any
GENCODE-annotated junction.

STAR reports intron donor/acceptor coordinates directly in `SJ.out.tab` and is not subject to
this offset. It does, however, leave the strand undefined for junctions whose orientation it
cannot infer; in those cases the strand is recovered from the intron-motif code (e.g. GT–AG and
GC–AG → plus strand; CT–AC and CT–GC → minus strand). Junctions with a non-canonical motif,
which cannot be strand-assigned, are dropped rather than carried forward with an ambiguous
strand: an incorrect strand would yield a reverse-complemented contig and translation in the
wrong frame at the contig-assembly stage.

---

## 3. Junction Classification

Junctions are classified into a funnel of categories:

```
all junctions
  └─ annotated              (present in GENCODE annotation)          → discard
  └─ unannotated            (absent from GENCODE annotation)
       ├─ normal_shared          (also in matched normal)            → excluded from prediction (retained in TSV)
       ├─ gtex_pantissue_shared  (in GTEx pan-tissue population normal) → excluded from prediction (retained in TSV)
       └─ tumor_exclusive        (absent in matched normal AND GTEx) → neoepitope prediction
```

**Annotated junctions** are filtered against the GENCODE v47 GTF reference junction set.

**normal_shared junctions** are identified by comparing tumor junctions against the
matched normal sample. A junction present in the normal at ≥ 2 reads is classified
as `normal_shared` (germline or tissue-specific splicing) and excluded from prediction.
These are retained in the output TSV for reporting.

**gtex_pantissue_shared junctions**
([#211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211)/[#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212))
are unannotated junctions that survive matched-normal subtraction but are present in a GTEx
pan-tissue population normal (a population-level splicing blacklist). This always-on stage
catches tissue-specific-but-non-tumor junctions that a matched normal misses - important when
the matched normal is absent or not tissue-matched. Membership is checked *after* matched-normal
subtraction, so the two filters do not double-count. These are retained in the output TSV for
reporting.

**tumor_exclusive junctions** - unannotated and absent from both the matched normal and the
GTEx pan-tissue normal - are carried forward to neoepitope prediction.

When no matched normal sample is available, unannotated junctions not caught by the GTEx filter
are labeled `tumor_exclusive` with a warning.

### Comparison to related approaches

Two recently published resources are relevant to this pipeline's junction classification approach:

- **splice2neo** (Lang et al., 2024, *Bioinform Adv*) is **variant-driven**: it starts from
  somatic SNVs/indels, predicts altered splicing via the deep-learning models SpliceAI
  (Jaganathan et al., 2019) and MMSplice (Cheng et al., 2019), and derives candidate
  neoantigens from the predicted altered transcripts. The pipeline described here is
  **junction-driven**: it starts from observed RNA-seq junctions and classifies them by
  origin (annotated, normal-shared, tumor-exclusive). The two approaches are complementary
  in coverage: the variant-driven approach captures splicing changes attributable to
  specific somatic variants and can detect low-expression events missed by RNA-seq read
  pile-up; the junction-driven approach captures any tumor-exclusive junction regardless
  of underlying cause (somatic variant, intronic mis-splicing, intron retention, or
  splicing-factor dysregulation), but requires sufficient read support at the junction.
  We are evaluating whether to add a parallel variant-driven evidence stream; the
  outcome is reported in the Discussion.
- **AlphaGenome** (Avsec et al., 2026, *Nature*) is a deep-learning model that predicts
  regulatory and splicing outcomes from sequence context, with a dedicated splice-junction
  output that returns donor–acceptor pair probabilities. We are evaluating AlphaGenome as
  a candidate **computational normal filter** — an alternative or complement to the
  matched-normal RNA-seq filter described above, that does not require a per-patient
  normal sample. Inputs are junction coordinates and the GRCh38 reference; outputs are
  per-junction tissue-specific splicing probabilities. The validation outcome is reported
  in the Discussion.

---

## 4. HLA Typing

Patient HLA-A, -B, and -C alleles are typed using OptiType (Szolek et al., 2014).
OptiType aligns RNA-seq reads to an HLA-specific reference (IMGT/HLA) using razers3
and solves an integer linear programme to call the most likely genotype.

For each locus the aggregator (`workflow/scripts/aggregate_hla_alleles.py`) applies a
**tumor-first cascade**:

1. **OptiType from Primary Tumor** (first preference) — reflects what the tumor cell
   actually presents, including any HLA loss-of-heterozygosity (LOH) events that may
   alter the surface allele set.
2. **Clinical serology** (Red Cross / WGS) — germline ground truth, more accurate than
   RNA-seq-based typing when available; read from `serology_X1`/`serology_X2` columns
   in the samples TSV.
3. **OptiType from Solid Tissue Normal / Blood Derived Normal** — used when no tumor
   call passes the read threshold and no serology is available.
4. **Config-defined fallback alleles** (`config.mhcflurry.fallback_alleles`) — used
   only when no sample call is available or confident (reads < `min_reads_per_locus`,
   default 30).

Null alleles (e.g. `A*01:11N`) are not expressed at the cell surface and are excluded
from the prediction allele list, but are preserved in the QC output for transparency.

The aggregated per-patient alleles (`alleles.tsv`) and a QC file flagging the source
of each call and any normal/tumor discrepancies (`hla_qc.tsv`) are written to
`results/hla_typing/{patient_id}/`.

---

## 5. Contig Assembly and Translation

For each tumor-specific junction, a nucleotide contig is assembled by joining flanking
sequences extracted from the GRCh38 reference genome (via `bedtools getfasta`). The flank
size is `3 × (max_peptide_length − 1)` nucleotides on each side; with peptide lengths
8–10 this gives 27 nt symmetric flanks (54 nt contigs).

Reference soft-masking (lower-case bases marking repeat regions) is normalized to upper-case; contigs are retained regardless of masking.

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

## 6. MHC Presentation Prediction

Scoring proceeds in two levels.

**Level 1 — molecular prediction (MHCflurry).**
Junction-spanning peptides are scored using MHCflurry 2.x `Class1PresentationPredictor`,
a composite model that integrates MHC class I binding affinity with antigen processing
predictions (proteasomal cleavage, TAP transport). The predictor is called once with the
full patient HLA genotype (all HLA-A/B/C alleles, ≤6) to obtain the best-allele
prediction per peptide, and then called independently for each allele (one-allele genotype)
to obtain per-allele `presentation_score` and `presentation_percentile`. This dual-call
design is required because the genotype API reports only the best allele and discards scores
for the remaining alleles; the per-allele calls recover those scores.

Each prediction row contains the following scores:

| Column | Description |
|---|---|
| `best_allele` | Allele with the highest presentation score in the patient genotype |
| `ic50_nM` | Binding affinity (nM); informational |
| `processing_score` | Predicted antigen processing efficiency (0–1) |
| `presentation_score` | Composite presentation probability (0–1) for the best allele |
| `presentation_percentile` | Percentile rank of presentation_score for the best allele |
| `{allele}_presentation_score` | Per-allele presentation probability (one column per allele) |
| `{allele}_presentation_percentile` | Per-allele percentile rank (one column per allele) |

Peptides are assigned a single classification label based on `presentation_percentile`
(lower = better, best allele):

- **`presentation_class`** — strong (≤ 0.5%), weak (≤ 2%), non (> 2%)

*Note on terminology.* The 0.5% (strong) and 2% (weak) cutoffs originate as the canonical
IEDB/NetMHCPan **binder** percentile thresholds, originally defined on binding-affinity
percentile. We apply them here to MHCflurry's `presentation_percentile` and refer to them
as strong/weak **presenter** thresholds throughout this manuscript. These are also
MHCflurry's documented default thresholds (O'Donnell et al. 2020, *Cell Systems*); the
0.5% is analogous to the conventional IC50 ≤ 50 nM cutoff.

**Level 2 — genotype presentation score.**
Per-allele scores are combined into a single genotype-level probability using a
complementary-probability formula:

$$\text{genotype\_presentation\_score} = 1 - \prod_{i} \left(1 - w_i \cdot p_i\right)$$

where $p_i$ is the `presentation_score` of allele $i$ and $w_i$ is a locus weight
($w_{\text{HLA-A}} = w_{\text{HLA-B}} = 1.0$; $w_{\text{HLA-C}} = 0.5$ by default,
reflecting ~50% lower surface density of HLA-C relative to HLA-A/B (van Bergen et al.,
2004)). The HLA-C weight is configurable via `config.mhcflurry.hla_c_weight`.

Two supporting statistics are also computed:

- **`n_strong_alleles`** — number of alleles with `presentation_percentile` ≤ 0.5%
- **`best_presentation_percentile`** — minimum `presentation_percentile` across all alleles

**Ranking and quality gate.**
Candidates are ranked by `genotype_presentation_score` (descending) as the primary signal,
with `n_strong_alleles` (descending) as a secondary tie-breaker and
`best_presentation_percentile` (ascending) as a tertiary tie-breaker. A candidate is
excluded from the top-candidates list if `best_presentation_percentile > 2%` (no allele
reaches weak-presenter threshold), providing a quality gate against high-GPS candidates
composed entirely of non-presenting alleles.

---

## 7. Structural Validation (TCRdock)

The top strong-presenting candidate is submitted to TCRdock (Bradley, *eLife* 2023) for structural
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
| `results/hla_typing/{patient_id}/alleles.tsv` | Patient HLA-A/B/C alleles (tumor-first cascade) |
| `results/hla_typing/{patient_id}/hla_qc.tsv` | Per-locus source, read counts, discrepancies |
| `results/junctions/{patient_id}/novel_junctions.tsv` | All unannotated junctions with origin labels |
| `results/peptides/{patient_id}/peptides.tsv` | Junction-spanning 8/9/10-mers (contig_key, start_nt, peptide) |
| `results/predictions/{patient_id}/mhc_presentation.tsv` | Per-peptide MHC presentation predictions: best_allele, ic50_nM, processing_score, presentation_score, presentation_percentile, presentation_class, per-allele presentation_score / presentation_percentile columns, genotype_presentation_score, n_strong_alleles, best_presentation_percentile |
| `results/predictions/{patient_id}/tcrdock/top_candidate.pdb` | Predicted TCR–peptide–MHC ternary complex (PDB) |
| `results/predictions/{patient_id}/tcrdock/docking_scores.tsv` | TCRdock geometry metrics |
| `results/reports/{patient_id}/report.html` | Summary HTML report with HLA QC and Mol\* viewer |

---

## Known Limitations and Future Work

- **Reference-based contig assembly:** flanking exon sequences are extracted from the GRCh38 reference genome rather than assembled from patient reads. Somatic SNVs or indels within the exon flanks are therefore not captured, and predicted peptides assume a wild-type exonic context. This is a standard simplification; the dominant neoepitope signal arises from the novel exon–exon junction sequence itself. In hypermutated tumors (e.g. MSI-high, POLE-mutant), read-based local assembly would improve accuracy.
- **No proteome filter:** peptides are not cross-referenced against the full human
  proteome beyond the junction-spanning filter. A BLAST or exact-match check against
  the reference proteome would catch any remaining false positives.
- **Genome index build:** the production STAR index is built once per VM from the GRCh38
  primary assembly and the GENCODE v47 annotation; the low-memory HISAT2 path instead
  downloads a pre-built `genome_tran` index (GRCh38 + GENCODE splice sites) from the HISAT2 S3
  mirror. The chr22 test configuration builds a small index from a local FASTA.
- **HLA typing from RNA-seq:** OptiType is run on RNA-seq reads, which may have lower
  HLA coverage than WES/WGS. Low read depth (< 30 reads per locus) triggers fallback
  to configured default alleles with a warning.
- **No tissue-matched normal for patient_002:** its only normal is a CD3+ T-cell PBMC, which
  subtracts some junctions but is not tissue-matched to the bone tumor, so its tumor-exclusive
  set is not tissue-specificity-controlled (see Results / Discussion). A tissue-appropriate
  normal, or an explicit population-only mode, is required for a defensible candidate set
  ([#940](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/940)).
