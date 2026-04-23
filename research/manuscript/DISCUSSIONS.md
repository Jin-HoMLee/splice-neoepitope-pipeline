# Discussion

*Living document — intended to be refined into a publication discussion section.*

---

## Junction-spanning filter: chimeric codons and the complete-codon rule

The junction-spanning filter (see `METHODS.md` §5) requires each retained 9-mer to
contain at least one **complete codon** from each side of the splice junction. This is
a conservative approximation that warrants discussion.

### Why chimeric codons are not checked individually

A chimeric codon — one whose three nucleotides span the junction (e.g. 2 upstream nt +
1 downstream nt) — could in principle encode a novel amino acid not present in the
normal protein. A more precise implementation would translate each chimeric codon and
compare it against the reference proteome, keeping the 9-mer only when the amino acid
genuinely differs.

In practice this check is omitted for two reasons:

1. **Codon degeneracy makes chimeric codons unreliable.** The downstream exon may
   contribute a nucleotide that, combined with the upstream nucleotides, still encodes
   the same amino acid. This was observed directly in the first gastric cancer
   production run: the chimeric last codon of `YLADLYHFV` still encoded valine (V),
   making the entire 9-mer identical to SH3BP1 residues 209–217. There is no guarantee
   that a chimeric codon is novel without explicitly checking.

2. **Peripheral amino acid changes have limited biological impact.** MHC class I
   binding affinity is governed primarily by **anchor positions** — typically positions
   2 and 9 for HLA-A\*02:01 (P2 and PΩ anchor motif). A single amino acid change at
   a peripheral position (1 or 8) arising from a chimeric codon is unlikely to
   meaningfully alter binding affinity or T-cell receptor (TCR) recognition. Such
   peptides are therefore weak neoepitope candidates even if the amino acid does differ.

### Conservative bias and its justification

The complete-codon rule discards some true positives near the junction boundary.
Peptides retained by the rule contain multiple novel downstream amino acids and are more
likely to be genuinely foreign to the immune system. This conservative bias is
preferable in a discovery context: the cost of a missed weak candidate is lower than the
cost of pursuing a false positive through expensive downstream validation.

### Future refinement

A hybrid approach could be considered: apply the complete-codon rule as the primary
filter, then optionally recover chimeric-codon 9-mers where the altered amino acid falls
at an anchor position (P2 or P9). These would represent a small, high-confidence set of
junction-boundary candidates worth investigating further.

---

## Contig upstream length: 26 nt vs. 24 nt

The current contig design uses 26 nt upstream + 24 nt downstream. With `upstream_nt=26`,
the spanning condition is `2 ≤ start ≤ 23`, which means the first valid 27 nt window
across all three reading frames starts at `start=2` (frame 2). The nucleotides at
positions 0 and 1 are never the start of any valid junction-spanning window — they
contribute only as interior nucleotides of later windows.

Reducing `upstream_nt` to 24 would make `min_start=0`, so frame 0 windows starting at
`start=0` would be valid and no upstream nucleotides are wasted. This is a cleaner
design. However, `upstream_nt=26` is deliberately retained for the following reason:

**Chimeric codon data preservation.** The 2 extra upstream nucleotides (nt 0 and 1) are
excluded from valid window starts under the current complete-codon rule, but they remain
in the contig sequence. This preserves the raw data needed to handle chimeric codons —
codons that straddle the splice junction with 1 or 2 nucleotides on one side — if the
pipeline is extended to translate them in a future release.

A more principled future change would be to extend the downstream flank symmetrically
from 24 to 26 nt (giving 26 + 26 = 52 nt contigs), providing chimeric codon coverage on
**both** sides of the junction. This would also simplify the config to a single
`flank_nt` value instead of separate `upstream_nt` / `downstream_nt`.

---

## Reading frame annotation: why translation is not restricted to the canonical frame

The pipeline annotates each tumor-exclusive junction with its canonical reading frame
(derived from the GENCODE CDS, see Methods §5), but translates junction contigs in all
three frames rather than restricting to the annotated one.

### The case for restriction

For a junction whose splice donor matches an annotated CDS exon end in a gene not
otherwise perturbed by somatic mutations, the CDS-derived frame is the most likely frame
being translated. Restricting translation to that frame would reduce the peptide candidate
set by up to two-thirds for those junctions, lowering the false-positive burden on
downstream MHC binding prediction.

### Why restriction is not implemented

Restricting to the canonical frame introduces **false negatives** — true neoepitopes
permanently removed from the candidate set — in scenarios that are common in the tumors
where junction-derived neoepitopes are most clinically relevant:

- **Upstream frameshift indels.** A somatic insertion or deletion in an upstream exon of
  the same gene shifts the reading frame for everything downstream. These events are
  frequent in hypermutated tumors (MSI-high, POLE-mutant) and cannot be identified from
  RNA-seq junction data alone; WGS or WES would be needed to account for them.
- **Structural variants.** A gene fusion or large rearrangement can place exons in a
  reading frame context that has no GENCODE counterpart.
- **Upstream novel junctions.** A second tumor-exclusive junction upstream in the same
  gene may shift the reading frame before reaching the junction of interest. Although the
  pipeline detects co-occurring novel junctions in the same patient, short-read RNA-seq
  cannot phase two junctions to the same transcript, making the propagated frame
  unknowable without transcript assembly.
- **Alternative reading frames.** Some loci encode multiple proteins in different frames
  (e.g. CDKN2A p16/p14ARF). The CDS annotation captures only the canonical frame per
  donor; ARF-frame peptides would be silently dropped by a hard restriction.

In all of these cases, the biologically active frame in the tumor differs from the
GENCODE-derived canonical frame. Crucially, the risk is highest in hypermutated tumors —
precisely those expected to harbour the most actionable neoepitopes overall.

### On six-frame (sense + antisense) translation

Translating both strands of each contig was considered. For strand-specific libraries
(e.g. dUTP second-strand marking, as confirmed for patient_001's gastric cancer samples:
KAPA RNA HyperPrep with RiboErase), only the first-strand cDNA is amplified. HISAT2
assigns the correct strand to all canonical GT-AG junctions via the XS auxiliary tag
(derived from splice-site dinucleotide sequence). Genuine antisense transcription then
appears as junctions on the opposite strand and is already translated in the correct
orientation. Antisense translation of a strand-corrected contig would correspond to the
non-transcribed DNA strand and has no established biological basis for MHC-I presentation.

For non-stranded RNA-seq libraries, strand assignment relies entirely on splice-site
sequence inference and may be incorrect for non-canonical splice sites. In that context,
six-frame translation would be more appropriate. The strandedness of patient_002's
osteosarcoma samples has not been verified; if those samples turn out to be non-stranded,
the strand annotation of their junctions should be treated with caution.

### Current approach

Translation proceeds in all three sense-strand frames. The `reading_frame` annotation in
`novel_junctions.tsv` is retained as metadata: it records the canonical CDS-derived frame
for biological interpretation and downstream stratification of candidates, but does not
gatekeep any peptide from analysis.

---

## Normal sample filtering: junction level vs. peptide level

Tumor-specific junctions are currently defined by absence in the matched normal sample
at the **junction level** (see `METHODS.md` §3). An alternative approach would be to
run MHC binding prediction on both tumor and normal samples and subtract normal
predictions from tumor predictions at the **peptide level**.

The junction-level approach was chosen because:

- It is computationally cheaper (normal predictions never need to be run).
- A junction present in normal tissue is not tumor-specific by definition, regardless
  of the peptide it produces.

The peptide-level approach would additionally catch cases where a tumor-specific junction
produces a peptide that happens to match a normal protein at a completely different
genomic locus. The junction-spanning filter addresses the most common version of this
(same-locus exonic sequence), but a cross-proteome BLAST check would be more thorough.
This remains an open improvement.

---

## HISAT2 vs. STAR for novel junction detection

HISAT2 is used for local development and testing (macOS M1, 8 GB RAM). Benchmarks
consistently show STAR to be more sensitive for novel/unannotated junction detection,
which is the critical step for this pipeline. STAR requires ~32 GB RAM for the full
GRCh38 index and is therefore unsuitable for local runs but appropriate for cloud
production runs.

A planned comparison (issue #17) will run both aligners on the same gastric cancer
samples and compare the number and quality of tumor-specific junctions detected. The
HISAT2 production run provides a baseline.

---

## Impact of missing matched normal: patient_002

Patient_002 (osteosarcoma IPISRC044) has no matched RNA-seq normal sample. Blood WGS
DNA is available but cannot substitute: `regtools junctions extract` requires reads with
spliced CIGAR operations (`N`), which are absent from DNA-seq alignments. Running the
pipeline without a normal labels all unannotated junctions `tumor_exclusive`.

The patient_002 T0 run completed with 347,046 raw tumor junctions; 291,131 were
annotated (GENCODE v47) and discarded, leaving 55,915 unannotated junctions.
BG003082_N0_WES (WES, DNA) was used as the normal input. Although HISAT2 produced
106,474 junctions from the WES alignment, only 3 overlapped with the tumor set —
consistent with WES-derived spliced reads being largely alignment artifacts rather
than biological splice events. Those 3 were discarded as normal-shared, leaving
55,912 tumor-exclusive candidates.

For reference, patient_001 had an 8.9% normal-shared rate among unannotated junctions
using a matched RNA-seq normal. The WES proxy provides minimal filtering (effectively
none), so a corresponding fraction of the 55,912 candidates likely represent
patient-specific but non-tumor splicing. Given the downstream MHCflurry and TCRdock
filtering steps, the impact on the final top candidates is expected to be limited, but
results should be interpreted with this caveat.

If a matched RNA-seq sample becomes available from the osteosarcoma dataset in the
future, the pipeline can be re-run with the normal to apply the full junction-level
filter.
