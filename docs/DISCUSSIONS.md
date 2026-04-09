# Discussions

Deeper reasoning behind design decisions, known tradeoffs, and open questions.
Intended as a living document to inform a future Discussion section in a publication.

---

## Junction-spanning filter: chimeric codons and the complete-codon rule

The junction-spanning filter (see `METHODS.md` §5.2) requires each retained 9-mer to
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
