# The functional-validation base for splice-junction neoantigens is critically thin and source-clustered

**Draft manuscript section** (destined for Results / Discussion) · Issue [#737](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/737) · leaf of [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680).
All numbers are computed by [`notebook.ipynb`](notebook.ipynb) from [`registry.tsv`](registry.tsv); figures are in [`outputs/`](outputs/).
This file is the canonical prose; the experiment deck ([`slides.qmd`](slides.qmd)) condenses it.

## Summary claim

Across the entire published literature we could assemble exactly **79 splice-junction-derived peptides** with a *measured* T-cell response and the minimum metadata to score them (sequence + HLA restriction + positive functional readout).
That base is not merely small - it is **concentrated**: two studies supply 67% of it, its mechanism spread is worth ~4 independent categories, and its HLA spread is worth **~1.3 alleles** - effectively a single-allele (HLA-A\*02:01) resource.
The matched negative side is scarcer still: **one** hard true-negative peptide exists field-wide, against the ~19-31 negatives a powered discrimination probe would require.
The scarcity is therefore not a temporary gap to be closed with more curation effort - it is a **structural property of the field's evidence base**, and it bounds what any benchmark built on this ground truth can claim.

## Methods (registry + dual-gate definition)

A peptide enters the registry only if it passes two independent gates.
**Gate 1 (splice origin):** the peptide is derived from an aberrant splice junction - alternative 3′/5′ splice-site use, exon skipping, intron retention, a minor-intron event, or an exonized transposable-element junction - and *not* from a single-nucleotide variant or indel.
**Gate 2 (measured immunogenicity):** its T-cell reactivity was assayed functionally (IFN-γ ELISpot, peptide-MHC tetramer, cytotoxicity / granzyme-B / CD107a, or 4-1BB upregulation), *not* predicted in silico.
Peptides that are MHC-presented by immunopeptidomics but never functionally assayed fail Gate 2 and are held out as a separate untested-decoy tier (see Negative set).

Of 96 registry rows, 92 pass both gates; **79 are scorable** (`tier == functional-scorable`: a published amino-acid sequence, an HLA restriction, and a positive readout).
Concentration is reported three ways per axis: the top contributor's share, the Herfindahl-Hirschman index (HHI = Σpᵢ²), and the **effective number** of contributors (inverse-Simpson, 1/HHI) - the count of *equally-weighted* sources the observed diversity is actually worth.

## Result 1 - the positive base is a few-study assembly

The 79 scorable positives come from **10 studies**, but the distribution is heavily skewed (Figure 1).
The single largest source - Bigot et al. 2021, an SF3B1-mutant uveal-melanoma A2:N tetramer panel - contributes 35 peptides (44%); together with Kim et al. 2025 (SF-mutant leukemia, 18 peptides) the **top two studies supply 67%** of the entire field-wide base.
The effective number of independent studies is **3.7** (HHI = 0.272): ten nominal sources carry the information content of fewer than four.
This source-clustering matters because per-study idiosyncrasies - the tumor type, the splicing lesion (SF3B1 hotspot vs. spliceosome-mutant leukemia vs. exon-TE), the assay, and the HLA background - are confounded with the peptides themselves.
A model that appears to "predict immunogenicity" on this set may instead be learning to recognize the two dominant studies' peptide chemistry.

## Result 2 - the positive base is an HLA-A\*02:01 monoculture

The HLA skew is the sharpest sparsity signal (Figure 2, left).
**70 of 79 scorable positives (89%) are restricted to HLA-A\*02:01**; the effective number of alleles is **1.27** (HHI = 0.789).
The remaining nine peptides are spread thinly across A\*11:01 (4), C\*04 (2), and one each of C\*08, A\*24:02, and A\*02:06.
Every numerically dominant source - the Bigot SF3B1-UM panel, the IR-CRC positives, the Merlotti exon-TE set, and the Kim leukemia fold - is an A2-restricted panel by construction, so the skew compounds rather than averages out.
Any allele-stratified analysis collapses to an A\*02:01 result; claims of cross-allele generalization are unsupported by the available ground truth.
The mechanism axis is comparatively healthier (Figure 2, right): aberrant 3′ splice-site use leads at 40% but the effective number of mechanisms is **4.2**, reflecting genuine representation of exon-TE junctions, frameshift neojunctions, and intron retention alongside the SF3B1-driven 3′SS events.

## Result 3 - the negative set is the binding constraint

Ranking peptides requires negatives, and here the field is at its thinnest (Figure 3).
Exactly **one** peptide qualifies as a hard true-negative - MHC-presented by immunopeptidomics **and** scored non-reactive in a functional assay (`VELEDHVML`, Fisher et al. 2026).
Eight further *soft* negatives (Manoharan IR-CRC) failed to prime in healthy-donor in-vitro sensitization - a categorically weaker claim than a measured non-response, because failure to prime naïve donor T cells does not establish that a patient's repertoire cannot respond.
Pooling the 13 Tier-2 presented-but-untested decoys lifts the *usable-decoy* ceiling to 22, but those peptides were never assayed and cannot anchor a specificity estimate.

To size the deficit, a Hanley-McNeil power calculation (detect a true AUC against the 0.5 null, balanced arms, two-sided α = 0.05, 80% power) requires **19 negatives to detect AUC = 0.75 and 31 to detect AUC = 0.70**.
The positive side clears these thresholds comfortably (79 ≫ 31); the negative side does not approach them.
With a single hard true-negative, the false-positive rate of any splice-neoantigen immunogenicity predictor is, at present, **unmeasurable** field-wide - a specificity cannot be estimated from n = 1.

## Interpretation - a probe, not a powered benchmark

These three results converge on one conclusion: a head-to-head benchmark on splice-neoantigen immunogenicity ground truth is today a **stress-test probe, not a powered evaluation**.
The positives suffice to *rank-order* a predictor's scores, but only within an A\*02:01, ~4-study, ~4-mechanism slice; the absence of validated negatives forecloses any specificity or calibrated-probability claim.
This bounds the #680 program's own headline metric and should be stated wherever a splice-immunogenicity AUC is reported: the denominator is small, single-allele, and source-clustered, and there is essentially no true-negative anchor.
The scarcity is itself a publishable finding - it quantifies *why* the field has no established splice-immunogenicity benchmark, and it identifies the two reagents whose scarcity is rate-limiting (non-A\*02:01 functionally-validated positives, tracked in [#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839); and measured true-negatives) so that future curation and assay effort can target them rather than adding more A2 positives to an already-saturated axis ([#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839) and [#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911) respectively).

## Caveats

- **Source-string hygiene.** Two studies (SNAF, IRIS) carry split source strings in `registry.tsv`; they were collapsed to one canonical label each for the study tally, so the distinct-study count is 10 (not the 12 raw strings).
A future source-key pass is tracked separately; it does not change the concentration conclusion (collapsing *reduces* the apparent study count, strengthening the clustering finding).
- **The 79 is a floor, not a census.** Sequence-blocked positives (e.g. Zhao 2025, GNAS/RPL22) are gate-passing but excluded pending sequence retrieval; folding them would raise the count modestly but does not relieve the A\*02:01 / negative-set constraints.
- **Soft ≠ tested-negative.** The eight IR-CRC soft negatives must never be pooled with the hard negative as equal-strength true-negatives; the tiering is enforced in [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §7 and any metric must declare which tier(s) it used.
- **Power calc is illustrative.** The Hanley-McNeil sizing assumes balanced arms and a single global AUC contrast; it sets the order of magnitude of the negative requirement, not an exact target for a stratified design.
