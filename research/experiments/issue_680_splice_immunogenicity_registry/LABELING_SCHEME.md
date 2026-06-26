# Labeling scheme - open splice-immunogenicity registry

**Issue:** [#735](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/735) (leaf of epic [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680))
**Date:** 2026-06-26 · **Author:** Scientist
**Status:** canonical; Tasks 2-5 of #735 and the #736 scoring harness conform to this document verbatim.

---

## 1. Purpose

This document is the citable, machine-checkable rule set that makes the registry's positive/negative labels reproducible for the [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) splice-immunogenicity benchmark.
It codifies the rules the adversarial curation already follows, adds a graded `evidence_strength` axis to de-conflate the existing `confidence` column, establishes a `junction_mapping_grade` ladder for peptide-to-junction traceability, and specifies how decoy negatives are constructed without pooling tiers of unequal claim strength.
No row's label is re-derived from scratch here: every previously adversarially-verified row is treated as ground truth, and only boundary cases where a written rule genuinely disagrees with the prior ad-hoc call are adjusted (each adjustment logged in `PROVENANCE.md`).
The scheme is intentionally narrow: it answers "positive or negative, by what assay, how confidently, and how well is the junction mapped?" - nothing more.

---

## 2. Positive labels (`label=positive`)

A row is `label=positive` if it passed both registry gates: (1) genuinely splice-junction-derived (not SNV / indel / fusion / canonical-TAA / proteasomal-cis-spliced) and (2) a functional T-cell assay was performed on that exact peptide.
Positives are further graded by a new `evidence_strength` axis that captures the *type* of assay, not its sample origin:

### `evidence_strength=strong`

A measured **effector** readout was obtained for the peptide.
Qualifying readouts: IFN-γ (ELISA, intracellular cytokine staining, or ELISpot), cytotoxicity / killing (LDH-release, Incucyte confluency, caspase activation), granzyme-B / GZMB secretion, CD107a degranulation, CD137 / 4-1BB upregulation, cloned or engineered-TCR (TCR-T) activation, or MHC-stabilization assay paired with a functional effector readout on the same peptide.
Mass spectrometry (`MS`) detection alone is presentation evidence, not T-cell function, and does NOT upgrade a row to `strong`.

### `evidence_strength=weak`

Antigen-specific T-cell **detection** only, with no effector readout on that peptide.
Qualifying readout: tetramer / dextramer / multimer binding.
A `weak` positive is still a registry-eligible positive; it is retained so the full field-wide validation base is visible and scorable, with the evidence hierarchy transparent.

---

## 3. Negative labels (`label=negative`)

A row is `label=negative` if it is splice-junction-derived (passes gate 1) AND was functionally tested (passes gate 2) AND produced no measurable T-cell response.
Two subtypes of experimental true-negatives are distinguished by `evidence_strength`:

### `evidence_strength=hard`

The peptide is splice-derived, MHC-presented (confirmed by immunopeptidomics or MHC-stabilization), and produced no T-cell response in a functional assay.
This is the gold-standard discrimination negative: both the presentation evidence and the functional-negative evidence are direct.
Example: `VELEDHVML` (Fisher 2026, MS-presented + ELISpot-negative).

### `evidence_strength=soft`

The peptide was tested negative in a healthy-donor in-vitro sensitization (IVS) protocol - i.e., it failed to prime T cells in a moDC-primed donor assay.
Soft negatives are categorically weaker than hard negatives: failed-to-prime in an IVS context does **not** imply the peptide is intrinsically non-immunogenic in a tumor-bearing patient.
The IVS context is precisely what makes the claim soft.
Example: the 8 IR-CRC 9-mers (Manoharan 2026, `candidate-negative` tier).

---

## 4. Excluded from the scorable set

These rows are **retained in `registry.tsv`** (never deleted) but are excluded from the benchmark's scored peptide set.
The `presentation-prevalence` and `negative-control-not-splice` tiers carry `evidence_strength=na` (not applicable); the `functional-nonscorable` and `candidate` tiers retain their real `strong` or `weak` strength from their assay readout.

### `label=untested` (`tier=presentation-prevalence`)

The peptide has MHC-presentation evidence and/or tumor-prevalence / survival signal, but **no functional T-cell assay has been performed**.
It fails gate 2 and therefore cannot be assigned a positive or negative label.
It is retained for completeness and to document that the field's splice-immunogenicity coverage is sparse.

### `tier=functional-nonscorable`

The peptide has been functionally validated (passes both gates), but the **exact amino-acid sequence was not published** in any accessible source.
Because the benchmark keys on sequence, this row cannot be scored and is excluded from the scorable set.
It is retained with its source provenance, and the `notes` column records the reason the sequence is unavailable (e.g., figure-locked MS spectrum, no supplementary peptide table).

### `tier=candidate`

A row in this tier was initially proposed as a `label=positive` but a **second adversarial-verify pass did not confirm** the label.
It is held below scorable until re-verification resolves the dispute.
Do not include `candidate` rows in any positive set.

### `tier=negative-control-not-splice`

The peptide **fails the splice gate** (gate 1) - it is a constitutive / canonical isoform control, not a splice-junction-derived neoepitope.
It belongs to no benchmark score; it is retained for methodological transparency.
`evidence_strength=na` for this tier.

---

## 5. Confidence de-conflation

The existing `confidence=high/medium` column conflates two orthogonal axes:

- **Assay type** - what kind of readout was obtained (effector vs detection).
This axis is now captured explicitly by `evidence_strength` (introduced in this issue, #735).
- **Provenance context** - whether the assay was run in a patient (ex-vivo), a healthy donor (in-vitro sensitization / IVS), or an engineered / cloned-TCR system.
This axis will be captured by a new `assay_context` column to be added in [#823](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/823).

The full decomposition is: `confidence = f(evidence_strength, assay_context)`.
The `confidence` column is **retained unchanged** until #823 lands; once #823 adds `assay_context`, the `confidence` column becomes fully derivable from the two new axes and can be validated for consistency.
This issue (#735) introduces `evidence_strength` only; `assay_context` formalization is explicitly out of scope here.

---

## 6. Junction mapping grade ladder

For each peptide, the `junction_mapping_grade` column records the best junction evidence the source actually published for that peptide.
The four valid values and their exact definitions are:

| Grade | Definition |
|---|---|
| `coords` | The source published genomic donor and acceptor coordinates for that peptide's junction, or published a coordinate that maps 1:1 to those positions without inference. |
| `event-id` | The source published a transcript or splice-event identifier as the unit (e.g., SF3B1 S8 junction, IRIS event ID) but did not publish genomic coordinates. |
| `gene-mechanism` | Only the gene name and canonical splice mechanism are recoverable from the source (e.g., "KRAS intron-4 retention"); no junction identifier or coordinates are available. |
| `none` | No junction information is recoverable; the reason is recorded in the row's `notes` column. |

The ladder is ordered by decreasing specificity: `coords` > `event-id` > `gene-mechanism` > `none`.
Assign the highest grade that the source material actually supports.

### No-inference rule

**Never reconstruct or infer a junction coordinate that the source did not publish.**
It is not permitted to translate a peptide sequence back to a splice-junction coordinate using a genome reference, a splice-site predictor, or any other computational tool not shown in the original source.
If the source published only a gene name, the grade is `gene-mechanism`, not `coords`.
If the source published only an event identifier, the grade is `event-id`, not `coords`.
This rule is strict: the `coords` grade is reserved for coordinates the source explicitly provides (or that map trivially and verifiably 1:1 to a published coordinate, such as a supplement table row whose coordinate and sequence both appear).
Coordinate sparsity - most rows landing at `gene-mechanism` or `none` - is a documented finding about the field's publication practices, not a gap to paper over.

---

## 7. Decoy-negative construction

The benchmark's negative set is built from three tiers of decreasing claim strength.
Each tier is defined here; Tier 1 and Tier 2 are materialized within this issue; Tier 3 generation is deferred to the scoring harness ([#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736)).

### Tier 1 - Experimental true-negatives (strongest claim, smallest n)

The 1 hard-negative and 8 soft-negative rows already in `registry.tsv` (sections 3 above).
These are the only rows with direct experimental evidence of non-immunogenicity for a splice-junction-derived peptide.
`VELEDHVML` (hard) and the 8 IR-CRC 9-mers (soft) are materialized here.

### Tier 2 - Presented decoys from #681

Splice-junction-derived peptides that are MHC-presented (confirmed by immunopeptidomics) but for which **no functional T-cell assay has been performed**.
These are real presented peptides with no measured response - not synthetic constructs.
They are sourced from the 13 peptides enumerated in the registry `README.md` (SNAF Supp Fig 4 + Supp Fig 7, coordinated via [#681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681)) and materialized as `decoy_negatives/presented_decoys_681.tsv`.
The absence of a measured response does not mean a response was sought and found absent; these are untested, not tested-negative.
That distinction makes Tier 2 weaker than Tier 1 as a discrimination signal.

### Tier 3 - Matched synthetic (weakest claim, abundant)

Length-matched and HLA-allele-matched shuffled or decoy-junction peptides, drawn from a non-immunogenic presented background.
These serve as the abundant ranking-denominator floor when the benchmark computes ranking metrics.
The construction algorithm is specified here (length + allele matching, shuffled sequence or decoy-junction source); actual generation is deferred to the #736 harness as a runtime operation, not a fixed artifact committed to this registry.

### Hard rule - tiers are never pooled as equal

**Tiers 1, 2, and 3 must never be pooled into a single undifferentiated negative set for scoring.**
A synthetic decoy (Tier 3) is not equivalent to an experimental hard-negative (Tier 1, `hard` subtype), and must not be counted alongside one as if it carries the same claim strength.
Any score, AUC, or ranked-recall metric that uses negatives from this registry **must** report which negative tier(s) it used, so readers can interpret the result at the appropriate claim level.
A result computed on Tier 1 only (9 rows) is a different claim from one computed on Tier 1+2 (22 rows) or Tier 1+2+3 (many rows); all three are valid but non-interchangeable.
