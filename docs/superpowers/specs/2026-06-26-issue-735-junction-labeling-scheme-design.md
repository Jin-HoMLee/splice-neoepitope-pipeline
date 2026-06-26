# Design — Junction mapping + documented pos/neg labeling scheme (Issue #735)

**Issue:** [#735](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/735) (leaf of epic [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680), carries original #680 AC#2)
**Arc:** `immunogenicity-benchmark` · **Milestone:** `i6 - S3 - Data Preparation` · **Size:** M
**Date:** 2026-06-26 · **Author:** Scientist

## Problem

The open splice-immunogenicity registry (`research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv`, 79 rows) was adversarially curated, but its positive/negative labeling lives in *informal* columns (`label`, `tier`, `confidence`, free-text `readout`) with no documented, citable rule set, and no peptide→junction mapping.
A benchmark needs both: every peptide tied to its originating junction (where the source published it) and a reproducible scheme stating which assay outcomes count as positive, how negatives and hard-negatives are defined, and how confidence is graded.

This issue delivers the **documented scheme**, the **junction annotation**, the **applied labeled set with per-peptide rationale**, and the **decoy-negative construction approach** — without re-litigating the verified curation and without inferring data a source did not publish.

## Scope guard (what this issue does NOT do)

- Does **not** generate synthetic decoys or build the FASTA+label benchmark file → that is the scoring harness, [#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736).
- Does **not** add the `assay_context` column (patient vs healthy-donor-IVS) → that is the data-quality pass, [#823](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/823). This design *defines* the axis and documents how it feeds confidence, but leaves the column to #823.
- Does **not** translate/reconstruct junction coordinates a source did not publish (the README's explicit Zhao rule; no inference).
- Does **not** re-derive labels from scratch — existing adversarially-verified rows are treated as ground truth; only genuine rule-violations among boundary cases are adjusted, each logged.

## Approach (decisions, with rationale)

Four design forks, resolved:

1. **Scheme basis — Hybrid (formalize-in-place + re-audit boundary cases).** Codify the rules the existing curation already follows; then re-examine only the ~14 boundary rows (2 `candidate`, 8 `candidate-negative`, the `presentation-prevalence` row, the 2 `functional-nonscorable`, the `negative-control-not-splice` row) against the written rule and adjust where the rule genuinely disagrees with the ad-hoc call. Bounds churn; the explicit scheme earns its keep by surfacing divergences rather than rubber-stamping.

2. **Positive boundary — graded positive.** Both effector-confirmed and tetramer-detection-only rows stay `label=positive`, separated by a new `evidence_strength` axis (`strong` = measured effector function; `weak` = antigen-specific detection only). Preserves all 68 positives (a strict effector-only rule would drop ~40, gutting an already field-tiny set) while encoding the real evidence hierarchy so a scorer can threshold or weight.

3. **Negatives — tiered, never pooled.** Three documented tiers of decreasing claim strength: (1) experimental true-negatives (1 hard + 8 soft IVS), (2) MS-presented-but-untested splice peptides from #681 as "presented decoys", (3) length+allele-matched synthetic decoys. The scheme forbids scoring a synthetic decoy as equivalent to an experimental hard-negative.

4. **Junction mapping — tiered annotation, no inference.** A `junction_mapping_grade` ladder (`coords` > `event-id` > `gene-mechanism` > `none`) captures the best junction evidence each source actually published. Coordinate-level coverage is expected to be thin; that thinness is a documented finding, not a gap to paper over.

## Deliverables

All under `research/experiments/issue_680_splice_immunogenicity_registry/`:

| Artifact | Type | Serves |
|---|---|---|
| `LABELING_SCHEME.md` | new | AC#2 (scheme) + AC#4 (decoy approach) |
| `registry.tsv` | updated (4 new columns + re-audited rows) | AC#1, AC#3 |
| `decoy_negatives/presented_decoys_681.tsv` | new (materialized Tier-2 seed, 13 rows) | AC#4 (concrete #681 coordination) |
| `README.md`, `PROVENANCE.md` | updated (new columns + re-audit changelog) | traceability |

### New `registry.tsv` columns

| Column | Values | Notes |
|---|---|---|
| `junction_id` | genomic `chr:donor-acceptor:strand`, OR transcript/event ID, OR blank | only what the source published |
| `junction_mapping_grade` | `coords` / `event-id` / `gene-mechanism` / `none` | honest coverage ladder |
| `evidence_strength` | positives `strong`/`weak`; negatives `hard`/`soft` | assay-type axis (de-conflated from `confidence`) |
| `label_rationale` | one line per peptide | why this label, per AC#3 |

Existing columns are preserved; `confidence` (`high`/`medium`) is retained and documented as `evidence_strength` × `assay_context` (the latter formalized later in #823).

## The labeling scheme (`LABELING_SCHEME.md` contents)

### Positives (`label=positive`)
- **`strong`** — a measured effector readout: IFN-γ, cytotoxicity/killing, ELISpot, granzyme-B/CD107a degranulation, engineered/cloned-TCR activation, or MHC-stabilization paired with function.
- **`weak`** — antigen-specific detection only: tetramer/dextramer/multimer binding with no effector readout on that peptide.

### Negatives (`label=negative`)
- **`hard`** — splice-derived, MHC-presented, functional assay performed, no T-cell response (e.g. `VELEDHVML`). The gold-standard discrimination negative.
- **`soft`** — tested-negative in healthy-donor in-vitro sensitization (failed to prime ≠ intrinsically non-immunogenic; e.g. the 8 IR-CRC 9-mers). Categorically weaker; the IVS context is what makes it soft.

### Excluded from the scorable set (retained in TSV, flagged)
- `presentation-prevalence` — presented + prevalence/survival evidence but **no functional assay** → `label=untested`.
- `functional-nonscorable` — validated but exact sequence unpublished → cannot be keyed.
- `candidate` — positive whose 2nd adversarial-verify pass did not confirm → held until re-verified.
- `negative-control-not-splice` — fails the splice gate (constitutive control) → not a splice benchmark row.

### Confidence de-conflation (documented, not yet columnar)
`confidence=high/medium` currently mixes two orthogonal axes:
- **assay type** → captured now as `evidence_strength`.
- **provenance context** (patient ex-vivo vs healthy-donor IVS vs engineered) → to be captured as `assay_context` in #823.
The scheme documents the mapping `confidence = f(evidence_strength, assay_context)` so the split is reproducible once #823 lands.

## Junction mapping (AC#1)

For each peptide, record the best junction evidence the source actually published, graded:

| Grade | Meaning |
|---|---|
| `coords` | source published genomic donor/acceptor (or recoverable 1:1 from a published coordinate) |
| `event-id` | source published a transcript/splice-event identifier as the unit (e.g. SF3B1 S8 junction, IRIS event ID) but not genomic coords |
| `gene-mechanism` | only gene + canonical mechanism available |
| `none` | unmappable; `notes` records the reason |

**No inference.** Never translate a peptide back to coordinates the source did not publish. Coordinate-level coverage is expected to be a minority; the per-grade tally is reported in the README as a finding (junction-level grounding mirrors the thin functional base).

## Decoy-negative construction (AC#4)

A documented 3-tier hierarchy; the scheme specifies *construction* for each; only Tiers 1–2 are materialized here (Tier 3 generation belongs to the harness #736).

1. **Experimental true-negatives** — the 1 hard + 8 soft rows already in the registry. Strongest claim, smallest n. (Already present; documented as Tier 1.)
2. **Presented decoys (from #681)** — the 13 splice peptides that are MHC-presented (immunopeptidomics) but have **no functional T-cell assay**. Real, presented, no measured response → materialized as `decoy_negatives/presented_decoys_681.tsv` (sequences already enumerated in the registry README, SNAF Supp Fig 4 + Supp Fig 7). This is the concrete #681 coordination.
3. **Matched synthetic** — length + allele-matched shuffled / decoy-junction peptides, drawn from a non-immunogenic presented background. Abundant ranking-denominator floor; weakest claim. **Construction algorithm specified here; generation deferred to #736** (it is a harness-time operation, not a fixed artifact).

**Hard rule:** tiers are never pooled as equal. Any score must report which negative tier(s) it used; a synthetic decoy is never counted alongside an experimental hard-negative as if equivalent.

## Acceptance-criteria mapping

| AC | Delivered by |
|---|---|
| #1 peptide→junction where annotated; unmappable flagged | `junction_id` + `junction_mapping_grade` columns + `none`-grade reasons |
| #2 documented labeling scheme committed | `LABELING_SCHEME.md` |
| #3 scheme applied → fixed labeled set + per-peptide rationale | updated `registry.tsv` + `label_rationale` column + re-audit changelog |
| #4 decoy-negative construction approach specified | `LABELING_SCHEME.md` §Decoy + materialized `presented_decoys_681.tsv` |

## Testing / verification

- **Schema integrity:** a small check that every row has a valid `junction_mapping_grade` and `evidence_strength` from the controlled vocabularies, and that `label`↔`tier`↔`evidence_strength` are mutually consistent per the scheme (e.g. no `strong` row whose `readout` is detection-only). Lives with the registry's existing validation if present; otherwise a short standalone script.
- **No-inference audit:** every `coords`/`event-id` row's junction value is traceable to a source location in `PROVENANCE.md`; spot-check that no coordinate was synthesized.
- **Re-audit log:** each boundary-row change records old → new label + the rule that drove it.
- **Decoy seed:** the 13 `presented_decoys_681.tsv` sequences match the README enumeration exactly (length, allele, gene).

## Risks

- **Re-audit scope creep** — re-examining boundary cases could spill into re-litigating verified rows. Mitigation: the re-audit is rule-driven and logged; a change requires a cited rule violation, not a fresh opinion.
- **`assay_context` coupling** — leaving the column to #823 means `confidence` stays partly conflated until then. Accepted: the scheme documents the decomposition so #823 is mechanical.
- **Coordinate sparsity reads as failure** — most rows landing at `gene-mechanism`/`none` could look like incomplete work. Mitigation: framed as a finding with a per-grade tally, consistent with the registry's "thin base" thesis.
