# Design note: two-resolution registry (junction + peptide) with a nullable peptide

**Status:** proposal (for board decision)
**Created by:** Scientist
**Relates to:** [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) (registry epic), [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) (detection benchmark), [#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911) (measured true-negatives), [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817) (triggering case: Zhao 2025)

## Motivation

The registry is peptide-keyed: every row needs an amino-acid sequence to exist.
That makes a whole class of otherwise-eligible sources unusable.

The triggering case is Zhao et al. 2025 (*Cell Research*, DOI `10.1038/s41422-025-01199-0`, [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)).
Its Supplementary Table S1 lists **139 candidate AS antigens by genomic coordinate + gene + per-allele HLA affinity**, with **zero peptide sequences** anywhere in the supplement.
The sequences are only in the paywalled main text or unpublished; recovery is gated behind an author request.
Peptide-keyed, Zhao contributes 0 rows today. Junction-keyed, it contributes 139 candidate rows immediately.

## The asymmetry (and why the naive framing is a red herring)

The intuitive argument is "peptide -> coordinate is easier than coordinate -> peptide, so retro-translate."
That is only partly true and is not the real lever.

- **Coordinate -> peptide is non-unique and is forbidden here.**
  A coordinate span does not determine the peptide: the exact junction breakpoint, the reading frame (in-frame vs frameshift), and which specific presented k-mer was synthesized are all lost.
  Zhao's own Table S1 proves it: several coordinates appear twice with different affinity vectors (one junction -> multiple distinct peptides).
  Deriving peptides from Table S1 coordinates is exactly the operation [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)'s hard constraint bans (fabrication / verify-against-source), and it mirrors the registry's existing no-inference rule (`junction_id` is never inferred from the peptide; `LABELING_SCHEME.md` section 6).
- **Peptide -> coordinate is tractable but does not help the blocked cases** (there is no peptide to retro-translate when the source is sequence-blocked).
- **The real lever:** the junction is the *natively published* unit of a splice-neoantigen study.
  These papers detect aberrant junctions from RNA-seq; the coordinate is their primary object.
  Keying on the junction unblocks every coordinate-first source, independent of the translation direction.

## What already exists (build on it, do not reinvent)

The schema half-anticipated this axis:

- `junction_id` (column 18): the most specific junction identifier the source published (coordinate pair, transcript/event id, or empty), never inferred.
- `junction_mapping_grade` (column 19): `coords | event-id | gene-mechanism | none`.
- A no-inference rule (`LABELING_SCHEME.md` section 6) already forbidding coordinate<->peptide fabrication in the other direction.

What is missing is (a) the junction as a *first-class key* rather than an optional attribute, (b) support for peptide-null rows, and (c) a junction-resolution evaluation path.

## Proposal

### 1. Registry entity = `(junction_id, peptide?, hla?, readout, provenance)`; junction becomes the stable key

The peptide becomes nullable.
A canonical, always-present `junction_id` (harmonized coordinate + genome build, or a canonical event id) becomes the row identity.
Dedup keys on `(junction_id, peptide, hla)` with `peptide` allowed null.

### 2. A missing peptide is a *typed* null, never a bare blank

A blank cell conflates at least three states that demand different follow-up:

| `peptide_status` | meaning | follow-up |
|---|---|---|
| `published-recovered` | sequence in hand from an authoritative source | none |
| `published-pending` | sequence exists, not yet obtained (e.g. Zhao, awaiting authors) | keep chasing |
| `unpublished-idonly` | source only ever gave coordinates / internal IDs | author contact only |
| `na-junction-level` | never a peptide (e.g. a normal-tissue true-negative junction) | none |

This extends an *existing* convention rather than inventing one: the long-read-UM rows (SEPTIN6, AMZ2P1) already carry an empty `hla` + `tier = functional-nonscorable` because the restriction was unpublished.
A null peptide is the same move one column over, with the additions that identity shifts to the junction and absence carries an explicit reason code.
A null peptide correctly nulls its derivatives (`length`, the scorable tier).

### 3. Two evaluation resolutions, one per pipeline stage

The pipeline has two evaluable stages living at different resolutions; the benchmark should not collapse them.

- **Peptide-level** (presentation / immunogenicity scorer, the [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) purpose): filter `peptide IS NOT NULL`.
  Immunogenicity is a property of the peptide-HLA pair and is irreducible; our MHCflurry scorer consumes peptides.
- **Junction-level** (detection + tumor-specific filtering + burden ranking, the [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) axis): use all rows.
  Validated tumor junctions are high-confidence positives for "did we detect and rank the right junction."

Switching the *scorer* benchmark target from peptide to coordinate would quietly measure a different thing (junction immunogenic-potential, not peptide immunogenicity) and could not score a peptide-level predictor without an aggregation rule. So this is additive, not a replacement.

## Payoffs

- **Unblocks Zhao today** at junction resolution (139 candidate rows), independent of the sequence recovery.
- **Attacks the scarcest reagent, measured true-negatives ([#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911)).**
  Junctions present in normal tissue (GTEx and friends) are coordinate-native and abundant; peptide-keyed, true-negatives are near-unsourceable.
- **Aligns with how the field's tools evaluate.**
  SNAF (`10.1126/scitranslmed.ade2886`), IRIS (`10.1073/pnas.2221116120`), and ASNEO (`10.18632/aging.103516`) are junction-discovery pipelines: the junction is the detection unit, the peptide is downstream (grounded in our Zotero shelf; a closer read of their exact eval metrics is a to-do if we want metric-level alignment).

## Costs and caveats

- **Coordinate harmonization is the new work.**
  Genome build (GRCh37/38 liftover), 0- vs 1-based conventions, and donor/acceptor vs anchor-outer semantics (cf. the regtools BED12 anchor-outer gotcha in `CLAUDE.md`).
  This is deterministic and scriptable (liftover + a canonical junction-ID scheme), unlike paywall-gated peptide recovery.
- **Junction-level label semantics are a compound claim.**
  A junction positive asserts detected + translated + presented + immunogenic; for a *detection* benchmark we only need "this is a real tumor-specific junction," so validated-immunogenic junctions are a high-confidence subset of real junctions.
- **The scorer still needs peptide labels**, so peptide-level rows remain the [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) core; this axis is complementary.

## Open questions (for the board)

1. Canonical `junction_id` scheme: coordinate normalization + build, or a standard junction-ID (e.g. `chr:donor-acceptor:strand` on GRCh38)?
2. Exact dedup / merge rule when two studies report the same junction with different peptides.
3. How a junction-level positive is defined for the [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) detection benchmark (any validated peptide -> junction positive?).
4. Scope + ownership: a new leaf sub-issue under [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680), cross-linked (native `blockedBy`) to [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) and [#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911).
