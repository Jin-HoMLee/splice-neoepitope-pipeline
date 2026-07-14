# Open splice-neoepitope → measured T-cell immunogenicity registry

**Issue:** [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) (seed); grown by [#733](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/733) (standing-watch + library-sweep additions) and [#734](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/734) (4-DB no-category audit + IEDB free-text recovery) · **Status:** seed (23 rows, 2026-06-11) + #733 pass (→ 44 rows): IR-CRC (Manoharan) + exon-TE (Merlotti) folded; Tier-2 triaged + #734 pass (→ 79 rows): Bigot 2021 SF3B1-UM folded (5 `high` + 30 `medium`); 4-DB no-splice-category finding documented in [`issue_734_db_audit/`](issue_734_db_audit/) + [#838](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/838) pass (→ 96 rows): POSTN sequence-corrected (`KTEGPTLTK`→`TVYTTKIITK`), GNAS/RPL22 promoted non-scorable→scorable, Kim mis-splicing (+13 scorable, coords-grade) + long-read UM (+4 non-scorable, HLA-unresolved) folded from the #734 IEDB-recovered candidates + [#904](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/904) pass (96 rows; scorable 79 → 81): 2 IRIS `candidate` rows (MAN2C1/AP2A1) promoted → `functional-scorable` (coords-grade) after locating the IRIS supplement on the Zotero **parent** item + [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) standing-watch pass (97 rows; scorable 81 -> 82): COL6A3-FLNV (`FLLDGSANV`, A\*02:01, Kim/Immatics 2022 tumor-stroma, `cloned_tcr`/`medium`) folded, 3-vote verified · **Snapshot:** 2026-07-07

> **📊 Reviewing this arc?** [`slides.qmd`](slides.qmd) is the 10-slide human review gate over the registry arc: what the registry holds, the A\*02:01 skew, the peptide-keying wall, the two-resolution schema that removed it, and the calls that need a human's judgment. Render with `quarto render slides.qmd`; figures regenerate from `registry.tsv` via [`figures/_regenerate_figures.py`](figures/_regenerate_figures.py), so the registry stays canonical and no slide number is hand-typed.

## What this is

The first openly-assembled, provenance-tracked registry of **splice-junction-derived (aberrant-splicing / neojunction / intron-retention) cancer neoantigens with experimentally MEASURED T-cell immunogenicity**. Each entry passed a two-gate test, each gate adversarially verified (3 independent skeptic votes, ≥2/3 required): (1) genuinely splice-derived - not SNV / indel / fusion / canonical-TAA / proteasomal-cis-spliced; (2) an actual functional T-cell assay was performed on that exact peptide (IFN-γ, tetramer/dextramer, cytotoxicity, ELISpot, or engineered-TCR readout).

## What this is NOT - a probe, not a powered benchmark

A feasibility gate established that only ~tens of such peptides exist field-wide (TESLA explicitly excluded splice isoforms; the four major immunogenicity DBs - CEDAR/IEDB/NEPdb/dbPepNeo2.0 - have **no splice category**, rigorously documented in [`issue_734_db_audit/`](issue_734_db_audit/)). **This is a small, transparently-underpowered stress-test probe**, not a TESLA-scale calibrated benchmark. Its value is being first to assemble an open, sequence-resolved splice-immunogenicity set with documented provenance - and the honest demonstration that the field's functional-validation base is critically thin. Source clustering is real: the entries come from ~a dozen studies, and **73/97 rows are HLA-A\*02:01** (see Caveats).

## Sparsity analysis (#737)

The scarcity is quantified rigorously in [`issue_737_sparsity/sparsity_writeup.md`](issue_737_sparsity/sparsity_writeup.md) (manuscript-ready prose), computed by [`issue_737_sparsity/notebook.ipynb`](issue_737_sparsity/notebook.ipynb) → [`issue_737_sparsity/outputs/`](issue_737_sparsity/outputs/), and condensed in the experiment deck [`issue_737_sparsity/slides.qmd`](issue_737_sparsity/slides.qmd). Headline numbers (on the **81 scorable positives**, computed on the n=81 snapshot 2026-06-30):

> ⚠️ The COL6A3-FLNV addition (row 82, A\*02:01, +1 study; 2026-07-07) is not yet folded into these effective-count figures. The qualitative conclusions (A\*02:01 monoculture, ~1 hard negative, few effective studies) are unchanged; a recompute of the #737 sparsity analysis is tracked as a follow-up.

- **Few-study assembly:** 10 studies, but top study = 43%, top two = 65%, **effective ≈ 3.8** independent studies (inverse-Simpson of the per-study shares).
- **A\*02:01 monoculture:** 72/81 (89%) are A\*02:01 → **effective ≈ 1.26 alleles**. Mechanism spread is healthier (effective ≈ 4.2).
- **Negatives are the binding constraint:** exactly **1** hard true-negative field-wide (+ 8 soft failed-to-prime); a powered AUC probe needs **19-31** negatives (Hanley-McNeil), so field-wide specificity is currently unmeasurable.
- The two rate-limiting reagents named by the analysis: non-A\*02:01 functionally-validated positives ([#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839)) and measured true-negatives ([#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911)).

**Category-absence ≠ data-absence (#734).** The no-splice-category finding holds at the *schema* level - you cannot retrieve splice-neoantigens as a class from any of the four DBs (even deposited ones, e.g. IRIS-CLASP1 and RCAN1-4 in IEDB, are filed as plain protein "isoforms"; the splice origin survives only in the free-text reference title). **But free-text/reference-title mining of IEDB's IQ-API *does* recover real splice-neoantigens the manual literature sweep missed** - most importantly the Bigot 2021 SF3B1-mutant uveal-melanoma panel (+35 rows, folded here). So the DBs hold recoverable splice-immunogenicity data; it is simply uncategorized. The thin-functional-base point still stands: most recovered rows are tetramer-detection-only (`medium`), not effector-confirmed.

## Tiers (`registry.tsv`)

| Tier | n | Meaning |
|---|---|---|
| `functional-scorable` | 82 | Splice-derived + functional T-cell assay + exact sequence & allele → directly usable to score a predictor. Confidence sub-graded `high` (32) / `medium` (50) (e.g. healthy-donor IVS, tetramer-detection-only, or `Positive-Low` ELISPOT → `medium`). Includes the #838 Kim mis-splicing fold (+13, A\*02:01, coords-grade junctions) the #904 IRIS promotion (+2, MAN2C1/AP2A1, coords-grade), and the #680 standing-watch COL6A3-FLNV fold (+1, A\*02:01, engineered-TCR / `cloned_tcr`, `medium`). |
| `candidate` | 0 | Was 2 (`LLLGIAKLLKV`/MAN2C1, `FLSELEPPA`/AP2A1, IRIS) - **both promoted → `functional-scorable` in [#904](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/904)** once the IRIS supplement was located on the Zotero parent item; the prior "2nd-verify did not confirm" hold was a supplement-access failure (parent-key trap), not negative evidence. |
| `candidate-negative` | 8 | IR-CRC (Manoharan) IR-derived 9-mers **tested-negative** in healthy-donor IVS → *soft* negatives (failed to prime ≠ intrinsically non-immunogenic). A decoy-negative seed, categorically weaker than the hard negative. |
| `functional-nonscorable` | 4 | Functionally validated but not scorable. The 4 long-read-UM peptides (SEPTIN6, AMZ2P1, MZT2B; #838) - splice-confirmed + ELISPOT-positive but **no per-peptide HLA published**, so no allele key. (GNAS/RPL22, previously here, were promoted to `functional-scorable` in #838 once their sequences were recovered from the Kwok IEDB deposit.) |
| `presentation-prevalence` | 1 | `TEFQTRRAM`/SLC45A2 - a real splice neoantigen with prevalence/survival evidence + predicted binding, but **no functional T-cell assay performed** (distinct from the 5 SNAF IFN-γ-tested peptides) |
| `hard-negative-true-splice` | 1 | `VELEDHVML` - splice-derived, MS-presented, experimentally NON-immunogenic. The gold-standard discrimination negative. |
| `negative-control-not-splice` | 1 | `VFVDGLCRAKF` - RCAN1-1 constitutive same-locus control (fails splice gate) |

**Two machine-readable columns added in #733 (AC#4):** `splice_mechanism_canonical` (controlled vocabulary - see below) and `hla_resolution` (`2-digit` for SNAF C\*04/C\*08 + the SNAF A\*02 prevalence row, `4-digit` otherwise).

### `splice_mechanism_canonical` controlled vocabulary

Free-text `splice_mechanism` is preserved; `splice_mechanism_canonical` normalizes it for grouping:

`alt_5p_ss` (alt 5′ splice site / donor) · `alt_3p_ss` (alt 3′ splice site / acceptor) · `alt_splice_junction` (alternative-splicing junction, sub-event unspecified - IRIS/Kim/Xiong/POSTN) · `intron_retention` · `minor_intron` (U12-type) · `poison_exon` (NMD-targeting cassette) · `neojunction_frameshift` · `neojunction_inframe` · **`exon_te_junction`** (exon–transposable-element junction or TE-induced frameshift - **registered in #733** for the Merlotti set) · `not_splice` (constitutive control) · `ambiguous` (splice-derived, mechanism unresolvable from source - Fisher unnamed-gene peptides).

### New columns added in #735

[`LABELING_SCHEME.md`](LABELING_SCHEME.md) is the authoritative rule set for all four columns; the entries below are cross-references.

- **`evidence_strength`** (`strong | weak | hard | soft | na`) - grades the type of experimental evidence for the row's label, de-conflating assay type from provenance context.
  - `strong` = measured effector readout (IFN-γ, cytotoxicity, granzyme-B, CD107, 4-1BB, TCR-T activation).
  - `weak` = antigen-specific detection only (tetramer / dextramer).
  - `hard` = splice-derived, MS-presented, and functionally tested negative (gold-standard discrimination negative).
  - `soft` = tested negative in healthy-donor in-vitro sensitization (IVS), categorically weaker than hard.
  - `na` = not applicable, specifically the `presentation-prevalence` tier (untested) and the `negative-control-not-splice` tier (fails the splice gate). The `functional-nonscorable` and `candidate` tiers are non-scorable too but keep their real `strong`/`weak` strength, not `na`.
  - Full rule definitions: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §2-4.
- **`label_rationale`** - one-line human-readable rationale for the `evidence_strength` assignment, citing the key readout keyword or the rule that applies.
- **`junction_id`** - the most specific junction identifier the source published for that peptide: a genomic coordinate pair, a transcript/event identifier, or empty when only gene + mechanism is recoverable.
Never inferred from peptide sequence or genome reference (no-inference rule: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §6).
- **`junction_mapping_grade`** (`coords | event-id | gene-mechanism | none`) - grades how specifically the peptide is traceable to its originating junction.
`coords` = explicit genomic donor/acceptor coordinates published by the source; `event-id` = transcript or splice-event identifier published (not full coordinates); `gene-mechanism` = only gene name and splice mechanism recoverable; `none` = no junction information recoverable (reason in `notes`).
Full ladder + no-inference rule: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §6.
Per-source grade rationale: [`junction_evidence_by_source.md`](junction_evidence_by_source.md).

### New column added in #823

- **`assay_context`** (`patient_exvivo | patient_til | healthy_donor_ivs | cloned_tcr | animal_syngeneic | prevalence_only | unspecified | na`) - **a T-CELL-SOURCE axis and nothing else** ([#1120](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1120)); the in-vivo/in-dish *setting* lives in the separate `in_vivo_model` column below - records *which immunological system* produced the functional readout, so a scoring run can weight rows by assay realism (patient ex-vivo detection > healthy-donor IVS priming > engineered-TCR readout) **programmatically**, instead of re-parsing free-text `notes`. Previously this distinction lived only in prose, yet it already drives the `high`/`medium` confidence split (Manoharan IVS → `medium`; patient tetramer → `high`/`medium`).
  - `patient_exvivo` (42) = patient PBMC/blood ex-vivo tetramer⁺ or functional. `patient_til` (3) = patient tumor-infiltrating (or draining-LN) lymphocytes. `healthy_donor_ivs` (11) = healthy-donor in-vitro-sensitized T cells (all Manoharan IR-CRC, positives + negatives). `cloned_tcr` (4) = engineered/cloned-TCR functional readout, no primary patient/donor detection (IRIS JPTCR clones). `prevalence_only` (1) = population prevalence/presentation, no per-peptide T-cell assay (the SNAF `TEFQTRRAM` row). `na` (1) = constitutive non-splice control.
  - **`unspecified` (17)** = the assay was reported but the held provenance does **not** determine the T-cell source (SNAF / Kim / Xiong / Fisher / POSTN / Kwok). Per the verify-against-source rule the patient-vs-donor system is **not guessed**; a scoring run should apply no context weight to these rows. Pinning them is a future read-the-methods pass, not a fabrication.
  - Assignment is **source-keyed** (the `source` column is the stable curation key; rationale per source in [`PROVENANCE.md`](PROVENANCE.md)), with Merlotti's ex-vivo-vs-TIL split read from the per-context detail in `notes`. Rule: [`derive_assay_context.py`](derive_assay_context.py); enforced (vocabulary + `healthy_donor_ivs`↔IVS-marker + `prevalence_only`↔tier cross-checks) by [`validate_registry.py`](validate_registry.py).

### New column added in #1001

- **`in_vivo_model`** (`none | xenograft | syngeneic | unspecified`) - records the **setting** of the functional readout (in a dish, or in a living animal), **orthogonal** to `assay_context`'s T-cell-source axis. Added in [#1120](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1120), which was originally filed asking for an "animal model" value *inside* `assay_context` - a **category error**: an animal is a *venue*, not a T-cell source. Both of our in-vivo rows (`IFSESETRAKF`/RCAN1-4, `FLLDGSANV`/COL6A3) are **human T cells inside an immunodeficient NSG mouse** (verified first-hand; see [`PROVENANCE.md`](PROVENANCE.md)), and **an NSG mouse cannot produce T cells at all** - so an "animal" T-cell source would be false for both, and on COL6A3 would have overwritten `cloned_tcr`. A row can be `cloned_tcr` **and** in-vivo-confirmed; one column cannot hold both facts.
  - **Precedence: an in-vivo setting never overrides the T-cell source.** Derived by [`derive_in_vivo_model.py`](derive_in_vivo_model.py), gated on the in-vivo marker in `readout` (so a source's non-in-vivo rows - e.g. Xiong's constitutive `VFVDGLCRAKF` control - are never stamped) and source-keyed for the model *type*, with a no-guess `unspecified`. [`validate_registry.py`](validate_registry.py) enforces the vocabulary, the `readout`-marker cross-check **in both directions**, and the coupling `animal_syngeneic` ⟹ `syngeneic` (an immunodeficient host has no T cells to be the source of).
  - `animal_syngeneic` (in `assay_context`) is the **only** sense in which an animal is a T-cell source - an immunocompetent host responding with its own T cells. **0 rows today**, deliberately forward-looking for a Burbage-2023-shaped mouse exon-TE fold ([#699](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/699)).
  - **Realism is not monotone in this column alone:** a xenograft result is strong evidence of effector function but runs in an immunodeficient host with transferred human T cells, so it says nothing about whether an intact immune system would respond. Read it *with* `assay_context`, never instead of it.

- **`venue_type`** (`journal | preprint | db_recovered`) - records the **publication venue class** of each row's source, so peer-reviewed evidence is **queryable, enforceable, and auditable at a glance** instead of living only in free-text `source` / PROVENANCE prose. Standard evidence-grading practice (GRADE, systematic-review methodology) down-weights non-peer-reviewed sources; marking venue lets a scoring run flag (and optionally down-weight or hold out) preprint-sourced rows in a sensitivity analysis.
  - **All 97 current rows are `journal`** (0 preprints), per the first-hand venue audit in [`PROVENANCE.md`](PROVENANCE.md) (2026-07-04). This is *not* a policy against preprints - a preprint row can be legitimately folded, it just must be marked. The three preprints encountered so far (Lin 2025, TEtrans/Li 2025, the SNAF exon-TE preprint) were each deferred for **sequence-unavailability**, not because they are preprints. `db_recovered` is reserved for a future row recovered from a public deposit (IEDB/CEDAR) whose underlying publication venue is ambiguous; none exist today.
  - Assignment is **source-keyed** (same convention as `assay_context`; the source map lives in [`labeling_constants.py`](labeling_constants.py), applied by [`derive_venue_type.py`](derive_venue_type.py)). **Two guards** in [`validate_registry.py`](validate_registry.py) catch a mis-marked venue: (1) a source matching **no** map key derives to an out-of-vocab sentinel (`unclassified`) that the validator **rejects** - so a genuinely new source cannot slip in venue-unmarked; and (2) because a study-substring key cannot self-distinguish a preprint from the journal version of an *already-mapped* study (a future `"SNAF ... bioRxiv"` source would first-match `"snaf"` → `journal`), the validator separately **fails** any row whose `source` names a preprint server (`PREPRINT_MARKERS`) but whose `venue_type` isn't `preprint`. The validator also prints the **preprint-row count** as an advisory line (reports, does not reject).
  - **Zotero cross-check (offline-optional).** Rather than *derive* venue from Zotero live (which would couple this committed column to a network service and break hermetic CI / offline runs), the value stays a repo-local map and [`validate_registry.py`](validate_registry.py) **cross-checks** it against Zotero's authoritative `itemType`: when `ZOTERO_*` creds are present it resolves each source by DOI in collection `Z38GTJNW` and fails on any `venue_type`↔`itemType` disagreement; with no creds the check skips cleanly (hermetic runs stay green). DOI-keyed, so it sidesteps the parent-vs-attachment-key trap. This is the [`annotate-canary`](../../../.github/workflows/annotate-canary.yml) pattern - an independent source of truth wired as a *verification*, not a dependency - and it is what caught the initial hand-entry errors in the venue-name audit.

### New column added in #1086 - the registry became two-resolution

- **`peptide_status`** (`published-recovered | published-pending | unpublished-idonly | na-junction-level`) - types the *absence* of a peptide, so a blank sequence cell no longer conflates "sequence exists, keep chasing the authors" with "the source only ever published coordinates".

The registry used to be **peptide-keyed**: every row needed an amino-acid sequence to exist, which locked out every coordinate-first source. `junction_id` and `peptide` are now **both nullable** under an at-least-one-non-null invariant, and row identity is the coalesced key `COALESCE(junction_id, peptide)` disambiguated by `hla`. Full rules + the harmonization with `tier` / `junction_mapping_grade`: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) section 7.

- **All 97 current rows are `published-recovered`** - every row in the registry today carries its sequence. The schema change is what makes a *future* coordinate-first fold possible, most importantly coordinate-native measured true-negatives ([#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911)), the scarcest reagent in the set.
- ⚠️ **Zhao 2025's Table S1 is NOT the first fold, contrary to what this README and the design doc used to say.** The first-hand audit ([#1089](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1089)) found its 139 rows carry only *predicted* binding, so every row is `label=untested` with **no legal tier**, and they resolve to **103 distinct junctions**, not 139 (65 rows share a coordinate). Nullable identity was necessary to admit a coordinate-first source but is not sufficient - the row still needs a legal `(label, tier)`. Whether a new predicted-only tier is created: [#1125](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1125). The *validated* Zhao antigens remain sequence-blocked and are additionally not joinable to Table S1's coordinates ([#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)). Full audit: [`PROVENANCE.md`](PROVENANCE.md).
- The grade→`junction_id` rule was **replaced, not relaxed**: a `coords`/`event-id` grade still requires a `junction_id`, but a `gene-mechanism`/`none` row may now carry one recovered later from an authoritative source. The old rule barred 64% of rows (the 62 `gene-mechanism`/`none` rows) from ever gaining a junction. The no-inference rule is untouched.
- Dedup keys on the **full `(junction_id, peptide, hla)` triple**, not on the coalesced identity: one junction legitimately yields several distinct peptides. Exactly one such case exists today - the Kim 2025 constitutive-intron event `ci@16:719606:720123:+|16:719606:719607:+` carries `FLWPGLGPS`, `FLWPGLGPSV`, and `ILGSLTWSC`. `registry_dedup.junction_view()` gives junction-level consumers that grouping.
- ⚠️ **`junction_id` is stored source-verbatim and no genome build is recorded**, so cross-source junction comparison is **not yet possible**. The 17 `coords` rows span four incompatible formats, two of which encode multi-junction *events*. Canonical scheme + build audit: [#1100](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1100).

## Coverage (97 rows, 12 studies, 75 genes)

- **Alleles (11 distinct `hla` values + 4 HLA-unresolved):** HLA-A\*02:01 (73), A\*11:01 (9), A\*24:02 (2), **C\*04 (2), C\*08** (HLA-C breadth from SNAF supplements), + A\*34:01, A\*31:01, A\*24:07, A\*02:07, A\*02:06 (allele *diversity* added by IR-CRC), + the 2-digit **A\*02** (the `TEFQTRRAM` prevalence row), + **4 blank** (the #838 long-read-UM rows with no per-peptide HLA published). **A\*02:01 skew deepened to 73/97** - the Bigot 2021 SF3B1-UM panel (+35), the IR-CRC positives, the Merlotti exon-TE set, and now the #838 Kim mis-splicing fold (+13) are all A\*02:01/A2 (see Caveats; per-allele table + scoring implication in "Per-allele coverage and the A\*02:01 scoring skew (#839)" below; rebalance tracked in [#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839)).
- **Splice mechanisms:** + `alt_3p_ss` (Bigot SF3B1 aberrant 3′SS, 30) and `neojunction_frameshift` (Bigot S8-junction-validated frameshift+NMD, 5) dominate alongside `exon_te_junction` (Merlotti, 10) and `intron_retention` (IR-CRC, 12); the #838 Kim fold adds 13 `alt_splice_junction`/`alt_5p_ss` rows (skipped-exon, constitutive-intron, mxe, a5ss events).
- **Sources (12 studies, 12 source strings - one string per study):** SNAF (Li 2024, Sci Transl Med), IRIS (Pan/Xing, PNAS), Kim 2025 (SF-mutant leukemia, expanded +13 via #838), Xiong 2025 (GBM), Fisher 2026 (CoREST), Kwok 2024 (Nature), POSTN-203 study, Manoharan 2026 (IR-CRC), Merlotti 2023 (NSCLC exon-TE), Bigot 2021 (SF3B1 uveal melanoma, +35 - #734 IEDB mine), **Long-read UM 2023 (Cancer Immunol Res, SEPTIN6/AMZ2P1/MZT2B, +4 non-scorable - #838)**, Kim GB 2022 STM (COL6A3 tumor-stroma, +1 - #680 standing watch). (The Kim string was canonicalized in #838; the remaining SNAF/IRIS/Xiong 2-string splits were canonicalized in [#1106](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1106), taking 15 strings -> 12. `source` is the study group-by key, so a split study was silently double-counted; `validate_registry.py` now **fails** when two `source` spellings resolve to the same substring key, so the 1:1 mapping cannot silently re-break. Rewrite recorded in [`PROVENANCE.md`](PROVENANCE.md).)
- **Genes (75 distinct strings):** + Bigot SF3B1-UM source genes - NF1, USP39, NET1, ATP8B2, MAPK8IP2 (the 5 S8-junction-validated), plus ZFYVE27, SRSF1, SF3A2, SEPSECS, UBA1, CTDNEP1, ARIH1, MRPS10, VPS51, MINDY4, ZDHHC16, SOAT1, OXA1L, HADHA, SLC3A2, RNF38, ATP5MC2, VARS2 (tetramer-only); + the #838 Kim genes USF1, AP4B1, EZH2, BIN3, CASP2, SLC20A1, RHOT2, NCOA7, LTBR, MAN2C1, SH3GL1 + the long-read-UM genes SEPTIN6, AMZ2P1, MZT2B + COL6A3 (Kim/Immatics tumor-stroma). (5 Bigot tetramer-only peptides have no source gene in IEDB → `not-named`.)

### Per-allele coverage and the A\*02:01 scoring skew (#839)

Per-row `hla` tallied across all 97 rows, split by scorability.
The **scoring-relevant** column is the 82 `functional-scorable` positives (the set a predictor is actually scored on); the total column includes the negative tiers and the 4 HLA-unresolved rows.
(Recompute from this directory (`conda activate snakemake` first). Both columns, including the load-bearing scorable split (`label=='positive' and tier=='functional-scorable'`, per [`LABELING_SCHEME.md`](LABELING_SCHEME.md)):
`python -c "import csv,collections as c; rows=list(csv.DictReader(open('registry.tsv'),delimiter=chr(9))); tot=c.Counter((r['hla'].strip() or '(unresolved)') for r in rows); sco=c.Counter(r['hla'].strip() for r in rows if r['label']=='positive' and r['tier']=='functional-scorable'); print('total',dict(tot)); print('scorable',dict(sco))"`
The second Counter reproduces the 73/82 A\*02:01 and 8/82 non-A\*02 figures the analysis pivots on; the first reproduces the Registry-rows column.)

| Allele | Registry rows | Scorable positives | Role of the non-scorable rows |
|---|---:|---:|---|
| HLA-A\*02:01 | 73 | 73 | - |
| HLA-A\*11:01 | 9 | 4 | 1 hard + 4 soft negatives (the registry's true-negative set) |
| HLA-A\*24:02 | 2 | 1 | 1 RCAN1 constitutive non-splice control |
| HLA-C\*04 | 2 | 2 | - (HLA-C breadth, SNAF) |
| HLA-A\*02:06 | 1 | 1 | - (A\*02 family) |
| HLA-C\*08 | 1 | 1 | - (HLA-C breadth, SNAF) |
| HLA-A\*02:07 | 1 | 0 | IR-CRC allele-diversity (non-scorable) |
| HLA-A\*02 (2-digit) | 1 | 0 | `TEFQTRRAM` prevalence-only row |
| HLA-A\*24:07 | 1 | 0 | IR-CRC allele-diversity (non-scorable) |
| HLA-A\*31:01 | 1 | 0 | IR-CRC allele-diversity (non-scorable) |
| HLA-A\*34:01 | 1 | 0 | IR-CRC allele-diversity (non-scorable) |
| (unresolved) | 4 | 0 | long-read UM, no per-peptide HLA published |
| **Total** | **97** | **82** | |

**Family rollup on the scoring set:** the whole A\*02 family is **74/82 scorable positives (90%)**; A\*02:01 alone is **73/82 (89%)**.
Only **8 of 82 scorable positives (10%) are non-A\*02**, and they span just **four** alleles - A\*11:01 (4), C\*04 (2), A\*24:02 (1), C\*08 (1).
A\*11:01 is the only non-A\*02 allele with more than two scorable rows.

**Scoring implication.**
Any AUC, ranked-recall, or sensitivity figure computed on this registry is, to first order, an **A\*02:01-specific** result - a cross-allele generalizability claim is unsupportable, because seven of the eight non-A\*02 scorable rows sit at n≤2 per allele (no allele but A\*02:01 has enough rows to estimate a per-allele operating point).
A **balanced** per-allele scoring split is therefore not achievable from the current set; the only defensible stratification today is **A\*02:01 vs. pooled-non-A\*02:01** (73 vs. 8), and even that pooled arm is too thin for a stable estimate.
This is a documentation-and-sourcing matter, not a data-quality defect: every A\*02:01 row is legitimately validated (see "What this is NOT" and Caveats).
The remedy is to **add** non-A\*02:01 functionally-validated splice-neoantigens (the allele-targeted sweep below), never to remove A\*02:01 rows.

### Allele-targeted rebalance sweep - outcome (#839, AC2/AC3)

An allele-targeted hunt for non-A\*02:01 functionally-validated splice-junction-derived neoantigens was run 2026-07-09.
**Outcome: no new foldable rows.**
The field's functionally-validated splice-immunogenicity is genuinely A\*02:01-dominated; this is a real limitation of the field's validation base, not a curation gap.

Sweep performed:
- Zotero `Z38GTJNW` re-scan (our curated splice-immunogenicity shelf).
- Re-check of the #734 IEDB free-text recovered-candidate list ([`issue_734_db_audit/recovered_candidates.tsv`](issue_734_db_audit/recovered_candidates.tsv)) - zero unfolded non-A\*02 functional rows.
- Live IEDB `tcell_search` free-text mine (reference title matching `splic`), with proteasome-*cis*-spliced-peptide papers excluded (also gate-1-excluded); schema confirmed, endpoint rate-limited before exhaustive allele-stratified enumeration, so the completed #734 audit remains the authoritative IEDB splice coverage.
- Recency web sweep (2025-2026) for non-A\*02 functional splice-neoantigen reports.

Every non-A\*02 splice-immunogenicity signal the sweep surfaced was already accounted for:
- **Already folded** (the current 8 non-A\*02 scorable rows): RCAN1-4 (Xiong 2025, GBM, A\*24:02), POSTN-203 (A\*11:01), and the SNAF HLA-C entries (C\*04, C\*08).
- **Sequence-blocked - the one live rebalance lever:** Zhao 2025 (HCC AS-vaccine, *Cell Research*) reports patient-TIL functional hits restricted to A\*02:01 / **A\*11:01 / A\*24:02**, but the exact peptide sequences are paywalled. Recovery is tracked in [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817); if the sequences land, this is the single largest realistic non-A\*02 rebalance available.
- **Gate-1 excluded (documented boundary case):** survivin-2B (`AYACNTSTL`, HLA-A\*24, functional CTL killing, oral/colorectal cancer). It is splice-variant-derived and non-A\*02, but survivin-2B is a *naturally-occurring physiological alternative isoform* of the canonical shared tumor-associated antigen survivin (an established peptide-vaccine target; per the literature its expression *decreases* in later tumor stages), not a tumor-specific aberrant-splicing neojunction. It fails **gate-1** (canonical-TAA splice isoform, excluded by design) and is therefore not folded.

**Conclusion (AC3 - documented null).**
A meaningful allele rebalance is **not achievable from currently-published functional data**.
The registry's A\*02:01 monoculture (89% of scorable positives) reflects where the field's functional-validation effort has actually gone, and must be carried as a registry limitation: **any AUC / ranked-recall / sensitivity claim in the manuscript must be scoped to A\*02:01**, not framed as cross-allele generalization.
The single concrete future lever is [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817) (Zhao HCC sequence recovery, A\*11:01 / A\*24:02); the standing splice-immunogenicity watch remains allele-alert for new non-A\*02 functional reports.

### Junction-mapping coverage (#735)

Per-row `junction_mapping_grade` tallied across all 97 registry rows (derivable with `research/.venv/bin/python -c "import pandas as pd; d=pd.read_csv('research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv',sep='\t'); print(d['junction_mapping_grade'].value_counts().to_string())"`):

| Grade | Count |
|---|---|
| `gene-mechanism` | 58 |
| `event-id` | 18 |
| `coords` | 17 |
| `none` | 4 |
| **Total** | **97** |

Coordinate-level junction grounding jumped from 2 to 15 with the #838 Kim fold - the Kim supp (mmc2) publishes genomic junction coordinates per peptide, the first sizable `coords`-grade source - and to **17** with the #904 IRIS promotion (Dataset S3 publishes per-peptide junction coordinates). It is still a minority (17/97); most rows rest at `event-id` or `gene-mechanism`, mirroring the field's sparse coordinate-publishing practice.
Per-source rationale for each grade assignment: [`junction_evidence_by_source.md`](junction_evidence_by_source.md).

### Labeling re-audit (#735)

All 14 boundary rows were audited against the documented labeling scheme ([`LABELING_SCHEME.md`](LABELING_SCHEME.md)) after `evidence_strength` was derived algorithmically.
The 14 rows cover every non-standard tier: the two `candidate` rows, the two `functional-nonscorable` rows, the one `presentation-prevalence` row, the one `negative-control-not-splice` row, and the eight `candidate-negative` rows.
**Zero label changes resulted.**
The scheme codifies the existing adversarially-verified curation exactly; all `evidence_strength` assignments and `label` values matched prior manual decisions for every row.

## Caveats

- **Sequence not published (GNAS, RPL22):** exact 9-mers are figure-locked (MS-spectra peak plots) or only in a generic method schematic; supplementary tables hold no peptide list. Both A*02:01 (no allele-diversity loss). Recoverable only via author contact / raw MS.
- **Medium-confidence rows:** MAN2C1 / AP2A1 (IRIS) failed a second verifier pass; Fisher's RSQGWLFLR / AEHAHRVPL / VELEDHVML lack named genes and have ambiguous restriction (loaded A*11:01 but predicted B*40:01).
- **Readout heterogeneity:** several Kim entries are tetramer-only (antigen-specific T-cell detection) vs IFN-γ/cytotoxicity for others - recorded per-row.
- **Allele granularity:** SNAF HLA-C entries are 2-digit (C*04 / C*08); others are 4-digit.
- **Negatives:** one true splice *hard* negative (`VELEDHVML`, MS-presented + ELISpot-negative) + **8 IR-CRC *soft* negatives** (`candidate-negative`) added by #733 - failed-to-prime in healthy-donor IVS, categorically weaker than the hard negative but a usable decoy-negative seed for #681.
- **IR-CRC (Manoharan) provenance = healthy-donor IVS, not patient:** the 3 positives (CHD7/FERMT3/TPCN2) and 8 negatives were measured in **healthy-donor moDC-primed T cells** (in-vitro sensitization), *not* patient T cells → positives graded `medium` (a peptide that *can* prime ≠ primed in a patient). FERMT3/RQDPAPQQV was functionally restricted on **A\*02:01** only (A\*24:07 predicted, not used).
- **Merlotti (exon-TE) confidence split:** 3 with cloned-TCR functional validation (`high`); 7 with **patient tetramer-DETECTION only** (`medium` - no per-peptide IFN-γ/killing). ~10 further single-replicate patient tetramer⁺ hits were **not** folded as scorable rows (recorded in PROVENANCE as candidates).
- **Bigot 2021 (SF3B1-UM) confidence split (#734):** 5 with effector function (IFN-γ / cytotoxicity / granzyme B / CD107 / TNF-α) **and** S8 alternative-mRNA-junction validation (frameshift+NMD) → `high` (NF1, USP39, NET1, ATP8B2, MAPK8IP2). The other 30 are **patient ex-vivo tetramer-DETECTION only** → `medium`; their splice origin rests on the A2:N panel design (all SF3B1-aberrant-splicing-predicted) + IEDB source-gene curation, not a per-peptide junction shown in the supplements we hold. Sequences are `direct` from local supplementary **Table S2**, cross-confirmed against IEDB (35/35) and - for the 5 `high` - JEM 2024 (PMC10986814). **The published Correction (PMID 35257149 / CD-22-0009) could not be obtained; its content is unread.** Sequence risk is mitigated by the two/three-way cross-confirmation, but the correction should be reviewed before citing.
- **A\*02:01 skew (worsened sharply by #734, again by #838):** **73/97 rows are A\*02:01** (73/82 = 89% of scorable positives) - the Bigot SF3B1-UM panel (35), the IR-CRC positives, the Merlotti exon-TE set, and the #838 Kim mis-splicing fold (13) are all A2. The HLA-C/allele-diversity story is a small fraction of the whole; weight scoring analyses accordingly. **The [#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839) allele-targeted rebalance sweep (2026-07-09) returned a documented null** - no new foldable non-A\*02 functional splice-neoantigens exist in currently-published data (full per-allele table + sweep outcome under "Per-allele coverage and the A\*02:01 scoring skew (#839)"). The one live lever is [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817) (Zhao HCC A\*11:01/A\*24:02, sequence-blocked); scope all generalizability claims to A\*02:01.

## Decoy-negative tiers (#735)

The benchmark's negative set is built from three tiers of decreasing claim strength.
Tiers must never be pooled as equal (full rule: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §7).

| Tier | n | Claim | Status |
|---|---|---|---|
| **Tier 1** - Experimental true-negatives | 9 | Splice-derived + functional T-cell assay + no response measured | Materialized in `registry.tsv` (1 `hard` + 8 `soft`) |
| **Tier 2** - Presented decoys | 13 | Splice-derived + MHC-presented (immunopeptidomics) + no functional assay performed | Materialized in [`decoy_negatives/presented_decoys_681.tsv`](decoy_negatives/presented_decoys_681.tsv) |
| **Tier 3** - Matched synthetic | TBD | Length- and allele-matched shuffled / decoy-junction peptides drawn from a non-immunogenic presented background | Generation deferred to the #736 scoring harness (runtime, not a committed artifact) |

**Tier 1** contains the 1 hard-negative (`VELEDHVML`, MS-presented + ELISpot-negative, `evidence_strength=hard`) and the 8 IR-CRC soft-negatives (failed to prime in healthy-donor IVS, `evidence_strength=soft`).
**Tier 2** is the 13-peptide SNAF MS-presented seed coordinated via [Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681).
These are real presented splice-junction-derived peptides with no measured response - not tested-negative, but untested.
**Tier 3** serves as the abundant ranking-denominator floor.
Its construction algorithm is specified in [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §7, and generation runs at harness time ([#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736)), not here.

Any score, AUC, or ranked-recall metric **must** report which tier(s) it used.
A Tier-1-only result (9 rows) is a materially different claim from a Tier-1+2 result (22 rows).

## Deferred - MS-presented / immunogenicity-untested tier (→ Issue #681)

13 distinct splice-junction-derived neoepitopes are MHC-presented (immunopeptidomics) but have **no functional T-cell assay** - they fail the functional gate and are NOT in this registry.
These 13 form the Tier-2 decoy-negative seed ([`decoy_negatives/presented_decoys_681.tsv`](decoy_negatives/presented_decoys_681.tsv), materialized in #735) and are coordinated via [Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681). They belong with #681 (public immunopeptidome mining) and are a natural feeder for the decoy negative set. Source: SNAF Supp Fig 4 (`NQDEDPLEV`/C6orf52, `KGPWYPLSL`/C20orf204, `VAPGEAKNL`/RASA3, `YALANIKWI`/DYNLT5, `KEKLDQLVY`/FBXO7, `TELQRTLSL`/NGLY1) + Supp Fig 7 (`SQTPKSRAL`/PSMF1, `RKLEAPYLL`/MCF2L, `LSWPRSTPM`/CPN1, `VSTGCAVVL`/SERPINE2, `RRLPNPPAV`/RGS12, `IVKRPRSEL`/EXO1, `AVPLLQTNR`/ETV4).

**#733 Tier-2 triage → #681.** Four library-sweep papers were read for a per-peptide functional T-cell assay and routed away from the functional registry (functional gate absent for Courcelles / Lin / Bathini; gate-PASS-but-sequence-blocked for Zhao):
- **Courcelles 2026** (CRC MSI/MSS, MCP) - splice-derived aeTSAs (intron retention, retroelement junctions) but immunogenicity is *in-silico only* (PRIME/ImmuneApp) → **#681**.
- **Lin 2025** (HCC, bioRxiv) - neoTSTs are 96.6% splicing-derived, but human evidence is MS + prediction; functional assays only in mice (mixed vaccine, no sequences) → **#681**.
- **Bathini/MHC1-TIP 2026** (Comm Biol) - immunopeptidomics *methods* paper; **no splice peptides and no functional assay** → out of scope (not even a #681 splice feeder).
- **Zhao 2025** (HCC AS-vaccine, *Cell Research*) - **gate-1 + gate-2 PASS**: human per-peptide ELISpot / tetramer / 4-1BB on HCC patient TILs (A\*02:01 / A\*11:01 / A\*24:02), splice (AS) source. **Registry-eligible but sequence-blocked** - exact peptides are referenced only by internal ID (pA02-28, pA11-10, …) in the 19 local supplements; amino-acid strings live in the main paper / a peptide table not in our Zotero set. Folding deferred to a sequence-retrieval pass (do **not** translate junction coords - that would be inference).
- **Later #832 standing-watch find (same sequence-blocked basis, not part of the #733 four):** **TEtrans / Li 2025** (gastric TE-derived transcripts, bioRxiv) - **gate-1 + gate-2 PASS** at transcript level (A\*02:01-restricted CD8+ killing + IFN-gamma; T4180/LTR30 MS-presented in STAD) but the 387 candidate neopeptides are **count-only, no sequence at source** (no supplement table, no data deposit). Documented as a transcript/antigen-level exon-TE source in [`PROVENANCE.md`](PROVENANCE.md); **no `registry.tsv` rows** ([#832](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/832)).

## How it was curated (provenance)

Per-peptide source locations + provenance grades (direct / inferred / agent-local / agent-web / unpublished) are in **`PROVENANCE.md`** and the `provenance_grade` / `provenance_ref` columns of `registry.tsv` - so anyone scoring against this set knows which sequences are rock-solid vs need a manual confirm.

1. Adversarially-verified novelty check (25/25 claims, 0 killed) confirmed no prior open splice-immunogenicity benchmark exists.
2. Multi-agent literature + database sweep → dual-gate verification.
3. Local Zotero full-text PDF mining (paywall-free) - surfaced sources the literature sweep missed (Kim leukemia, Zhao PDAC, Ji, Fisher) and **doubled** the set; added intron-retention biology.
4. Human + machine read of SNAF + Kwok supplementary files to recover exact sequences (HLA-C peptides resolved; GNAS/RPL22 confirmed unpublished).

## Sources

SNAF: Li et al. 2024, *Sci Transl Med* (eade2886). IRIS: Pan et al., *PNAS*. RCAN1-4: Xiong et al. 2025 (GBM). Kim et al. 2025 (mis-splicing neoantigens in SF-mutant leukemias). Fisher et al. 2026 (CoREST). Kwok et al. 2024/2025, *Nature* (s41586-024-08552-0). POSTN-203 study. Manoharan et al. 2026, *Sci Rep* (s41598-026-43687-2; IR-CRC), Zotero `ZAT8678F`. Merlotti et al. 2023, *Sci Immunol* (abm6359; NSCLC exon-TE junctions), Zotero `5ZADNVCB`. **Bigot et al. 2021, *Cancer Discov* 11(8):1938-51 (10.1158/2159-8290.CD-20-0555; SF3B1-mutant uveal melanoma), PMID 33811047; Correction CD-22-0009 / PMID 35257149 (unread); Zotero `DJAS2BJ2` (+ supp tables S1–S9 local). Recovered via the #734 IEDB IQ-API free-text mine (ref 1040678).**

## Next steps

- ~~Board-commit #680 → branch → land this folder.~~ **Done.** ~~#733 data-quality follow-ups (canonicalize `splice_mechanism`; HLA-resolution flag; FMDDYIFV length; split IRIS `candidate`).~~ **Done (#733 AC#4).** ~~IR-CRC + Merlotti folding.~~ **Done (#733 AC#1/#2).**
- **Zhao 2025 sequence-retrieval** ([Issue #817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)) - fetch the main paper PDF / peptide-sequence table to recover the amino-acid strings for the human per-peptide functional hits (pA02-/pA11- IDs), then fold as `functional-scorable` (A\*02:01/A\*11:01/A\*24:02). Sequence-blocked, not gate-blocked.
- **Merlotti single-replicate candidates** ([Issue #818](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/818)) - ~10 peptides recorded in PROVENANCE; re-assess if a stricter or looser inclusion bar is adopted.
- ~~**Data-quality pass** ([Issue #823](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/823)).~~ **Done** - the machine-readable `assay_context` column (see above) is added + populated for all rows, and the 11 Manoharan IR-CRC rows were upgraded `agent-web` → `direct` after a third independent first-hand read of PMC13096189 Table 1 confirmed all 11 sequences verbatim (11/11). Pinning the 34 `unspecified` rows to specific contexts (a read-the-methods pass for SNAF / Kim / Xiong / Fisher; the #838 Kim + UM rows landed `unspecified` pending a T-cell-source read) remains future work.
- [Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681) immunopeptidome mining → MS-presented tier (+ the #733 Tier-2 routes, now routed via a comment there) + decoy negative construction.
- ~~CEDAR/IEDB free-text mine.~~ **Done (#734)** - see [`issue_734_db_audit/`](issue_734_db_audit/) for the 4-DB no-category audit + the recovered-candidate list. Bigot 2021 folded here; the remaining ~20 candidates (long-read UM, POSTN 2nd peptide, Kim extras, Kwok GNAS/RPL22 sequence recovery) + the **A\*02:01 rebalance** are in #734 follow-ups.
- ~~**Sparsity quantification + writeup + experiment deck** ([Issue #737](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/737), #680 AC#4).~~ **Done** - [`issue_737_sparsity/sparsity_writeup.md`](issue_737_sparsity/sparsity_writeup.md) + [`issue_737_sparsity/notebook.ipynb`](issue_737_sparsity/notebook.ipynb) + [`issue_737_sparsity/slides.qmd`](issue_737_sparsity/slides.qmd); see "Sparsity analysis (#737)" above.
- Bigot **Correction (PMID 35257149)** retrieval + review - confirm no peptide/gene change before publication.
- Pre-publication: systematic completeness sweep.
