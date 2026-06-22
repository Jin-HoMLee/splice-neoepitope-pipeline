# Open splice-neoepitope → measured T-cell immunogenicity registry

**Issue:** [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) (seed); grown by [#733](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/733) (standing-watch + library-sweep additions) · **Status:** seed (23 rows, 2026-06-11) + #733 pass (→ 44 rows): IR-CRC (Manoharan) + exon-TE (Merlotti) folded; Tier-2 triaged · **Snapshot:** 2026-06-21

## What this is

The first openly-assembled, provenance-tracked registry of **splice-junction-derived (aberrant-splicing / neojunction / intron-retention) cancer neoantigens with experimentally MEASURED T-cell immunogenicity**. Each entry passed a two-gate test, each gate adversarially verified (3 independent skeptic votes, ≥2/3 required): (1) genuinely splice-derived — not SNV / indel / fusion / canonical-TAA / proteasomal-cis-spliced; (2) an actual functional T-cell assay was performed on that exact peptide (IFN-γ, tetramer/dextramer, cytotoxicity, ELISpot, or engineered-TCR readout).

## What this is NOT — a probe, not a powered benchmark

A feasibility gate established that only ~tens of such peptides exist field-wide (TESLA explicitly excluded splice isoforms; the four major immunogenicity DBs — CEDAR/IEDB/NEPdb/dbPepNeo2.0 — have no splice category). **This is a small, transparently-underpowered stress-test probe**, not a TESLA-scale calibrated benchmark. Its value is being first to assemble an open, sequence-resolved splice-immunogenicity set with documented provenance — and the honest demonstration that the field's functional-validation base is critically thin. Source clustering is real: the entries come from ~6 studies.

## Tiers (`registry.tsv`)

| Tier | n | Meaning |
|---|---|---|
| `functional-scorable` | 29 | Splice-derived + functional T-cell assay + exact sequence & allele → directly usable to score a predictor. Confidence sub-graded `high`/`medium` (e.g. healthy-donor IVS or tetramer-detection-only → `medium`). |
| `candidate` | 2 | `LLLGIAKLLKV`/MAN2C1, `FLSELEPPA`/AP2A1 (IRIS) — positives whose 2nd adversarial-verify pass did not confirm → held below scorable until re-verified (#733 AC#4). |
| `candidate-negative` | 8 | IR-CRC (Manoharan) IR-derived 9-mers **tested-negative** in healthy-donor IVS → *soft* negatives (failed to prime ≠ intrinsically non-immunogenic). A decoy-negative seed, categorically weaker than the hard negative. |
| `functional-nonscorable` | 2 | Functionally validated (GNAS, RPL22) but exact 9-mer not published anywhere accessible → cannot be keyed |
| `presentation-prevalence` | 1 | `TEFQTRRAM`/SLC45A2 — a real splice neoantigen with prevalence/survival evidence + predicted binding, but **no functional T-cell assay performed** (distinct from the 5 SNAF IFN-γ-tested peptides) |
| `hard-negative-true-splice` | 1 | `VELEDHVML` — splice-derived, MS-presented, experimentally NON-immunogenic. The gold-standard discrimination negative. |
| `negative-control-not-splice` | 1 | `VFVDGLCRAKF` — RCAN1-1 constitutive same-locus control (fails splice gate) |

**Two machine-readable columns added in #733 (AC#4):** `splice_mechanism_canonical` (controlled vocabulary — see below) and `hla_resolution` (`2-digit` for SNAF C\*04/C\*08 + the SNAF A\*02 prevalence row, `4-digit` otherwise).

### `splice_mechanism_canonical` controlled vocabulary

Free-text `splice_mechanism` is preserved; `splice_mechanism_canonical` normalizes it for grouping:

`alt_5p_ss` (alt 5′ splice site / donor) · `alt_3p_ss` (alt 3′ splice site / acceptor) · `alt_splice_junction` (alternative-splicing junction, sub-event unspecified — IRIS/Kim/Xiong/POSTN) · `intron_retention` · `minor_intron` (U12-type) · `poison_exon` (NMD-targeting cassette) · `neojunction_frameshift` · `neojunction_inframe` · **`exon_te_junction`** (exon–transposable-element junction or TE-induced frameshift — **registered in #733** for the Merlotti set) · `not_splice` (constitutive control) · `ambiguous` (splice-derived, mechanism unresolvable from source — Fisher unnamed-gene peptides).

## Coverage (44 rows, 9 sources, 39 genes)

- **Alleles (11 distinct):** HLA-A\*02:01 (24), A\*11:01 (9), A\*24:02, **C\*04, C\*08** (HLA-C breadth from SNAF supplements), + A\*34:01, A\*31:01, A\*24:07, A\*02:07, A\*02:06 (allele *diversity* added by IR-CRC). **A\*02:01 skew deepened** to 24/44 — the IR-CRC positives + the entire Merlotti exon-TE set are A\*02:01/A2-predicted (see Caveats).
- **Splice mechanisms:** + `exon_te_junction` (Merlotti, 10) and `intron_retention` (IR-CRC, +8) on top of the seed's alt-5′/alt-3′/minor-intron/poison-exon/neojunction set.
- **Sources (9):** SNAF (Li 2024), IRIS (Pan), Kim 2025, Xiong 2025 (RCAN1-4), Fisher 2026, Kwok 2024, POSTN-203, **Manoharan 2026 (IR-CRC, +11)**, **Merlotti 2023 (NSCLC exon-TE, +10)**.
- **Genes (39):** seed ~16 + IR-CRC (CHD7, FERMT3, TPCN2, LRRC8B, OPLAH, FANCB, PLEKHG5, WDR19, RNF123, ITGAM, CDH24) + Merlotti (SRSF7, FRG1B, SCGB3A2, SBNO2, DCBLD2, KYNU, TENM1, ATAD3B, SLC39A11, PDE4D).

## Caveats

- **Sequence not published (GNAS, RPL22):** exact 9-mers are figure-locked (MS-spectra peak plots) or only in a generic method schematic; supplementary tables hold no peptide list. Both A*02:01 (no allele-diversity loss). Recoverable only via author contact / raw MS.
- **Medium-confidence rows:** MAN2C1 / AP2A1 (IRIS) failed a second verifier pass; Fisher's RSQGWLFLR / AEHAHRVPL / VELEDHVML lack named genes and have ambiguous restriction (loaded A*11:01 but predicted B*40:01).
- **Readout heterogeneity:** several Kim entries are tetramer-only (antigen-specific T-cell detection) vs IFN-γ/cytotoxicity for others — recorded per-row.
- **Allele granularity:** SNAF HLA-C entries are 2-digit (C*04 / C*08); others are 4-digit.
- **Negatives:** one true splice *hard* negative (`VELEDHVML`, MS-presented + ELISpot-negative) + **8 IR-CRC *soft* negatives** (`candidate-negative`) added by #733 — failed-to-prime in healthy-donor IVS, categorically weaker than the hard negative but a usable decoy-negative seed for #681.
- **IR-CRC (Manoharan) provenance = healthy-donor IVS, not patient:** the 3 positives (CHD7/FERMT3/TPCN2) and 8 negatives were measured in **healthy-donor moDC-primed T cells** (in-vitro sensitization), *not* patient T cells → positives graded `medium` (a peptide that *can* prime ≠ primed in a patient). FERMT3/RQDPAPQQV was functionally restricted on **A\*02:01** only (A\*24:07 predicted, not used).
- **Merlotti (exon-TE) confidence split:** 3 with cloned-TCR functional validation (`high`); 7 with **patient tetramer-DETECTION only** (`medium` — no per-peptide IFN-γ/killing). ~10 further single-replicate patient tetramer⁺ hits were **not** folded as scorable rows (recorded in PROVENANCE as candidates).
- **A\*02:01 skew (worsened by #733):** 24/44 rows are A\*02:01 — both the IR-CRC functional positives and the entire Merlotti exon-TE set are A2. The HLA-C/allele-diversity story is now a smaller fraction of the whole; weight scoring analyses accordingly.

## Deferred — MS-presented / immunogenicity-untested tier (→ Issue #681)

13 distinct splice peptides are MHC-presented (immunopeptidomics) but have **no functional T-cell assay** — they fail the functional gate and are NOT in this registry. They belong with #681 (public immunopeptidome mining) and are a natural feeder for the decoy negative set. Source: SNAF Supp Fig 4 (`NQDEDPLEV`/C6orf52, `KGPWYPLSL`/C20orf204, `VAPGEAKNL`/RASA3, `YALANIKWI`/DYNLT5, `KEKLDQLVY`/FBXO7, `TELQRTLSL`/NGLY1) + Supp Fig 7 (`SQTPKSRAL`/PSMF1, `RKLEAPYLL`/MCF2L, `LSWPRSTPM`/CPN1, `VSTGCAVVL`/SERPINE2, `RRLPNPPAV`/RGS12, `IVKRPRSEL`/EXO1, `AVPLLQTNR`/ETV4).

**#733 Tier-2 triage → #681 (functional gate confirmed absent).** Four library-sweep papers were read for a per-peptide functional T-cell assay; all routed away from the functional registry:
- **Courcelles 2026** (CRC MSI/MSS, MCP) — splice-derived aeTSAs (intron retention, retroelement junctions) but immunogenicity is *in-silico only* (PRIME/ImmuneApp) → **#681**.
- **Lin 2025** (HCC, bioRxiv) — neoTSTs are 96.6% splicing-derived, but human evidence is MS + prediction; functional assays only in mice (mixed vaccine, no sequences) → **#681**.
- **Bathini/MHC1-TIP 2026** (Comm Biol) — immunopeptidomics *methods* paper; **no splice peptides and no functional assay** → out of scope (not even a #681 splice feeder).
- **Zhao 2025** (HCC AS-vaccine, *Cell Research*) — **gate-1 + gate-2 PASS**: human per-peptide ELISpot / tetramer / 4-1BB on HCC patient TILs (A\*02:01 / A\*11:01 / A\*24:02), splice (AS) source. **Registry-eligible but sequence-blocked** — exact peptides are referenced only by internal ID (pA02-28, pA11-10, …) in the 19 local supplements; amino-acid strings live in the main paper / a peptide table not in our Zotero set. Folding deferred to a sequence-retrieval pass (do **not** translate junction coords — that would be inference).

## How it was curated (provenance)

Per-peptide source locations + provenance grades (direct / inferred / agent-local / agent-web / unpublished) are in **`PROVENANCE.md`** and the `provenance_grade` / `provenance_ref` columns of `registry.tsv` — so anyone scoring against this set knows which sequences are rock-solid vs need a manual confirm.

1. Adversarially-verified novelty check (25/25 claims, 0 killed) confirmed no prior open splice-immunogenicity benchmark exists.
2. Multi-agent literature + database sweep → dual-gate verification.
3. Local Zotero full-text PDF mining (paywall-free) — surfaced sources the literature sweep missed (Kim leukemia, Zhao PDAC, Ji, Fisher) and **doubled** the set; added intron-retention biology.
4. Human + machine read of SNAF + Kwok supplementary files to recover exact sequences (HLA-C peptides resolved; GNAS/RPL22 confirmed unpublished).

## Sources

SNAF: Li et al. 2024, *Sci Transl Med* (eade2886). IRIS: Pan et al., *PNAS*. RCAN1-4: Xiong et al. 2025 (GBM). Kim et al. 2025 (mis-splicing neoantigens in SF-mutant leukemias). Fisher et al. 2026 (CoREST). Kwok et al. 2024/2025, *Nature* (s41586-024-08552-0). POSTN-203 study. **Manoharan et al. 2026, *Sci Rep* (s41598-026-43687-2; IR-CRC), Zotero `ZAT8678F`. Merlotti et al. 2023, *Sci Immunol* (abm6359; NSCLC exon-TE junctions), Zotero `5ZADNVCB`.**

## Next steps

- ~~Board-commit #680 → branch → land this folder.~~ **Done.** ~~#733 data-quality follow-ups (canonicalize `splice_mechanism`; HLA-resolution flag; FMDDYIFV length; split IRIS `candidate`).~~ **Done (#733 AC#4).** ~~IR-CRC + Merlotti folding.~~ **Done (#733 AC#1/#2).**
- **Zhao 2025 sequence-retrieval** ([Issue #817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)) — fetch the main paper PDF / peptide-sequence table to recover the amino-acid strings for the human per-peptide functional hits (pA02-/pA11- IDs), then fold as `functional-scorable` (A\*02:01/A\*11:01/A\*24:02). Sequence-blocked, not gate-blocked.
- **Merlotti single-replicate candidates** ([Issue #818](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/818)) — ~10 peptides recorded in PROVENANCE; re-assess if a stricter or looser inclusion bar is adopted.
- [Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681) immunopeptidome mining → MS-presented tier (+ the #733 Tier-2 routes, now routed via a comment there) + decoy negative construction.
- Pre-publication: systematic completeness sweep + the CEDAR/IEDB free-text mine (deferred).
