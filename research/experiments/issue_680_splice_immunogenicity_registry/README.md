# Open splice-neoepitope → measured T-cell immunogenicity registry (DRAFT)

**Issue:** [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) · **Status:** draft for review (not yet board-committed or branched) · **Snapshot:** 2026-06-11

## What this is

The first openly-assembled, provenance-tracked registry of **splice-junction-derived (aberrant-splicing / neojunction / intron-retention) cancer neoantigens with experimentally MEASURED T-cell immunogenicity**. Each entry passed a two-gate test, each gate adversarially verified (3 independent skeptic votes, ≥2/3 required): (1) genuinely splice-derived — not SNV / indel / fusion / canonical-TAA / proteasomal-cis-spliced; (2) an actual functional T-cell assay was performed on that exact peptide (IFN-γ, tetramer/dextramer, cytotoxicity, ELISpot, or engineered-TCR readout).

## What this is NOT — a probe, not a powered benchmark

A feasibility gate established that only ~tens of such peptides exist field-wide (TESLA explicitly excluded splice isoforms; the four major immunogenicity DBs — CEDAR/IEDB/NEPdb/dbPepNeo2.0 — have no splice category). **This is a small, transparently-underpowered stress-test probe**, not a TESLA-scale calibrated benchmark. Its value is being first to assemble an open, sequence-resolved splice-immunogenicity set with documented provenance — and the honest demonstration that the field's functional-validation base is critically thin. Source clustering is real: the entries come from ~6 studies.

## Tiers (`registry.tsv`)

| Tier | n | Meaning |
|---|---|---|
| `functional-scorable` | 18 | Splice-derived + functional T-cell assay + exact sequence & allele → directly usable to score a predictor |
| `functional-nonscorable` | 2 | Functionally validated (GNAS, RPL22) but exact 9-mer not published anywhere accessible → cannot be keyed |
| `presentation-prevalence` | 1 | `TEFQTRRAM`/SLC45A2 — a real splice neoantigen with prevalence/survival evidence + predicted binding, but **no functional T-cell assay performed** (distinct from the 5 SNAF IFN-γ-tested peptides) |
| `hard-negative-true-splice` | 1 | `VELEDHVML` — splice-derived, MS-presented, experimentally NON-immunogenic. The gold-standard discrimination negative. |
| `negative-control-not-splice` | 1 | `VFVDGLCRAKF` — RCAN1-1 constitutive same-locus control (fails splice gate) |

## Coverage

- **Alleles:** HLA-A*02:01, A*24:02, A*11:01, **C*04, C*08** (HLA-C breadth recovered from SNAF supplements — breaks the otherwise-severe A*02:01 skew).
- **Splice mechanisms:** exon-junction-derived, intron-retention, minor-intron, poison-exon, alt-5′/alt-3′.
- **Genes (~16):** SLC45A2, CDH19, PMEL, CLASP1, SCAMP3, MAN2C1, AP2A1, CLK3, RHOT2, c16orf70, ATG3, MYO1F, RCAN1, POSTN, GNAS, RPL22 (+ Fisher's two unnamed-gene peptides).

## Caveats

- **Sequence not published (GNAS, RPL22):** exact 9-mers are figure-locked (MS-spectra peak plots) or only in a generic method schematic; supplementary tables hold no peptide list. Both A*02:01 (no allele-diversity loss). Recoverable only via author contact / raw MS.
- **Medium-confidence rows:** MAN2C1 / AP2A1 (IRIS) failed a second verifier pass; Fisher's RSQGWLFLR / AEHAHRVPL / VELEDHVML lack named genes and have ambiguous restriction (loaded A*11:01 but predicted B*40:01).
- **Readout heterogeneity:** several Kim entries are tetramer-only (antigen-specific T-cell detection) vs IFN-γ/cytotoxicity for others — recorded per-row.
- **Allele granularity:** SNAF HLA-C entries are 2-digit (C*04 / C*08); others are 4-digit.
- **Negatives are scarce:** only one true splice hard negative (`VELEDHVML`). A usable discrimination probe needs a constructed decoy negative set — see deferred work.

## Deferred — MS-presented / immunogenicity-untested tier (→ Issue #681)

13 distinct splice peptides are MHC-presented (immunopeptidomics) but have **no functional T-cell assay** — they fail the functional gate and are NOT in this registry. They belong with #681 (public immunopeptidome mining) and are a natural feeder for the decoy negative set. Source: SNAF Supp Fig 4 (`NQDEDPLEV`/C6orf52, `KGPWYPLSL`/C20orf204, `VAPGEAKNL`/RASA3, `YALANIKWI`/DYNLT5, `KEKLDQLVY`/FBXO7, `TELQRTLSL`/NGLY1) + Supp Fig 7 (`SQTPKSRAL`/PSMF1, `RKLEAPYLL`/MCF2L, `LSWPRSTPM`/CPN1, `VSTGCAVVL`/SERPINE2, `RRLPNPPAV`/RGS12, `IVKRPRSEL`/EXO1, `AVPLLQTNR`/ETV4).

## How it was curated (provenance)

Per-peptide source locations + provenance grades (direct / inferred / agent-local / agent-web / unpublished) are in **`PROVENANCE.md`** and the `provenance_grade` / `provenance_ref` columns of `registry.tsv` — so anyone scoring against this set knows which sequences are rock-solid vs need a manual confirm.

1. Adversarially-verified novelty check (25/25 claims, 0 killed) confirmed no prior open splice-immunogenicity benchmark exists.
2. Multi-agent literature + database sweep → dual-gate verification.
3. Local Zotero full-text PDF mining (paywall-free) — surfaced sources the literature sweep missed (Kim leukemia, Zhao PDAC, Ji, Fisher) and **doubled** the set; added intron-retention biology.
4. Human + machine read of SNAF + Kwok supplementary files to recover exact sequences (HLA-C peptides resolved; GNAS/RPL22 confirmed unpublished).

## Sources

SNAF: Li et al. 2024, *Sci Transl Med* (eade2886). IRIS: Pan et al., *PNAS*. RCAN1-4: Xiong et al. 2025 (GBM). Kim et al. 2025 (mis-splicing neoantigens in SF-mutant leukemias). Fisher et al. 2026 (CoREST). Kwok et al. 2024/2025, *Nature* (s41586-024-08552-0). POSTN-203 study.

## Next steps

- Board-commit #680 (Backlog→Ready, PM-coordinated) → branch → land this folder.
- #681 immunopeptidome mining → MS-presented tier + decoy negative construction.
- Pre-publication: systematic completeness sweep + the CEDAR/IEDB free-text mine (deferred).
