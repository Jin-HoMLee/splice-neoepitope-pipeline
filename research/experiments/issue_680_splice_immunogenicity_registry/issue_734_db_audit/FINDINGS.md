# 4-DB splice-category audit + IEDB free-text recovery (#734)

**Parent:** [#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) · **Leaf:** [#734](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/734) (carries the original #680 AC#1) · **Date:** 2026-06-22

**Headline:** None of the four major immunogenicity databases (CEDAR, IEDB, NEPdb, dbPepNeo2.0) exposes a splice category — splice-derived neoantigens are **not separable by schema**. *But* free-text/reference-title mining of IEDB's IQ-API **recovers real, functionally-validated splice-neoantigens** that a manual literature sweep missed. **Category-absence ≠ data-absence.**

---

## 1. The null: no DB has a splice category (AC#1)

Each database classifies neoantigens along axes that have **no place for a generation mechanism** (splice junction / intron retention / aberrant 3′SS / neojunction). Confirmed per-DB with the exact evidence below.

| DB | How it classifies | Splice category? | Evidence |
|---|---|---|---|
| **IEDB** | `structure_type` (linear/discontinuous), source antigen (protein/organism), `e_modification` (PTM), disease/host/MHC | **No** | IQ-API (`query-api.iedb.org`) free-text: every splice-*mechanism* keyword in `antigen_description` returns **0** rows — `intron retention`, `neojunction`, `aberrant splic`, `cryptic exon`, `exon skipping`, `retained intron`, `frameshift`. No schema field encodes generation mechanism. |
| **CEDAR** | 3 mutually-exclusive categories: **Neoantigen / Viral / Germline-Self-Host** | **No** | "Using CEDAR" (PMC12543298): mutation-specific search is explicitly *future work* ("we are working on including the option to search for neoantigens encoded by a specific mutation"); **no API** ("it is not possible to … We plan to provide an API"). `query-api.cedar.iedb.org` does not resolve. |
| **NEPdb** | Response / tumor type / HLA / gene / year; scope = **"non-synonymous mutations"** | **No** | Front. Immunol. 2021 (10.3389/fimmu.2021.644637) + live `nep.whu.edu.cn`: no mutation-type facet at all; the "non-synonymous mutations" scope is a point-mutation framing that structurally excludes splice. 173 immunogenic positives. |
| **dbPepNeo2.0** | confidence tiers (LC/MC/HC); source types **SNV / INDEL / fusion / non-coding** | **No** | Front. Immunol. 2022 (10.3389/fimmu.2022.855976): no occurrence of splice/splicing/intron-retention anywhere in tiers, source-type list, or search fields. Closest non-canonical sources = fusion + non-coding. 746 HC class-I. (Live site down at audit; paper-confirmed.) |

**The sharpest demonstration (IEDB, row-level):** splice-neoantigens that *are* deposited are filed as ordinary protein peptides. IRIS-CLASP1 `SLDGTTTKA` → "Cytoplasmic linker associated protein 1"; RCAN1-4 `IFSESETRAKF` → "calcipressin-1 **isoform c** / Isoform 2 of Calcipressin-1". The splice origin survives **only in the free-text reference title** ("…alternative splicing of RCAN1…", "IRIS: … pre-mRNA alternative splicing"), never in a structured field. So even when the data is present, you cannot retrieve splice-neoantigens *as a class*.

## 2. Free-text recovery works (AC#2)

Despite the no-category finding, mining IEDB's IQ-API by **reference title** (`reference_search.reference_title ILIKE *splic*` → 35 references) and by free-text antigen fields surfaces the splice-neoantigen literature. Triaging the 35 splice-titled references:

- **Already in our registry:** RCAN1/Xiong, Kim leukemia, POSTN, IRIS.
- **Proteasome cis/trans-spliced peptides** (≈12 refs: MHC-I-spliced immunopeptidome, KRAS-G12V spliced epitope, etc.) — a *different* phenomenon (post-translational, not RNA-splicing); correctly **excluded** (cf. the splice-junction-terminology collision).
- **Non-cancer** (T1D insulin, ALS-FTD TDP-43, lupus spliceosome, HBV, HMSD minor-H antigen, CD20/desmoglein isoforms) — out of scope.
- **Net-new cancer splice-neoantigens with functional T-cell assays** (the recovery): see `recovered_candidates.tsv`.

**55 net-new functionally-validated candidate peptides** were recovered (gate-2 confirmed by IEDB assay records — tetramer, IFNγ, cytotoxicity, granzyme B, degranulation):

| Source | PMID | Net-new positives |
|---|---|---|
| Bigot 2021 — SF3B1-mutant uveal melanoma (*Cancer Discov*) | 33811047 | 35 |
| Kim mis-splicing leukemias — additional positives beyond our 5 | — | 13 |
| Long-read alt-splicing, uveal melanoma (*Cancer Immunol Res* 2023) | 37756564 | 4 |
| Kwok/Nature — GNAS/RPL22 (registry rows currently `functional-nonscorable`/sequence-unpublished) | 39972144 | 2 |
| POSTN splicing-junction (*Genes Immun* 2025) — 2nd peptide | 40181162 | 1 |

## 3. Dedup + fold (AC#3)

Recovered peptides were deduped against `registry.tsv` (exact sequence). **Bigot 2021 folded this pass** (35 rows: 5 `high` + 30 `medium`) — dual-gate verified (gate-1 = supp. Table S8 alternative-mRNA-junctions + panel design + supp. methods; gate-2 = IEDB effector/tetramer; sequences `direct` from local supp. Table S2, cross-confirmed 35/35 vs IEDB). See `PROVENANCE.md` → "Bigot et al. 2021". The remaining ~20 candidates are deferred to #734 follow-ups (each needs its own primary-source gate-1 pass).

## 4. Why this is a citable result

1. **A documented null with method evidence**, not an assertion: the four canonical immunogenicity DBs cannot surface splice-neoantigens as a class — a real obstacle for anyone trying to assemble a splice-immunogenicity benchmark, and the empirical justification for why this registry had to be hand-built from primary literature.
2. **A method correction:** the obstacle is *schema*, not *data*. Reference-title / free-text mining recovers the data the categories hide — and caught Bigot 2021, a functionally-validated splice-neoantigen paper the Zotero-centric manual sweep had missed. Any future splice-immunogenicity curation should run this DB-text mine, not rely on category browse.

## Reproduce

- IEDB IQ-API (PostgREST): `https://query-api.iedb.org/` — tables `tcell_search`, `epitope_search`, `reference_search`, `antigen_search`. No auth.
- Mechanism-keyword null: `tcell_search?antigen_description=ilike.*<kw>*` with `Prefer: count=exact` (read `Content-Range`).
- Reference-title recall: `reference_search?reference_title=ilike.*splic*`; then `tcell_search?reference_id=eq.<id>` per reference.
- Candidate extraction + dedup: see `recovered_candidates.tsv` (peptide, gene, hla, measure, assay, in_registry, source_ref, pmid).

---

## 5. #838 follow-up dispositions (per-source gate-1, 2026-06-30)

Folding the remaining ~20 recovered candidates ([#838](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/838)), each with its own primary-source gate-1 pass. Dispositions diverge from the §2 table framing where the primary source warranted it.

### POSTN — `TVYTTKIITK` is a **correction, not an addition** (§2 framed it as a "2nd peptide")

Primary-source verification (Liu et al., *Genes Immun* 2025, PMID 40181162, Zotero `HMXA22SW`; supp MOESM1–5 + IEDB) overturns the §2 "2nd peptide" framing:

- The supp **functional** assays (S5/S6: IFN-γ ELISA + flow, caspase-3, LDH cytotoxicity) validate a **single** epitope, "POSTN-203_A11".
- Supp prediction tables (S2/S3, MOESM3/4) tile **both** `KTEGPTLTK` (pos 35–43) and `TVYTTKIITK` (pos 8–17) as *predicted binders* from the **same** POSTN-203 (ENST00000379747) neojunction region — they are not two independent epitopes.
- **IEDB is decisive:** `tcell_search?linear_sequence=eq.TVYTTKIITK` → **4 Positive A11 rows** (ref1046253); `KTEGPTLTK` → **0 rows** (never assayed).
- The pre-existing registry row keyed the validated epitope as `KTEGPTLTK` with `provenance_grade=agent-web` and note `…VERIFY` — an **agent-web mis-key** that recorded the top *predicted* A*11:01 binder instead of the *assayed* peptide.

**Action:** corrected the existing functional-scorable row in place (not a new row): `KTEGPTLTK`→`TVYTTKIITK`, `length` 9→10, `provenance_grade` agent-web→`direct`, refreshed provenance (PMID/DOI/Zotero/IEDB). **HLA resolution downgraded 4-digit→2-digit**: the paper types the line only as `A*11+` (SF10281) and IEDB curates `A11`; the 4-digit `A*11:01` appears solely in the binding-prediction tables. `hla` kept at `A*11:01` (predominant A11 allele, needed for scoring) with the caveat noted in the row. **Net registry count unchanged** (correction, not addition). Stakes: this is a `functional-scorable` ground-truth label feeding the #736 benchmark — the wrong sequence would have poisoned it.

### Kwok GNAS + RPL22 — `functional-nonscorable` → `functional-scorable` (sequences recovered + confirmed)

Both rows existed as `NOT-PUBLISHED` / `functional-nonscorable` ("exact 9-mer not published, figures only"). The Kwok IEDB deposit (ref1045680, PMID 39972144) carries the sequences; primary-source verification confirms both gene/junction assignments:

- **GNAS** ← `SLLLPSFHL` (9-mer, A*02:01). IEDB `parent_source_antigen` = "Guanine nucleotide-binding protein G(s) subunit alpha" → **GNAS-confirmed**; gate-1 = NEJ^GNAS frameshift (Kwok). Sequence not legible in the paper figures → recovered from the authors' IEDB deposit. `provenance_grade` unpublished→`direct`.
- **RPL22** ← `GIMDAANFFL` (**10-mer**, A*02:01). IEDB does **not** structurally annotate its gene (`parent_source_antigen=None`) — so the gene was *not* assignable from the mine alone. Confirmed independently: Kwok names exactly two public NEJ neoantigens (NeoA^GNAS + NeoA^RPL22); GNAS is taken by `SLLLPSFHL`, and `GIMDAANFFL` maps cleanly onto **RPL22 (UniProt P35268)** WT `…G·IMDAANFE·QFLQER…` with the paper's **in-frame 6-nt (EQ) loss** → `GIMDAANF·FLQER` → A*02:01 `GIMDAANFFL`. This **supersedes the prior illustrative `~LALDVLQGYSL`** placeholder. `provenance_grade` unpublished→`direct`.

### Kim mis-splicing leukemias — 13 net-new → `functional-scorable` (the "screening hits" framing was wrong)

§2 / the #838 body cautioned the ~13 remaining Kim candidates "may be screening hits the original curation deliberately excluded." Reading the **Kim supplementary tables** (Zotero item `XB3CPX5P`, not the PDF-attachment key `LX6DMXTL` the registry cites — same dual-key trap as Kwok) overturns that:

- **mmc2 (Supp Table S2)** gives every candidate a specific AS event + parent gene + **genomic junction coordinates** (`se`=skipped exon, `ci`=constitutive intron, `mxe`=mutually exclusive exons, `a5ss`=alt 5′SS per the table legend). → gate-1 confirmed per peptide; `junction_mapping_grade=coords` (the **top** grade — *better*-provenanced than the original 5 rows, which sit at `gene-mechanism`).
- **mmc5 (Supp Table S5)** shows each has an **A*02:01 dextramer** panel; IEDB records IFN-γ ELISPOT positivity. Per §2.3 **ELISPOT is a `strong` effector readout**, so these are legitimate scorable positives, not weak screening hits.
- The original curation folded only 5 because only those 5 were legibly named in the *figures*; the other 13 lived in the supp tables, which weren't read at the time.

**Folded 13 as `functional-scorable`** (A*02:01, `evidence_strength=strong`, `provenance=direct`, `junction_mapping_grade=coords`): USF1, AP4B1, EZH2, BIN3, CASP2, SLC20A1 (6× IEDB-`Positive`, `confidence=high`); RHOT2 ×3, NCOA7, LTBR, MAN2C1, SH3GL1 (7× IEDB-`Positive-Low`, `confidence=medium`).

### Long-read UM — 4 → `functional-nonscorable` (HLA-unresolved, *not* "A*02:01")

The #838 body called these "~4 A*02:01 positives." The supp (Table 13, item `6P5JCCIB`) confirms gate-1 — `MADAGAMAA`/`RQVGEGCRT`→SEPTIN6, `KNILNGSRA`→AMZ2P1 (archaemetzincin-2 *pseudogene*), `RGAGLGRAL`→MZT2B (all alt-splicing neojunctions) — **but publishes no per-peptide HLA** (IEDB = `human`; ELISPOT on bulk PBMC). Without an allele key these cannot be scored → folded as `functional-nonscorable` (retained for the validation base, excluded from the scored set). `RGAGLGRAL`/`RQVGEGCRT` are IEDB mixed Negative/Positive(-Low).

### Net effect of #838 on the registry

| Change | Count |
|---|---|
| POSTN sequence-corrected (`KTEGPTLTK`→`TVYTTKIITK`) | 1 row fixed |
| GNAS + RPL22 promoted non-scorable → scorable | 2 rows |
| Kim net-new folded scorable | +13 rows |
| UM net-new folded non-scorable | +4 rows |

`functional-scorable`: 64 → **79**. Total registry rows: 79 → **96**. Validator: `PASS: 96 rows valid`.
