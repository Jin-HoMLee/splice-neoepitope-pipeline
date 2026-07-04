# Per-peptide provenance — splice-immunogenicity registry (#680)

Each sequence is graded by **how it was obtained**, so it can be re-verified manually. Local files are under `~/Zotero/storage/<KEY>/`.

> **⚠️ Zotero key convention — cite the PARENT item key, not the PDF-attachment key.** A paper's supplementary tables hang off the **parent item**, while its main PDF often has its own attachment key. Several rows here cite the *attachment* key (Kwok `EA2NMFNZ` → parent `5ZT8KC8X`; Kim `LX6DMXTL` → parent `XB3CPX5P`), which hides the supplements: a `children` query on the attachment key returns **nothing**, so the supp tables look absent when they're present on the parent. This trap cost real folds twice (Kwok/Kim, recovered in #838) and a stale "NOT in Zotero" note once (IRIS `QTGXGQZM`, [#904](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/904)). **When a source's supp looks missing, resolve the parent item first** (`gh`/curl the attachment key → `.data.parentItem` → query *its* children).

## Provenance grades

| Grade | Meaning |
|---|---|
| `direct` | Exact sequence read first-hand at the cited figure/table/sheet |
| `inferred` | Sequence derived by matching a figure's 3-letter prefix to the unique catalog entry — **verify against the synthesized-peptide list** |
| `agent-local` | Extracted by an analysis agent from the **local Zotero PDF**; location recorded, not every panel personally eyeballed |
| `agent-web` | Extracted by a web-reading agent; source **not in Zotero**; exact figure/table not captured — **verify most carefully** |
| `unpublished` | Exact sequence not published in any legible form |

---

## Venue audit (#1001) - the `venue_type` column

Publication venue class per source, resolved first-hand (2026-07-04) and recorded machine-readably in the `venue_type` column (`journal | preprint | db_recovered`; controlled vocab in [`labeling_constants.py`](labeling_constants.py), derivation in [`derive_venue_type.py`](derive_venue_type.py)). Venue is source-keyed; each source's venue is already named in its per-source section below. **Result: all 14 folded sources (96 rows) are peer-reviewed journal articles - 0 preprints, 0 db-recovered.**

| Source (key) | Venue | `venue_type` |
|---|---|---|
| Bigot 2021 (SF3B1 UM) | *Cancer Discovery* | `journal` |
| Kim 2025 (SF-mutant leukemia) | *Cell* | `journal` |
| Manoharan 2026 (IR-CRC) | *Scientific Reports* | `journal` |
| Merlotti 2023 (NSCLC exon-TE) | *Science Immunology* | `journal` |
| SNAF / Li 2024 | *Science Translational Medicine* | `journal` |
| Long-read UM 2023 (Yao) | *Cancer Immunology Research* | `journal` |
| IRIS (Pan/Xing) | *PNAS* | `journal` |
| Fisher 2026 (CoREST) | *JCI Insight* | `journal` |
| Kwok 2024 | *Nature* | `journal` |
| Xiong 2025 (GBM) | *Cellular & Molecular Immunology* | `journal` |
| POSTN-203 (Genes & Immunity) | *Genes & Immunity* | `journal` |

*(Venue names verified against Zotero `publicationTitle` for collection `Z38GTJNW`, 2026-07-04.)*

**Not a policy against preprints.** A preprint-sourced row can be folded; it just must be **marked** `preprint`. The three preprints encountered to date (Lin 2025 bioRxiv, TEtrans/Li 2025 bioRxiv, the SNAF exon-TE preprint) were each deferred for **sequence-unavailability**, not because they are preprints. The guard against a *silent* preprint fold is structural: a source absent from the venue map derives to the out-of-vocab sentinel `unclassified`, which `validate_registry.py` rejects, so any newly-folded source must be venue-classified before it validates.

---

## SNAF — Li et al. 2024, *Sci Transl Med* 16, eade2886
Files: `A46YTSYY` (main PDF), `7SLW6B96` (sm.pdf), `WK4DHT6M` (Data S1–S15 .xlsx)

- **`RLLGTEFQT`** SLC45A2 / A\*02 — `direct`. Main PDF **p.10, Fig 5A legend**, verbatim: *"loaded with FLU and HCMV control peptides and RLLGTEFQT (RLL) and FQTRRAMTL (FQT) peptide neoantigens"*. IFN-γ response Fig 5C/D (`A*02_RLL`).
- **`FQTRRAMTL`** SLC45A2 / A\*02 — `direct`. Same Fig 5A/B legend (p.10); `A*02_FQT` MHC-stabilization.
- **`IPDSQGNDI`** SLC45A2 / C\*04 — `inferred`. Fig 5C shows only label `C*04_IPD (SLC45A2)` (41% IFN-γ⁺ CD8). Matched `IPD`→the unique SLC45A2 `IPD*` row in **Data S1–TCGA sheet** of `WK4DHT6M`: *"IPDSQGNDI | ENSG00000164175:E3.2_33963931-E4.2 | chr5:33954504-33963931(-) | SLC45A2"*. **Verify**: confirm against the 36-synthesized-peptide list (sm.pdf Methods, image-only).
- **`IIDNQEPVF`** CDH19 / C\*04 — `direct`. sm.pdf **Supp Fig 7B**, verbatim *"IIDNQEPVF"* + *"CDH19:E12.1-E13.2_66509195"*; functional in main Fig 5C `C*04_IID (CDH19)` (28%).
- **`STESITATL`** PMEL / C\*08 — `inferred`. Fig 5C `C*08_STE (PMEL)` (29%). Matched `STE`→Data S1 row *"STESITATL | ENSG00000185664:E11.9-E12.2_55956189 | chr12:55956189-55956949(-) | PMEL"*. Same verify caveat as IPDSQGNDI.
- **`TEFQTRRAM`** SLC45A2 / A\*02 *(presentation/prevalence — NOT functionally tested)* — `direct`. Main text verbatim: *"SLC45A2 (TEFQTRRAM) was detected in 212 of 472 patients from the TCGA SKCM cohort and was correlated with poor overall survival"*; sm.pdf Supp Fig 5B label *"SLC45A2 TEFQTRRAM"*.

## Kim et al. 2025 — mis-splicing neoantigens, SF-mutant leukemias
Parent item `XB3CPX5P` (DOI 10.1016/j.cell.2025.03.047); the registry/figures cite the **PDF-attachment key `LX6DMXTL`** — the supp tables `mmc1–5` hang off the **parent**, the dual-key trap that hid them on the first pass (same as Kwok `EA2NMFNZ`/`5ZT8KC8X`).

The original 5 (figures, `agent-local`):
- **`RLWGTWVKA`** CLK3 / A\*02:01 — **Fig S1B** *"CLK3 #1 (RLWGTWVKA)"*; MS Fig 3G / S4A; Fig S3A.
- **`CLLPPALFL`** RHOT2 / A\*02:01 — **Fig S1B** *"RHOT2 #5 (CLLPPALFL)"*; MS spectrum Fig 2C; Fig 7E.
- **`RLLAAVLEA`** c16orf70 / A\*02:01 — **Fig 2H** (MS + label *"RLLAAVLEA"*); Fig S1B *"C16orf70 (RLLAAVLEA)"*.
- **`LVAKWLLTI`** ATG3 / A\*02:01 — **Fig S2F** *"LVAKWLLTI"*; Fig S1B *"ATG3 (LVAKWLLTI)"*.
- **`FMDDYIFV`** MYO1F / A\*02:01 — **8-mer, confirmed (#733 AC#4).** Fig S2I legend (verbatim, from ft-cache): *"…the MYO1F peptide (FMDDYIFV). (J) FACS plot…"* — the functionally-assayed peptide is the 8-mer; the `RFMDDYIFV` in Fig S1B/S3B was the source-region label, not the tested peptide.

### #838 fold — 13 net-new from the supp tables (`direct`, **coords**-grade)
Gate-1 + gene + **genomic junction coordinates** from **mmc2 (Supp Table S2)** ("Candidate immunogenic peptides … predicted to bind HLA-A\*02:01"); A\*02:01 **dextramer** panel from **mmc5 (Supp Table S5)**; IFN-γ ELISPOT positivity from IEDB **ref1046844**. All A\*02:01, `evidence_strength=strong` (ELISPOT is effector per LABELING_SCHEME §2.3), `junction_mapping_grade=coords`. Event codes per the mmc2 legend: `se`=skipped exon, `ci`=constitutive intron, `mxe`=mutually exclusive exons, `a5ss`=alt 5′SS. This corrects the #734-body assumption that these were "screening hits deliberately excluded" — the original curation simply hadn't read the supp tables (only the figures, which named 5).
- **`high`** (IEDB Positive): `AVIQGQFFV`/USF1 (a5ss), `GLRDKASYV`/AP4B1 (mxe), `LLPRVSVWL`/EZH2 (se), `SLIRAQEEA`/BIN3 (se), `VLCDQTAQV`/CASP2 (ci), `YLSKYLNVNL`/SLC20A1 (ci).
- **`medium`** (IEDB Positive-Low): `FLWPGLGPS` + `FLWPGLGPSV` + `ILGSLTWSC`/RHOT2 (ci), `FTDGKVYSNV`/NCOA7 (ci), `NLQALLYQA`/LTBR (se), `SVAWPWPAV`/MAN2C1 (se), `YLQPNPGDAL`/SH3GL1 (se).

## Xiong et al. 2025 — RCAN1-4, glioblastoma
File: `JJZL6D84` · `agent-local`

- **`IFSESETRAKF`** RCAN1-4 / A\*24:02 — **Fig 3B** peptide table (*"RCAN1-4_22-32, %Rank EL 0.488, Strong"*); Fig 3A structure (*"…NSDIFSESETR–AKFESLFRTY"*, SJ exon4/exon5); functional Fig 3D–K, 5, 6.
- **`VFVDGLCRAKF`** RCAN1-1 control / A\*24:02 *(negative)* — RCAN1-1_77-87 constitutive control; **verify in Fig 3 control panels**.

## Fisher et al. 2026 — CoREST
File: `WLSZLMA4` · `agent-local` · genes not named in source

- **`RSQGWLFLR`**, **`AEHAHRVPL`**, **`VELEDHVML`** *(VELEDHVML is the negative)* — all **Fig 7I** (peptide table), **Fig 7J/K** (ELISpot), **p.15**; context text p.12. Loaded on A\*11:01 APCs (`AEHAHRVPL`/`VELEDHVML` predicted B\*40:01 — restricting allele ambiguous).

## IRIS — Pan et al., *PNAS* 2023 (10.1073/pnas.2221116120; Zotero parent `QTGXGQZM`; supp `pnas.2221116120.sd01-03.xlsx` + appendix `sapp.pdf`)
The supplement **is** in Zotero — it hangs off the **parent** item `QTGXGQZM`, not the cited PDF-attachment key (the parent-key trap; see the callout at top). The prior "⚠️ NOT in Zotero — verify in supplement" note was stale and is corrected here ([#904](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/904)). Functional readout for all four = **IFN-γ ELISA of TCR-transduced PBMCs co-cultured with K562-A2 + exogenous peptide** (appendix Fig S6); TCR clones isolated + screened in Jurkat-NFAT-GFP.

- **`SLDGTTTKA`** CLASP1 / A\*02:01 — TCR JPTCR-238 (IFN-γ ELISA + Incucyte cytotoxicity). `agent-web` (not re-graded in #904; supp now located for a future direct upgrade).
- **`STMYYLWML`** SCAMP3 / A\*02:01 — TCRs JPTCR-45/47/50/56 (IFN-γ). `agent-web` (as above).
- **`LLLGIAKLLKV`** MAN2C1 / A\*02:01 — TCR JPTCR-13. **Promoted `candidate` → `functional-scorable` (#904), `direct`.** Dataset S3a gives the AS event + junction coords `chr15:-:75655550:75655631:75655089:75656828` (gate-1); S3b/S3c show JPTCR-13 **specifically** recognizes the single epitope (gate-2). `junction_mapping_grade` upgraded `gene-mechanism` → `coords`.
- **`FLSELEPPA`** AP2A1 / A\*02:01 — TCR JPTCR-52. **Promoted `candidate` → `functional-scorable` (#904), `direct`.** Junction coords `chr19:+:50305790:50305856:50305398:50306205` (gate-1); S3c shows JPTCR-52 specifically recognizes the single epitope (gate-2). `coords`-grade. (The earlier "2nd verifier did not confirm" hold was a *supplement-access* failure, not negative evidence — the parent-key trap; resolved.)

## POSTN-203 study (Liu et al., *Genes Immun* 2025; PMID 40181162; Zotero `HMXA22SW`) — `direct` (corrected #838)
- **`TVYTTKIITK`** POSTN-203 (ENST00000379747) / A\*11:01 — the functionally validated **POSTN-203_A11** epitope (supp S5/S6: IFN-γ ELISA + caspase-3/LDH cytotoxicity; IEDB ref1046253, 4× Positive). **Corrected in [#838](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/838)** from `KTEGPTLTK`, an `agent-web` mis-key that had recorded the top *predicted* A\*11:01 binder (supp prediction table MOESM4) instead of the assayed peptide — `KTEGPTLTK` is **absent from IEDB**. HLA typed at 2-digit (paper: A\*11+ line SF10281; IEDB: `A11`); `A*11:01` retained as the predominant A11 allele for scoring. (Prior note: "⚠️ NOT in Zotero / verify" — resolved; the paper is Zotero `HMXA22SW`.)

## Kwok et al. 2025, *Nature* (s41586-024-08552-0; PMID 39972144; Zotero item `5ZT8KC8X`, PDF `EA2NMFNZ`) — `direct` (recovered #838)

- **GNAS `SLLLPSFHL`** & **RPL22 `GIMDAANFFL`** / A\*02:01 — the two public NEJ-derived neoantigens (NeoA-GNAS frameshift + NeoA-RPL22 in-frame). Functionally validated across **Fig 4b–h, Fig 5d (MS), Extended Data Fig 6e/f + 9b** — exact 9/10-mers **not legible in the paper figures**, but recovered from the authors' **IEDB deposit ref1045680** in [#838](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/838): `SLLLPSFHL` is IEDB-annotated as GNAS; `GIMDAANFFL` (gene-unannotated in IEDB) confirmed = RPL22 by mapping onto UniProt **P35268** WT `GIMDAANF·EQ·FLQER` with the paper's in-frame 6-nt (EQ) loss → `GIMDAANFFL`. **Promoted non-scorable → `functional-scorable`.** Supersedes the prior illustrative `~LALDVLQGYSL` placeholder (do not use).

## Manoharan et al. 2026 — intron-retention neoantigens, colorectal cancer (#733 AC#1)
Zotero `ZAT8678F` (no local PDF) · open access **PMC13096189** · all `direct` (upgraded from `agent-web` in [#823](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/823): a third independent first-hand read of Table 1 "Selected 9-mer peptides for in vitro validation" at the open-access PMC source confirmed all 11 sequences/genes/alleles verbatim, **11/11**, agreeing with the two prior PMC reads — the exact table is now captured, discharging the `agent-web` "verify most carefully" caveat). The article is permanently citable at PMC13096189; a local Zotero PDF copy remains optional (open-access source is stable).
Provenance note: all functional readouts are **healthy-donor moDC-primed in-vitro-sensitized T cells, NOT patient T cells** → positives `medium` (IVS-grade); negatives are *soft* (`candidate-negative`). Captured in the #823 `assay_context` column as `healthy_donor_ivs` (all 11 rows).

Table 1 "Selected 9-mer peptides for in vitro validation" (IR1–IR11). **3 immunogenic** (IR3/IR7/IR9), **8 tested-negative**:
- **`ILAFIAPLK`** CHD7 / A\*11:01 — IR3, immunogenic (ELISpot + tetramer across 2 donors).
- **`RQDPAPQQV`** FERMT3 / A\*02:01 — IR7, immunogenic (tetramer 2.56% + GZMB reporter + CD137). **Restricted on A\*02:01** (functional validation on A\*02:01 aAPCs only; A\*24:07 listed in Table 1 but not functionally used).
- **`AVLHGRLFL`** TPCN2 / A\*02:06 — IR9, immunogenic (IFN-γ ELISpot + GZMB + CD137).
- Tested-negative (no T-cell response, healthy-donor IVS): `LSCSPMMRK`/LRRC8B (A\*11:01, IR1), `QLLGPWVFK`/OPLAH (A\*11:01, IR2), `NLFQVMHIK`/FANCB (A\*11:01, IR4), `AVVAMVTLM`/PLEKHG5 (A\*34:01, IR5), `MWSFIHNNL`/WDR19 (A\*24:07, IR6), `SLPPGLRGT`/RNF123 (A\*02:07, IR8), `TGRGVLCLR`/ITGAM (A\*31:01, IR10), `SLLASSPAR`/CDH24 (A\*11:01; also predicted A\*31:01, IR11).

## Merlotti et al. 2023 — exon-TE (noncanonical splice) junctions, NSCLC (#733 AC#2)
Zotero `5ZADNVCB` (local main PDF + supp S1–S8) · all `direct` (sequences read first-hand from **Data file S7**; patient reactivity from **Data file S8** Fig 5D/5G + Fig S5C/D; cloned-TCR functional from Fig S4I/J). All **HLA-A\*02:01** (S7 = predicted-A2 in-vitro-assay set), validated in A2⁺ **patient** TILs/PBMC ex vivo.

Sequence source = S7 "Characteristics of predicted pJETs used for in-vitro assays" (ID → sequence). Reactivity tiering = patient tetramer⁺ CD8 across Fig 5D (ex-vivo), 5G (Day-20 TIL), S5C/D (juxta-tumor / draining-LN).

**Folded — TCR-validated (`high`):** `LLGETKVYV`/SRSF7 (SVA_D junction, ID 3228-1), `LLDRFGYHV`/FRG1B (HERVL74 junction, ID 321), `YLWTTFFPL`/SCGB3A2 (HERVH frameshift, ID 704) — each has patient tetramer⁺ **and** cloned-TCR Jurkat-TPR NFAT/NF-κB/AP-1 activation (Fig S4I/J).

**Folded — multi-context patient tetramer⁺ (`medium`, tetramer-DETECTION only):** `SLMQSGSPV`/SBNO2 (299), `AILPKANTV`/DCBLD2 (4111), `YLPYFLKSL`/KYNU (2242), `ILSGYGPCV`/TENM1 (2930-1), `ILANLPPAL`/ATAD3B (1825), `RLLHLESFL`/SLC39A11 (1375), `VLMWTMAHL`/PDE4D (4188).

**NOT folded — single-replicate patient tetramer⁺ candidates (recorded only here; below the scorable bar):** `ILTASITSI`/ERO1LB (3134-2), `AMDGKELSL`/MAST4 (1146), `CLIDEMPEA`/HDGF (3479), `ILPKANTVV`/DCBLD2 (W13; register-shifted vs 4111), `LLHLESFLV`/SLC39A11 (W2; vs 1375), `MLMKTVWQA`/RLN3 (3443), `GLLNISHTA`/CUL4B (2368), `ILHSLVTGV`/ERO1LB (3134-1), `YLQGLPLPL`/FUT8 (1000), `FLGTRVTRV`/COL28A1 (2914). Re-assess on a bar change. (Fig 4D donor-1–6 reactivity is **healthy-donor** screening, not patient — not used for eligibility.)

## Zhao et al. 2025 — AS-derived mRNA-vaccine neoantigens, HCC (#733 AC#3 — registry-eligible, sequence-blocked)
Zotero `B2MJ776X` (19 local supplements; **main paper PDF not local**). Gate-1 (alternative splicing) + gate-2 (human per-peptide IFN-γ ELISpot / tetramer / 4-1BB on HCC patient TILs, A\*02:01/A\*11:01/A\*24:02) **both pass** — but the 19 supplements reference peptides only by internal ID (pA02-28, pA11-10, …); the amino-acid strings are not present locally (MOESM19 is an ssGSEA gene list, not the peptide table). **No rows folded** — sequences must be recovered from the main text / a peptide-sequence table before entry. Do NOT translate Table-S1 junction coordinates (inference). Tracked in README "Next steps".

## Bigot et al. 2021 — SF3B1-mutant uveal melanoma (#734, IEDB free-text recovery)
Zotero `DJAS2BJ2` (supp tables S1–S9 + methods local; **main-text PDF not local**; Correction `EMSRGMQZ` / PMID 35257149 **has no attachment — unread**). Recovered via the #734 IEDB IQ-API free-text mine (reference 1040678). All **HLA-A\*02:01**.

Sequence source = **local Supplementary Table S2** ("peptides" sheet, the `A2:N` synthesized panel) — read first-hand → `direct`. All 35 IEDB-positive sequences matched S2 verbatim (35/35), so the unread Correction cannot have silently altered them. Source gene = IEDB `r_object_source_molecule_name` (30/35; 5 panel peptides have no gene in IEDB → `not-named`). Splice mechanism = SF3B1-mutant aberrant 3′ splice-site selection (the paper's central finding; supp. methods confirm "endogenous aberrantly spliced" transcripts). Confidence split is **dual-source-grounded**: the 5 with effector function are also the 5 with explicit alternative-mRNA-junction validation in **Supplementary Table S8** (frameshift + predicted NMD) — full table in [`issue_734_db_audit/`](issue_734_db_audit/).

**Folded — effector + S8 junction-validated (`high`, 5):** `LLIRWQHFL`/NF1 (A2:14, frameshift +14 NMD+; TNF-α+activation), `AALPILFQV`/ATP8B2 (A2:17, +20 NMD+; degranulation), `ALLLQLFTL`/USP39 (A2:18, +40 NMD+; IFN-γ+cytotoxicity+granzyme B+CD107), `ALLPGLPAA`/NET1 (A2:26, +22; IFN-γ+cytotoxicity+degranulation), `RLPGVLPRA`/MAPK8IP2 (A2:37, +16; IFN-γ+cytotoxicity+degranulation). Sequences cross-confirmed in JEM 2024 (PMC10986814) by `UM-A2-NN` ID.

**Folded — patient ex-vivo tetramer⁺ only (`medium`, 30):** the rest of the A2:N positives — `FLSSVALAL`/ZDHHC16, `FLPPGGAPV`/VPS51, `MLAAPISGL`/SLC3A2, `GLVETFYFT`/HADHA, `FQVQSLPAA`/OXA1L, `CLFSPQLLV`/not-named, `FMALDDPVI` & `ALDDPVIFL`/SEPSECS, `FLSRKLSSL`/MINDY4, `KLLEDDLLI` & `LLYLGIIKL`/NF1, `TMDVVIPGL`/VARS2, `VLFLTLNEV` & `FLTLNEVFL` & `TLNEVFLIL` & `FLILQVHGT`/ZFYVE27, `FLWCFPFSL`/SOAT1, `SLLPYVFVI`/RNF38, `FLFIGMSQM`/ARIH1, `SLSPALPGA`/SF3A2, `FMDDAKILF`/not-named, `FLFLDTIRT`/MRPS10, `ALDNVDAPL`/UBA1, `RLFPISPET`/NET1, `SLTFLPLYV`/ATP5MC2, `AVAEATAGV`/not-named, `KLMGPEVQV`/SRSF1, `FILKPILFL`/CTDNEP1, `RLGEVRHPV`/not-named, `GLNTFALGL`/not-named. Tetramer-DETECTION only (no per-peptide effector or junction shown in the supplements we hold); splice origin via panel design + IEDB gene.

**NOT folded:** the 4 IEDB-negative panel peptides (`AILPVIADL`/SRSF1, `FLPLYVIPA`, `GVLPRAPGV`, `RLPQGVYPV`) and the 3 controls (Melan-A `ELAGIGILTV`, CMV `NLVPMVATV`, Flu `GILGFVFTL`).

## Long-read uveal melanoma — Cancer Immunol Res 2023 (#838, IEDB free-text recovery)
PMID 37756564; DOI 10.1158/2326-6066.cir-23-0083; Zotero `6P5JCCIB` (supp figures `GBXC9V2P`, supp tables `SNAB6XA5`). IEDB ref1043193. `direct` but **`functional-nonscorable`**.

- **`MADAGAMAA`** & **`RQVGEGCRT`** / SEPTIN6, **`KNILNGSRA`** / AMZ2P1 (archaemetzincin-2 *pseudogene 1*), **`RGAGLGRAL`** / MZT2B — alt-splicing neojunctions; gene + neojunction-ORF context from **supp Table 13** ("Sequence information for SEPTIN6, AMZ2P1 and MZT2B"; the long ORFs `MADAGAMAATDIARQVGEGCRT` / `…SVKNILNGSRATVKHISIA` / `…CPRGAGLGRALHGAAPRMGQRL`). Gate-1 confirmed. **No per-peptide HLA published** (IEDB = `human`; IFN-γ ELISPOT on bulk PBMC) → no allele key → non-scorable (retained for the validation base). `RGAGLGRAL` and `RQVGEGCRT` are IEDB mixed Negative/Positive(-Low).

## TEtrans / Li et al. 2025 - transposable-element-derived transcript neoantigens, gastric cancer (#832, transcript-level source; no peptide in source)
bioRxiv 2025.03.01.640928; Zotero `IJ55K4WH` (full 55-page PDF read 2026-06-22). Surfaced via the splice-immunogenicity standing watch. **Gate-1 (exon-TE / TE-junction source) and gate-2 (human functional T-cell) both pass, but no peptide sequence is resolvable at source -> no rows folded.** Represented here as a transcript/antigen-level exon-TE source annotation only.

- **Passes at transcript/antigen level:** T4180 (LTR30) and T2432 TEtrans transcripts elicit **HLA-A\*02:01-restricted CD8+ killing + IFN-gamma**; T4180 is **MS-presented in gastric (STAD) tumor tissue** (STAD = our chr22 test tumor type); validated in 2 gastric-cancer patients plus an in vivo mRNA-vaccine model.
- **Fails the peptide gate (why no row):** all functional validation is **transcript-level** (mRNA transfection / K562-expressing-TEtrans), never minimal-peptide pulsing. The 387 candidate neopeptides are reported **count-only** (Fig 6A funnel) - **no sequence table** in the embedded supplement, **no data deposit** (Data Availability lists only the ASJA tool repo). A peptide-keyed `registry.tsv` row is not constructable without fabricating a sequence, so none was added (no-inference rule: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) section 6). Comparable to the Zhao 2025 (sequence-blocked) and long-read UM (functional-nonscorable) entries above.
- Zotero note `NIGQEC6J` "vs our pipeline" line was refreshed at filing (2026-06-22) to reflect transcript-level-only / no peptide in source (previously stale: "row needs minimal-epitope check").

---

## Junction mapping grades (#735)

The `junction_mapping_grade` column added in [#735](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/735) grades how specifically each peptide is traceable to its originating junction (`coords > event-id > gene-mechanism > none`).
Grades are assigned strictly from evidence quoted in this file (PROVENANCE.md); no junction coordinate or identifier is ever inferred from the peptide sequence or a genome reference (no-inference rule: [`LABELING_SCHEME.md`](LABELING_SCHEME.md) §6).
Per-source grade rationale: [`junction_evidence_by_source.md`](junction_evidence_by_source.md).

## Curation chain (how the registry was built)

1. Adversarially-verified novelty check (25/25 claims held) — confirmed no prior open splice-immunogenicity benchmark exists.
2. Multi-agent literature + 4-database sweep → 3-vote dual-gate verification.
3. Local Zotero full-text PDF mining (paywall-free) — surfaced Kim/Zhao/Ji/Fisher (missed by the web sweep), doubled the set, added intron-retention biology.
4. Human + machine read of SNAF + Kwok supplementary files — recovered the HLA-C sequences; confirmed GNAS/RPL22 unpublished.
