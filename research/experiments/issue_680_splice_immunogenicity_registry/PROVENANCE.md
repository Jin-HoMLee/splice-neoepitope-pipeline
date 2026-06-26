# Per-peptide provenance — splice-immunogenicity registry (#680)

Each sequence is graded by **how it was obtained**, so it can be re-verified manually. Local files are under `~/Zotero/storage/<KEY>/`.

## Provenance grades

| Grade | Meaning |
|---|---|
| `direct` | Exact sequence read first-hand at the cited figure/table/sheet |
| `inferred` | Sequence derived by matching a figure's 3-letter prefix to the unique catalog entry — **verify against the synthesized-peptide list** |
| `agent-local` | Extracted by an analysis agent from the **local Zotero PDF**; location recorded, not every panel personally eyeballed |
| `agent-web` | Extracted by a web-reading agent; source **not in Zotero**; exact figure/table not captured — **verify most carefully** |
| `unpublished` | Exact sequence not published in any legible form |

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
File: `LX6DMXTL` · all `agent-local`

- **`RLWGTWVKA`** CLK3 / A\*02:01 — **Fig S1B** *"CLK3 #1 (RLWGTWVKA)"*; MS Fig 3G / S4A; Fig S3A.
- **`CLLPPALFL`** RHOT2 / A\*02:01 — **Fig S1B** *"RHOT2 #5 (CLLPPALFL)"*; MS spectrum Fig 2C; Fig 7E.
- **`RLLAAVLEA`** c16orf70 / A\*02:01 — **Fig 2H** (MS + label *"RLLAAVLEA"*); Fig S1B *"C16orf70 (RLLAAVLEA)"*.
- **`LVAKWLLTI`** ATG3 / A\*02:01 — **Fig S2F** *"LVAKWLLTI"*; Fig S1B *"ATG3 (LVAKWLLTI)"*.
- **`FMDDYIFV`** MYO1F / A\*02:01 — **8-mer, confirmed (#733 AC#4).** Fig S2I legend (verbatim, from ft-cache): *"…the MYO1F peptide (FMDDYIFV). (J) FACS plot…"* — the functionally-assayed peptide is the 8-mer; the `RFMDDYIFV` in Fig S1B/S3B was the source-region label, not the tested peptide.

## Xiong et al. 2025 — RCAN1-4, glioblastoma
File: `JJZL6D84` · `agent-local`

- **`IFSESETRAKF`** RCAN1-4 / A\*24:02 — **Fig 3B** peptide table (*"RCAN1-4_22-32, %Rank EL 0.488, Strong"*); Fig 3A structure (*"…NSDIFSESETR–AKFESLFRTY"*, SJ exon4/exon5); functional Fig 3D–K, 5, 6.
- **`VFVDGLCRAKF`** RCAN1-1 control / A\*24:02 *(negative)* — RCAN1-1_77-87 constitutive control; **verify in Fig 3 control panels**.

## Fisher et al. 2026 — CoREST
File: `WLSZLMA4` · `agent-local` · genes not named in source

- **`RSQGWLFLR`**, **`AEHAHRVPL`**, **`VELEDHVML`** *(VELEDHVML is the negative)* — all **Fig 7I** (peptide table), **Fig 7J/K** (ELISpot), **p.15**; context text p.12. Loaded on A\*11:01 APCs (`AEHAHRVPL`/`VELEDHVML` predicted B\*40:01 — restricting allele ambiguous).

## IRIS — Pan et al., *PNAS* ⚠️ NOT in Zotero — verify in supplement
all `agent-web` (web read of the paper + **Dataset S3c**)

- **`SLDGTTTKA`** CLASP1 / A\*02:01 — TCR JPTCR-238 (IFN-γ ELISA + Incucyte cytotoxicity).
- **`STMYYLWML`** SCAMP3 / A\*02:01 — TCRs JPTCR-45/47/50/56 (IFN-γ).
- **`LLLGIAKLLKV`** MAN2C1 / A\*02:01 — TCR JPTCR-13. **Medium** (2nd verifier did not confirm).
- **`FLSELEPPA`** AP2A1 / A\*02:01 — TCR JPTCR-52. **Medium** (2nd verifier did not confirm).

## POSTN-203 study ⚠️ NOT in Zotero — verify
- **`KTEGPTLTK`** POSTN-203 (ENST00000379747) / A\*11:01 — `agent-web`; extraction agent flagged `partial`. Lowest certainty; verify against the source directly.

## Kwok et al. 2024/2025, *Nature* (s41586-024-08552-0) — `unpublished`
File: `EA2NMFNZ` (main + Extended Data)

- **GNAS** & **RPL22** / A\*02:01 — functionally validated across **Fig 4b–h, Fig 5d (MS), Extended Data Fig 6e/f (MS) + 6i (NJ aa-sequence track, generic schematic), Extended Data Fig 9b (alanine scan)** — but exact 9-mers are **not legible** (MS peak plots; 6i illustrative). RPL22 region ≈ `…LALDVLQGYSL…` (illustrative, do not log). Recoverable only via author contact / raw MS. Recorded as validated, non-scorable.

## Manoharan et al. 2026 — intron-retention neoantigens, colorectal cancer (#733 AC#1)
Zotero `ZAT8678F` (no local PDF) · open access **PMC13096189** · all `agent-web` (two independent PMC reads of Table 1 + Results/Fig 3 agreed on all 11 seqs/genes/alleles)
Provenance note: all functional readouts are **healthy-donor moDC-primed in-vitro-sensitized T cells, NOT patient T cells** → positives `medium` (IVS-grade); negatives are *soft* (`candidate-negative`).

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

Sequence source = **local Supplementary Table S2** ("peptides" sheet, the `A2:N` synthesized panel) — read first-hand → `direct`. All 35 IEDB-positive sequences matched S2 verbatim (35/35), so the unread Correction cannot have silently altered them. Source gene = IEDB `r_object_source_molecule_name` (30/35; 5 panel peptides have no gene in IEDB → `not-named`). Splice mechanism = SF3B1-mutant aberrant 3′ splice-site selection (the paper's central finding; supp. methods confirm "endogenous aberrantly spliced" transcripts). Confidence split is **dual-source-grounded**: the 5 with effector function are also the 5 with explicit alternative-mRNA-junction validation in **Supplementary Table S8** (frameshift + predicted NMD) — full table in [`db_audit_734/`](db_audit_734/).

**Folded — effector + S8 junction-validated (`high`, 5):** `LLIRWQHFL`/NF1 (A2:14, frameshift +14 NMD+; TNF-α+activation), `AALPILFQV`/ATP8B2 (A2:17, +20 NMD+; degranulation), `ALLLQLFTL`/USP39 (A2:18, +40 NMD+; IFN-γ+cytotoxicity+granzyme B+CD107), `ALLPGLPAA`/NET1 (A2:26, +22; IFN-γ+cytotoxicity+degranulation), `RLPGVLPRA`/MAPK8IP2 (A2:37, +16; IFN-γ+cytotoxicity+degranulation). Sequences cross-confirmed in JEM 2024 (PMC10986814) by `UM-A2-NN` ID.

**Folded — patient ex-vivo tetramer⁺ only (`medium`, 30):** the rest of the A2:N positives — `FLSSVALAL`/ZDHHC16, `FLPPGGAPV`/VPS51, `MLAAPISGL`/SLC3A2, `GLVETFYFT`/HADHA, `FQVQSLPAA`/OXA1L, `CLFSPQLLV`/not-named, `FMALDDPVI` & `ALDDPVIFL`/SEPSECS, `FLSRKLSSL`/MINDY4, `KLLEDDLLI` & `LLYLGIIKL`/NF1, `TMDVVIPGL`/VARS2, `VLFLTLNEV` & `FLTLNEVFL` & `TLNEVFLIL` & `FLILQVHGT`/ZFYVE27, `FLWCFPFSL`/SOAT1, `SLLPYVFVI`/RNF38, `FLFIGMSQM`/ARIH1, `SLSPALPGA`/SF3A2, `FMDDAKILF`/not-named, `FLFLDTIRT`/MRPS10, `ALDNVDAPL`/UBA1, `RLFPISPET`/NET1, `SLTFLPLYV`/ATP5MC2, `AVAEATAGV`/not-named, `KLMGPEVQV`/SRSF1, `FILKPILFL`/CTDNEP1, `RLGEVRHPV`/not-named, `GLNTFALGL`/not-named. Tetramer-DETECTION only (no per-peptide effector or junction shown in the supplements we hold); splice origin via panel design + IEDB gene.

**NOT folded:** the 4 IEDB-negative panel peptides (`AILPVIADL`/SRSF1, `FLPLYVIPA`, `GVLPRAPGV`, `RLPQGVYPV`) and the 3 controls (Melan-A `ELAGIGILTV`, CMV `NLVPMVATV`, Flu `GILGFVFTL`).

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
