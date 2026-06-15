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
- **`FMDDYIFV`** MYO1F / A\*02:01 — **Fig S2H** *"FMDDYIFV/FMDDYIFVP"*; Fig S1B/S3B context *"MYO1F (RFMDDYIFV)"*. ⚠️ confirm exact tested length (8-mer vs R-extended).

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

---

## Curation chain (how the registry was built)

1. Adversarially-verified novelty check (25/25 claims held) — confirmed no prior open splice-immunogenicity benchmark exists.
2. Multi-agent literature + 4-database sweep → 3-vote dual-gate verification.
3. Local Zotero full-text PDF mining (paywall-free) — surfaced Kim/Zhao/Ji/Fisher (missed by the web sweep), doubled the set, added intron-retention biology.
4. Human + machine read of SNAF + Kwok supplementary files — recovered the HLA-C sequences; confirmed GNAS/RPL22 unpublished.
