# Issue #1207 - snowball registry-additions sweep

Backward-citation (snowball) sweep of three human splice-neoantigen catalogs surfaced from the reference lists of two 2026-07-16 Zotero adds.
Method: two-gate verification (Gate 1 = genuinely splice-junction-derived; Gate 2 = per-peptide minimal-peptide functional T-cell assay) + 3-vote adversarial refutation, per the registry method.

## Outcome: zero functional-registry additions (registry stays at 97 rows)

| Catalog | Zotero | Determination |
|---|---|---|
| Exitron pan-cancer (Wang et al., Mol Cell 2021) | `HTUR4MKR` | 54 exitron-derived splice neoantigens (Tables S4+S5), **all MS-presentation-only, no functional T-cell assay** -> Gate 1 PASS, Gate 2 FAIL. Zero functional rows; all 54 route to the MS-presented tier ([Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681)). |
| U2AF1 neoantigens (Biernacki et al., JITC 2023) | `GUUF7E7H` | **Out of scope - Gate 1 FAIL.** The epitopes are the U2AF1 **Q157R missense** mutant-protein sequence (plus ScanProsite self-peptide controls), not aberrant-junction peptides. Functionally validated (cloned TCR, tetramer) but not splice-derived. |
| Proteogenomics-ovarian (Zhao et al., Cancer Immunol Res 2020) | `3Z5UHGFQ` | 103 MS-eluted TSAs (Table S5: 20 mTSA + 83 aeTSA). A non-canonical / aberrant-expression repertoire, **not splice-junction-derived**. Per the main text (source-read), aeTSAs arise from *aberrantly expressed* nonexonic regions (intronic 29%, intergenic 22%) regulated by copy number + DNA methylation - **transcriptional dysregulation, not splicing**. The 28 Intronic aeTSAs are therefore **low-confidence splice-adjacent candidates only** (routed to [Issue #681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681) flagged low-confidence, not confirmed splice); the other 75 are out of splice scope. Gate-2 MS-only is **source-confirmed** (main text: "mass spectrometry pipeline", no T-cell assay) -> zero functional-registry rows. |

## Artifact

`exitron_ms_presented_S4_S5.tsv` - the 54 exitron MS-presented splice neoantigens recovered from the paywalled supplement (`mmc1.pdf`, Document S1, Tables S4+S5). Columns: peptide, gene, chr/start/end (exitron locus), sample, cohort/MHC-type, table, mhc_class, gate flags, route, provenance. 33 class-I, 21 class-II. Neither table gives the specific HLA allele (S4: none; S5: MHC class only), so these enter #681 at class resolution.

## AC-4 payoff cross-reference (against our n=2156 candidate presenters)

Both catalogs' surviving MS-presented peptides were intersected against our n=2156 splice presentation scores (patient_001 n=395 + patient_002 n=1761, fetched from R2):

- **Exitron** - 33 class-I peptides -> **zero overlap** (exact + substring).
- **Ovarian** - all 103 TSAs (class-I, allele-resolved) -> **zero overlap** (exact + substring).

Verified as a real negative via a ceiling control (presenters span 8/9/10-mers, matching the catalogs' class-I lengths; positive control fires). Biologically expected: independent cohorts, disjoint HLA (our patients carry A\*31:01 / A\*26:01 / B\*18:01 / B\*15:63 / C\*07:01 / C\*03:03), no shared recurrent target.

Reproduce: `xref_exitron_vs_presenters.py` (R2 creds from project-root `.env`, `research/.venv/bin/python`).
