# Deep-research run #1 - splice-immunogenicity registry enrichment

**Date:** 2026-06-30 (completed ~23:30 UTC) · **Run ID:** wf_870575e5-9a3 · **Task:** wxq6ep3w3
**Raw output (durable):** `tasks/wxq6ep3w3.output` (session f1a362d9-5aad-4b03-ac35-0ff07b9e3ae4)
**Scale:** 36 agents, ~3.0M subagent tokens, 547 tool uses, ~48 min.
**Pipeline:** 6 discovery lenses -> consolidate/dedup -> adversarial dual-gate verify -> synthesize.
**Counts:** 56 raw -> 28 consolidated -> **16 CONFIRM_ADD**, 4 NEEDS_SOURCE_ACCESS, 10 REJECT. (Funnel is approximate: the three verify buckets sum to 30, not the 28 consolidated - a 2-row discrepancy at the consolidation/verify boundary; treat the counts as indicative, not exact.)

## Honest headline

- 🎯 **Hard true-negatives: ZERO new.** The prize constraint did NOT move. Field still holds exactly one (`VELEDHVML`/Fisher 2026). Every "negative" recovered here is SOFT (Tg-mouse vaccination priming-failure, no MS-presentation), categorically not a hard negative.
- 🌐 **A\*02:01 monoculture: moved, meaningfully.** First sizeable non-A2 cluster = **5 HLA-A\*24:02 functional positives** from Oka 2021 (Genome Biology, HLA-A24 transgenic-mouse peptide-vaccination IFN-g ELISpot).
- ⚠️ **Caveat tempering the win:** all 12 Oka rows (5 pos + 7 soft neg) come from ONE study, ONE allele, a SURROGATE species (Tg-mouse de-novo vaccination, not patient ex-vivo / not human cloned-TCR). Diversifies the HLA axis but concentrates provenance + assay-context risk.
- The verifier earned its keep: caught a **label-flip** (PKM `SFLGFSILL` proposed negative -> actually a weak positive on a 400-dpi figure read) and **rejected 4 Oka frameshift-DNA peptides** as gate-1 fails (not splice).

## The 16 confirmed registry-ready rows (NOT yet folded into registry.tsv)

Column order = registry.tsv schema. These await a fold PR (#838/#904-style) + a decision on the Tg-mouse surrogate caveat.

### Non-A2 positives (5) - HLA-A*24:02, Oka 2021 (Zotero S8ZVY4YB; DOI 10.1186/s13059-020-02240-8)
- `TYTTIKINF` HOOK2 - alt last exon - strong - functional-scorable
- `FLLPAPFPF` TUFM - intron retention - strong - functional-scorable
- `VYFTSDFKV` PKM - alt last exon - strong - functional-scorable (label-conflict resolved to positive)
- `SYFETIAAL` SELENBP1 - intron retention - soft - functional-scorable (one-mouse-dominant)
- `SFPLVFLFF` CEACAM6 - alt last exon - strong - functional-scorable (weakest of the 8)

### Soft negatives (7) - HLA-A*24:02, Oka 2021, candidate-negative tier (NOT hard; no MS-presentation)
- `KFSPEPSQF` CEACAM6, `NYFNLGMVV` ERO1A, `AWPKHLDLM` MCEE, `KYEEVAQLY` UQCRB, `QYSLATAFL` SCGB3A2, `RFQPHGDGW` TMC4, `TWLTSGPHL` TMEM45A
- All: splice-junction-derived, peptide-pulsed IFN-g ELISpot NO response in Tg mice, NetMHC-selected (not MS-presented) -> soft.

### Presentation-prevalence decoys (4) - #681 feeders, untested
- `LYIPALAAL` FGGY (HLA-A*24:02) - Zhao & Andersen iScience 2026; MS-presented, authors disclaim T-cell validation. **Non-A2 decoy.**
- `FIQENMVMMVA` PDIA3 (A*02:01 or A*11:01) - Zhao 2026 PDAC bioRxiv; MS-presented, transcript-level functional only.
- `IVLPPWPPK` unnamed neoTST (predicted A*02:01/A*11:01) - Zhao 2026 PDAC bioRxiv; MS-presented, mixed-mRNA assay only.
- `QANSFPLTF` unnamed neoTST (predicted A*02:01) - Zhao 2026 PDAC bioRxiv; lowest-confidence, verify junction+presentation before any fold.

(Full verbatim TSV rows with provenance_ref + label_rationale are in the raw output file `registry_ready_rows`.)

## NEEDS HUMAN ACCESS (4) - the route to the field's 2nd hard negative is here

These are the proactive "please fetch" asks (per feedback_ask_human_for_files):

1. **Zhao et al. PDAC neoTST preprint** - bioRxiv 10.64898/2026.02.10.705024 (Zotero XT37BWFZ). Supp tables: (a) the 10 neoTST AA sequences, (b) per-neoTST/allele IFN-g ELISpot status (the NEGATIVE subset behind Fig 6D-G / S6C-D), (c) the LC-MS/MS immunopeptidome list, (d) per-neoTST splice annotation. **PRIORITY hard-negative feeder.** CAVEAT: main assay is mRNA-transfection -> confirm whether any MINIMAL-PEPTIDE (synthetic) ELISpot-negative exists, else negatives stay transcript-level (soft).
2. **Lin/Zhao et al. HCC neoTST preprint** - bioRxiv 10.64898/2025.12.07.692877 (Zotero T3UNU59Y). Supp Table S1 (14 neoTST fragments + HLA) + per-neoTST ELISpot behind Fig 2J-K (7/14 non-responders) + immunopeptidome table. Same Fudan group, same logic + same transcript-level caveat.
3. **Courcelles et al. 2026 MCP** CRC immunopeptidome (Zotero GWAES24H). MS-eluted HLA-B\*44:03 + A\*11:01 splice ligands. NOT a current candidate (in-silico immunogenicity only) but the strongest UNTAPPED goldmine: running ELISpot/tetramer on these would generate non-A2 hard negatives - moving BOTH constraints at once. (Requires wet-lab, not just a fetch.)
4. **Zhao 2025 HCC AS-vaccine** Cell Research (Zotero B2MJ776X / 9XTA5B2L; Issue #817). Main-text peptide table mapping internal IDs (pA02-*/pA11-*/pA24-*) to AA sequences - sequence-blocked across all 18 local supp PDFs + MOESM19 xlsx; main PDF not in Zotero.

## Completeness critique + single highest-value next probe

> Obtain the two neoTST preprint supplementary tables (PDAC + HCC) AND explicitly determine whether a minimal-peptide (synthetic, not mRNA) ELISpot/tetramer negative exists for any MS-presented splice neoTST. **That one check decides whether the field's hard-negative count can move from 1 to 2.**

Other unprobed: (2) the 12 Oka adds rest on one 2021 figure read, uncorroborated in IEDB - a high-res source-data request to the Oka authors would harden them. (3) No PRIDE/MassIVE raw immunopeptidome reanalysis was attempted - where uncatalogued splice-junction MS evidence lives.

## Rejects (10), by reason
1. **Gate-1 fail (frameshift-DNA, not splice):** Oka `TYMAGSEAW`+`VWPTAWASL` (NDST1), `RYCRGVMLL` (ARL8B), `SFTDISIYL` (SENP2). 3 of 4 were ELISpot-POSITIVE but ineligible (not splice).
2. **Label flip:** PKM `SFLGFSILL` proposed negative -> actually weak positive (Fig 6c asterisk; reconciles paper's "8 of 17").
3. **Out-of-scope species:** Ji 2025 MC38 LSTV peptides pass both gates but in MOUSE H2-Kb/Db (registry is human-HLA-scoped; routes to #681).
4. **Gate-2 fail (predicted/MS-only/transcript-level):** Wickland 2024 SPLICE-neo, Courcelles CRC, mRCC AS, NeoGuider, Smart 2018, ~20 method papers.
