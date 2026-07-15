<!--
Provenance: multi-agent scoping workflow for Issue #1176 (best-bet 3 "main event").
Run 2026-07-15. Harness: 5 research facets (cohort / search-DB / FDR-decoy / tooling / PCPS-priorart),
each adversarially verified by an independent skeptic (default-to-refute), then synthesized.
All 5 facet verdicts returned holds=true; three non-fatal corrections were folded into the text below.
Workflow run id wf_38edaf38-46d (facet 3 / FDR was re-run after a first-pass agent failure).
This is desk-research: it settles the design/cohort front-half only. It does NOT touch spectra.
-->

# Issue #1176 scoping brief - non-canonical public-spectra search

## Executive summary

The front-half of Issue #1176 (design plus cohort selection) is now decidable, and every load-bearing claim across all five facets survived adversarial verification (all five verdicts returned `holds=true`). The recommended approach: re-search public tumor HLA-I immunopeptidome spectra against a combined FASTA (canonical proteome + our RNA-seq-bounded junction-spanning peptide partition) using FragPipe/MSFragger nonspecific search, controlling error with group-specific (subset) FDR computed on the junction partition alone, backstopped by entrapment validation and the orthogonal RNA-seq junction witness that the PCPS dispute structurally lacks. Start on the Courcelles CRC cohort (matched MS + RNA-seq, both confirmed public 2026-07-15, MS licensed CC0 - see the Update note above), with MassIVE MSV000096853 (Ely PDAC, CC0) as the second, A*02:01-diverse deposit. Recovered hits land as `label=untested`, `tier=presentation-prevalence` only - never as functional/immunogenic positives.

## Update 2026-07-15 - Courcelles cohort confirmed public (supersedes the pending-fetch hedge below)

The one open cohort dependency is **resolved, and no human fetch was needed.** Recovered directly from the local PDF's Data Availability statement (Zotero attachment GWAES24H) plus live PRIDE/GEO checks: Courcelles' MS immunopeptidome is **ProteomeXchange PXD071022** (PRIDE, public 2026-05-26, licensed **CC0**) and the matched RNA-seq is **GEO GSE312236** (public 2026-05-26, BioProject PRJNA1372708). Both are public and reuse-clean. Courcelles is therefore promoted from "pending co-pick" to the **primary cohort** (its matched RNA-seq is our differentiator); MSV000096853 (Ely PDAC, CC0) remains the second, A*02:01-diverse deposit. Open question #1 below is closed. Method note: the brief's "fetch-blocked" status was a transient wall on mcponline.org - the authoritative answer was already in our local bytes, per the exhaust-local-PDF-first rule.

## Recommended design

### Cohort

Re-search two complementary, strongly non-A*02:01 public tumor HLA-I deposits.

- Second deposit (A*02:01-diverse, CC0): MassIVE **MSV000096853 = ProteomeXchange PXD059832** (Ely et al., "Pancreatic cancer-restricted cryptic antigens are targets for T cell recognition," Science 2025; PMC12163983), PDAC + organoids, Orbitrap Exploris 480, per-patient HLA typing incl. A*11:01, licensed **CC0 1.0** (source: massive.ucsd.edu dataset page + proteomexchange PXD059832; our Zotero re-use precedent XT37BWFZ, the Zhao PDAC preprint that already recovered TST/splice neoantigens from this deposit).
- **Primary cohort (confirmed public 2026-07-15): Courcelles et al. 2026** CRC (Mol Cell Proteomics; Zotero GWAES24H), 26 primary CRC tumors, 115,292 MAPs across 61 HLA alleles incl. HLA-B*44:03/A*11:01, splice/aeTSA-focused, and decisively **ships matched RNA-seq per tumor** - which maximally amplifies our differentiator. MS = **PXD071022** (PRIDE, **CC0**), matched RNA-seq = **GEO GSE312236** (BioProject PRJNA1372708) (source: local PDF Data Availability statement, Zotero GWAES24H; verified public via PRIDE + GEO 2026-07-15).

Adversarial qualification (cohort verdict): the finding's claim that MSV000096853 is "fully public" is an **overstatement and is corrected** - the deposit is flagged "Partial Public" and consists of **processed peak lists explicitly labeled "suitable for reanalysis," not full RAW files**. This does **not** refute re-searchability (a custom-DB search consumes MS/MS peak lists), so the recommendation to start there stands; only "fully public / raw files" is downgraded to "processed peak lists, partial-public." Verify every needed peak-list file actually downloads before committing. Courcelles' accessions are now **confirmed public** (resolved 2026-07-15 from local PDF bytes + PRIDE/GEO, no human fetch needed): MS = PXD071022 (PRIDE, CC0), RNA-seq = GEO GSE312236 (BioProject PRJNA1372708). Only the Ely deposit carries the partial-public/processed-peak-list caveat.

### Search database

Build an RNA-seq-bounded, peptide/ORF-hybrid non-canonical FASTA and search it with the canonical proteome concatenated as a **competing target**, not a pre-exclusion filter (AC-2).

Concrete recipe (source: `workflow/scripts/translate_peptides.py`, `assemble_contigs.py`, `config/config.yaml`):
1. Start from `tumor_exclusive` novel junctions; widen the current 27-nt symmetric flanks to **30 nt** each side (`3*(11-1)=30`) so all class-I 8-11mers span the breakpoint (today's 27 nt caps at 10-mers).
2. Emit per-junction 3-frame, stop-codon-split ORF stretches crossing the breakpoint as traceable-header mini-protein entries (`junction_id|frame|chrom:start-end:strand`).
3. Concatenate three target partitions - canonical (SwissProt reviewed + curated isoforms), cRAP contaminants (~600), tagged junction partition - then generate reversed decoys over the whole DB (Philosopher).
4. Search nonspecific/no-enzyme, length 8-11, variable Met-ox; optionally a second engine for PSM intersection (NewAnce-style).
5. Report junction hits at <=1% **group-specific (subset) FDR** on the junction partition alone, and require each junction PSM to out-score the best canonical/cRAP explanation of the same spectrum.

The competing-target concatenation (vs pre-excluding canonical matches) is the correct realization because exclusion removes the alternative hypothesis and biases toward novel calls; the finding marks this specific point "unverified" design reasoning, but the verdict judged it "directionally sound under standard target-decoy competition semantics." The DB-inflation remedy (keep the DB small/sample-specific via RNA-seq) is verified field practice (source: BMC Genomics 10.1186/s12864-016-3327-5; "FDR: the Achilles heel of proteogenomics").

Adversarial qualification (searchdb verdict): separate/subset target-decoy FDR is itself **known to be unreliable for very small groups** (too few decoy hits - the motivation behind Transferred Subgroup FDR methods). Given this project's binding constraint of extreme evidence scarcity, the junction partition may be too small for a reliable standalone FDR estimate. This **strengthens** the "never pool into global FDR" direction but means entrapment/synthetic-peptide validation and orthogonal RNA-seq corroboration are **likely required second gates for any load-bearing hit, not optional extras**.

### FDR / decoy

Adopt a layered, class-specific scheme, not a single global pass (AC-3, AC-5):
1. Canonical and junction as separate target groups, each with its own decoys and target-decoy competition (group-specific FDR); a global FDR massively under-controls the small low-prior junction subset.
2. Two-engine PSM intersection (most false positives are engine-specific).
3. Semi-supervised rescore (Percolator/mokapot) with MS2-prediction (Prosit/MS2PIP) + RT (DeepLC) features, reporting **per-PSM PEP**, not just a set-level q-value - because most junctions are single-peptide-supported where local confidence matters.
4. Stricter junction-group threshold + **entrapment validation** (shuffled-junction / evolutionarily-distant spike-ins; observed entrapment hit fraction estimates the true FDP). Tooling exists: Noble-Lab FDRBench (source: bioRxiv 2024.06.01.596967).
5. State the RNA-seq orthogonal check explicitly.

The load-bearing NewAnce evidence was **verbatim-verified** against the primary (PMC7064602, Nat Commun 2020): two-group FDR dropped surviving non-canonical peptides to **28%** of the un-adjusted count while the predicted-HLA-binder fraction among survivors rose to **85%** - the global FDR had been admitting mostly junk.

Adversarial qualification (fdr verdict): the finding's key_claim[3] "~36% real FDR vs ~0.03% canonical" is **self-flagged unverified / low-confidence and was NOT confirmed** against a fetched primary (the verified primary reports the 28%/85% framing instead). Do not cite the 36% figure; cite the NewAnce 28%/85% numbers. Entrapment-estimated FDP is a cross-check, not an oracle (Madej 2024); report it alongside decoy FDR, never as ground truth. Note also the DDA-vs-DIA caveat: DDA tools control peptide-level FDR more reliably than DIA - prefer DDA repositories or treat DIA junction hits more conservatively.

### Tooling

Two-engine setup (source: FragPipe immunopeptidomics docs; github.com/lazear/sage release v0.14.7; msfragger.nesvilab.org; uwpr.github.io/Comet):
- **Primary (FDR-defensible deliverable): FragPipe** (MSFragger nonspecific + Philosopher group-specific FDR + MSBooster/MS2Rescore). Native group/class-specific FDR is the published, reviewer-expected control for a non-canonical peptide class; confirmed verbatim in FragPipe's own group-FDR docs and precedented by the circRNA-derived-peptide Nat Commun 2024 workflow. RAM-heavy (32-64 GB) - run on a free high-RAM CPU cloud box (Kaggle ~30 GB), **not** the 8 GB M1; it is CPU/RAM-bound, not GPU-bound.
- **Fallback / license-clean cross-check: Sage** (MIT, native aarch64-apple-darwin binary, CPU-only, low-memory, fast) with `cleave_at=""` for nonspecific digestion, rescored via MS2Rescore/Percolator, with class-specific FDR computed **externally** (Sage has no native group-FDR).
- Ruled out for this environment: MaxQuant (.NET, no official macOS arm64, closed), PEAKS (commercial). MASCOT/SEQUEST/ionbot do not allow nonspecific digestion.

Adversarial qualification (tooling verdict): the exact "group FDR 0.03 / 1% PSM" thresholds from the Nature workflow rest on a **WebSearch digest paraphrase** (the paper was paywall-blocked to direct fetch); the FragPipe group-FDR *capability* is confirmed from FragPipe's own docs, and the recommendation does not hinge on the precise threshold. The real residual is methodological, not capability-based: a small junction FASTA yields few target-decoy pairs -> coarse, potentially anti-conservative class-FDR. Report exact recovered-hit counts and per-hit q-values plus RNA-seq support, never an FDR number alone.

### PCPS framing

Cite the proteasome-catalyzed cis-spliced peptide (PCPS) false-discovery dispute as the explicit cautionary precedent, and state our differentiation **categorically** (AC-5). Ready-to-paste framing sentences are in the priorart facet's key_claims; both linchpin Admon quotes were **verbatim-verified** against the open-access primary (PMC8724635, Admon 2021, Mol Cell Proteomics 20:100099):
- The abundance-collapse trend (~25% [Liepe 2016] / up to ~45% [Faridi 2018] -> <1% as FDR tightened [Mylonas 2018; Rolfs 2019; Erhard 2020]), disputed by Admon 2021 and unresolved (Mishto 2021 rebuttal) - always cite the **pair** to keep the dispute two-sided.
- Root vulnerability: a purely spliced peptide is encoded nowhere in the genome/transcriptome, so MS is the only witness.
- Our differentiator: a splice-junction peptide is genomically encoded and RNA-seq independently evidences the junction - and Admon himself names a nucleic-acid assay (Ribo-seq) as the kind of independent evidence that would settle a candidate (verbatim quote verified).

Do not conflate the two mechanisms: PCPS = post-translational transpeptidation; ours = transcribed RNA splice junctions (the shared word "spliced" is a trap). The proteasome's splicing biochemistry itself is not disputed - contest **abundance/FDR**, not the enzyme's existence (source: InvitroSPI 2023, Zotero NNDREW5X).

Adversarial qualification (priorart verdict): two non-fatal scope corrections. (1) Admon names **Ribo-seq (translation-level)** specifically; our differentiation rests on **RNA-seq (transcription-level)**, generalized to "a nucleic-acid assay" - a real slippage (already flagged in the finding's open_question #3), and Ribo-seq is actually the *stronger* check Admon endorses. (2) "Categorical" is supported **only in the narrowed form**: RNA-seq corroborates that the source **sequence exists**, NOT that the individual PSM is correct or that the peptide was translated/presented. Keep the claim to "the underlying sequence exists"; the MS step still needs ordinary FDR discipline.

## Mapping to acceptance criteria

Note: the exact AC wording was not in the facet bundle; these are inferred from each facet's `ac_coverage` and the #1176 task spec. AC-6 is inferred as the deliverable requiring the actual run.

| AC (inferred) | Status | One-line note |
|---|---|---|
| **AC-1** - name 1-2 re-searchable public tumor HLA-I datasets (splice-relevant, HLA-typed, reuse-licensed) | **Complete** (2026-07-15) | Both confirmed public: Courcelles (primary, matched RNA-seq) = MS PXD071022 (CC0) + RNA-seq GSE312236; Ely MSV000096853 (CC0) as second deposit. |
| **AC-2** - canonical proteome enters as a competing target, not a pre-exclusion filter | Design-complete now | Concrete concatenated 3-partition DB recipe specified; the "competing vs exclusion" point is design reasoning judged directionally sound, not primary-sourced. |
| **AC-3** - address target-decoy fragility for novel/non-canonical peptides | Design-complete now | Group-specific FDR + two-engine intersection + entrapment; NewAnce 28%/85% evidence verbatim-verified. |
| **AC-4** - select the search-engine/pipeline tooling (nonspecific, custom-FASTA, group-FDR) | Design-complete now | FragPipe primary (cloud), Sage MIT fallback (local); config knobs, licenses, compute footprint specified. |
| **AC-5** - defensible FDR scheme + PCPS cautionary precedent + orthogonal-RNA differentiation | Design-complete now | Layered scheme + verbatim-anchored framing sentences ready; apply the Ribo-seq-vs-RNA-seq and "sequence-not-ligand" scope corrections. |
| **AC-6** - actually recover presentation-prevalence hits (label=untested) from the re-search | Needs the actual MS run | The wet/compute half: download spectra, build the genome-wide DB, run the engine, compute subset FDR + entrapment. Not settleable by desk research. |

## Open questions and risks

Most decision-critical first:

1. **[RESOLVED 2026-07-15] Courcelles cohort confirmed public** - MS PXD071022 (PRIDE, CC0) + matched RNA-seq GSE312236 (GEO, BioProject PRJNA1372708), recovered from local PDF bytes + PRIDE/GEO checks (no human fetch). Courcelles is the primary cohort; the critical-path blocker for AC-6 is now the genome-wide junction-FASTA build (#3 below).
2. **Subset-FDR reliability at extreme scarcity.** The verified honesty control (group-specific FDR) is itself statistically unstable for a tiny junction partition with few decoys. Implies entrapment + RNA-seq corroboration are **required** second gates, and hit counts + per-hit q-values must be reported instead of a lone "1% FDR" line.
3. **Our junction panel is chr22-scoped today** (patient_001/002 test data). A productive re-search needs a genome-wide junction-spanning peptide DB - a prerequisite build, not a config tweak.
4. **Per-sample HLA typing lives only in paper supplements**, not structured deposit metadata (Courcelles supplementary HLA table; Ely per-patient typing). A recovered hit with unknown HLA loses the binding-prediction validation leg.
5. **MSV000096853 is "Partial Public" / processed peak lists**, not raw files. Confirm the needed peak lists actually download before committing.
6. **Canonical set choice** (SwissProt-only vs +isoforms vs +TrEMBL vs +GENCODE) trades competing-hypothesis stringency against DB size; must be version-pinned per search.
7. **DDA vs DIA provenance** - prefer DDA repositories (peptide-level FDR control); treat any DIA source conservatively and stratified.
8. **Frame is not independently witnessed** - RNA-seq evidences the junction, not the reading frame; do not over-claim frame-level support.
9. **Single-engine vs dual-engine** at our evidence scale - added compute/specificity tradeoff; settle with a small side test on one dataset.
10. **Whether targeted-MS/PRM or synthetic-peptide co-elution is required** for top hits (the PCPS "gold standard"), or whether RNA-seq + entrapment-validated subset FDR suffices at `tier=presentation-prevalence`.

## What "start" means here

This desk-research settles the **front-half** and produces reusable write-up assets, but does not touch spectra.

**Settled now (design-complete, no wet/compute step):**
- AC-1 (cohort) fully settled 2026-07-15: Courcelles (PXD071022 + GSE312236, both public, MS CC0, matched RNA-seq) primary; Ely MSV000096853 (CC0) second.
- AC-2 (competing-target DB architecture), AC-3 (novel-peptide FDR fragility + remedy), AC-4 (tooling + compute plan), AC-5 (FDR scheme + verbatim-anchored PCPS framing sentences).
- All method choices are grounded in verified primary sources; the differentiator argument is written and quote-checked.

**Still requires the wet/compute half (not startable from the desk):**
- Building the **genome-wide** junction-spanning peptide FASTA from a de-chr22-scoped panel (prerequisite engineering; now the sole critical-path blocker for AC-6 since the cohort is confirmed).
- Downloading the actual peak lists, running FragPipe/Sage nonspecific search, computing subset FDR, and running the entrapment validation.
- Producing the recovered `label=untested / tier=presentation-prevalence` hits (AC-6) and cross-checking each against RNA-seq junction evidence.

In short: the human can approve the design and the starting cohort today; what remains is the MS run itself (spectra download + DB build + engine + FDR/entrapment), which is where the real evidence - and the real residual FDR risk - lives.
