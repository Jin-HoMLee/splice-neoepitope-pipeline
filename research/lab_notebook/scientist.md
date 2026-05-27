# Lab Notebook — Scientist

Per-role lab notebook for Scientist sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-27

### 13:55 UTC — Editor: Scientist

#### Phase 1 — Morning briefing → 2 Zotero adds (HCC AS off-the-shelf mRNA vaccine + Prime-Target neoantigen vaccination)

Two papers surfaced from current-year search across cancer neoepitope / splice immunology / vaccinology; neither already in Zotero collection `Z38GTJNW`.

**[Zotero B2MJ776X — Zhao et al., *Cell Research* 35:970–986, 2025-11-28](https://www.nature.com/articles/s41422-025-01199-0)** — *Harnessing alternative splicing for off-the-shelf mRNA neoantigen vaccines in hepatocellular carcinoma*. Reframes AS-derived neoantigens as **shared** rather than patient-specific: AS events occur **>59× more frequently** than somatic mutations across the HCC cohort, yielding 50.94% population coverage vs 4.40% for mutation-derived neoantigens. Stringent filter reduces to **34 prioritized AS neoantigens**; mouse-model proof-of-concept mRNA vaccine encoding the panel produces tumor regression + neoantigen-reactive TIL infiltration. **vs. our pipeline:** same source biology (splice-junction-spanning translated peptides), opposite end of the public-vs-personalized axis (cf. Kwok et al. 2025 / [Sub-Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) DISCUSSION subsection). Strong manuscript anchor for the **population-coverage** framing alongside SNAF / SpliceMutr / Tumour-wide RNA splicing aberrations references.

**[Zotero XNH627GQ — Prime-Target preprint, bioRxiv 2026-01-13](https://www.biorxiv.org/content/10.64898/2026.01.13.699214v1)** — *Prime-Target neoantigen vaccination unleashes unprecedented T cell immunity within "cold" immunosuppressive tumors*. Heterologous prime-boost regimen: mRNA prime + peptide target boost of the same antigen produces stronger systemic T cell immunity inside "cold" tumors than either modality alone in mouse models. **vs. our pipeline:** delivery-strategy paper, not a prediction tool — no overlap with junction-calling stack. Manuscript DISCUSSION hook on **downstream use** of predicted AS-neoepitopes; contrast point against single-modality regimens (autogene cevumeran, GNOS-PV01).

#### [Issue #503](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/503) — glossary additions: Prime-Target / heterologous prime-boost terminology closed via [PR #504](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/504)

Three new entries in [research/glossary.md](research/glossary.md) for terms surfaced by the Prime-Target paper that will land in the manuscript DISCUSSION on downstream use of predicted AS-neoepitopes:

- **Prime-Target** (P section, between PoN and PSI) — the regimen itself
- **mRNA prime + peptide-target boost** (M section, between MHC1-TIP and MSI) — the delivery sequence inside it
- **single-modality regimen** (S section, between SGE and SLURM) — the contrast term covering mRNA-only / peptide-only / DNA-only

Cross-references between entries follow existing glossary convention (plain-prose term references, not bracketed shorthand — the original AC mentioning `[Prime-Target]` was reworded to match the file's actual style). +6 lines total (3 entries × 2 lines each including blank separator). Bot review via `@claude review` returned **LGTM, no issues found**; CI all green. `bash scripts/audit_and_merge.sh 504` shipped after Issue #503 backfilled with a **Priority rationale** line (script gate caught its absence on first invocation).

#### [Issue #478](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/478) — RESULTS.md stale `junctions.tsv` path reference closed via [PR #505](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/505)

Single-line catch-up to [PR #477](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/477)'s rename (`alignment/<sample>/junctions.tsv` → `alignment/<sample>/raw_junctions.tsv`) at [research/manuscript/RESULTS.md:54](research/manuscript/RESULTS.md#L54). PR #477 deliberately carved this out as Scientist territory and filed Issue #478 as the follow-up.

**Scope verification before merge** (in response to user pushback "shouldn't we check the other manuscript files too?"): `grep -rEn "junctions\.tsv|junctions\.bed" research/manuscript/` across all 6 manuscript files (CONCLUSIONS / DISCUSSIONS / INTRODUCTION / METHODS / REFERENCES / RESULTS) returned exactly 4 hits — the 1 fixed in this PR plus 3 surviving `novel_junctions.tsv` references in METHODS.md:154/261 and DISCUSSIONS.md:138, all on patient-level paths (`results/junctions/{patient_id}/...`). `novel_junctions.tsv` is a **different file** (the downstream-labeled tumor-exclusive output, not the raw per-sample HISAT2/regtools extraction) and was intentionally untouched by [PR #477](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/477). `\bjunctions.tsv` (word-boundary) returned zero hits — confirming no bare stale references remain. INTRODUCTION.md, CONCLUSIONS.md, REFERENCES.md have zero references of any kind.

Bot review via `@claude review` returned **LGTM, no issues found** (landed 2 minutes before my Monitor armed — `since=$last` filter missed it; the post-mortem on Monitor cold-start is a minor follow-up). `bash scripts/audit_and_merge.sh 505` shipped after Issue #478 ACs were ticked.

#### Phase 2 — Standup hygiene: replied to PM [2026-05-26 12:05 UTC] re-raise + archived own [2026-05-21 10:24 UTC] Sci post

Two halves of the standup hygiene run from the morning routine:

- **Reply posted** to PM's [2026-05-26 12:05 UTC] re-raise at [2026-05-27 12:38 UTC] — confirmed all 5 items from PM's [2026-05-22 08:21 UTC] substantive reply caught up (ImmSET added to [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) verdict list, sizing nudge ack'd, due_on curiosity investigated, [Issue #433](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/433) Dev-handoff flag taken, parent-Status mechanism staying in PM court). Invited PM to flip both [2026-05-22] and [2026-05-26] posts to Done.
- **Own [2026-05-21 10:24 UTC] post archived** to `team_standup_archive/2026-05.md` (chronological append, sender-owned); `_index.md` count bumped 94 → 95. Post had been Done since 2026-05-22 (≥3 days) per the standup hygiene rule.

#### Reflection — closure-audit bot caught the lab-notebook gap this session

Both [Issue #503](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/503) and [Issue #478](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/478) closures triggered `closure-audit bot` comments flagging the missing `## 2026-05-27` block in this file. Individually each Issue close is borderline routine (single-PR-closes-single-Issue), but the session as a whole crosses both Issues plus Phase 1 Zotero adds plus Phase 2 standup hygiene — that bundle of cross-Issue work is **the non-routine criterion** in [shared/feedback_lab_notebook.md](.claude/memory/shared/feedback_lab_notebook.md). The closure-audit bot's heuristic is a fair backstop; this entry backfills via [Issue #511](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/511). Next time, the entry should ride the *first* closing PR of a non-routine session rather than chase a backfill PR after the fact.

---

## 2026-05-26

### 11:55 UTC — Editor: Scientist

#### [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) — HERMES evaluation closed as (a) integrate; sub-issue [#492](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/492) filed under milestone 29

[HERMES — Visani et al., *PNAS* 2025-10-21 (DOI 10.1073/pnas.2504783122)](https://doi.org/10.1073/pnas.2504783122). 3D-equivariant ML model pre-trained on the protein universe (~10k CASP12 chains); predicts amino-acid propensities at sites from local 3D structural environments. **Zero-shot on TCR-pMHC** (no domain-specific training data) yet achieves 0.72 correlation with experimental binding affinities + up to 50% T cell activation rate for de novo peptide designs across 3 TCR-MHC systems. Feb 2026 follow-up ([bioRxiv 707744](https://www.biorxiv.org/content/10.64898/2026.02.24.707744v1.full)) applies HERMES to AF-predicted TCR-pMHC structures and discriminates interacting from non-interacting complexes — exactly our use case (TCRdock outputs AF-Multimer PDBs).

**Decision: (a) integrate as post-TCRdock structure-based confidence cross-check.** Three reasons:

1. **No integration blockers.** Open-source ([StatPhysBio/hermes](https://github.com/StatPhysBio/hermes), MIT-licensed), pretrained weights bundled + on Zenodo, Python 3.9 + PyTorch + e3nn==0.5.0 — drops into our `workflow/envs/` pattern cleanly. CLI: `python run_hermes_on_pdbfiles.py -pd pdbs -m hermes_bp_000 -o results.csv`.
2. **Orthogonal signal to AF confidence.** TCRdock outputs AF-Multimer ipTM + pLDDT + PAE; HERMES outputs per-site amino-acid propensities from the 3D environment. The Feb 2026 follow-up demonstrates HERMES discriminates interacting from non-interacting predicted TCR-pMHC complexes — concrete evidence the signal complements AF's own confidence on the same input.
3. **Niche complementarity with NetTCR-struc** ([Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) / [Issue #433](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/433)). NetTCR-struc predicts DockQ (docking quality); HERMES predicts amino-acid propensity at the interface (binding affinity proxy). Together they cover the two structural-QC dimensions for TCR-pMHC.

**Modeling sub-issue filed: [Issue #492](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/492)** under [`milestone 29 (i5 - S5 - Modeling - TCR-pMHC Scorer Integration)`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/29). Title: `feat(scoring): integrate HERMES as post-TCRdock structure-based confidence cross-check`. `role:developer`, P2. Gated by TCRdock end-to-end pipeline functional (same prerequisite as [Issue #433](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/433) NetTCR-struc sibling).

**Manuscript carry-forward to [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432):** HERMES + NetTCR-struc as a complementary structural-QC pair (binding-affinity + docking-quality dimensions). Concrete (a)-integrate verdict to pair with today's ImmSET (b)-decline ([Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201)) — balanced manuscript example.

**Verdict progress on [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) per-scorer AC:** 2 of 5 scorers resolved today (ImmSET → b, HERMES → a). Remaining: Boltz-2 ([Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188)), t2pmhc/TCRLens ([Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236)), [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316) outcome.

---

### 11:48 UTC — Editor: Scientist

#### [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) — ImmSET evaluation closed as (c) not-relevant — (b) decline comment posted; manuscript-citation kept

[ImmSET — Garcia Noceda et al., arXiv 2603.26994 (2026-03-27)](https://arxiv.org/abs/2603.26994), Adaptive Biotechnologies. Transformer-based "Immune Synapse Encoding" predictor for TCR-pMHC specificity; outperforms AF2/AF3-based pipelines on A\*02:01 with sufficient training data.

**Decision: (c) not relevant** for our pipeline integration → posted (b) decline-with-rationale comment on [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) per [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432)'s per-scorer verdict AC. Three hard blockers:

1. **Code + model weights not released** — paper explicitly states *"The code from this study will not be made available."* No inference path exists.
2. **Training data is proprietary** — Adaptive Biotechnologies' MIRA + pairSEQ. No VDJdb / IEDB / McPAS path; reproduction requires Adaptive's data.
3. **OOD weakness on novel peptides + non-A\*02:01 alleles** — IMMREP25 B\*40:01 per-peptide AUROC range 0.336–0.832; scaling exponents α=0.078 (peptides) / β=0.036 (TCRs/peptide) sublinear. Splice-junction neoepitopes are by construction OOD for any sequence-based model trained on public TCR-pMHC pairs.

**Carry-forward to [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) manuscript subsection:**

- **Shortcut-learning critique.** ImmSET authors identify TCR-motif memorization as a failure mode that inflated prior sequence-based methods' reported performance; mitigated by enforcing Levenshtein-4+ peptide separation between train and test. Citable quality-bar consideration for any future sequence-based TCR-pMHC scorer.
- **Sequence-vs-structure data point.** Sequence wins on well-trained alleles + sufficient data; structural maintains advantage on undertrained alleles + when interpretability matters. Concrete head-to-head reference for the DISCUSSION subsection.

**Zotero:** D6P6JXSH in collection Z38GTJNW (preprint type, manuscript-citation-candidate tag, three-section note attached).

**Verdict progress on [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) per-scorer AC:** 1 of 5 scorers resolved (ImmSET → b). Remaining: HERMES ([Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218)), Boltz-2 ([Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188)), t2pmhc/TCRLens ([Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236)), [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316) outcome.

---

## 2026-05-25

### 19:24 UTC — Editor: Scientist

#### [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439) — patient_002 WGS-vs-WES doc correction note landed on [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) body

Closeout for [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439) (P4 doc fix; [PR #434](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/434) audit follow-up) via path 1: top-of-body correction note prepended to [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) (patient_002 onboarding, closed 2026-04-21). The misleading "Normal Sample — Blood WGS (DNA only)" section on [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) described a 115 GB Jun 2024 UCLA file that was never onboarded; what actually got integrated (via [Issue #67](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/67) B2 migration) was a Dec 2022 BostonGene WES normal (~9.9 GB). Root cause: two physically separate files on osteosarc.com/data conflated by the original [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) author.

**Why path 1 over path 2.** [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439) offered an alternative: don't edit [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37), rely on canonical metadata files (`config/samples/patient_002.tsv`, RESULTS.md, DISCUSSIONS.md) which already say WES correctly. Rejected because future readers landing on [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) directly (e.g. via cross-links in onboarding docs) would not encounter the correction at all. Top-of-body note co-locates the correction with the misleading text.

**Provenance pointers in the correction note.** [Issue #67](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/67) (B2 migration — where the WES file landed), [PR #434](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/434) (VCF audit — full investigation), [Issue #277](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/277) (WES → PBMC scRNA-seq switch — current state; WES no longer used), and [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439) (this follow-up).

**Closure status.** [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439) AC ticked (2 of 3 boxes; AC #3 — note-only path — removed as obsolete per closure-ritual "remove if obsolete"). Closeout comment posted on [Issue #439](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/439); final `gh issue close` deferred until this lab notebook PR lands so the closure-audit bot sees both events in the same merge-day date block.

---

### 18:47 UTC — Editor: Scientist

#### [PR #472](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/472) ([Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271)) merged — AlphaGenome DISCUSSION shipped; closure-audit follow-up filed as [Issue #473](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/473)

[PR #472](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/472) (AlphaGenome predicted-normal filter NO-GO DISCUSSION subsection) merged at 18:46:59 UTC, squashed as [`70a2b9a`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/70a2b9a) on `main`; [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) auto-closed at 18:47:00 UTC. @-claude review pass returned an approve verdict with 4 pre-submission style notes; 1 applied in-cycle ([`4a8c71d`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/4a8c71d), "invert" → "AlphaGenome-unique junctions … emerge" precision fix), 3 deferred to the pre-submission proofread window alongside REFERENCES.md Open-items #2 and #3.

**Closure-audit follow-up.** The closure-audit bot flagged that the 18:16 UTC entry below doesn't reference [PR #472](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/472) — root cause: the 18:16 entry was written before the PR was opened, and lab notebook entries are immutable once committed. This 18:47 entry closes the gap forward (not by edit). Filed [Issue #473](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/473) as the tracking carrier so the follow-up branch could use `gh issue develop` per the ALWAYS-use rule, rather than `git checkout -b`.

**Generalisable lesson.** When a lab notebook entry is committed *before* the corresponding PR opens, the closure-audit bot (which scans the merge-day date block for the PR number) will flag a gap even though the audit-trail intent is preserved across the entry + commit history. Two patterns avoid this: (i) defer the lab notebook entry until *after* the PR is opened so the PR number is in scope from the first commit (preferred when the work is small enough that the PR is imminent), or (ii) write the entry first and add a post-merge supplementary entry referencing the PR (this PR's path). Not codifying as a feedback memory until a second incident calibrates the preference.

---

### 18:16 UTC — Editor: Scientist

#### [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) — AlphaGenome DISCUSSION subsection drafted and inserted into `research/manuscript/DISCUSSIONS.md`

Drafted the 5-point AlphaGenome DISCUSSION subsection per the concrete scope appended to [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) body in the 12:22 UTC close-out. Placed as a `###` under "Normal sample filtering: junction level vs. peptide level" (line ~206), structurally parallel to the existing `### GTEx pan-tissue filter` subsection — both are "extensions of normal filtering," one population-level (works), one sequence-based (doesn't). Title: "AlphaGenome as a predicted-normal filter: foundation-model evaluation (NO-GO)". Five paragraphs cover: (1) three-experiment validation outcome with verbatim Exp 1 F1 = 0.300 / Exp 3 catch counts MN 91 / GTEx 483 / AG 124 of 1,872 / union 503 / AG-unique-vs-GTEx 0.0% / NO-GO over-determined; (2) tissue-prior-vs-patient-specific-predictor framing, generalised as a foundation-model failure mode (model trained on reference tissue without per-individual signal collapses to a smoothed tissue prior); (3) two-axis production filter stack (MN + GTEx) — AG dropped, no third axis to add; (4) three niche-use angles preserved as principle-only (confidence proxy for GTEx hits via intersection scoring; sensitivity-tuned alternative; full-genome robustness check) so the manuscript doesn't read as "AG is useless"; (5) one-sentence cross-link to splice2neo (Lang et al., *Bioinform Adv* 2024) as a convergent design choice from an independent tool.

**Number verification.** All decision-rule numbers cross-checked against [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body Exp 3 Result subsection verbatim: F1=0.300 at τ=3.16 (P=0.238, R=0.405) ✓; n=1,872 tumor ✓; MN 91 (4.9%), GTEx 483 (25.8%), AG 124 (6.6%), union 503 (26.9%) ✓; AG-unique-vs-GTEx 0.0% ✓; NO-GO clauses (F1 < 0.5 + 5% fallback failure) ✓; Sub-Issue #381 deferral noted, not load-bearing ✓.

**REFERENCES.md updates.** Cross-reference summary table line for AlphaGenome (`UZWZ5QEB`) extended from "METHODS.md §3" to "METHODS.md §3; DISCUSSIONS.md". Avsec et al. 2026 entry's "Cited in" appendix updated to add the DISCUSSIONS context (predicted-normal NO-GO per [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)). Open-items item 4 (AlphaGenome validation strategy entry pending [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271)) struck-through and resolved: no new in-text citations beyond the existing Avsec 2026 + Lang 2024 entries are required, since splice2neo already has a manuscript citation in METHODS §2-3 and is being cross-linked rather than newly introduced.

**Cross-link policy.** Used `[Sub-Issue #381](url)` (Exp 2 deferred) and `[Issue #211](url)` (full-genome robustness check) inline in the DISCUSSION prose — following the existing precedent in DISCUSSIONS.md line 217 (`(issue #17)` for HISAT2 vs STAR planned comparison). These are forward-pointing operational links to tracked follow-up work; deletion at the pre-submission proofread is straightforward if the journal house style rejects GitHub references.

**Next:** tick [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) AC, commit + PR (separate steps per the chain rule).

---

### 12:22 UTC — Editor: Scientist

#### [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (AG predicted-normal filter) — epic CLOSED; close-out carrier [Issue #470](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/470); manuscript scope into [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271)

Post-merge close-out for [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452). Sequence: ticked [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Sub-issues A ([Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) carrier), B ([Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) Exp 1 + [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) Exp 2 carve-out), C ([Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) Exp 3 chr22); appended Exp 3 Result subsection into [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body; updated [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) body with 5-subsection DISCUSSION scope (unblock note + concrete subsection breakdown); filed [Issue #470](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/470) as the decision-memo + close-out carrier; ticked the last [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) task with [Issue #470](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/470) + [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) carrier refs; posted summary comment + closed [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) at 12:22 UTC.

**Decision recap.** AG NOT adopted as additional normal-junction filter — over-determined NO-GO (F1=0.300 trips `<0.5` clause; 0% AG-unique-vs-GTEx fails fallback tier's `≥5%`; F1=0.300 is 2.7× below strongest tier's `≥0.8`). [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) (Exp 2, WGS-blocked) is not load-bearing and remains independently open as a sanity check. Production filter stack stays two-source: MN (patient-specific, when available) + GTEx pan-tissue (always-on via [Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126) / [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) / [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)). AG dropped from `workflow/scripts/filter_junctions.py` integration plan.

**Niche-uses nuance preserved.** User raised the framing "AG would be only interesting if we need a filter less restrictive than GTEx" — correct, with caveats. Three potential niche uses NOT ruled out (folded into [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) DISCUSSION subsection 4 so the manuscript doesn't read as "AG is useless"): (i) confidence proxy for GTEx hits via intersection scoring, (ii) sensitivity-tuned alternative if downstream candidate counts bottleneck, (iii) full-genome scale-up robustness check on the `AG ⊂ GTEx` subset relationship.

**Process note — decision-memo location.** Considered creating a new `research/decisions/YYYY-MM-DD-<short>.md` convention dir for go/no-go scientific decisions; chose instead to put the memo inline in [Issue #470](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/470) body. Rationale: the verdict is already captured in 4 durable artifacts ([Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Exp 3 Result subsection, [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) body, [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) merged commit, this lab notebook entry); adding a 5th in a new convention dir would be premature design without a 2nd consumer. Issue body is citable, dated, GitHub-discoverable, and durable. The new convention can be promoted later if a 2nd go/no-go decision warrants it (mirroring the `_shared/` lazy-promotion rule for cross-experiment data in CLAUDE.md).

**Outstanding follow-ups (not blocking close-out):** [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) Exp 2 (WGS-blocked sanity check); [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) production GTEx panel → full-genome AG re-run as robustness check; [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) DISCUSSION prose writing (now unblocked).

---

### 11:43 UTC — Editor: Scientist

#### [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) — second @-claude review pass; 3 minor items addressed in `bf87d96`; pre-merge audit clean

Re-requested @-claude review on [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) because 5 commits had landed after the first review pass (2026-05-21 18:31 UTC): the entire slide-deck scope shipped 2026-05-22 wasn't yet reviewed. Second pass returned 3 minor items, all addressed in commit [`bf87d96`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/bf87d96):

1. **`workflow/envs/alphagenome.yaml` — `scikit-learn` undeclared as direct dep.** Both `notebook.ipynb` §2(b) and `figures/_regenerate_figures.py` line 27 import `sklearn.metrics`; the package was only available transitively via the `alphagenome` pip closure. Added `scikit-learn >=1.3` to the conda-forge deps so a future `alphagenome` version drop won't silently break either consumer.
2. **`figures/_regenerate_figures.py:155` — print label mismatch.** `f"Outputs dir: {FIGURES_DIR.relative_to(REPO_ROOT)}"` printed the `figures/` path under the wrong label (FIGURES_DIR is not OUTPUTS_DIR). Relabelled to "Figures dir:".
3. **README Outputs section — slide deck + deck-only figures missing.** The 3 notebook artifacts (chr22_gtex_panel.parquet, filter_overlap_table.tsv, filter_venn_chr22.png) were listed but `slides.qmd`/`slides.html` and `figures/pr_curve.png`/`caught_bar.png` were not — discoverability gap for a future reader using the README as a deliverables map. Added 2 lines.

All 3 items were independently verified against the actual code before applying (reviewer claim → grep/Read check → apply); no items pushed back. CI re-ran green on `bf87d96` at 11:39:11 UTC (3/3 checks: pipeline-snakemake-dry-run, pipeline-pytest, ci-tools-pytest). Pre-merge audit clean: [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) Test plan + [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) AC both 0 unticked. Mergeable: MERGEABLE / CLEAN. Closure-ritual gate predicted to pass.

**Process lesson — declare transitive deps explicitly.** Item 1 (scikit-learn) is a generalisable pattern: when a notebook or analysis script imports from a package that's only present transitively via another env-yaml entry, the dep is invisible to `conda list` and silently breaks on upstream version drops. Rule: every package directly imported by a notebook or script in the experiment dir should appear as an explicit entry in `workflow/envs/*.yaml`. Future-applicable when adding any new analysis dep.

**Pending user OK to merge** via `bash scripts/audit_and_merge.sh 452 --squash --delete-branch`. Post-merge: update [parent Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body Exp 3 row with the decision-rule outcome (NO-GO; F1=0.300; % AG-unique vs GTEx = 0.0%) — explicit post-merge step flagged in the [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) body.

---

## 2026-05-22

### 13:08 UTC — Editor: Scientist

#### [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (filter strength deck) — slide deck added to [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) + slides-co-location convention documented

User asked for the [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) slide deck to ship in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) (vs deferring to [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455)) on the grounds that notebook + deck "belong together and are necessary for human review". Agreed; the slides-co-location convention agreed informally in the 2026-05-21 19:38 UTC entry (below) is now formalised in CLAUDE.md as part of this PR.

**Scope shipped here:**

- `research/experiments/issue_225_*/slides.qmd` (12 slides, lab-seminar-quality, reveal.js HTML render)
- `figures/_regenerate_figures.py` + `figures/pr_curve.png` + `figures/caught_bar.png`
- `outputs/filter_venn_chr22.png` re-used from the notebook (notebook stays canonical per CLAUDE.md)
- `refs.bib` — 7 entries (5 reused from issue_393 deck, 2 new: Wilks-Snaptron + GTEx v8)
- `slides.html` committed for [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452)-review ergonomics (divergence from `research/slides/README.md`'s gitignore rule — tracked as a follow-up under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455))
- CLAUDE.md "Slide decks for experiment Issues" section updated with the co-location rule + rationale; also fixed a stale inverted cross-reference in the "Experiment notebooks" section
- [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body descoped: removed the CLAUDE.md convention-docs item and the "[Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) deck follow-up" Out-of-scope line (both shipped here)
- [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) Test plan amended with 6 deck-related items

**Process note.** Followed the brainstorming → writing-plans → subagent-driven-development superpowers flow. Spec at `docs/superpowers/specs/2026-05-22-issue-225-slide-deck-design.md`; plan at `docs/superpowers/plans/2026-05-22-issue-225-slide-deck.md`. The plan structure paid off — small bite-sized tasks gave clean per-task commits, and the spec's "Divergence from prior art" section flagged the committed-`slides.html` decision up front instead of discovering it mid-render. One in-flight plan correction: the Snaptron DOI initially provided (`10.1093/bioinformatics/bty025`) was wrong — caught at Task 1, corrected to `10.1093/bioinformatics/btx547` (PMID 28968689); plan + refs.bib updated accordingly.

**Headline numbers unchanged from notebook re-run on 2026-05-21:** F1 = 0.300, % AG-unique vs GTEx = 0.0%, decision = NO-GO. The deck just re-renders the same story in lab-seminar form.

**Ready to merge after this commit lands.** Closure-ritual gate (test plan + [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) ACs + [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) carrier step) will re-check before invoking `bash scripts/audit_and_merge.sh 452`.

---

## 2026-05-21

### 19:38 UTC — Editor: Scientist

#### [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) — @claude review iteration; slides-co-location convention agreed; migration carrier [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) filed

@claude review on [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) returned 6 items; all addressed in this turn before the merge gate. None blocking, two recommended ("README status" + "inline GTEx coord assertion"), four polish/defensive. Listed below verbatim by reviewer ordering so future readers can cross-reference:

1. **README status field** — flipped from `in progress` to `complete — NO-GO verdict (2026-05-21)`. Reviewer's point: first thing a future reader checks when this becomes a frozen reference.
2. **Inline GTEx coord-convention assertion** — added a notebook-internal `259/259 ground_truth introns present ✓` guard at the end of §2(c). The lab notebook (18:14 UTC entry) documented this empirically but the notebook itself only had non-empty + chr22-only asserts. Now self-contained.
3. **`snaptron_to_key_set` simplification** — collapsed the row-loop to a set comprehension, dropping the redundant `int()` casts inside the 880K-row hot loop. Functional behaviour unchanged.
4. **`gencode_introns_chr22` inverted-intron guard** — added `if donor >= acceptor: continue`. GENCODE v47 doesn't have this pathology in practice (intron count unchanged at 7,731 after the guard) but the function is now self-documenting + robust to non-GENCODE annotation swaps.
5. **Cross-experiment path coupling** — already explicitly listed under #455's "Cross-experiment dependency fix-ups" section before the review arrived. No additional action; tracked.
6. **CLAUDE.md slide deck section — stale figure-source path** — pre-existing reference to `research/notebooks/<exp>_outputs/*.parquet` updated to `research/experiments/issue_NNN_<short>/outputs/*.parquet` matching the new convention introduced in this PR.

All numbers unchanged after re-execution: F1=0.300, % AG-unique vs GTEx = 0.0%, decision = NO-GO. New `259/259 ✓` coord-validation message added to §2(c) output; intron count + ground-truth count both stable.

**Slides-co-location convention agreed mid-session.** Original CLAUDE.md draft kept `research/slides/issue_NNN_<short>/slides.qmd` as a parallel top-level folder. Discussion concluded co-located form (`research/experiments/issue_NNN_<short>/slides.qmd`) is preferable: notebook + outputs + deck rename / archive / migrate as a unit, deck figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`, fewer broken-link footguns. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/`. Convention documentation deferred to #455 (so the same migration PR can move both notebooks AND slides + update CLAUDE.md in one shape).

**Migration carrier filed: [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).** Bundled scope: notebook migrations (`issue_224_*` from `research/notebooks/`, `issue_299_*` from same) + slide migrations (`issue_393_*` from `research/slides/`) + cross-experiment path fix-ups in [#225's notebook](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/main/research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb) once #224 moves + CLAUDE.md update to document the slides-co-location rule. Per-patient stable notebooks (`patient_001_results.ipynb`, `patient_002_results.ipynb`) explicitly out of scope — they're long-lived manuscript-supporting, not per-Issue experiment work.

**Process lesson.** Reviewer caught item 2 (inline coord assertion missing) precisely because I had documented the validation in the **lab notebook** but not in the **notebook itself**. The lab notebook is for narrative + rationale; the notebook is the canonical artifact that must be self-verifying when re-executed. Reverse-direction rule: whenever a coord-convention or schema-shape gets validated empirically during writing, the validation belongs as an `assert` in the notebook, not as prose in the lab notebook. Future-applicable when adding any external-data-source loader.

**Ready to merge after this commit lands.** Closure ritual re-check pending; will confirm before invoking `bash scripts/audit_and_merge.sh 452`.

---

### 18:14 UTC — Editor: Scientist

#### [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (research: normal-junction filter strength on patient_001 chr22) — Exp 3 ran end-to-end; verdict **NO-GO** for AG-as-3rd-filter

Picked up [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (sub-issue of parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Experiment 3) and ran the three-way filter-strength comparison on patient_001's chr22 tumor junctions. New notebook at `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`, established under a new `research/experiments/issue_NNN_<short>/` convention documented into CLAUDE.md as part of the PR.

**Headline numbers.** On chr22 (test-config harness; tumor = SRR9143066, matched-normal = SRR9143065, AG predictions = #224's cached parquet @ F1-max τ, GTEx = Snaptron hg38 GTEx v2 endpoint ≥1 sample):

- **Exp 1 F1 (universe-restricted, MN ∩ GENCODE positives = 259 / universe = 7,731 GENCODE introns):** 0.3000 at τ = 3.16 (P=0.238, R=0.405). Matches the universe shape #224 §5 reported exactly (positives=259, negatives=7,472, AG-scored=5,728 / 74.1%).
- **Tumor (n=1,872) caught by each filter:** MN 91 (4.9%), GTEx 483 (25.8%), AG 124 (6.6%). Caught-by-any (union) = 503 (26.9%).
- **% AG-unique vs GTEx: 0.0%.** Every chr22 tumor junction AG catches is also caught by GTEx — AG is fully subsumed at the F1-max threshold.

**Decision-rule outcome for [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Exp 3 row: NO-GO — treat as tissue prior.** F1=0.30 trips the <0.5 NO-GO clause directly; even at a higher F1, the 0% AG-unique-vs-GTEx number would block the "fallback" tier (which requires ≥5% unique). Exp 2 (germline-aware AG) deferred to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS — but Exp 1's signal alone is conclusive for the no-go branch.

**Process: plan deviation in §2(b).** The implementation plan called for a per-unique-threshold set-op loop. On 770K unique scores across 2.6M predictions that's intractable (hours of O(N²) work). Switched to `sklearn.metrics.precision_recall_curve` over universe-restricted (universe = GENCODE chr22 introns) score/label vectors — same math (TP / FP / FN computed against the same universe), O(N log N), few seconds. Documented inline in the cell + in the §7 caveats. The universe choice (rather than the plan's unrestricted ground_truth) is what makes the F1 number comparable to #224's reported F1 and to the τ-thresholds in the #203 decision rule (which were calibrated against universe-restricted F1).

**Snaptron coord-convention validation.** Snaptron GTEx v2 `start` is 1-based inclusive intron donor; `end` is 1-based inclusive acceptor (= 0-based exclusive). Validated empirically: 259/259 (100%) of ground_truth (MN ∩ GENCODE) introns appear in the Snaptron panel under the `start-1, end` normalisation; both off-by-one alternatives give 0% overlap. Cached as `outputs/chr22_gtex_panel.parquet` (880,769 rows, 6.3 MB — fits the <10 MB checked-in band per the new size-guidance rule in CLAUDE.md).

**Convention established.** Per the design spec, this PR introduces `research/experiments/issue_NNN_<short>/` (with `README.md`, `notebook.ipynb`, `outputs/`) as the per-Issue notebook convention — mirrors the `research/slides/issue_NNN/` deck convention. Cross-experiment data sharing rules (default-own / `_shared/` lazy-promotion at 2nd consumer / explicit path reference for single consumer / production `resources/` promotion) and size bands (<10 MB / 10–100 MB + regenerator / >100 MB GCS manifest) all documented into CLAUDE.md. Migration of existing per-Issue work (`issue_224_*` / `issue_299_*`) from `research/notebooks/` punted to a follow-up Issue filed post-merge.

**Follow-ups:**
- Update parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body Exp 3 row with these numbers after this PR merges (Task 20 of the implementation plan).
- File a migration Issue: move `research/notebooks/issue_224_*` and `issue_299_*` under the new `research/experiments/` convention.
- Re-run §2(c) against the production GTEx panel once [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) lands — the Snaptron proxy is good enough for the no-go decision but the production panel will tighten the % AG-unique number.
- Optional: persist `best_threshold` from #224's notebook so #225 doesn't have to recompute it. Small QoL refactor; not blocking.

**Outputs:** `outputs/chr22_gtex_panel.parquet` (6.3 MB), `outputs/filter_overlap_table.tsv` (248 B), `outputs/filter_venn_chr22.png` (60 KB).

---

## 2026-05-20

### 21:29 UTC — Editor: Scientist

#### [Issue #384](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/384) (research(data): audit somatic VCF availability) — full cohort audit; verdict **Implement** for [Issue #416](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/416) (gated on [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413))

Picked up [Issue #384](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/384) from Backlog after the 14:51 UTC parent-Status hygiene work cleared the board. PM's skeletal scope ("Sci to fill in") explicitly left the technical scope open — what counts as "available", whether RNA-seq-derived germline calls qualify, which patients to include. The audit-process arc that landed today re-defined "proper audit" against my own initial reflex to infer from samplesheet alone.

**Process: three pushback-driven re-scopes.** First proposed verdict (defer, "no DNA-seq for either patient") rested on `config/samples/patient_001.tsv` and `patient_002.tsv` only — both RNA-seq-only. User pushback "we should do a proper audit" rejected the samplesheet-inference shortcut; second pushback "you didn't do any web search to look for more data?" surfaced that the audit needs external research (GDC cohort donor counts, hybrid-calling literature, vendor-deliverable inventories). Third pushback ("https://osteosarc.com/data/ mentions VCF files, can you check?") cracked the patient_002 picture wide open — the verdict flipped from defer → research-feasible-but-engineering-heavy → fully Implement-tier in three iterations. The 14:51 UTC parent-Status hygiene lesson repeats: **work-completion state lives in places memory doesn't index automatically; surface every relevant data source explicitly before drawing conclusions.**

**patient_001 inventory — full BioProject PRJNA545281 / Moisseev et al. 2020 *Biomedicines* 8(3):67 ([DOI 10.3390/biomedicines8030067](https://doi.org/10.3390/biomedicines8030067), PMID [32210001](https://pubmed.ncbi.nlm.nih.gov/32210001/)).** Single 80-year-old female gastric cancer patient case study by Moisseev et al. (Sechenov University + MGH; Buzdin NOT on the paper despite the Shemyakin IBCh data-deposit affiliation suggesting the Oncobox group). Patient has stomach + esophagus tumor involvement. The BioProject contains **12 sequencing runs across 5 anatomically distinct samples**; we currently use 2 of 12.

| Alias | Tissue | Description | RNA-seq run | WXS run | In pipeline? |
|---|---|---|---|---|---|
| **LST** | stomach | Tumor, surgical section | SRR9143066 (1.96 GB) | **SRR9143067 (34 GB)** | RNA only |
| **LSN** | stomach | Normal | SRR9143065 (1.61 GB) | **SRR9143064 (40 GB)** | RNA only |
| LET | esophagus | Tumor, surgical section | SRR9143070 (1.98 GB) | SRR9143073 (34 GB) | – |
| LEN | esophagus | Normal | SRR9143063 (2.39 GB) | SRR9143068 (20 GB) | – |
| GC_1 | stomach | Primary biopsy | SRR9143072 (1.62 GB) | SRR9143069 (31 GB) + 2× Targeted-Cap | – |

**Pre-computed somatic VCFs for patient_001: none.** Moisseev 2020 deposited raw FASTQ to ENA but published only Supplementary Tables (gene-list / spreadsheet form, not standards-compliant variant files): Supp 1 platform comparison, Supp 2 = 502 germline+somatic from tumor-only, Supp 3 = 386 mutations common across tumor samples, Supp 5 = 137 off-target FoundationOne genes. None are SpliceAI/MMSplice-consumable. Producing a usable VCF for the variant-driven prong requires re-calling from the ENA FASTQ (BWA-MEM2 + Mutect2 or Strelka against GRCh38) — engineering, not data acquisition.

**patient_002 inventory — osteosarc.com Sijbrandij longitudinal cohort.** Patient is IPISRC044 (one patient, 4 clinical time points T0 Nov 2022 → T3 Apr 2025). Sample provenance and somatic VCF deliverables enumerated from the 23,571-line public `https://b2.osteosarc.com/manifest.txt` (Backblaze B2, no auth, HTTPS). **Pre-computed Sarek 3.5.1 somatic VCFs exist across the full longitudinal arc** (cell counts = VCF files per caller × timepoint, including raw + VEP- + snpEff-annotated copies of each somatic call set):

| Timepoint | Platform | Source | Mutect2 | Strelka | Manta | CNVkit | FreeBayes |
|---|---|---|---|---|---|---|---|
| **T0** Nov 2022 | WES | BostonGene | ✅ 4 | ✅ 10 | ✅ 9 | ✅ 2 | – |
| **T0** Nov 2022 | WGS | Personalis | ✅ 4 | ✅ 10 | ✅ 9 | ✅ 2 | – |
| **T1** Jun 2024 | WGS | UCLA / Sarek 3.5.1 | ✅ 8 | ✅ 10 | ✅ 9 | ✅ 5 | ✅ 4 |
| **T2** Jan 2025 | WGS | UCLA / Sarek 3.5.1 | ✅ 6 | ✅ 10 | ✅ 9 | – | – |

Every Mutect2/Strelka somatic call has VEP- and snpEff-annotated forms. T2 uses cross-timepoint pairing (T2 tumor vs T1 blood normal). Crucially, the Sarek 3.5.1 deliverables match the input format SpliceAI/MMSplice/splice2neo natively consume — **zero re-calling effort for patient_002**.

**[Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) WGS/WES discrepancy resolved.** [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) body says patient_002's blood-derived normal is "Blood WGS (DNA only) ~115 GB at `gs://osteosarc-genomics/genomics/...`". What actually got integrated (via [Issue #67](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/67) B2 migration) is `BG003082_WES-normal_1/2.fastq.gz` from BostonGene's Dec 2022 WES — ~9.9 GB total, not 115 GB. The osteosarc.com data portal makes the source of the confusion clear: there ARE both files, the **Dec 2022 BostonGene WES** (~9.9 GB normal, ~19.5 GB tumor) AND a separate **Jun 2024 UCLA Blood Normal WGS** (~115 GB) — and [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37)'s author conflated the two. Doc-fix Issue filed separately; not load-bearing on this audit.

**Hybrid RNA-tumor + WES-normal somatic calling: rejected as a path to the variant-driven splice prong.** Surfaced during patient_002's mid-audit exploration when only the WES germline normal was confirmed (before the matched WES tumor was found). SpliceAI scores variants in the ±50 bp window around exon boundaries — intronic positions. RNA-seq variant callers (RNA-MuTect, RNA-SSNV, Mutect2-on-RNA) are blindest exactly there: spliced reads don't span intronic positions, so SpliceAI's scoring window is invisible to RNA-derived calls. The literature search backs this — splice2neo ([Lang et al., Bioinformatics Advances 2024](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae080/7684965)), NeoSplice ([Bioinformatics Advances 2022](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac032/6581739)), and the NAR Cancer 2023 cryptic-splice paper all use **WES/WGS-derived somatic calls** as SpliceAI/MMSplice input; no 2024–2026 splice-neoantigen paper feeds RNA-derived calls. Documented for future-readers; the path may resurface for other questions (HLA-LOH calling, mutational burden) but is structurally mismatched to this prong.

**Future-cohort options.** GDC API queries (release 45.0, queried today): TCGA-STAD has 415 donors with matched RNA-Seq + WXS (434 open-access masked MAFs, raw VCFs dbGaP-controlled under phs000178); TARGET-OS has 81 donors with matched RNA-Seq + WXS (167 open-access masked MAFs, raw VCFs dbGaP-controlled under phs000468). Critical nuance: **GDC's open-access "Masked Somatic Mutation" MAFs filter out intronic/splice-region variants by construction** — they're coding-region-only and useless for SpliceAI's scoring window. Full raw VCFs require dbGaP DAR. osteosarc.com's Sijbrandij cohort offers the lowest-friction expansion path: same B2 access envelope as patient_002, additional IPISRC patients with native somatic VCFs.

**Decision feeding [Issue #416](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/416): Implement.** Original conditional from [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365): *"If [Issue #384](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/384) shows somatic VCFs available for ≥1 patient AND [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413) shows SpliceAI/MMSplice installable → Implement."* Audit answer: ✅ pre-computed somatic VCFs for patient_002 across 4 timepoint × 2-platform combinations, all publicly accessible. [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413) (Dev-side install) is still open — the final Implement gate is the second condition, not the first. patient_001 has matched-WXS inputs available but no pre-computed VCFs (re-calling needed); audit-tier verdict for patient_001 is "feasible with bounded engineering, not data-blocked".

**Five follow-up Issues filed (separately, not in this PR):**

1. **P1** — `feat(samples): integrate patient_002 pre-computed somatic VCFs (Sarek 3.5.1, osteosarc.com)` — sub of [Issue #416](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/416); the variant-driven-prong PoC. Lightest engineering since Sarek already did the calling.
2. **P2** — `feat(pipeline): patient_001 somatic-calling sub-pipeline (BWA-MEM2 + Mutect2/Strelka)` — sub of [Issue #416](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/416); unblocks patient_001 + any future ENA-sourced patient.
3. **P3** — `docs(samples): patient_001 LET/LEN esophagus + GC_1 biopsy — cohort expansion option` — captures the 8 unused PRJNA545281 runs; schema implications (anatomical-site column) deferred until spatial heterogeneity becomes a manuscript question.
4. **P4** — `docs(scientist): clarify patient_002 normal sample WGS vs WES (correct [Issue #37](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/37) body)` — small doc-fix.
5. **P5** — `feat(cohort): osteosarc.com Sijbrandij cohort onboarding — additional IPISRC patients` — long-horizon cohort expansion.

**Patient_001 cohort-scope recommendation (audit-tier).** Just matched-pair WXS (LST + LSN — SRR9143067 + SRR9143064) for the variant-driven prong question; LET/LEN/GC_1 documented as known expansion options but punted to P3 above. Three reasons: (i) the splice-prong question is binary feasibility, matched-pair suffices; (ii) LET (esophagus) lacks clinical staging context in the ENA metadata to interpret a 2-site comparison cleanly; (iii) adding multiple anatomical sites requires a samplesheet schema change that's its own design decision, not load-bearing on splice-prong feasibility.

**Zotero added.** Moisseev et al. 2020 added to collection Z38GTJNW with tags `Issue-384 patient-001 gastric-cancer Oncobox cohort-source`; three-section HTML note (Findings / Methods / vs. our pipeline) per convention. Used as canonical reference for patient_001's clinical context + variant-data-deposit-form caveat going forward.

**Structural finding worth keeping.** The audit-process arc (samplesheet-inference → external research → vendor-deliverable inventory) is the same shape as the [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) "patient_001 has no WGS — only RNA-seq" finding from 2026-05-16, which I cited as load-bearing context for the initial defer-verdict here. That earlier finding was reached by the same shortcut: read samplesheet, infer DNA-seq absence, write down conclusion. The shortcut was wrong then too — patient_001 has matched WXS at ENA, never integrated. **Lesson: when a prior audit/note reads as "patient has no X data", treat it as a hypothesis to re-verify against the original data source, not a fact to inherit.** The samplesheet records what the pipeline currently uses; the BioProject records what's actually available. The two are not the same.

---

### 14:51 UTC — Editor: Scientist

#### Parent-Status drift resolutions — [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) closed as Done (Exp 1 shipped via chr22 PoC); [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) flipped to Backlog (all unblocked sub-issues done, #271 gated on #203)

Resolving the 2 parent-Status drifts PM surfaced in [2026-05-20 12:37 UTC standup ask](.claude/memory/shared/team_standup.md) — completion drift on [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) and forward drift on [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232). Both are board-hygiene fixes; neither involves new science.

**[Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) (research(filter): patient_001 AlphaGenome predicted-normal validation, Exp 1+2) — close as Done.** The chr22 PoC ([Sub-Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393), closed via [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398) on 2026-05-18) is the operational Exp 1 — [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398)'s body explicitly states *"Executes Experiment 1 end-to-end on chr22 PoC scope"*. Headline metrics shipped: AP 0.214, F1 0.300 at τ=3.16, recall 0.405 (95% bootstrap CI [0.258, 0.333]); decision call **GREEN-with-caveats** (viable secondary stream, not standalone matched-normal replacement). Public-facing [Quarto deck](research/slides/issue_401_alphagenome_chr22_poc/slides.qmd) presented the same result via [PR #403](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/403) ([Issue #401](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/401)) on 2026-05-19.

**Why close — chr22 is the only feasible Exp 1 scope.** patient_001 has no whole-genome data; the chr22 test pipeline (`config/test_config.yaml`, `results/patient_001_test/...`) IS the canonical patient_001 dataset for this pipeline. The Exp 1 question — *does AlphaGenome predict known tissue-expressed splicing from reference alone?* — has been answered against the patient_001 matched-normal ground truth at the only scale the patient has data for. A whole-genome scale-up would need a non-chr22 test dataset and would be a fresh engineering issue, not a continuation of #224's scope. Exp 2 (germline-aware) remains tracked at [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) (blocked on WGS acquisition); independent of this closure.

**Task back-ticking before close.** The 6 unticked body tasks were executed inside [Sub-Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393) but never back-ticked on [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) — stale, not aspirational. Ticked each against [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398) / [PR #403](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/403) evidence (notebook §-references for §1 ground-truth load, §3 universe construction, §4 client + sweep, §5/§6 metrics + decision). Closure-ritual gate cleared. Close action lands as a separate post-merge step (not in this PR), per the convention that issue closures don't ride in code/docs PRs.

**[Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) (docs(manuscript): 2024–2026 splice neoantigen tooling landscape) — flip parent to Backlog; [Sub-Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) stays Backlog.** 7 of 8 ACs ticked on the body (CNNeoPP via [Sub-Issue #267](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/267), Rojas 2023 via [Sub-Issue #268](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/268), METHODS splice2neo+AlphaGenome via [Sub-Issue #269](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/269), DISCUSSION SpliceMutr+ENEO via [Sub-Issue #270](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/270), refs via [Sub-Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272), Pan et al. via [Sub-Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280), PSR_GTEx re-derivation via [Sub-Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299)). Remaining AC — DISCUSSION AlphaGenome ([Sub-Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271)) — is explicitly gated on [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) closing, which is itself gated on [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (Exp 3, comparative filter strength — still open). No active drafting on this epic; parent's *In progress* overclaims. **Backlog** honestly reflects "all unblocked work done, remainder gated on #203". When #203 closes, both [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) + [Sub-Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) lift to Ready.

**Structural finding worth keeping.** Two complementary failure modes on the same parent-Status field: forward drift ([Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) — parent claims In progress while subs are all Backlog) and completion drift ([Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) — parent claims In progress while sub work has already shipped). Both pathologies share the same root cause: **work completion happens at the sub-issue level, but the parent's Status only flips if someone actively edits it.** Yesterday's [Sub-Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393) close (via [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398)) should have prompted a [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) parent-Status review as a same-session reflex — same lesson as yesterday's [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) carve entry (10:22 UTC), and same lesson PM caught at epic-scale in [PR #415](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/415) (3 epics auto-fixed). The recheck-hook mechanism in [PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) catches drift on a fixed schedule; the missing piece is sub-issue-close-time prompts on the parent. Surfacing to PM in the standup follow-up as a potential mechanism-tier follow-up.

---

### 12:42 UTC — Editor: Scientist

#### [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) (research(tcr): NetTCR-struc hybrid structural-QC eval) — classified mode (c), recommended **post-TCRdock structural-QC filter**; first eval session on the freshly-carved issue

NetTCR-struc ([Deleuran & Nielsen, *Frontiers in Immunology* 2025](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1616328/full), DOI [10.3389/fimmu.2025.1616328](https://doi.org/10.3389/fimmu.2025.1616328)) — the **fifth category** in the TCR-pMHC scoring landscape (hybrid structural-QC), orthogonal to the four prior eval issues ([Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET, [Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2, [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES, [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) t2pmhc + TCRLens).

**AC verification on [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422)** (all 5 met):

| AC | State | Evidence |
|---|---|---|
| 1. Methods read + inference-mode classification | ✅ Done | Mode **(c)** — *"Structures of TCR-pMHC complexes were modeled using an AF-Multimer version 2.3 based pipeline"* (Methods §2.1). Input is the AF-Multimer-predicted PDB; not sequence-only. |
| 2. Pipeline-fit recommendation | ✅ Done | **Post-TCRdock structural-QC filter.** Pre-TCRdock ruled out by construction (input is the structure itself); no-fit ruled out (the predicted-DockQ signal closes a real gap). |
| 3. Zotero entry with 3-section note | ✅ Done | Key `FQ58DCAE`; tags `Issue-422 NetTCR-struc tcr-pmhc structural-validation alphafold` |
| 4. Glossary entries for missing terms | ✅ Met | All 9 relevant terms already in [research/glossary.md](research/glossary.md) (GVP-GNN, ESM-IF1, DockQ, AF-Multimer, AF_confidence, CAPRI, CDR3, pLDDT, RMSD); GVP-GNN + ESM-IF1 definitions both already reference NetTCR-struc inline. No edits needed. |
| 5. Lab notebook entry before close | ✅ Met by this entry |

**Classification rationale.** NetTCR-struc fits mode (c) of the [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) scheme — full AF2/TCRdock structure upstream, then scoring. The paper's own framing: *quality-assessment tool for AF-Multimer output*, not an end-to-end binding predictor. Per-residue node features (φ/ψ/ω dihedrals, Cβ direction vectors, edge Cα distances) all derive from the input 3D coordinates — no inference path bypasses the upstream structure requirement.

**Pipeline-fit verdict: post-TCRdock structural-QC filter.** Three reasons:

1. **Mechanical.** Mode (c) input requirement (PDB from AF-Multimer) leaves no "pre-TCRdock" position available — the structure must already exist for NetTCR-struc to run.
2. **Information role.** Current TCRdock-fronted pipeline relies on `AF_confidence` alone for structural quality (Spearman ~0.68 with DockQ per the paper). NetTCR-struc closes that gap: a post-TCRdock gate on `GNN-AF score ≥ τ` filters geometrically-bad predictions before binding analysis sees them.
3. **Stacking with HERMES ([Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218)).** HERMES asks *"does this TCR bind this pMHC well?"* (biology). NetTCR-struc asks *"is this AF-Multimer prediction geometrically accurate?"* (structural QC). Orthogonal questions on the same input — clean stack is NetTCR-struc filtering structural failures → HERMES scoring binding on survivors. Neither replaces the other.

**Headline numbers.** Spearman 0.681 → 0.855 on DockQ prediction (+25% over `AF_confidence` baseline); GNN-AF top-1 selection *"completely avoids selection of 'Incorrect' candidates"* per CAPRI classification on the 25-target post-AF-M-training-cutoff benchmark. Per-target mean DockQ for top-1 selections rises 0.615 → 0.673 (+9.4%). Training set: 80 PDB TCR-pMHC class I crystal structures → 60,000 raw AF-M decoys (5 model weights × 30 candidates × 5 perturbation seeds per target) → 12k–16.5k retained models after filtering, depending on homology-reduction regime (12,057 = full DockQ range + homology-reduced via Hobohm 1; 16,541 = DockQ ≥ 0.5 + no homology reduction; paper §2.3). Bot-review math fix on [PR #426](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/426) — original draft conflated the raw-generation product with the post-filter retained set.

**Trial blocker (intentional defer).** No follow-up trial sub-issue scoped from this eval. The natural trial pairs NetTCR-struc with HERMES under a shared post-TCRdock benchmarking harness, and HERMES itself hasn't been trial-scoped yet ([Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) still pre-trial). The empirical questions worth answering (does NetTCR-struc + HERMES reject different bad predictions? what's the false-positive rate of an `AF_confidence`-only gate?) also need an evaluation dataset that's gated on the AlphaGenome track outputs ([Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) / [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)). Trial sub-issue filed when [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) reaches the same eval-completion state and the AlphaGenome track lands.

**Structural finding worth keeping.** This eval session took ~30 min from issue read to lab notebook — distinctly faster than the original [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) cycle. The acceleration comes from the carve-pattern: [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) was filed *with* classification-scheme link, the 5-niche landscape table, and the stack-with-HERMES hypothesis pre-loaded in the body — so the eval pass was verifying claims rather than discovering them. Lesson for future carves: front-load the eval scaffolding into the issue body (not just the AC list); it converts a multi-session research drift into a single-session check.

---

### 10:22 UTC — Editor: Scientist

#### [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) (research(tcr): t2pmhc + TCRLens hybrid TCR-pMHC scoring eval) — closing as Done; [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) carved for NetTCR-struc structural-QC niche

PM's [2026-05-19 11:56 UTC standup ask](.claude/memory/shared/team_standup.md) board sweep flagged [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) as 14d quiet since the 2026-05-05 NetTCR-struc scope-expansion proposal. Three-option frame: (1) still active, (2) Backlog, (3) close + carve. Pick: **Option 3** — close [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) as Done for the original t2pmhc + TCRLens scope; file separate eval issue for NetTCR-struc.

**AC verification on [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236)** (all 5 met or honestly-deferred):

| AC | State | Evidence |
|---|---|---|
| 1. Candidates classified (a)/(b)/(c) with citation | ✅ Done | [2026-05-02 16:59 UTC comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) — t2pmhc = mode (c), TCRdock cite; TCRLens = mode (b), tFold-TCR cite |
| 2. Pipeline-fit recommendation per tool | ✅ Done | Same comment — t2pmhc → confidence cross-check; TCRLens → triage + cross-check |
| 3. "Trial in pipeline" follow-up sub-issue scoped | ⏸️ Deferred | Both got "yes-trial" recs but benchmarking dataset gates on the AlphaGenome track ([Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) / [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) / [Issue #384](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/384)). Trial sub-issues scoped when those outputs land. |
| 4. "No-fit" rationale | N/A | Both got "yes-trial" recs; AC 4 condition didn't trigger |
| 5. Lab notebook entry before close | ✅ Met by this entry |

**Carve rationale for [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422).** NetTCR-struc occupies a DIFFERENT niche than the original t2pmhc + TCRLens — **hybrid structural-QC** (predicts how good an AF-Multimer-predicted TCR-pMHC structure IS — predicted DockQ) rather than hybrid binding-score (predicts whether a TCR binds its pMHC given a structure). Folding NetTCR-struc into [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236)'s scope would have widened the eval beyond its original first-deliverable framing; carving keeps the landscape decomposition honest — one niche per eval issue, mirroring the existing [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET / [Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2 / [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES / [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) t2pmhc+TCRLens split. NetTCR-struc is the **fifth category** in this landscape.

**Priority on the carve: P3.** Lower than [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236)'s P2 because the 14d quiet indicates the NetTCR-struc proposal is scope-creep that hasn't generated active reading energy — informational eval to scope a fifth category, not on the critical path. User can flip to P2 (parity with the eval-batch siblings) if they prefer; surfaced both readings in [standup follow-up](.claude/memory/shared/team_standup.md).

**Sequencing.** [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) was filed before this lab notebook entry so the close-summary comment on [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) can forward-reference [Issue #422](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/422) (no dangling "to be filed later"). Close action on [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) happens in the follow-up post-merge step — issue closures don't ride in code PRs.

**Structural finding worth keeping.** The 2026-05-02 comment on [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) said *"first deliverable complete"* but the board Status never flipped from In progress, and the 2026-05-05 NetTCR-struc proposal kept the visual "still active" reading even though no further work happened. The lesson: when a comment signals deliverable-completion, the Status flip should be a same-session reflex — not "we'll see if more work happens, then decide". The 14d drift here is the single-issue analog of the parent-vs-children Status drift PM caught yesterday in 3 epic flips ([2026-05-19 11:56 UTC FYI block](.claude/memory/shared/team_standup.md)) — same pathology at a different scale.

---

## 2026-05-19

### 14:55 UTC — Editor: Scientist

#### [PR #414](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/414) closes [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) (SpliceAI/MMSplice Sci-tier evaluation — literature, METHODS framing, cross-confirmation strategy) — merged

Sci-tier deliverables for the variant-driven splice-prediction evaluation (descendant of [parent Issue #222](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/222), splice2neo eval). Three artifacts:

1. **Zotero adds**: SpliceAI primary (Jaganathan 2019, *Cell*, key `P348JRG9`) and MMSplice primary (Cheng 2019, *Genome Biol*, key `THP7T3Q7`); splice2neo (`Z4FAE6QM`) verified. Each carries the 3-section HTML note (Findings / Methods / vs. our pipeline) per the Zotero note convention. Cleaned the auto-injected `morning-reading` default tag from both new entries (they're not morning-news finds, they're targeted additions for this eval).

2. **METHODS.md framing extension** ([L64-88](research/manuscript/METHODS.md#L64-L88)): the splice2neo bullet now names SpliceAI + MMSplice as the underlying deep-learning predictors and articulates "complementary in coverage" framing — *"the variant-driven approach captures splicing changes attributable to specific somatic variants and can detect low-expression events missed by RNA-seq read pile-up; the junction-driven approach ... requires sufficient read support at the junction."* Mirrors the AlphaGenome bullet's "We are evaluating ... outcome reported in the Discussion" pattern. The `@claude` review caught two wordsmithing nits ("prong" → "approach", "porting" → restructured to "add a parallel variant-driven evidence stream"); accepted after technical evaluation, addressed in follow-up commit `275c19f`.

3. **Cross-confirmation strategy** (recorded inline in [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) body): `evidence_tier` column on `tumor_exclusive` output with values `both` / `rna_only` / `variant_only`. Default treatment: `both` ranks above `rna_only`; `variant_only` suppressed unless `--include-variant-only` flag set. Rationale rests on splice2neo's FDR 0.04–0.07 being conditioned on cross-confirmation, plus the "variant_only = no surface evidence in tumor RNA" reasoning. Three open implementation questions deferred to [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413) (variant→junction conversion, DAG position, cohort-permutation FDR).

**Structural carve-out (two splits in one session).** Original [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) was cross-role (`role:scientist` + `role:developer`) and held both Sci-tier framing AND the final implement/defer/won't-do decision. First split early in the session: carved Dev install/sketch tasks into [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413), kept [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) cross-role for the final decision. Second split after the Sci-tier work landed: discovered that `gh issue develop` had auto-linked the branch as a closing-PR for [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365), which would force [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) to close on merge — but [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365)'s final decision was still gated on [Issue #413](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/413) + [Issue #384](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/384) and shouldn't close yet. Resolved by carving the final decision into [Issue #416](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/416) and re-scoping [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) to "Sci-tier evaluation only" so the merge could close [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) legitimately. 1 issue → 3 issues over the session.

**Structural finding worth keeping.** Decision issues with deferred conclusions are structurally incompatible with `gh issue develop`-style PR branches that auto-close on merge — GitHub's public GraphQL has no mutation to remove a `ConnectedEvent` post-creation (verified via schema introspection: `disconnectFromIssue` / `removePullRequestFromClosingIssues` do not exist as public mutations; only `createLinkedBranch` / `deleteLinkedBranch` exist, but the latter requires a `LinkedBranch` node that doesn't get created when only a `ConnectedEvent` is in play). Pre-emptive split into Sci-tier scope + Dev-tier sub-task + final-decision follow-up (three issues from day one) is the clean pattern. Carved-after-the-fact is the recovery pattern (what I did this session); works but adds restructuring noise.

**Closure-audit bot flagged two gaps post-close.** (1) Missing **Priority rationale** line on [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) body — fixed inline via `gh issue edit` (the bot's check is a literal `"priority rationale"` substring; the existing inline `**Priority:** P3 — exploratory; not blocking current pipeline.` didn't include "rationale" as a word). (2) Missing **lab notebook entry** for `developer.md` referencing [PR #414](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/414) — false-positive in the bot's role-attribution: it picks only the alphabetically-first `role:*` label, so a cross-role issue with both `role:scientist` + `role:developer` triggers a `developer.md` check ("developer" sorts before "scientist"). The actual responsible role was Scientist; the actual gap is the missing `scientist.md` entry (this one). This PR fixes that.

---

### 13:00 UTC — Editor: Scientist

#### [PR #403](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/403) closes [Issue #401](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/401) (Quarto slide convention + AlphaGenome chr22 PoC pilot deck) — `@claude review` round-2 fixes, merged

Reviewer's second pass flagged that three issues from review 1 were still unresolved on the branch (commit `0933df8` had only addressed slide-layout/content, not the technical items), plus one new issue. Followed up with commit `5a7cde0`:

**Blockers:**
1. Decision slide emoji vs label — 🟡 + "GREEN-with-caveats" contradicted; resolved to "🟡 AMBER — real signal, not standalone" (matches the honest-headline conclusion of 40% recall = viable secondary stream, not standalone).
2. Remote CSL fetch — `csl: https://www.zotero.org/styles/nature` replaced with locally-cached `research/slides/nature.csl`, referenced as `../nature.csl` from each deck. Offline-render safe; no CDN round-trip per render. This becomes the convention default — README documents the rationale.
3. `quarto=1.9` conda pin in README — dropped; let conda-forge resolve current stable.
4. **GENCODE citation mismatch (new)** — `frankish2021gencode` (NAR 2021, describes v38) was cited but the experiment used v47. Replaced with `mudge2025gencode` (NAR 2025, 53:D966–D975, DOI `10.1093/nar/gkae1078`), which is the canonical citation describing the v47 release. Two in-prose cite keys swapped accordingly.

**Reviewer-correction worth flagging:** Reviewer recommended Cell Reports Methods for the RegTools published version, but PubMed lookup (36949070) verified the actual journal is **Nature Communications** (14:1589, DOI `10.1038/s41467-023-37266-6`). Used the verified citation, not the reviewer's suggestion. Reminder to always verify reviewer-suggested citations rather than blindly applying — even careful reviewers can misremember a journal.

**Nits also addressed:** RegTools entry → published Nat Comm version (vs preprint); `incremental: false` added to pilot deck YAML to match template; `_` → `tx_id` in `_regenerate_figures.py` groupby unpacking for readability.

**Verification:** `quarto render slides.qmd --to revealjs` runs clean; bibliography auto-populates with all 6 entries correctly resolved (Mudge 2025 + Cotto 2023 Nat Comm + Avsec 2026 + Kim 2019 + Pedregosa 2011 + O'Donnell 2020); decision slide reads "🟡 AMBER — real signal, not standalone"; no pandoc warnings. Re-ran the figure regen script post-merge — produces byte-identical PNGs from the cached parquet (AP=0.2136, F1=0.3000 at τ=3.1633, CI=[0.2577, 0.3333]), confirming Test plan #2 reproducibility.

**Closure-ritual reminder.** The `audit_and_merge.sh` gate counts ticked boxes but does not verify them. I relied on the script's clean pass without re-running each Test-plan claim against the current branch tip — items #2 (regen reproducibility), #3 (slide overflow), #4 (Mermaid SVG render) were verified at PR-open time but I didn't re-check after the round-2 fix commit. The closure-audit bot then flagged the missing priority rationale on the Issue + the missing lab notebook entry (this entry). Lesson: closure audit = visually verify each box against current tip, not just rely on the gate's tick-counting.

---

## 2026-05-18

### 11:33 UTC — Editor: Scientist

#### Project convention: research-notebook outputs live under `research/notebooks/<notebook>_outputs/`, not `results/`

User pulled on the conceptual split between `results/` and `research/`. Settled on: **research-notebook-produced artefacts belong under `research/`** because the existing project split is `results/` = Snakemake-rule outputs (pipeline), `research/` = exploratory + manuscript stuff (notebooks, lab notebook, manuscript). Mixing notebook caches into `results/` muddles the `rm -rf results/` reset and obscures where to look for notebook artefacts.

Convention adopted: `research/notebooks/<notebook_name>_outputs/` is the sibling-output dir for any notebook. Concrete change in this PR — moved `results/alphagenome/issue_224_exp1/chr22_stomach_predicted_junctions.parquet` → `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet`. `.gitignore` extended with `research/notebooks/*_outputs/` so the 16 MB parquet stays out of git. Earlier lab notebook entries (this morning's 10:41 UTC + this afternoon's 11:27 UTC) reference the old `results/alphagenome/...` path — they're historically accurate at the time of writing and left as-is per the immutable-entries rule. The notebook setup cell + §6 operational notes carry the new path going forward.

**Why it matters as a convention.** Sets a clear rule for the next research notebook (and Exp 3 in [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)): its outputs go to `research/notebooks/issue_225_*_outputs/`. Discoverable next to the notebook, easy to gitignore, doesn't pollute the production pipeline cache.

---

### 11:27 UTC — Editor: Scientist

#### [PR #398](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/398) review response — AP convention + dense threshold grid bumped headline F1 0.282→0.300, recall 0.351→0.405

`@claude review` on the morning's PR caught three issues + one CLAUDE.md nit. Web-source cross-check (user prompted) was the difference between accepting the bot's framing and verifying its actual impact — turned out two of the bot's items were materially larger than first read.

**Issues triaged + their verified impact:**

1. **Bootstrap CI documentation gap** (bot called it minor). Apply — added explicit "positives-only resampling → conservative" note in §5 markdown + §6 caveat 6.
2. **AUC-PR `distinct` mask first-of-tied + trapezoid integration** (bot called it minor). Verified against [sklearn docs](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html): sklearn explicitly warns trapezoidal AUC "uses linear interpolation and can be too optimistic." Implemented `average_precision_numpy` matching sklearn's AP formula `Σₙ (Rₙ−Rₙ₋₁) Pₙ` exactly — **delta vs `sklearn.metrics.average_precision_score` = 0.000000**. Tie-handling fix (first→last of tied) had near-zero impact on our data (biggest tie cluster at score=0, n=2003, sits at bottom of ranking). Net AP shift: 0.210 → 0.214 (the trapezoid number was actually slightly LOW because of missing recall=0 anchor — not the universal "trapezoid > AP" the literature describes; depends on curve shape).
3. **Threshold grid 41-quantile points is coarse** (bot called it cosmetic). **NOT cosmetic.** Dense grid (`np.unique(scores)`, 5,729 points) shifted best F1 from 0.282 at τ=3.50 to **0.300 at τ=3.16** — a real ~6% improvement that the coarse grid hid. Recall jumped 0.351 → **0.405**. The bot under-rated this one and I would have under-applied if I hadn't run the dense-grid comparison.
4. **CLAUDE.md `--` workaround ordering**: applied; `--` is now the canonical first option (order-independent vs the order-dependent "target before --configfile" alternative).

**Bonus design improvement:** Cached-parquet load branch added to §4 sweep cell. Re-execution is now idempotent — skips the 49 × 1 Mb API calls if `results/alphagenome/issue_224_exp1/chr22_stomach_predicted_junctions.parquet` exists. Made headless re-execution via `jupyter nbconvert --execute` viable without burning more API quota; verified the post-fix notebook outputs match my standalone sklearn cross-check exactly.

**Corrected headline numbers (sklearn-verified):**

| Metric | Morning value (committed in `5bad904`) | Corrected value (this revision) |
|---|---|---|
| AP / AUC-PR | 0.210 (trapezoid + first-of-tied) | **0.214** (sklearn AP exact) |
| Best F1 | 0.282 at τ=3.50 | **0.300** at τ=3.16 |
| Recall at best F1 | 0.351 | **0.405** |
| Precision at best F1 | 0.235 | 0.238 (≈ same) |
| F1 95% bootstrap CI | [0.238, 0.321] | [0.258, 0.333] |
| Best-F1 confusion matrix | TP=91 / FP=296 / FN=168 / TN=7176 | TP=105 / FP=336 / FN=154 / TN=7136 |

**Decision call unchanged.** GREEN-with-caveats holds: recall 40% still means AG is not a standalone matched-normal replacement; it remains a viable secondary evidence stream for the multi-filter design in [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203). The improved numbers tighten the case rather than overturn it.

**Process meta-note.** User's "double-check with web sources" instinct on my acceptance of the bot review was a real save. My initial read of issues #2 and #3 was "small documentation tweaks" — pulling the actual sklearn docs (which explicitly warn against trapezoid AUC-PR) and running the dense-grid comparison showed both were materially larger. Worth a Scientist-facing memory: when a bot review touches a methodology that has a canonical reference implementation (sklearn, scipy, etc.), verify against the reference before deciding "small."

**Scope of this PR revision:** notebook §4 sweep cell (cached-parquet branch) + §5 metrics (dense grid) + cb450049 (sklearn-style AP + last-of-tied) + §5 markdown (AP convention) + §5 plot (titles + symlog x-axis) + §6 markdown (corrected numbers + caveats 5-6 + operational notes). Plus `CLAUDE.md` `--` workaround reorder. Notebook headless-re-executed via `jupyter nbconvert` (sweep loaded from parquet — no API quota burned).

---

### 10:41 UTC — Editor: Scientist

#### [Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393) (AlphaGenome chr22 PoC) — Framing 1 executed end-to-end; AUC-PR 0.210, recall 0.35 → GREEN-with-caveats for scale-up

Morning session, Framing-decision-through-execution arc. F1 framing locked in after recommendation review (annotated-chr22-introns universe; 7,731 introns / 259 positives / 7,472 negatives — production-filter use case mirror). Notebook §4 wired from scratch (3 cells: env + smoke + sweep), §5 refactored to universe-restricted F1 + AUC-PR + bootstrap CI + 3-panel plot. End-to-end run completed.

**Local env setup (one-time).** `conda env create -f workflow/envs/alphagenome.yaml` (~3 min); kernel `splice-neoepitope-alphagenome (Python 3.12.13)` auto-discovered by VSCode. Note: user's normal pattern is pyenv + .venv from `research/requirements.txt` (3.14.4), but `alphagenome 0.6.1` supports ≤3.13 only — conda env used here as a per-notebook isolation. Avoid the option of downgrading the project-wide Python pin for one notebook's SDK constraint.

**AG API operational learnings:**
- **Per-minute MB quota** caps free tier at ~15 MB/min — first sweep hit `RESOURCE_EXHAUSTED` at tile 16 (`gdmscience.googleapis.com`). Mitigation: `time.sleep(5)` between calls (~12 MB/min) + exponential backoff retry on `grpc.StatusCode.RESOURCE_EXHAUSTED` (60/120/240 s).
- **Sweep range restriction** — first naive sweep tiled chr22:0-50.8M (49 tiles); tiles 1-10 (0-10.5 Mb) burned ~10 MB on the acrocentric N region for 0 returned junctions. Annotated-range restriction (`[annotated.donor.min(), annotated.acceptor.max()]`) cut to 39 tiles, all in productive range. Total sweep: 443 s wall, 2,702,077 raw junctions → 2,632,983 unique after dedup.
- **Single stomach GTEx track** in AG inventory (index 212, `biosample_name='stomach'`, `polyA plus RNA-seq`) — max-over-stomach aggregation reduces to a single-column lookup. Coverage asymmetry vs. liver (4 tracks) flagged in §6 caveats.

**Results (Framing 1, chr22 PoC):**

| Metric | Value | Read |
|---|---|---|
| AUC-PR | 0.210 | 6.3× baseline (0.034); the rule-of-thumb 0.7 referred to AUC-ROC, AUC-PR is harder on this 3.4%-prevalence class imbalance |
| Best F1 | 0.282 at τ=3.50 | TP=91, FP=296, FN=168 — the operating point we'd ship if AG-only |
| Recall at best F1 | **0.351** | **The honest headline number** — depth confounder doesn't bias recall. AG identifies ~⅓ of confirmed tissue-expressed introns |
| Precision at best F1 | 0.235 | Lower bound — many of the 296 "FPs" are plausibly tissue-expressed introns the 500K-read matched-normal missed |
| F1 95% bootstrap CI | [0.238, 0.321] | 1000-iter resample of 259-positive bag at best τ — wide as expected |
| AG-predicted introns | 5,728 / 7,731 (74.1%) | Fraction of annotated chr22 introns with non-zero stomach signal |

**Decision call: GREEN-with-caveats for full-genome scale-up.** AG is not a standalone matched-normal replacement at this scale, but a viable **secondary evidence stream** to stack with GTEx pan-tissue ([Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)) + matched-normal-where-available. The Exp 3 comparative-strength experiment ([Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)) is the natural next step. Recall-as-primary-readout matches the production decision shape: "what fraction of normal-expressed candidates does AG correctly call out, so we don't falsely flag them as tumor-specific?"

**Caveats baked in to every headline number** (per §6): chr22 PoC scope; 500K-read matched-normal depth confounder (precision is a lower bound); single stomach GTEx track; annotated-only ground truth (no novel-splicing signal); n=259 positives (wide CI).

**Scope of the PR.** Notebook §4 wiring (3 cells: env loader + smoke + sweep with rate-limit handling) + §5 refactor (universe-restricted F1 + AUC-PR + bootstrap CI + 3-panel plot) + §6 outcome population (this metrics table + decision call + operational notes) + AlphaGenome predictions parquet at `results/alphagenome/issue_224_exp1/chr22_stomach_predicted_junctions.parquet` (2.6M rows). Closes [Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393); parent [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) and [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) stay open (scale-up + decision-rule integration remain).

---

## 2026-05-17

### 19:28 UTC — Editor: Scientist

#### [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) (AlphaGenome Exp 1) — notebook scaffold + chr22 matched-normal BED produced; metric-framing question surfaced

Afternoon session, 1h time-box. Standup clean for Scientist (PM→Dev message about dev-i1 milestone close, not addressed). Work-route view picked [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) Task 1 (matched-normal junction BED via chr22 test config) as the bounded, foundational pick — without the BED, no ground truth, no §4-§5 metrics. Branch `research/scientist/issue-224-alphagenome-exp1-notebook` (created 2026-05-16 via `gh issue develop`) was at origin/main with no own commits yet.

**Local toolchain gap caught upfront: macOS `samtools` not installed.** `which samtools` returned not-found. Per [CLAUDE.md](../../CLAUDE.md) "hisat2.yaml — samtools omitted (libdeflate conflict)", the pipeline relies on system `samtools` (apt-get on Ubuntu cloud VMs). On macOS that's `brew install samtools` — one-time install, ran in background (~3 min, samtools 1.23.1 at `/opt/homebrew/bin/samtools`). Worth noting for future Scientist sessions: the chr22 local-dev path has a hidden brew prerequisite not captured in `prepare_test_data.sh`.

**Notebook scaffold shipped: [`research/notebooks/issue_224_alphagenome_exp1_patient_001.ipynb`](research/notebooks/issue_224_alphagenome_exp1_patient_001.ipynb).** Sixteen cells across six sections — §1 matched-normal load (coord-aware: TSV format is asymmetric `<chrom>:<donor_1based>:<acceptor_0based_exclusive>:<strand>\t<reads>`), §2 GENCODE chr22 intron derivation from exons, §3 intersection ground truth, §4 AlphaGenome predict (stub — actual client wiring deferred until env installed locally), §5 P/R/F1 sweep, §6 decision-rule outcome + caveats. The §1 load cell explicitly normalizes to 0-based half-open intron coords (`donor_0based = donor_1based - 1`; `acceptor_0based_excl` unchanged) so set-ops vs. GENCODE annotated introns are coord-aligned.

**chr22 pipeline run produced the BED.** `snakemake --cores 4 --use-conda results/.../junctions.tsv --configfile config/test_config.yaml` — index build (chr22 only, ~7s) + align + regtools junctions extract + `bed12_to_junctions.py` — total **~1 min wall time** end-to-end. Output `results/patient_001_test/alignment/SRR9143065_test/junctions.tsv` (47,994 bytes, 1,714 junctions, all on chr22). The hisat2 conda env was first-time build during this run (~3 min within the 1 min number — Snakemake reports rule time, not env build time).

**Snakemake 8 argparse gotcha hit + CLAUDE.md extended.** First run attempt failed instantly with `FileNotFoundError: results/.../junctions.tsv` — Snakemake was *interpreting the positional target as a second config file* because of `nargs="+"` on `--configfile`. The existing CLAUDE.md gotcha note covers the "two separate `--configfile` flags" form but not the "target follows configfile" form. Extended the note with the companion form + three workarounds (target before `--configfile`, or `--` separator, or use `-n` which accidentally breaks the nargs sequence). Bundled into this PR per user agreement when surfaced.

**Smoke test of §1-§3 cells against real chr22 data — empirical findings.** Ran `load + intersect` end-to-end in a Bash python invocation under the snakemake env. Numbers:

- Matched-normal junctions: **1,714** (median 1 read, max 23, mean 1.2 — sparse, expected for 500K reads)
- GENCODE chr22 annotated unique introns: **7,731** (from 21,211 exon records across 3,364 transcripts)
- Ground truth = mn ∩ ann: **259** (15.1% of matched-normal, **3.4% of annotated chr22 introns**)

The 3.4% recovery is a structural floor at 500K-read depth — most chr22 annotated introns belong to transcripts not expressed (or under-detected) in normal stomach at this depth. Recorded in the notebook §6 markdown alongside the metric-framing question (below) so the empirical context is captured even before §4-§5 are wired.

**Open metric-framing question surfaced + documented in notebook §6.** N=259 positives is small; AG predictions evaluated against this set will give noisy P/R unless we restrict the evaluable universe. Three framings sketched:

1. Restrict evaluable universe to **annotated chr22 introns** (7,731): pos = ground truth (259), neg = annotated minus matched-normal (7,472). AG evaluated on annotated set only; P/R measure tissue-expression classification ability. ← *cleanest framing for chr22 PoC*; scaffold defaults to this.
2. Whole-genome AG predictions vs. matched-normal ∩ annotated: pessimistic, conflates "should predict" with "predicts at all".
3. Scale up first — re-run on full FASTQs on production GPU so the matched-normal sample is denser; only then evaluate.

Decision needed before §4-§5 execution. Flagging as next-session topic — not blocking the scaffold + BED PR shipping now.

**Scope of this PR.** Notebook scaffold + chr22 matched-normal BED (Task 1 of [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)) + CLAUDE.md gotcha extension. **Does NOT close [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)** — Tasks 4-7 (AlphaGenome client wiring, prediction run, P/R/F1 figures, decision-rule outcome, final lab notebook entry on findings) remain open. Task 3 (env pin) was already ticked via [PR #386](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/386). The chr22 BED also closes Task 2 (ground truth = mn ∩ ann) since the smoke test produced the intersection numbers; ticking that AC.

---

### 14:50 UTC — Editor: Scientist

#### [Sub-Issue #385](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/385) (env pin) → [PR #386](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/386) shipped — Plan B chosen after auto-close foot-gun audit

Continuation of [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) work (Exp 1 patient_001 notebook). Morning routine itself was a no-op (Yu et al. bioRxiv 2026-02-24 already in Zotero + news_log + [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218)).

**Plan choice for scaffold work: 1 sub-issue for env yaml, rest on the existing notebook branch.** Auto-close audit caught a foot-gun: `gh issue develop`-linked branches auto-close their parent issue on PR merge regardless of body content. The user pushed back when I initially claimed only PR-body closing keywords could auto-close — re-read of `feedback_github_auto_close.md` corrected the picture (specific 2026-05-04 incident [PR #260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/260) → [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232)). Mitigation: file sub-issues upfront so each PR closes only what it should. Net: env yaml ships under [Sub-Issue #385](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/385); the existing notebook feature branch (created via `gh issue develop 224` yesterday) is reserved for the final scaffold + execution PR that legitimately closes [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224).

**[Sub-Issue #385](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/385) administration.** `gh issue create` with full body template (Created by + Parent + Scope + ACs + Priority rationale), linked to parent [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) via REST `POST /repos/.../issues/224/sub_issues`. Status field ID lookup needed live query (`PVTSSF_lAHOB17eGc4BSomPzhAHFf8`) — wasn't in memory. Set Size XS + Priority P2 + Status In progress on board.

**[PR #386](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/386) (`workflow/envs/alphagenome.yaml`, 16 lines).** Pinned `alphagenome>=0.6,<0.7` (PyPI 0.6.1 latest, pip-only; verified via `https://pypi.org/pypi/alphagenome/json`). Notebook-side deps on conda-forge (numpy, pandas, matplotlib, scipy, seaborn, pyarrow, ipykernel). Python `>=3.11,<3.13` matching `python.yaml`. Smoke tests: mamba dry-run resolved 160 conda packages cleanly; pip dry-run resolved alphagenome-0.6.1 + 21 transitive deps cleanly (macOS arm64). CI 3/3 green. Pre-merge `closingIssuesReferences` audit: #385 only, #224 NOT in list (correct). Merged 14:49:39Z; #385 auto-closed 14:49:40Z; #224 stayed open. Board auto-set #385 + PR #386 Status to Done.

**Two process misses caught by user.** (1) Shipped [PR #386](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/386) without a lab notebook entry preceding the merge — literal violation of `feedback_closure_ritual.md` "Entry precedes the ship action — never follows it"; recovery is this entry on its own docs branch (separate-entry option per user choice). (2) Yesterday's scope-reduction on [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) was added as a *comment*, but the body's 9-task list was never rewritten to the 7 Exp-1-only tasks — the body has been stale for 24h. Body rewrite + Task 3 tick lands as part of today's recovery alongside this entry.

---

## 2026-05-16

### 20:26 UTC — Editor: Scientist

#### [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) — design audit caught WGS gap; Exp 2 carved off into [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) before branching

Afternoon session, no science news (already shipped at 14:08 UTC the day before). Standup clean for Scientist. Queue audit surfaced [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) (patient_001 notebook for AlphaGenome predicted-normal validation, Exp 1+2) as the natural next pick — Backlog, P2, S, parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (In progress, P1).

**Design audit before branching.** User asked to re-read both #224 and parent #203 (latest revision 2026-05-13 [comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203#issuecomment-4443079527)) to verify scope still holds. Read [Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) closure (API access confirmed, junction connectivity verified via [`scripts/alphagenome_spike.py`](scripts/alphagenome_spike.py), [`docs/alphagenome_primer.md`](docs/alphagenome_primer.md) rich + current), cross-checked patient_001 inputs, conda envs, and recent pipeline state.

**Critical finding: patient_001 has no WGS — only RNA-seq.** [`config/samples/patient_001.tsv`](config/samples/patient_001.tsv) is RNA-seq-only (SRR9143066 tumor + SRR9143065 normal). The pipeline references "WGS" only in HLA serotyping context (germline typing source), never as actual input data. Exp 2 explicitly requires "AlphaGenome on patient_001 WGS with germline variants applied" — not runnable as designed.

**Secondary findings:** (a) [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) regtools BED12 coord-fix merged 2026-05-15 18:38 UTC invalidates any pre-fix patient_001 normal BED — affects ground-truth construction in Exp 1; (b) AlphaGenome SDK not pinned in any [workflow/envs/](workflow/envs/) — needs `alphagenome.yaml` or extension to `python.yaml`; (c) no local `results/` dir — patient_001 outputs live on GCS or require local re-run on chr22 test config.

**Resolution.** User chose (a) scope-reduction over (b) GATK RNA-seq SNV proxy or (d) TCGA-STAD donor substitution. Plus chr22 local re-run over GCS fetch — clean coords, smaller region, avoids the pre-fix bug. Issue admin executed:

1. **Comment on [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)** ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224#issuecomment-4468006508)) defers Exp 2, tightens scope to Exp 1 only, replaces task list with 7 Exp-1-only tasks.
2. **Created [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381)** under parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) for Exp 2 — carries the WGS-acquisition precursor (decision (a) RNA-seq SNV calls vs (b) TCGA-STAD donor substitution deferred into the sub-issue's scoping). Milestone i2-S4 (same as #224, due 2026-05-22) — flagged that Exp 2 likely spills past milestone given the WGS lift.
3. **Linked [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) to parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)** via REST `POST /repos/.../issues/203/sub_issues` — parent sub-issue count now 4 (was 3: closed [Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) + open [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) + open [Sub-Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225); now plus [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381)).
4. **[Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) Status → In progress** (option `47fc9ee4`) via `updateProjectV2ItemFieldValue` GraphQL mutation.

**Workspace prep.** Fast-forwarded `workspace/scientist` from `5eac679` to `f40edb4` (2 commits behind origin/main — [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) regtools fix + [PR #379](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/379) developer lab notebook). Branch `research/scientist/issue-224-alphagenome-exp1-notebook` created via `gh issue develop 224 --base main --checkout`. No code or scaffold yet — that lands next session.

**Why this matters for the manuscript.** [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)'s "publishable either way" framing holds either outcome of Exp 1 alone (tissue-prior vs. per-individual signal) — the carve-out preserves the publishability argument. Exp 2 was always the *patient-specificity* discriminator; deferring it to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) means the AlphaGenome filter validation arc lands in two stages instead of one, but Exp 1 alone yields a real result (with documented WGS-data-availability caveat).

**Process catch.** User's "re-read the issue before branching" instinct surfaced the WGS gap that would otherwise have been caught only after the first failed Exp 2 task. The original [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) body and the parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (refreshed 2026-05-01, revised 2026-05-13) both assumed "patient WGS" without ever cross-checking the sample TSV — a 2-week-old design assumption that aged poorly. Worth a project-level note about pre-branch input-availability checks on multi-experiment sub-issues.

---

## 2026-05-15

### 14:08 UTC — Editor: Scientist

#### Glossary — DDA + DIA-MS entries surfaced from Pepyrus news re-search

Morning Phase 1 (Science news) initially surfaced Prélot et al. 2025 (splice neoepitope methodological-divergence benchmark) — paper turned out to already be in Zotero `ZXAUQAJL` (added 2026-05-08 by Dev session during the `zotero_add.py` bioRxiv bug-hunt that produced [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307)) and [Issue #304](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/304) already cites it. **Dedup miss:** I grep'd `news_log.md` (clean) but not the Zotero collection — morning-news dedup currently walks news_log only. Worth adding a Zotero `/items/top` DOI query to the dedup sweep; the rule `shared/feedback_zotero_dedup_check.md` already exists for `zotero_add.py` calls but doesn't cover the morning-news pre-surface check.

**Fresh re-search** dumped the full 29-item Zotero collection as a denylist, then searched unexplored angles (ERAP1/proteasome, single-cell immunopeptidomics). Surfaced [Manakongtreecheep et al. — *Sensitive detection of cancer antigens enabled by user-defined peptide libraries* (Pepyrus)](https://www.nature.com/articles/s41587-026-03003-9) (Nat Biotech 2026, doi:10.1038/s41587-026-03003-9): user-defined *E. coli* peptide libraries paired with DIA-MS enable ID of low-abundance neoantigens that classical DDA misses, including patient-specific splice-derived ones. >75% recovery from >10K-entry libraries per single injection; 0.1 fmol detection in complex background. Direct manuscript hook: addresses the empirical-validation gap left by Prélot (cross-pipeline disagreement → which candidates are real?) and closes the loop for low-stoichiometry splice neoepitopes specifically.

**This PR ships glossary only.** User asked for DIA-MS to be defined to help understand the Pepyrus paper; DDA bundled since DIA-MS only makes sense in contrast. Two entries inserted alphabetically in section **D** (DDA after DAG, DIA-MS after DDA, both before DockQ). Pepyrus + Prélot light news_log entries queued for the next batch (await user go-ahead on Pepyrus + Zotero add via paywall-deferred manual fetch per `feedback_zotero_defer_inaccessible.md`).

**Why DIA-MS matters in one paragraph.** DDA fragments only the top-N most intense MS1 precursors per cycle — stochastic and intensity-biased, low-abundance neoantigens routinely missed across replicates. DIA-MS fragments every precursor in a stepped m/z window regardless of intensity — comprehensive coverage but multiplexed MS2 spectra require a reference library for deconvolution. Public DIA libraries cover canonical proteomes but not patient-specific neoepitopes, so Pepyrus generates the custom library from *E. coli* expression of predicted candidate sequences before running DIA on the patient sample. This is the first practical MS-validation path for splice-derived neoepitopes at our pipeline's output stoichiometry.

---

## 2026-05-13

### 14:35 UTC — Editor: Scientist

#### [Issue #222](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/222) — splice2neo evaluation (close as benchmark; spun off [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365))

Picked up [Issue #222](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/222) (P2 tool-landscape eval, Scientist-led, decision-only) after `/compact` from earlier in this session. Closest published methodological analog to our pipeline; already framed at a high level in [METHODS.md L64-83](research/manuscript/METHODS.md#L64-L83) from [PR #274](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/274) (variant-driven vs junction-driven METHODS subsection).

**Tool baseline.** [splice2neo](https://github.com/TRON-Bioinformatics/splice2neo) is an actively maintained R package (TRON-Bioinformatics, v0.6.14 Oct 2025). Inputs: somatic VCFs (parsed via SpliceAI v1.3.1 / MMSplice v2.1.1 outputs) + RNA-seq junctions (parsed via RegTools, LeafCutter, SplAdder, others). Filtering: GENCODE canonical exclusion + GTEx 1,740-sample / 53-tissue population-level filter + detection rule `SpliceAI ≥ 0.35 OR MMSplice ≤ −0.35` with ≥3 junction reads. Outputs: peptide candidates with junction position + frameshift annotation. HLA prediction external (NeoFox in paper). Paper cohort: 85+27 melanoma; reports 1.7 candidates/patient at FDR 0.04–0.07.

**Comparison vs our pipeline.** Overlap is the *upstream* peptide-generation step (both accept regtools-style RNA-seq junctions and produce junction-spanning peptides). Differentiation runs both directions:

| Axis | splice2neo | This pipeline |
|---|---|---|
| Evidence streams | Variant-driven (SpliceAI/MMSplice) + RNA-seq detection, cross-confirmed | RNA-seq detection alone |
| Tumor-specificity filter | GTEx population-level (1740 samples, 53 tissues) | Matched normal (per-patient; clinically conservative) |
| FDR estimate | Permutation across cohort (FDR 0.04–0.07) | None |
| HLA presentation | External (NeoFox) | Integrated MHCflurry 2.x `Class1PresentationPredictor` (genotype-aware) |
| Structural validation | None | TCRdock (AlphaFold-based MHC-pep-TCR complex) |
| Genome | GRCh37 (hg19) per paper; hg38 support in repo undocumented | GRCh38 (UCSC hg38) |
| Stack | R / Bioconductor (component) | Python + Snakemake (end-to-end automation) |

**Decision: (b) benchmark, not (a) replace, not (c) irrelevant.** Wrapping an R component mid-Snakemake-DAG would be a significant adapter layer for marginal gain — splice2neo's *peptide-level* output overlaps with our pipeline's mid-stage; the downstream HLA/structural value-add is ours alone. The methodologically valuable piece is the **variant-driven evidence stream** (SpliceAI/MMSplice on somatic VCF), but porting it natively in Python is closer in spirit to [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (AlphaGenome computational normal filter — sibling on the upstream evidence axis) than to "integrate splice2neo".

**Follow-up.** Filed [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) (`role:developer`, P3 exploratory) for Dev-side eval of porting the variant-driven prong as a parallel Python step. Gating question: do current/planned cohort patients have matched DNA-seq + somatic VCFs? If no, defer until inputs available. No DISCUSSIONS hook added — [METHODS.md L64-83](research/manuscript/METHODS.md#L64-L83) already frames the comparison.

Zotero entry `Z4FAE6QM` (splice2neo paper) note retro-updated to convention (Findings / Methods / vs. our pipeline, ~63 words, `<ul><li>` bold-keyword pattern). Initial draft used wrong format (`<p>` prose) and was re-pushed in correct format after user correction.

---

### 13:30 UTC — Editor: Scientist

#### [Issue #361](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/361) — INTRODUCTION.md Sahin 2026 follow-up to [PR #359](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/359)

Held-back follow-up from the [10:14 UTC entry below](research/lab_notebook/scientist.md), picked up post-merge after the [PR #359](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/359) closure-audit / at-mention cleanup wrapped.

**Edit.** [`INTRODUCTION.md:144`](research/manuscript/INTRODUCTION.md#L144) — appended `Sahin et al., *Nature* 2026` to the existing vaccine-slot citation paren supporting the "typically 10–20 candidates in current clinical trials" claim. Citation list grows: `Sahin 2017; Ott 2017` → `Sahin 2017; Ott 2017; Sahin 2026`. Single-line diff; kept the original line-break position (after "primary") rather than reflowing — L144 grows from 75 → 105 chars (longest in file by 8 over L105's 97), accepted as a wrap-stretch since the alternative reflow cascaded into 3-4 surrounding lines for no semantic gain. No other prose change. REFERENCES.md untouched — Sahin 2026 entry at L116-117 + tracking-row at L208 already present from [PR #359](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/359).

**Scope correction (archaeology).** The 10:14 UTC entry below (and the [Issue #351 closure-audit comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/351#issuecomment-4441387311), and the [PR #359](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/359) body) all named "L144 and L192" as the held-back INTRODUCTION locations. Pre-edit grep showed only L143-144 has a Sahin citation (one combined paren with Ott 2017); L192 in current `main` is the Patient_002 osteosarcoma block — no Sahin reference. The L192 ref appears to have been a confabulated line number that propagated through three artifacts before today's correction. [Issue #361](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/361) body documents the correction; the prior artifacts (immutable PR body + audit comments) remain uncorrected — relying on this notebook entry as the canonical correction.

**Argument strength.** The 2017 citations support the slot-count claim with two papers from a single year. Adding Sahin 2026 demonstrates the 10–20 range has been platform-stable across a 9-year span and multiple tumor types (melanoma → PDAC → TNBC). Strengthens the empirical anchor for the "limited number of peptide slots" argument, which is load-bearing for the GPS-as-primary-ranker rationale in this paragraph.

---

### 10:14 UTC — Editor: Scientist

#### [Issue #351](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/351) — Sahin et al. 2026 (TNBC mRNA vaccine) citation in DISCUSSIONS clinical-translation

Picked up [Issue #351](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/351) from Backlog after [PR #349](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349) merged at 09:48 UTC — manuscript context (DISCUSSIONS §"Clinical translation") was hot, and the Zotero entry (key `XUERZBR8`) was already attached from the morning briefing.

**DISCUSSIONS.md edit.** Added a single sentence to the autogene cevumeran / Rojas-Sethna paragraph (DISCUSSIONS.md:459-465), positioning Sahin 2026 as *convergent evidence* for cross-tumor durability rather than starting a new paragraph: same BioNTech personalized mRNA platform applied to TNBC (different tumor type), 11/14 post-surgical patients relapse-free at up to 6 years, vaccine-induced T-cell responses functional for several years. The "convergent" framing avoids redundancy with the Rojas/Sethna detail above (which carries the 3.2-yr / 7.7-yr T-cell-lifespan numbers) and emphasizes what Sahin 2026 adds: tumor-type generalization of the durability finding.

**REFERENCES.md edits.** Added Sahin et al., *Nature* 2026 entry under `## S` between Sahin 2017 and Sethna 2025 (alphabetical by year within author). Full DOI metadata (`10.1038/s41586-025-10004-2`), with `→ Cited in DISCUSSIONS.md (clinical translation section; cross-tumor durability — 11/14 TNBC patients relapse-free at 6y)` annotation. Added matching row to *All current in-text citations* cross-reference table between Sahin 2017 and Sethna 2025. Not marked with ★ (not in the original 8-paper inventory of [parent Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) / [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272)).

**Scope discipline.** INTRODUCTION.md also has two Sahin 2017 citations (L144, L192) in the vaccine-slot rationale ("10-20 candidates per trial") where Sahin 2026 — which used "up to 20 neoantigens per dose" — would be a fresh datapoint. Per [feedback_scope_discipline](shared/feedback_scope_discipline.md), kept Issue #351 scoped to DISCUSSIONS only; INTRODUCTION expansion can be filed as a follow-up if useful.

---

### 09:34 UTC — Editor: Scientist

#### [PR #349](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349) pre-merge grep — one more drift fix at INTRODUCTION:147; section header retained

Final pre-merge `grep -nE "binder|binding"` across `research/manuscript/*.md` surfaced one more caveat-coherence drift hit the bot's two passes missed:

- [`INTRODUCTION.md:147`](research/manuscript/INTRODUCTION.md) — `at least one allele must reach strong or weak binder threshold` → `strong or weak presenter threshold` (describes our pipeline's quality gate, same category as METHODS §6 fix).

**Section-header question closed without rename.** Grep also flagged DISCUSSIONS.md:251 section header `## MHC binding prediction: composite presentation score over affinity-only` and made the REFERENCES.md "Cited in" tweak from [87be2db](research/manuscript/REFERENCES.md) look like it pointed at a non-existent literal header. Inspected: the section's title structure is `Topic: our choice`, where "MHC binding prediction" names the broader field category and "composite presentation score over affinity-only" names our position within it. Renaming `binding` → `presentation` in the header would (a) collapse the rhetorical contrast the section is built around, and (b) create the awkward `"presentation prediction: composite presentation score..."` tautology. The REFERENCES tweak remains valid as an internal-vocabulary description of where O'Donnell is cited; it doesn't require literal header-match because the cross-ref isn't an anchor link.

All other `binder/binding` grep hits classify as biological reality (binding affinity, peptide-binding grooves), Yewdell-quote retention, broader step-name references (the pipeline step is conventionally "MHC binding prediction" across the field), or the §6 caveat itself naming the field convention.

---

### 09:25 UTC — Editor: Scientist

#### [PR #349](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349) bot 2nd-pass follow-up — 2 more drift fixes in Calibration note + L421 parenthetical drop

[@claude review 2nd pass](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349#issuecomment-4439100502) confirmed the DISCUSSIONS apposition from the [08:45 UTC commit](research/lab_notebook/scientist.md) is appropriately explanatory (not redundant) and surfaced 3 new hits in the *Calibration note: presentation_score vs. presentation_percentile* subsection that weren't in scope of the first review.

**Fixes (caveat-coherence):**
- [`DISCUSSIONS.md:411`](research/manuscript/DISCUSSIONS.md) — `weak-binder threshold` → `weak-presenter threshold` (describes pipeline gate behavior).
- [`DISCUSSIONS.md:423`](research/manuscript/DISCUSSIONS.md) — `weak-binder percentile threshold` → `weak-presenter percentile threshold` (same).

**L421 parenthetical drop — scientifically motivated, not just drift:**

Bot called L421 (`(non-binder territory)`) "defensible either way" — covered by the §6 caveat as field-convention orientation. User pushed back with a sharper critique: the phrase *sounds* like quantitative inference of poor binding from a 3–5% percentile, but `presentation_percentile` is **rank-relative to a random peptide background for that allele**, not absolute affinity. On a promiscuous allele like HLA-A\*02:01, a peptide at 3–5% percentile can still have a competitive IC50 (~200 nM) — it's just outranked by many similarly-scoring peptides in the calibration pool. The parenthetical undercuts the very scale-mismatch argument the surrounding paragraph is trying to motivate.

Resolution: drop the parenthetical entirely on L421. L423 immediately below already explicitly states the gate-failure consequence (`"while all alleles remain below the weak-presenter percentile threshold"`), so dropping is genuinely lossless and cleaner scientifically.

**Style call retained as deferred:** bot's first-pass observation that the parenthetical `(no allele reaches weak-presenter threshold)` on METHODS.md:230 is redundant with `> 2%` — still open, separate style call.

---

### 08:45 UTC — Editor: Scientist

#### [PR #349](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349) bot-review follow-up — 3 residual presenter-drift fixes + REFERENCES "Cited in" tweak

[@claude review](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/349#issuecomment-4434913139) on yesterday's manuscript §6 hygiene PR flagged 3 residual binder-vocab hits I'd intentionally left untouched (they describe our pipeline's own quality-gate output, originally classified borderline). Bot's argument was sharp: the terminology caveat I added asserts "throughout this manuscript", so leaving these contradicts the scope claim. Fixing is the cheaper path than narrowing the note's scope.

**Fixes.**
- [`METHODS.md`](research/manuscript/METHODS.md) §6 ranking-and-quality-gate paragraph: `weak-binder threshold` → `weak-presenter threshold`; `non-binding alleles` → `non-presenting alleles`.
- [`DISCUSSIONS.md`](research/manuscript/DISCUSSIONS.md) §"MHC presentation prediction": rewrote `prioritizing candidates that are both strongly bound and well processed` → `prioritizing strong presenters — candidates that combine high MHC affinity with efficient antigen processing` (avoids the "presented + processed" double-counting that MHCflurry 2.x's combined `presentation_score` would otherwise create).
- [`REFERENCES.md`](research/manuscript/REFERENCES.md) O'Donnell "Cited in" line: `MHC binding prediction section` → `MHC presentation prediction section`. Internal tracking should also follow the manuscript's chosen vocabulary, per the bot's minor catch.

**Deferred.** Bot's secondary observation that the parenthetical `(no allele reaches weak-presenter threshold)` is redundant with the `best_presentation_percentile > 2%` condition — left as a separate style call.

---

## 2026-05-12

### 21:10 UTC — Editor: Scientist

#### [Issue #347](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/347) — manuscript §6 hygiene: presenter terminology + phantom Jiang citation cleanup

Terminology-drift scan on `research/manuscript/*` surfaced two bundled issues.

**Issue 1 — phantom citation.** Three places (METHODS.md:202, DISCUSSIONS.md:280-282, REFERENCES.md:56) cited *"Jiang et al. (2024, Communications Biology)"* as the precedent for the 0.5%/2% percentile threshold on MHCflurry's `presentation_percentile`. REFERENCES.md entry had placeholder title "TBD" and the paper cannot be located via Zotero, PubMed, or web search — confirmed unverifiable. Pivoted to O'Donnell et al. 2020 (MHCflurry 2.0, *Cell Systems*, DOI `10.1016/j.cels.2020.06.010`) — the actual source of MHCflurry's documented default cutoffs, already in REFERENCES.md (line 87) but missing from Zotero. Added O'Donnell to Zotero (key `SBQGHWRP`) with three-section note.

**Issue 2 — presenter terminology drift.** Two lines described our pipeline's own outputs using legacy "binder" vocab: METHODS.md:233 ("top strong-binding candidate" → "top strong-presenting candidate") and INTRODUCTION.md:119 ("qualifies as a strong binder" → "…strong presenter"). The remaining 7 "X-binder threshold" hits across the manuscript reference the IEDB/NetMHCPan binder-percentile convention — rather than rewrite all 7, added an explicit *Note on terminology* caveat at METHODS §6 that acknowledges the adaptation and names the cutoffs as "strong/weak presenter thresholds" throughout. INTRODUCTION.md:125 "strongest-binding epitope" legitimately retained (cites Yewdell & Bennink 1999 immunodominance framework — sourced biological vocabulary).

**Edits.** 8 surgical edits across REFERENCES.md (4 — Jiang entry removal + O'Donnell "Cited in" update + tracking-table row + title-TBD list), METHODS.md (2 — caveat + line 233), DISCUSSIONS.md (1 — line 280-282 O'Donnell pivot), INTRODUCTION.md (1 — line 119).

### 12:30 UTC — Editor: Scientist

#### [Issue #342](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/342) — RESULTS.md prose: duplicate-peptide dedup misrepresentation

User scrutiny pass on the patient_001 RESULTS.md tables surfaced a count discrepancy: *Peptide Translation* reports `Total = 1,286,492` / `Unique sequences = 1,260,074`, while *MHC Presentation Predictions* reports `Total predictions = 1,286,492` — suggesting MHCflurry is called on non-deduped peptides (~26k wasted predictions). Traced [`run_mhcflurry.py`](workflow/scripts/run_mhcflurry.py): dedup happens at [line 427](workflow/scripts/run_mhcflurry.py#L427) (`peptides_df["peptide"].unique()`), prediction runs on the 1,260,074 unique sequences only, then the predictions are merge-joined back to all 1,286,492 source rows at [line 475](workflow/scripts/run_mhcflurry.py#L475) to preserve `contig_key`/`start_nt` traceability.

**Verdict:** compute is already efficient; the prose just misrepresents the design. Current wording *"one prediction per input peptide; duplicate sequences from different junctions or reading frames are each predicted separately"* is factually wrong — duplicates are predicted **once** and **replicated across source rows**, not predicted independently.

**Fix.** Rewrote `RESULTS.md` lines 72–76 to explicitly state: (a) one row per input peptide position in the output; (b) MHCflurry runs only on the 1,260,074 unique sequences internally; (c) predictions are joined back per source row for `contig_key`/`start_nt` traceability. Added inline code links to `run_mhcflurry.py:427` + `:475` so a future reader can verify directly. Patient_002 section grepped clean — no equivalent misleading prose to fix.

### 12:10 UTC — Editor: Scientist

#### [PR #340](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/340) review pass — citation-form fix on [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334)

Back from lunch. [@claude review](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/340#issuecomment-4429566946) flagged one substantive item on [PR #340](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/340): the in-text citation `(Onkar et al., *Cell Rep Med* 2026; Bhardwaj lab)` reads as a two-work citation because a semicolon inside `()` is the standard separator for multiple sources (e.g. `(Smith 2020; Jones 2021)`). The bot's catch is correct — lab attributions belong in prose, not citation brackets.

**Fix.** Moved Bhardwaj lab attribution into prose: `a recent field synthesis from the Bhardwaj lab (Onkar et al., *Cell Rep Med* 2026) identifies…`. Re-wrapped the paragraph to the file's ~78-char convention (line 438 had drifted to 97 chars post-edit). Diff scope: 7 lines changed in [`research/manuscript/DISCUSSIONS.md`](research/manuscript/DISCUSSIONS.md), no other files.

All other review items (alphabetical placement, cross-reference table, lab notebook ordering, no-★ marking, DOI/volume year discrepancy) cleared as correct — no further action needed.

### 10:16 UTC — Editor: Scientist

#### [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) — Onkar/Bhardwaj 2026 DISCUSSION citation

Picked up [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) from Backlog as a short ship before lunch. Work is the natural follow-on to the morning news_log entry; Zotero entry was already attached this morning (key `P6K7QWR4`) so this PR is purely manuscript-side.

**DISCUSSIONS.md edit.** Added a new opening paragraph to the *Clinical translation: durability of personalized neoantigen vaccines in low-mutation tumors* section (before the existing "Among multiple active personalized neoantigen vaccine trials..." paragraph). The opener uses Onkar et al. as the analytical-synthesis framing — neoantigen + ICI convergence + off-the-shelf shared-NA vaccines as the emerging strategy + manufacturing-time / biomarker as unsolved field-wide problems. Naturally bookends with the existing trial-pipeline reference (Iamukova & Alferova 2026) and sets up the Kwok subsection downstream where the personalized-vs-shared axis is concretely treated for splice neoantigens.

**REFERENCES.md.** Added Onkar et al. entry under `## O` (alphabetically between O'Donnell and Ott), with full DOI metadata (`10.1016/j.xcrm.2025.102575`) and `→ Cited in DISCUSSIONS.md (clinical translation section; field synthesis)` annotation. Added matching row to the *All current in-text citations* cross-reference table. Not marked with ★ since Onkar/Bhardwaj is not in the original 8-paper inventory of [parent Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) / [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272); it's a discretionary morning-news addition.

**Issue body typo.** AC 1 references `DISCUSSION.md` (singular) but the actual file is `DISCUSSIONS.md` (plural) — corrected in the Issue body when ticking the AC. Pre-merge body edit; not a hidden change.

### 09:06 UTC — Editor: Scientist

#### Morning routine — news_log + Zotero + [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) (Onkar et al. 2026 cancer-vaccine field synthesis)

Surfaced Onkar et al. (Bhardwaj lab) review in *Cell Reports Medicine* during morning routine. Synthesizes decades of cancer vaccine work, names neoantigen + ICI combinations as the convergent field direction, and off-the-shelf shared-neoantigen vaccines as the emerging strategy to bypass personalized manufacturing time/cost.

- **news_log**: bullet under `## 2026-05-12 / 09:06 UTC — Editor: Scientist`, linked to [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334).
- **Zotero**: added to collection `Z38GTJNW` (key `P6K7QWR4`), three-section HTML note (Findings / Methods / vs. our pipeline). Pre-checked DOI for dedup — no existing entry.
- **[Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) (parked in Backlog)**: AC checkboxes, Priority rationale (P2 — review reference; pairs with Kwok + Kim shared-splice-NA axis), Created-by attribution. Concrete DISCUSSION hook: clinical-translation section as field-state opener. No board Status/Priority flip — `gh project` auth scope missing this session; left for whoever picks up.

No manuscript files touched — citation work waits for [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) to be picked up.

**Closure-ritual note for self.** Briefly overshot: created the [#334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) linked branch (`docs/scientist/issue-334-bhardwaj-discussion`) via `gh issue develop` and started drafting the manuscript edit before user re-grounded the scope. Rolled back: empty branch deleted locally + remotely. Lesson: "create Issue for X" ≠ "start work on X" — Backlog parking is the default unless explicitly asked to pick up.

---

## 2026-05-10

### 17:04 UTC — Editor: Scientist

#### [PR #322](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/322) follow-up — Bradley/Alam TCRdock citation resolved (open item #1)

While waiting for the bot review pass on [PR #322](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/322), used the time to resolve the deferred open item #1 from the REFERENCES.md draft (orphan-branch flag: *"Bradley et al. (METHODS §7) ↔ Alam et al. Science 2023 — almost certainly the same paper; reconciliation deferred"*).

**Verification.** Three independent confirmations that TCRdock = Bradley P., *eLife* 2023 (single-author paper):

1. WebSearch returned the canonical citation: *Bradley P. Structure-based prediction of T cell receptor:peptide-MHC interactions. eLife 2023;12:e82813.* doi:10.7554/eLife.82813
2. The project's own infrastructure cites Bradley directly: [`docker/Dockerfile.pipeline:15`](docker/Dockerfile.pipeline#L15) clones from `phbradley/TCRdock`, [`workflow/scripts/run_tcrdock.py:12`](workflow/scripts/run_tcrdock.py#L12) references the same repo.
3. WebSearch for "Alam et al. Science 2023 TCR" returned no matching paper — the Alam attribution in the orphan draft appears to be a hallucination from the previous bot session.

**Fix applied across the manuscript:**

- [INTRODUCTION.md:68](research/manuscript/INTRODUCTION.md#L68) toolchain table: `(Alam et al., *Science* 2023)` → `(Bradley, *eLife* 2023)`
- [DISCUSSIONS.md:704](research/manuscript/DISCUSSIONS.md#L704) structural validation paragraph: `(Alam et al., *Science* 2023)` → `(Bradley, *eLife* 2023)`
- [METHODS.md:233](research/manuscript/METHODS.md#L233): `TCRdock (Bradley et al.)` → `TCRdock (Bradley, *eLife* 2023)` (was correct author, missing year/journal)
- [REFERENCES.md](research/manuscript/REFERENCES.md): dropped the speculative `Alam et al., *Science* 2023` entry, replaced the placeholder `Bradley et al.` with full `Bradley, *eLife* 2023` metadata. Cross-reference table updated; open item #1 removed (3 remain).

**Note on author form:** TCRdock is a single-author paper (Philip Bradley, Fred Hutch), so the citation is `(Bradley, *eLife* 2023)` not `(Bradley et al., …)`. Flagged in the REFERENCES.md entry for future reference.

Folded into [PR #322](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/322) rather than opened as a separate follow-up: same files, naturally bundles, resolves an open item the PR itself surfaced. Bot review will re-trigger on the new commit.

### 16:26 UTC — Editor: Scientist

#### [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) — REFERENCES.md cherry-pick + refresh against today's manuscript state

**Background.** [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) is the citation-finalisation umbrella for [parent Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) — three closed sibling issues ([#269](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/269), [#270](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/270), [#311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) closed today) explicitly defer their citation work here. A 2026-05-05 orphan branch (`claude/issue-272-20260505-1025`) has a 214-line REFERENCES.md draft from a previous bot run; PM agreed 2026-05-06 16:52 UTC to *(a) cherry-pick + refresh* against today's state.

**Survey before branching.** Cross-checked the 8-paper inventory against today's manuscript:

| Paper (Zotero key) | Cite present | Section |
|---|---|---|
| CNNeoPP (`6RWWUDPC`) | ✅ | INTRODUCTION.md |
| splice2neo (`Z4FAE6QM`) | ✅ | METHODS.md §2-3 |
| AlphaGenome (`UZWZ5QEB`) | ✅ | METHODS.md §3 |
| ENEO (`9T3C58HQ`) | ✅ | DISCUSSIONS.md |
| SpliceMutr (`VQMU6JWH`) | ✅ | DISCUSSIONS.md |
| Rojas et al. 2023 (`UIN9DIUP`) | ✅ | INTRODUCTION.md + DISCUSSIONS.md |
| Sethna et al. 2025 (`N2QF8MC6`) | ✅ | DISCUSSIONS.md |
| Iamukova & Alferova 2026 (`II24UIUZ`) | ✅ | DISCUSSIONS.md |

Plus two new papers added to the manuscript since the orphan draft (and not in its inventory): **Kwok et al. 2025** (`5ZT8KC8X`, *Nature*) added by [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) closing [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280); **Kim et al. 2025** (`XB3CPX5P`, *Cell*) added today by [PR #315](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/315) commit [`70591c9`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/70591c9) (closing [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311)).

**Sibling state.** 6 of 7 sub-issues closed; only [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271) (AlphaGenome validation strategy) still open — but it's deferred until [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) closes and is about adding a substantive DISCUSSION subsection rather than new citations. Proceeding now is consistent with the issue body's "depends on sibling sub-issues closing first" gate (effectively cleared).

**File written: `research/manuscript/REFERENCES.md` (~210 lines).** Cherry-picked the orphan's structure (alphabetical sections, ★-flagged inventory entries, cross-reference summary table at the bottom, open-items list). Refreshed:

- All 5 splice-tooling papers had "TBD title" in the orphan; pulled full author lists + titles + journal + DOI from Zotero (collection `Z38GTJNW`) for each.
- Sethna et al. and APJCO 2026 entries upgraded from "title TBD" to full metadata.
- **Kwok year fix:** orphan listed *Nature* 2024; correct is 2025 (DOI `10.1038/s41586-024-08552-0`, online 2026-02-19, in print 2025-03-13). Lab notebook caught this on 2026-05-07.
- **Kim et al. entry added** (sister paper to Kwok; cited in DISCUSSIONS.md line 661 via [`70591c9`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/70591c9)).

**Cross-reference verification.** Grepped all `(Author et al., *Journal* YEAR)`-pattern citations in INTRODUCTION/METHODS/DISCUSSIONS; every in-text cite has a matching reference entry, and every reference entry has at least one in-text use. Style discrepancy flagged in cross-ref table: orphan's house style is abbreviated journal names (`*Nat Methods*`, `*Cell Syst*`, etc.), but several recent in-text citations use full forms (`*Nature Methods*`, `*Frontiers in Immunology*`, `*NAR Genomics and Bioinformatics*`, `*Cancer Research Communications*`) — listed as open item #3 for pre-submission proofread sweep, non-blocking.

**Deferred to follow-up:** the orphan flagged METHODS.md §7's `Bradley et al.` (TCRdock; no year/journal in-text) as "almost certainly the same paper" as DISCUSSIONS's `Alam et al., *Science* 2023`. Claim is unverified; preserved both entries with cross-link warning, deferred reconciliation to a future PR with TCRdock provenance verification (open item #1). Also deferred: AlphaGenome validation strategy citations (pending [Issue #271](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/271)/[Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) closure).

#### Earlier today — [PR #315](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/315) review-fix follow-ups (merged) + [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) manual closure

**[PR #315](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/315) bot review** surfaced 3 bugs + 3 observations on yesterday's commit. Fixed across 4 follow-up commits, all merged via squash to main as [`d173ed2`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/d173ed2):

1. **[`e26baec`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/e26baec)** — bug fixes in `research/scripts/issue_299/`: signed-offset matching in `find_a3_candidates()` (replaced `abs(offset) == loss_nt` with `(acceptor − a) == sign·loss_nt` matching the notebook), endpoint-based canonical filter in aggregation loop (replaced `j["annotated"] != "1"` — Snaptron's `annotated` is a count of DBs, not a boolean — with `is_donor_annotated() and is_acceptor_annotated()`), and `Path(__file__).parent`-anchored output in `snaptron_query.py`. Empirical outcome unchanged (RPL22 1/9,662 REJECT; GNAS undetectable); canonical pool roughly doubled after the annotated-filter fix (GNAS 34→68, RPL22 28→113).
2. **[`70591c9`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/70591c9)** — DISCUSSIONS.md: replaced "discussed elsewhere" with `(Kim et al., *Cell* 2025)` inline citation + HTML TODO comment pointing at [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) for the eventual SF-mutation paragraph. Verified Kim et al. covers SRSF2/SF3B1/U2AF1-mutant myeloid malignancies via 2026-05-08 lab notebook entry before citing.
3. **[`fcbb41f`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/fcbb41f)** — DISCUSSIONS.md: tightened `~789 public NJ pool` to `789 characterized public NJs (the glioma-focused set taken through peptide-presentation validation)`. Verified against the Kwok PDF (downloaded via Zotero API): `789` appears 5+ times in the paper (Fig. 2 caption, Fig. 3 prose, line 697 peptide-presentation section, Extended Data Fig. 2). Bot's provenance concern was unfounded but the framing genuinely conflated the glioma-focused set with the pan-cancer pool (Kwok cites "94 public NJs per TCGA tumour type" on average across 12 cancer types, not enumerated in aggregate).
4. **[`d6df95b`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/d6df95b)** — `identify_nej_candidates.py` docstring: added "Prereq: run `snaptron_query.py` first" since the input TSVs are gitignored due to size.

Bot re-review confirmed *"This PR is ready to merge"* with no remaining blockers; merged via squash + branch deleted; project board flipped to Done automatically.

**[Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) manually closed** — work landed via [PR #315](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/315) commit [`70591c9`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/70591c9) as a side-effect of bot review obs #4. PR body had `Closes #299`, not `Closes #311`, so the auto-link didn't fire. Closed manually with a pointer comment; ACs #1, #2, #4 ticked; AC #3 (manuscript reference list) deferred to [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) — picked up in this same PR.

#### Process notes from this session

- **`shared/MEMORY.md` post-compaction gap.** I chained `git commit && git push` for 4 commits today (PR #315 follow-ups + the [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) closure body edit). User caught it. Root cause: the `Commit, push, merge — three separate steps` rule lives in `shared/MEMORY.md` Always-in-effect (line 31), but `shared/MEMORY.md` was loaded at session start and dropped during a `/compact` mid-session — leaving only `scientist/MEMORY.md` (which didn't have the rule inlined). Fix: copied the rule inline to `scientist/MEMORY.md` Always-in-effect (now line 13), with the post-compaction root cause documented in the bullet itself. Same pattern as the archive-don't-delete rule promoted by PM 2026-05-09 11:50 UTC.
- **Standup hygiene.** Flipped own [2026-05-09 11:10 UTC] post (archive-don't-delete proposal) to Done since PM actioned it; archived 3 stale Sci-authored Done >3 days messages (2026-05-06 dates) into `team_standup_archive/2026-05.md` (51→54 messages, `_index.md` count bumped); standup file shrunk 413→360 lines.
- **[Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) closure audit cleared.** PM's 2026-05-09 11:09 UTC nudge: 6/6 ACs unticked + missing Priority rationale. Verified each AC against [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) (DISCUSSIONS + news_log + lab notebook); all met. Body edited to tick all 6 + add `**Priority rationale:**` line; reply comment posted with one annotation: AC #5's "news_log line 23" actually landed at line 37 because line numbers shifted between issue creation and merge — same fix, different line.

---

## 2026-05-09

### 21:09 UTC — Editor: Scientist

#### [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) — PSR_GTEx validation tightening + DISCUSSION edit landing

**Notebook audit & tightening.** Yesterday's [`36771df`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/36771df) PSR_GTEx validation notebook produced the headline result (NEJ_RPL22 detected at offset −6 in 1/9,662 GTEx samples, NEJ_GNAS undetectable). A user-driven audit session this evening surfaced four logical gaps; fixed in commit [`c0f3965`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/c0f3965). Empirical outcome unchanged.

1. **Outcome states (§1) — exposed Kwok's relative-frequency inner gate.** Original two-row outcome table conflated `PSR_GTEx = 0%` with "no NEJ reads in any sample" and concluded "filter not the discriminator." But Kwok's PSR has an inner gate (sample counts only if `NEJ/canonical ≥ 1%` at same donor) — so PSR=0% can also mean "NEJ reads exist but all sub-1% relative to canonical." The RPL22 result (1 sample at ratio 0.92%, sub-1% gate) lands in the "discriminator IS active" row that the original table missed. Replaced with four-state table mapping each per-sample state to PSR and our absolute-presence filter independently.

2. **Signed A3 matching (§5) — `abs(offset)` → strict signed offset.** Original `abs(offset) == loss_nt` matching admitted both A3 *loss* and A3 *gain* offsets ("sign can go either way"). WebSearch verification against [rMATS-turbo README](https://github.com/Xinglab/rmats-turbo/blob/master/README.md) and a [worked + strand example](https://groups.google.com/g/rmats-user-group/c/LWWvruwr-pg) confirmed the sign is biologically determined: A3 "loss of N nt" shifts the acceptor downstream-in-mRNA by exactly N, mapping to **+N on + strand** and **−N on − strand** + strand coords. Tightened to `(acceptor − a) == sign·loss_nt`. Re-run via `research/.venv` `jupyter nbconvert --execute` reproduced same hits (RPL22 at −6, GNAS zero).

3. **Mechanism scope caveat (§9) — alt-splicing on WT DNA vs somatic indel.** Added a new caveat: the predicate assumes alt-splicing on reference-sequence DNA (cryptic acceptor exists in the genome and the spliceosome chooses it). Matches Kwok's public-NEJ premise (cross-patient recurrence on diverse mutational backgrounds requires WT splice sites) and our use of GTEx healthy tissue. **Somatic indel at the splice site** is a different mechanism: canonical AG destroyed, nearby cryptic AG forced into use, aligner-reported junctions may not satisfy the signed-offset rule. Indel-driven splice neoantigens are typically *private* and map to the personalized end of the public-vs-personalized DISCUSSION axis.

4. **Sample-level QC check (§9) — single positive sample 51959 metadata.** User flagged a real confounder: the discriminator finding rests on **one** sample (`rail_id 51959`, `GTEX-SNOS-1126-SM-4DM67`) with **one** NEJ read; if 51959 is in the 5% Kwok dropped to reach `n=9,166`, the comparison isn't apples-to-apples. Queried Snaptron's `/samples?ids=51959` endpoint: RIN 7.1, mapping rate 93.6%, **`SMAFRZE = "USE ME"`** (in GTEx's standard analysis freeze), Hardy 1, 54-min ischemic time, 33,902 genes detected — no red flags, likely in Kwok's subset under any reasonable QC. Tissue is **testis** (splicing-permissive); the single read is consistent with testis-specific splicing noise. Cohort-size caveat walked back from "immaterial" — the 9,662 vs 9,166 difference IS load-bearing for our absolute-presence filter (only PSR_GTEx itself is robust to denominator choice).

**DISCUSSION edit landed in [research/manuscript/DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md).** Three edits to the existing Kwok subsection (lines 628–674):
- Resolved the `<!-- TODO(#299) -->` HTML comment by replacing the "is not resolvable from the published manuscript" sentence with the empirical PSR finding (1/9,662 RPL22 at PSR ≈ 0%, hedged with sample-51959 inclusion contingency + testis context).
- Inserted a mechanism-scope sentence after the public-vs-personalized axis paragraph: both frameworks take the molecular event to be alt-splicing on (mostly) WT DNA, since Kwok's public criterion only holds when the cryptic acceptor exists in the reference genome. Indel-driven splice neoantigens map to the personalized end of the axis. Splicing-factor mutations (SF3B1/SRSF2 — Kim et al. *Cell* 2025, to be folded in via [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311)) flagged as discussed elsewhere.
- Removed the resolved TODO marker.

**Sparring note.** Productive user-driven audit. Each gap surfaced by a user question requiring WebSearch (rMATS conventions), the Snaptron samples endpoint (51959 QC), or careful biology reasoning. The user's pushback on "intron extends" language led to a precision improvement — the canonical intron doesn't *extend* when an alternative acceptor is used; the spliceosome *selects a different intron* with different boundaries (donor stays fixed, acceptor moves). Notebook §4 and §5 now use this more careful framing.

---

## 2026-05-08

### 15:03 UTC — Editor: Scientist

#### [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) follow-up — PSR=0% framing fix + spin-off [PR #302](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/302) glossary additions

**Kim et al. *Cell* 2025 — sister paper on the public-neoantigen axis.** Morning routine surfaced [Kim et al. *Cell* 2025](https://www.cell.com/cell/fulltext/S0092-8674(25)00399-X) (DOI `10.1016/j.cell.2025.03.047`): mis-splicing-derived neoantigens in SRSF2/SF3B1/U2AF1-mutant myeloid malignancies. Same off-the-shelf TCR-from-healthy-donor strategy as Kwok et al., but **mechanism-driven** (mutant splicing factor → stereotyped mis-splicing → CLK3 / RHOT2 neoantigens) rather than mechanism-agnostic recurrent NJs. MHCflurry 2.0 in their stack (alongside NetMHCpan 4.0) — partial overlap with our `Class1PresentationPredictor`. Concept-relevant for the public-axis end of [DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md) but not target-relevant (myeloid hotspots vs our PDAC/melanoma cohorts). Zotero add deferred — note drafted for user review, action q agreed but not yet pushed through `zotero_add.py`. One-line DISCUSSION addition agreed in principle; landing deferred to next manuscript-touching PR (option (c) of three considered) rather than scope-creeping PR #300 or opening a tiny standalone PR.

**PSR definition verified against Kwok PDF.** While drafting glossary entries, the question came up: what does PSR actually expand to? Read [Kwok et al. *Nature* 2025](https://www.nature.com/articles/s41586-024-08552-0) page 1, "Characterization of public, pan-cancer NJs" paragraph — definition is *"the percentage of samples in a cohort that express the NJ with a read frequency of ≥1% relative to the canonical splicing junction"*. Expansion: **Positive Sample Rate** (not "Percent Sample Recurrence" / "Percent Splicing Recurrence" as initially guessed). Two-tier metric: per-sample inner threshold (≥1% relative read frequency) AND outer cohort fraction. Kwok cites Nejo et al. (*Nat. Med.* 2023, citation 14) for the upstream NJ nomenclature.

**Caught DISCUSSIONS.md inaccuracy via the verified definition.** Yesterday's [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) text on line 656 said *"this pipeline applies `min_read_count: 1`, equivalent to `PSR_GTEx = 0%`"* — strictly NOT equivalent. Kwok's `PSR_GTEx = 0%` means no sample crosses the ≥1%-of-canonical inner threshold, so junctions with reads at <1% relative read frequency in any number of normal samples are still admitted by Kwok. Our `min_read_count: 1` excludes any normal-sample read regardless of relative frequency — strictly stricter. Worked example: a junction with 1 read in 50 GTEx samples each at 0.5% relative-frequency has `PSR_GTEx = 0%` (Kwok would include it) but is excluded by our pipeline. Fix landed as commit [`f01499f`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300/commits/f01499f) on the existing PR #300 branch — replaces "equivalent to" with "stricter than" + adds the per-sample threshold explanation. Self-initiated, not bot-requested. CI (`pipeline-snakemake-dry-run` + `pipeline-pytest`) green on the new HEAD; the [yesterday's bot re-review](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300#issuecomment-4399584920) "ready to merge" verdict carries forward.

**Spin-off — glossary additions in stand-alone [PR #302](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/302).** 6 new entries to [research/glossary.md](research/glossary.md): **CAR-T**, **HLA-LOH**, **NJ** (covering Kwok's **NEJ** subset notation), **PSR**, **TCR-T**, **TIL**. Seeded by the Kim et al. read (TCR-T needed glossary attention), then expanded to nearby gaps spotted while reading the file (NJ/NEJ for the new Kwok subsection vocabulary, CAR-T to disambiguate from TCR-T cleanly, TIL for the third major adoptive-cell-therapy pillar, HLA-LOH for the immune-escape mechanism). Per user instruction this session, treated as journal-style PR — branch `docs/scientist/glossary-2026-05-08` directly off `origin/main`, no parent issue, mirror of the news_log workflow ([PR #298](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/298) precedent). PSR entry lifted directly from the verified Kwok PDF text. NJ entry attributes the notation to Nejo et al. *Nat. Med.* 2023 per Kwok's citation 14. Skipped the `@claude review` offer per the journal-style convention.

**Memory cleanup — `git -C <CWD>` redundancy.** During glossary PR execution, was reflexively prefixing every git command with `git -C /Users/.../splice-neoepitope-pipeline-scientist/` even though CWD already matched. User flagged: *"why are you using git -C all the time? Just cd once into the repo and stay there."* Existing rule in [shared/feedback_no_cd.md](shared/feedback_no_cd.md) covered cross-repo discipline (*"use `git -C <path>` for git on OTHER repos"*) but didn't spell out the contrapositive (bare git for CWD). Extended `feedback_no_cd.md` with an explicit bullet + caught-example section, and tightened the inline rule on `shared/MEMORY.md` line 12 to make the "for CWD, bare git" path unambiguous. Per the skip-for-minor exception in `feedback_team_standup.md` line 72 (*"additional references that don't change behavior"*), the change is a clarification of an existing rule rather than a new behavioral rule — skipped the `team_memory_broadcasts.md` post (close call; user agreed with skip).

---

## 2026-05-07

### 17:18 UTC — Editor: Scientist

#### [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) review-feedback fixups — DISCUSSIONS.md year + spelling sweep

Bot review on [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) (1m 58s pass) flagged two issues. Both addressed:

1. **Year fix on line 487** — `(Kwok et al., *Nature* 2024)` → `(Kwok et al., *Nature* 2025)`. Paper published online 2026-02-19; the 2024 likely captured pre-publication when the DOI was assigned. This was a pre-existing inconsistency, surfaced cleanly by the new subsection citing the same author with the correct year.

2. **British → American spelling sweep across `research/manuscript/DISCUSSIONS.md`** — user established the new American convention at this point (see [feedback_american_spelling.md](feedback_american_spelling.md) in role memory; index entry in `MEMORY.md`). 32 word changes: `tumour`/`Tumour`/`intratumoural` → `tumor`/`Tumor`/`intratumoral` (~20), `personalised`/`Personalised`/`personalising` → `personalized`/`Personalized`/`personalizing` (~9), `recognise(s)` → `recognize(s)`, `prioritising`/`prioritisation` → `prioritizing`/`prioritization`, `behaviour` → `behavior`, `signalling` → `signaling` (3 occurrences), `favour` → `favor`, `Modelling` → `Modeling`. The other manuscript files (`INTRODUCTION.md`, `METHODS.md`, `CONCLUSIONS.md`, `RESULTS.md`) were already mostly neutral so no sweep needed there in this PR; future American-style edits will harmonize naturally.

**Advisory notes from the bot (not fixed in this PR — flagged for follow-up):**

- `~789 public NJ pool` — no parenthetical source attribution. Could optionally soften to "hundreds of" if the exact count is from `dakwok/SSNIP` rather than the published paper. Not blocking.
- GTEx denominator `9,166` (Kwok et al.'s usage) vs the commonly-cited `~9,662` for GTEx v8. Sanity-check during [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) when querying GTEx v8 directly.
- Junction coordinates for NEJ<sub>GNAS</sub> / NEJ<sub>RPL22</sub> — bot suggests recording them in [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) body so the Snaptron query is self-contained without re-derivation. Will fold in when starting #299.

---

### 16:19 UTC — Editor: Scientist

#### [Sub-Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) DISCUSSION — public-vs-personalised splice neoantigen axis (Kwok et al.)

**Author-name correction.** Yesterday's frozen [`lab_notebook.md`](research/lab_notebook.md) entry at lines 23 + 25 attributes [the *Nature* 2025 paper](https://www.nature.com/articles/s41586-024-08552-0) to "Pan et al." — the correct first author is **Darwin W. Kwok**. Same paper, same DOI (`10.1038/s41586-024-08552-0`), Zotero key `5ZT8KC8X`. The stale attribution propagated into [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280)'s title + body and `research/news_log.md` line 23 before catching it during today's verification work. Per the lab notebook immutability rule (`shared/feedback_lab_notebook.md`), the frozen-snapshot entries are not editable; this note serves as the corrected reference for future lookups. Fixed today: [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) title + body (Pan → Kwok, plus expanded ACs covering the threshold-tradeoff scope), and `research/news_log.md` line 37 (rolled into this PR per yesterday's option-2 decision).

**Verification path that surfaced the correction.** The user pushed back on yesterday's first-pass framing — *"biologically I don't know if junctions that are tumor-recurrent across patients AND show up sporadically in normal tissue are really relevant as neoepitopes? Do Pan et al. show that such junctions result in potent neoepitopes?"* — driving a structured re-read of the paper. WebFetch on PubMed by DOI (39972144) returned first-author Darwin W. Kwok; Zotero search by `Kwok` returned the existing entry already in the collection. The correction triggered the immediate news_log + Issue body fixes today.

**Placement.** New 3rd subsection `### Kwok et al.: public neoepitopes from recurrent splicing — the off-the-shelf end of the axis` inserted into [DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md) `## Comparison to related neoantigen-prediction tools` after the existing SpliceMutr (line 557) and ENEO (line 585) subsections. Three paragraphs: (1) what Kwok et al. did — recurrent *GNAS* / *RPL22* splicing across glioma / mesothelioma / prostate / liver, spatially conserved *GNAS* neojunction across multi-region biopsies, TCRs isolated from healthy donors with HLA-dependent peptide-dose-dependent killing of endogenously expressing tumor cell lines (Fig 5, Extended Data Fig 9-10); (2) public-vs-personalised axis — Kwok-class for off-the-shelf shared TCR targets across patients with recurrent NJ + matching HLA, our pipeline for patient-private candidates covering the residual majority; (3) GTEx threshold tradeoff — Kwok `PSR_GTEx < 1%` vs our `min_read_count: 1` (= `PSR_GTEx = 0%`) at the population level, with verifiability of NeoA<sub>GNAS</sub> / NeoA<sub>RPL22</sub> specifically deferred to follow-up [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299).

**Threshold-tradeoff framing — option (a), population-level only.** The user asked the sharper question — *"Were the GNAS or RPL22 neoepitopes derived from splice junctions with non-zero (but <1%) GTEx normals? Because only then can we maybe start to justify the <1% GTEx normals threshold, right?"* — and we exhausted public sources looking for the answer. None of the published material (main text, Extended Data figures 1-10, Springer supplementary MOESM1-5, [`dakwok/SSNIP`](https://github.com/dakwok/SSNIP) GitHub repo, peer review rebuttal MOESM5) reports the specific PSR<sub>GTEx</sub> values for the validated NEJ<sub>GNAS</sub> and NEJ<sub>RPL22</sub>. Extended Data Fig. 1D scatter plots show the public NJ population (colored dots) clustered near `PSR_GTEx = 0` but at the 0–1 axis scale visually indistinguishable from `PSR_GTEx ≈ 0.005`. The DISCUSSION subsection therefore states the tradeoff at the population level (*"at least some of Kwok et al.'s ~789 public NJ pool — which by their inclusion criterion spans `0% ≤ PSR_GTEx < 1%` — would be excluded by this pipeline's filter"*) without claiming the specific GNAS/RPL22 verdict. Verification deferred to follow-up [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) (lightweight Snaptron query against GTEx v8, ~1–2 hours focused work, no pipeline rerun needed).

**Tooling note.** Added `openpyxl>=3.1` to [research/requirements.txt](research/requirements.txt) for the supplementary-table scan (verifying MOESM1/3/4 contents — none had the per-NJ table). Rolled into this PR.

**HTML-comment TODO marker.** Added `<!-- TODO(#299): re-derive PSR_GTEx ... -->` after the new subsection in [DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md) so [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299)'s outcome has a clear insertion point in the manuscript prose without cluttering the rendered output.

---

## 2026-05-06

### 13:30 UTC — Editor: Scientist

#### [Sub-Issue #268](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/268) INTRODUCTION — Rojas 2023 foundational personalized neoantigen vaccine reference

**Placement.** Single sentence appended to the existing clinical-translation context block in section `## Cancer Neoepitopes` of [INTRODUCTION.md](research/manuscript/INTRODUCTION.md), inserted after the line "...personalized cancer vaccines and adoptive T-cell therapies." That line was the natural anchor — the sole pre-existing clinical-translation statement in INTRODUCTION, previously uncited. Considered Section 4 (HLA Genotype) where Sahin et al. 2017 + Ott et al. 2017 already cluster, and rejected — those are cited there for the *vaccine slot count* argument (10–20 candidates per formulation), a different rhetorical use; piling Rojas in would muddle the load-bearing function of the Section 4 citation. Considered a new paragraph and rejected — one foundational reference doesn't justify its own paragraph; INTRODUCTION should anchor and DISCUSSION should develop.

**Framing.** Rojas et al., *Nature* 2023, 618:144-150 (Zotero `UIN9DIUP`) — Phase 1 autogene cevumeran trial in PDAC, 8/16 vaccine-induced T-cell responders, mRNA-LNP synthesis with up to 20 patient-specific neoantigens. The cited fact ("8 of 16 patients developed vaccine-induced T-cell responses") is the single cleanest summary of the trial's primary endpoint; no need to also cite the recurrence-free survival trend (covered in DISCUSSION via [PR #260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/260)) or the 6-year follow-up (AACR 2026, surfaced [2026-05-03 14:07 UTC] in news_log; deferred to a future DISCUSSION update if/when published as a full paper).

**Cross-section consistency.** INTRODUCTION now cites Rojas 2023 as foundational; DISCUSSION already develops the autogene cevumeran clinical context and adds the 6-year/long-lived T-cell durability follow-ups. Sahin et al. 2017 + Ott et al. 2017 remain in Section 4 only, doing different work (vaccine slot count). No duplication.

**Citation style.** Author-year inline `(Rojas et al., *Nature* 2023)` matches the existing INTRODUCTION convention used by Yewdell & Bennink 1999, Sahin et al. 2017, Ott et al. 2017, Cai et al. 2026 (added in [PR #282](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/282)), Alam et al. 2023. Reference-list finalisation deferred to [Sub-Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272).

---
