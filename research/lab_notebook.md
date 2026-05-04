# Lab Notebook

---

## 2026-05-04

### 18:08 UTC — Editor: Scientist

#### Issue #232 — clinical-translation DISCUSSION section landed ([PR #260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/260))

**Sunday read → Monday landing.** Sunday afternoon's deep read of [Sethna et al., *Nature* 2025](https://doi.org/10.1038/s41586-024-08508-4) (PDF supplied directly via conversation since bioRxiv-style sources block WebFetch — the autogene cevumeran 2025 follow-up) produced an initial draft of a new manuscript-DISCUSSION section. User caught a **polyvalency-framing overstatement** in the first draft: I had characterised the paper as "arguing for polyvalent vaccines," but the paper's actual data finding is **clonal pruning at recurrence** — the polyvalent-vs-high-potency suggestion is offered as *speculative future directions*, not a tested hypothesis. Corrected the framing (now explicitly attributes the suggestion to authors, presents both alternatives, and uses Option B for the closing claim — flags the splice-clonality question as an open empirical hypothesis rather than overclaiming pipeline-pipeline fit).

**Final section content.** Three paragraphs landed in `research/manuscript/DISCUSSIONS.md` between *Allele breadth and immunodominance* and *Immune-pathway gene neoepitopes: the presentation paradox*: (1) field-landscape opener (Iamukova & Alferova, APJCO 2026; Sethna 2025; Rojas 2023) → durability data (RFS, 7.7-yr clone lifespans, *de novo* priming); (2) two pipeline-relevant clinical observations (low-mutation tumor indication; clonal pruning + author-proposed future directions, with inline gloss for "subclinical clones"); (3) honest pipeline-fit framing (pool expansion is solid; clonal-shared-vs-subclonal for splice antigens is open empirical question).

**PR #260 opened with `Refs #232` (not `Closes`)** — closes one AC of #232 only; the Issue stays open as the manuscript-integration umbrella with 6 remaining ACs. [Multi-PR coordination comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232#issuecomment-4373317088) posted on #232 announcing the lazy-sibling-Issue plan for the rest. Scope-expansion comment also posted earlier ([here](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232#issuecomment-4370397677)) folding the clinical-translation arc into #232's body.

#### Memory + workflow rules captured

**Comment-first then body-edit for scope changes** (saved to `shared/feedback_scope_discipline.md`). Body-alone is silent (GitHub edit history is buried); comment-alone is stale (body misleads). Comment-then-body lets reviewers intervene before any PR ships, then keeps the body as current state of truth. Asymmetry note: body-edit alone is fine for minor wording / clarification changes that don't change scope — comment-first is specifically for *scope* changes.

**AC-checkbox heuristic for upfront sub-issue decomposition** (saved to `shared/feedback_branch_planning.md`). When an Issue's acceptance criteria touch distinct sections / files / themes, decompose into sub-issues *before* opening any branch on the parent — AC checkboxes are ready-to-use sub-issue scopes. The mid-flight restructuring we hit on #232 today (1 AC landing, 6 ACs pending, branch already created against parent) is the exact failure mode this rule prevents. For #232 itself, Option Y (single Issue + multi-PR with `Refs`) is the cleanest path forward only because the Issue is already in flight; for new Issues going forward, sub-issues from the start are the right pattern.

### 11:00 UTC — Editor: Developer

#### Morning routine — three threads

**Standup slip caught by user.** I claimed "no Pending messages to Developer requiring new action" based on a truncated system-reminder snapshot of `team_standup.md`. User pushed back ("did you actually check?"); I re-read the file explicitly and found PM had posted a **2026-05-04 10:06 UTC** message asking me to action the closure ritual on [Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223), plus a **10:12 UTC** new shared rule. The slip surfaced a real failure mode — system-reminder snapshots are partial diffs, not full state, and I'd been pattern-matching them as authoritative. PM codified the fix as `shared/feedback_read_before_claiming.md` (promoted to MEMORY.md Always-in-effect), so the rule will load via `/cerebrum` next session.

**Closure ritual on [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223).** Edited the body to tick 5 of 6 acceptance-criteria boxes; box 2 (patient_001 variant context) **comment-deferred to [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)** because the [PR #254](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/254) spike used a synthetic A→G SNV at the chr22:42M region midpoint — chosen for API plumbing — rather than an actual patient_001 germline variant. The literal spec wasn't met, but the gating purpose (verify `predict_variant` API path works end-to-end) was. #224's Experiment 2 covers the real test naturally via patient FASTA + `predict_sequence`. Full reasoning in the [comment on #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223#issuecomment-4370281435).

**NeoGuider eval issue opened.** Engineer news briefing surfaced [NeoGuider (XuegongLab)](https://github.com/XuegongLab/neoguider) — end-to-end ML neoepitope ranking pipeline that handles splice variants alongside SNV/indel/fusion. Direct architectural peer to ours. 459 commits, 5 releases — active. Opened [Issue #258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/258) for Scientist eval, joining the existing eval family (#218 HERMES, #201 ImmSET, #236 hybrid models, #188 Boltz-2, #222 splice2neo). Logged in `news_log.md` (this PR) so it doesn't re-surface in future briefings. Pipeline-fit modes pre-framed (triage / replacement / cross-check / component reuse) using the same vocabulary as #236, but specific scoping deferred to Scientist.

**Process notes:**

- News-log branch follows the PM-mandated convention from yesterday: `docs/developer/news-log-YYYY-MM-DD-HHMM` (today: `docs/developer/news-log-2026-05-04-1059`). The local "morning news" step name doesn't leak into the cross-role branch name.
- PR Status will be flipped to "Ready for review" immediately after `gh pr create` per yesterday's PM ask now in MEMORY.md Always-in-effect — first PR where the rule applies.

---

## 2026-05-03

### 18:12 UTC — Editor: Developer

#### Issue #223 / PR #254 — review-cycle fixes

`@claude` review on the AlphaGenome primer flagged three actionable items. All addressed in this commit:

1. **Medium — "powers of 2" phrasing was misleading.** The primer (line 230) and the spike script (line 70) both said the four required input lengths were "powers of 2 ≤ 1Mb." Reviewer correctly pointed out that 2^15=32768, 2^16=65536, and 2^18=262144 are also powers of 2 but **not accepted** — only the four specific values 16384, 131072, 524288, 1048576 work. A reader taking the primer literally would expect the intermediate sizes to work and get a `ValueError`. Fixed both call sites to say "these four specific values" and explicitly list the rejected powers of 2.

2. **Low — lenient-threshold caveat.** The primer described the spike's mean |ref − alt| of ~0.00014 as evidence the variant-affected-junction filter is "very lenient." Reviewer pointed out the spike SNV was a synthetic A→G chosen for API plumbing — not designed to disrupt splicing. So the tiny mean delta likely reflects micro-perturbations from a single nucleotide change, not characterisation of how the threshold behaves for biologically meaningful splice-disruptive variants (where deltas would be orders of magnitude larger). Added the caveat explicitly. Also added a one-line acknowledgment that geographic windowing wasn't ruled out by our test (the non-zero-delta rate just makes effect-based filtering the more parsimonious explanation).

3. **Low — dead branch in `summarise_output()`.** The function used `getattr(val, "shape", None) is None` as a "this is a DataFrame, fall through to column listing" heuristic. But pandas DataFrames also have `.shape` (returns `(n_rows, n_cols)` tuple), so the column-listing else-branch was unreachable for `metadata` — instead the script printed `"metadata.shape: (367, 8)"`, which is less informative than the column inventory the comment promised. Fixed by handling `metadata` via its own explicit code path (always print columns) and looping the shape-printing only over `values + interval`. Verified with a mock-DataFrame smoke test.

**Why this matters now:** the spike script is intended as a reusable probe for [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)'s AlphaGenome validation notebook. Shipping with dead branches and misleading comments would be a paper-cut for whoever picks that issue up.

**Verified:** smoke test confirms `summarise_output()` now prints the column list for DataFrames; primer + spike consistent on the input-length constraint. No tests changed (script is a one-off probe, not under pytest coverage). Closes [Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) when this PR merges.

### 17:49 UTC — Editor: PM

#### Issue #245 — memory duplicate-check at write time (first pm-i1 same-day ship after #246)

Implemented as `shared/feedback_memory_duplicate_check.md` (not `pm/...` per the issue's original AC — shared/ is correct because the rule applies to all roles writing memories). Defines the overlap heuristic (3 signals: name field, description field, keyword/domain), the prompt format (append/replace/proceed-anyway, default append), and exemptions (memories from current session, MEMORY.md itself, reference memories pointing at external systems).

**Meta-test passed**: ran the to-be-codified rule on its own write — scanned shared/MEMORY.md for prior overlap with the proposed memory; closest candidate was `feedback_memory_escalation.md` (escalating an existing rule on repeat correction) which shares the *memory hygiene* domain but has different trigger and action. No overlap, proceeded.

**Smoke test passed**: hypothetical near-duplicate `feedback_cerebrum_project_separation.md` (description: "Cerebrum framework vs project — keep them distinct in scoping decisions") would hit 3-of-3 signals against existing `feedback_cerebrum_vs_project.md`. Heuristic correctly flags.

**Indexed in shared `MEMORY.md` Always-in-effect** as a one-liner pointing at the full rule. Sister rule to `feedback_closure_ritual.md` — together they bracket memory hygiene at write (this) and at close (ritual).

**Sequencing note for the rest of pm-i1**: with #245 + #246 both shipped, the remaining 3 (#244 ask-for-help, #247 capacity recheck, #243 Rulesets) can proceed in any order. P1s (#244, #247) before P2 (#243).

### 15:26 UTC — Editor: Developer

#### Issue #214 / PR #240 — round-3 fix: close the report.tsv funnel

`@claude` round-2 review verified all four round-1 items as closed and approved correctness. One non-blocking design observation: `report.tsv` consumers (RESULTS.md authors) saw a non-closed funnel because the three patient-level totals (`junctions_extracted_total`, `junctions_annotated_discarded`, `junctions_unannotated_total`) didn't include `mean_reads_filtered` — that category only lived in `junction_filter_stats.tsv`, so anyone reading just `report.tsv` would see `extracted_total ≠ annotated_discarded + unannotated_total` with an unexplained gap.

**Fix:** added `junctions_mean_reads_filtered` as a fourth `junction_filtering` row in `_build_report_tsv` (notes: "all tumor samples"). The four rows now reconcile arithmetically: `extracted_total = mean_reads_filtered + annotated_discarded + unannotated_total`. New test `test_junction_funnel_totals_reconcile_in_report_tsv` enforces this invariant on `report.tsv` (separate from the equivalent invariant on `junction_filter_stats.tsv` from round 1). Updated existing tests to include `mean_reads_filtered` rows in their fixture stats data and assert on the new row's value.

**Why not done in round 1/2:** Issue #214's spec named exactly three metrics (extracted/annotated/unannotated). Round 1 implemented exactly what the spec said; round 2 fixed reconciliation in the upstream stats artefact. Round 3 closes the gap between "internal artefact reconciles" and "report.tsv consumers see a closed funnel" — reviewer flagged this in round 2 as non-blocking but worth a conscious decision.

**Verified:** 98/98 tests pass; Snakemake dry-run clean; rule graph unchanged.

### 14:05 UTC — Editor: Developer

#### Morning routine — standup cleanup + news briefing

Standup cleanup: deleted 4 own Done messages from 2026-04-29 (>3 days, per the standup file rule). One Issue #79 scope-expansion notification + three follow-ups to PM. The durable record stays in commits / lab notebook / project board; standup is for active conversation only.

News briefing surfaced three pipeline-relevant items, all logged in [`research/news_log.md`](news_log.md) so they don't re-surface in future morning briefings (gap caught when I re-surfaced the `googlebatch` executor today even though it was already tracked in [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66)):

1. **PyTorch 2.7 + CUDA 12.6 = last P100-supporting combo.** PyTorch's [dev-discuss thread](https://dev-discuss.pytorch.org/t/cuda-toolkit-version-and-architecture-support-update-maxwell-and-pascal-architecture-support-removed-in-cuda-12-8-and-12-9-builds/3128) describes Maxwell/Pascal/Volta as "feature-complete with no further enhancements planned" — 2.8 dropped Pascal kernels in cu128/cu129 builds. We already pin `torch>=2.0,<2.5` in `python.yaml`; this confirms the pin is **permanent on this hardware**, not a temporary workaround. Updated CLAUDE.md's `python.yaml — PyTorch SM 6.0 / P100 compatibility` section with this detail (this PR).
2. **Snakemake `googlebatch` executor plugin** — already covered by [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) (the issue body literally links the plugin catalog page). User correctly flagged the re-surface; logged it in `news_log.md` so we don't trip on it again.
3. **Snakemake 9.x deprecates `--use-conda`** in favour of `--software-deployment-method conda`. Affects every snakemake call site in CLAUDE.md, `run_cloud_gpu.sh`, `setup_local.sh`, and possibly internal docs. Posted a [comment on Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200#issuecomment-4366348544) listing affected files — bundle into the 9.x migration PR rather than fixing separately, since 8.x still accepts both forms.

Memory rule reinforced (didn't have to add): morning-routine news scan must check `news_log.md` first — items already logged shouldn't re-surface unless there's a meaningful update. The googlebatch slip happened because the original Issue #66 was opened before `news_log.md` existed (introduced 2026-05-01); now-logged → won't repeat.

---

## 2026-05-02

### 20:32 UTC — Editor: Developer

#### Issue #223 — live AlphaGenome API spike + primer doc

User registered for the AlphaGenome API key in the same session, so the live test moved up from "next session" to now. Two-call probe in [`scripts/alphagenome_spike.py`](../scripts/alphagenome_spike.py) (committed under `scripts/` for [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) reuse): one `predict_interval` and one `predict_variant` against a 131kb chr22 region. SDK is `pip install alphagenome` from PyPI; key loaded from `.env` as `ALPHAGENOME_API_KEY` matching the `zotero_add.py` `load_env` pattern.

**Empirical findings, beyond the morning recon:**

- **Latency**: 0.75–2.6s per call for a 131kb interval. Within the documented "<1s prediction" range (the rest is round-trip).
- **Required input lengths are constrained to {16384, 131072, 524288, 1048576}** — not arbitrary. First call with 100k failed with a `ValueError`. Plumbing for [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) needs to handle this. Documented in the primer.
- **367 splice-junction tracks** total: 313 ENCODE + 54 GTEx; 195 total RNA-seq + 172 polyA-plus RNA-seq. The 54 GTEx tracks are the relevant subset for the predicted-normal filter ([#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)) since GTEx profiled healthy human tissues; ENCODE includes cancer-derived cell lines (HeLa/K562) that are poor "normal patient tissue" proxies.
- **Output shape is `(n_junctions, n_tracks)`** — a matrix where junctions are SHARED across tracks and tissue specificity lives in the score distribution within each column. This was a major mental-model fix for me; I'd been imagining one set of junctions per track. Documented in the primer as a load-bearing concept.
- **`predict_variant` filters output to variant-affected junctions**, empirically verified: `predict_interval` returned 6084 junctions, `predict_variant` returned 503; 503/503 of those have ref ≠ alt in at least one track (none have identical values across all 367). User correctly hypothesised this was an effect-based filter (not the geographic windowing I'd guessed). Inclusion threshold appears very lenient — mean |ref − alt| of ~0.00014 still qualifies. So the 503 number is "any model-detectable variant response" rather than "biologically meaningful effect" — meaningful signal for [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) lives in delta magnitudes, not the count.
- **`predict_variants` (plural) is batched single-variant queries, not a multi-variant compositor** — the docstring confirms "Variant outputs for each DNA interval and variant pair." Each variant is processed independently, returning N separate `VariantOutput`s with `max_workers` parallelism. So patient_001's combined germline-variant effect (Experiment 2) needs externally-built patient FASTA + `predict_sequence`, not `predict_variants`. The latter is the right tool for per-variant attribution (Experiment 3).

**Documentation: [`docs/alphagenome_primer.md`](../docs/alphagenome_primer.md)** ships in this PR. ~290 lines covering tracks (the central concept), output shape, API methods as orthogonal query shapes (not different assays), `data_source` as training-time provenance, per-track interpretation, operational details, and the predicted-normal filter sketch. Iterated through several mental-model corrections with the user — the doc preserves the accurate framing rather than my initial half-right intuitions. Originally drafted with a "per-experiment API mapping" section; user pointed out [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)'s Experiment 1 design is Sci's call and shouldn't be encoded in a Dev primer — section dropped; suggestion to refine Exp 1 (use GRCh38 + annotated-only ground truth instead of patient WGS + observed normal) will be posted separately on [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) for Scientist.

**Aggregation choice for the predicted-normal filter:** `max(axis=1)` over the 54 GTEx tracks — captures "if ANY healthy tissue shows signal, treat as normal" — aligned with Scientist's vaccine-CTL safety framing in [#211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211)/[#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212). Mean would smooth tissue-specific normal signals toward zero and miss them — wrong direction for safety. Max is sensitive to single-track noise; top-k mean / quantile threshold are noise-robust alternatives Sci can pick if needed. Documented.

**Closes [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223)** when this PR merges. [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) and [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) unblock for Scientist.

### 19:53 UTC — Editor: PM

#### Issue #246 — implemented closure audit step in morning warmup + tightened "acceptance criteria visibly met"

Picked the coolest pm-i1 P1 from this morning's slate to implement same-day. Added "Step 0.5 — Closure audit" to `pm/feedback_morning_warmup.md`: per-issue audit checklist (milestone, acceptance criteria, priority rationale consistency, sub-issue completeness, lab notebook entry), cross-role consistency checks, comment-not-reopen feedback loop. New section header convention `## 🔍 Closure audit` (rendered only when issues actually closed in the past 24h).

**Smoke test on 9 recent closures surfaced 2 immediate fixes** — applied: backfilled milestone for #216 (closed 2026-04-30 with no milestone) and moved #234 from `i2 - S1 - Tool Landscape Evaluations` to `pm-i2 - PM Self-Improvement Tooling` (mis-milestoned pre-`pm-i<N>`-axis-convention).

**Bigger finding came from user pushback** — they asked why acceptance-criteria checkboxes are so often unchecked. Re-running the audit with tightened criteria (operationalising "acceptance criteria visibly met" as `- [x]` boxes ticked, not just "PR exists") flipped the verdict on 6 of the 9 closures: 23 unchecked boxes across #79, #190, #191, #197, #219, #221. Even today's #235 (closed by me this morning) had plain bullets, not checkboxes — so the boxes weren't even tickable. Backfilled #235's body: converted bullets to `- [x]` checkboxes with closure-attestation note.

**Codified two related rules:**

1. **`shared/feedback_closure_ritual.md`** — before closing any issue, tick `- [x]` every met acceptance-criteria box, OR comment-defer with a follow-up link. Closing with unchecked boxes ≠ closing. Applies to all roles, not just PM. PM enforces via morning audit; each role self-enforces.
2. **Companion rule in same memory** — every new issue's acceptance-criteria section MUST use `- [ ]` checkbox format from creation, not plain `-` bullets. Plain bullets make the ritual unenforceable.

**Why this lands deeper than expected:** the original closure-audit checklist used "linked PR merged" as the proxy for "criteria met". That's a weak signal — every closed issue has a PR, so the check always passes. The tightened version asks the strong question: "are the boxes ticked?" — which immediately exposed that 6/9 recent closures skipped the ritual entirely. The smoke test caught a real systematic gap on its first run, which is exactly what the audit is for.

**Pre-ritual closures (before 2026-05-02) are intentionally not retroactively ticked** — that requires verifying each criterion was actually met, which is dangerous to fake. They're the baseline. The ritual enforces forward.

**Process meta-note:** this was the first PM-meta-work issue completed via the new `pm-i<N>` axis. Ships entirely as memory edits (no project-repo PR needed) plus the lab notebook entry. Smoke-tested same-day; full morning routine using the new step happens tomorrow.

---

### 18:46 UTC — Editor: Scientist

#### Issue #236 — both candidates characterised, first deliverable complete

Picked up [#236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) and worked through inference-mode classification per candidate.

**TCRLens — mode (b)** via WebFetch on the [paper's Methods section](https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag066/8496266). Requires 3D structures generated via the tFold-TCR pipeline. tFold-TCR is **5000× faster than AlphaFold-Multimer** (eliminates MSA via ESM-PPI-TCR), so it qualifies as a genuine fast-structure-predictor-inline. TCRLens is structure-source-agnostic — viable as both **triage** (using tFold-TCR upstream of MHCflurry → TCRdock) AND **cross-check** (reusing our TCRdock outputs). The only candidate so far that can sit upstream of TCRdock.

**t2pmhc — mode (c)**, corrected from yesterday's abstract-only claim to a properly-cited methods-section read. WebFetch blocked on bioRxiv (403), so user supplied the PDF directly in conversation. Methods § Structure Prediction confirms: *"Structures of all TCR-pMHC complexes were predicted using TCRdock (v2.0.0)."* **t2pmhc uses TCRdock specifically** — the same tool we already run — so cross-check is essentially free in our pipeline (only the GNN scoring step adds). Pipeline-fit: cross-check only (mode c rules out triage and replacement). The discussion section also flags structure-prediction quality as the accuracy ceiling, which we'd inherit equally with them.

Consolidated [comment on #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236#issuecomment-4364300411) with full pipeline-fit table. Second deliverable (pipeline-fit recommendation with scoped follow-ups) deferred until [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) lands patient_001 outputs — the natural shared test bed.

#### Workflow revision — WebFetch envelope and PDF policy

User pushback clarified the actual WebFetch envelope: body text yes, table content yes, figure captions yes, but figure images, layout-dependent content, and supplementaries — no. For deep analysis on eval candidates (benchmarking plots, structural diagrams, supplementary tables), text-only WebFetch is genuinely insufficient and I'd undersold this earlier. Extended `shared/feedback_zotero_defer_inaccessible.md` with an explicit *"always ask for PDF when deeper analysis is needed"* rule. Future flow: surface the WebFetch limitation, ask user to attach the PDF to the Zotero entry, then I fetch via `/items/{key}/file` or use the user-attached PDF in conversation.

#### Zotero hygiene — bioRxiv replacements + HERMES correction

User replaced the three bioRxiv entries that had been added yesterday/this morning via direct API workaround. New keys: `78BQ23IV` (Benchmarking foundation models for splice site and exon annotation), `STQYEVAQ` (t2pmhc, replacing deleted `E3WRMAAH`), `4N2J7SIH` (AI predicted TCR-pMHC structures — first time in Zotero). All three now have PDFs attached. I re-applied project tags via API: `manuscript-METHODS` + `manuscript-DISCUSSION` for `STQYEVAQ`; `manuscript-DISCUSSION` for the other two.

**HERMES correction:** I claimed twice today that HERMES was missing from Zotero based on a title-substring search returning zero hits. Wrong both times. HERMES is in Zotero as `MWZFINV6` — the paper title is *"T cell receptor specificity landscape revealed through de novo peptide design"* (Visani et al., *PNAS* Oct 2025); **HERMES is the method name** introduced in this paper, but the title doesn't contain the word "HERMES" so the title-based search missed it. [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218)'s body correctly references `MWZFINV6` — the bug was in my recall, not the issue. Tagged `manuscript-DISCUSSION` to position it alongside t2pmhc/TCRLens in the cross-check role.

#### Standup — Developer cleared the AlphaGenome gate

[Pinged Developer at 16:20 UTC](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) asking them to prioritise [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) (gates P1 Sci work in milestone #17, due 2026-05-22). Developer replied at 18:18 UTC: AlphaGenome is **$0 for non-commercial use** → Low Cost branch hit per #203's decision tree → [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) and [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) unblock as soon as Developer finishes the live API test next session. Cohort-expansion bonus is on the table once validation results come in.

### 18:18 UTC — Editor: Developer

#### Issue #223 — AlphaGenome API spike: recon-only pass

Picked up Scientist's [16:20 UTC standup ask](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues) (gates [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) + [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) in milestone #17, due 2026-05-22). Recon-only pass without an API key — answered the gating questions for downstream Sci work, deferred live API test to the next session. Full findings on [#223 comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223#issuecomment-4364440234); summary here:

**Cost: $0** for non-commercial use. Cost decision branch in [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) → **Low cost** → sub-issues unblock; cohort-expansion bonus is on the table when Scientist sees validation results.

**Junction connectivity confirmed.** [`dna_client.py`](https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/models/dna_client.py) exposes `OUTPUT_TYPE_SPLICE_JUNCTIONS` as a first-class output type — distinct from `OUTPUT_TYPE_SPLICE_SITES` (per-position probabilities) and `OUTPUT_TYPE_SPLICE_SITE_USAGE`. So the API returns explicit donor → acceptor pairs per tissue, not just per-position scores. This was [#223 task #3](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) — the make-or-break verification — and it's a clean ✅. Without it the whole AlphaGenome arc wouldn't have made sense.

**Rate limits intentionally undocumented and dynamic** — the team's posture per [community FAQ](https://www.alphagenomecommunity.com/t/what-are-the-alphagenome-api-limits/673) is "increase parallel workers until `RESOURCE_EXHAUSTED`, then back off." Default 5 workers, < 1s per prediction, 1Mb max sequence per call, single variant per request. Stated comfort zone is "thousands of predictions" — millions is where it breaks down. patient_001 across Experiments 1–3 lands solidly in the well-suited band.

**SDK is open-source** ([github.com/google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome), since 2026-01-28). Standard pattern is `model = dna_client.create(API_KEY); model.predict_variant(interval=..., variant=...)` with `genome.Interval` and `genome.Variant` value objects. Tissue specified via UBERON ontology terms; concrete track inventory not enumerated in the README — needs a key to query.

**Why this is enough to unblock #224/#225:** the two gating questions Scientist's design depended on were "does this cost real money?" and "does it return junction connectivity?" Both have hard answers now. The remaining recon items (live latency, exact tissue track count, output schema details) inform implementation but don't change the experiment design.

**Next step:** user is registering for an API key now. Next session does the live test (one chr22 region, GRCh38 vs patient_001 variant), confirms output schema, then closes [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223). Hand-off to Scientist sub-issues happens at close.

**Process note:** entry ships on its own time-suffixed branch (`docs/developer/lab-notebook-2026-05-02-1818`) per the multi-session pattern in `shared/feedback_lab_notebook.md` — [#223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) stays open across sessions, so the notebook entry is independent of issue lifecycle.

### 16:48 UTC — Editor: Developer

#### Issue #214 / PR #240 — review-cycle fixes

`@claude review` on [PR #240](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/240) caught a real arithmetic bug in this morning's funnel design: `junctions_raw` was captured **before** the per-sample mean-reads filter, but `annotated_discarded` / `normal_shared` / `tumor_exclusive` were counted **after** — so the three patient-level totals didn't reconcile (`extracted_total ≠ annotated + unannotated`). The gap is exactly the noise junctions removed by `reads > mean`. My morning test even demonstrated it (4 raw → 3 classified, gap of 1) and I missed it.

**Fix (Option A from reviewer):** added a 4th intermediate category `mean_reads_filtered` to `junction_filter_stats.tsv` so it now records the full 5-step funnel (`junctions_raw → mean_reads_filtered → {annotated_discarded, normal_shared, tumor_exclusive}`). The 3 patient-level rows in `report.tsv` are unchanged in name but the underlying schema now reconciles arithmetically. Added an explicit test `test_funnel_reconciles_arithmetically` asserting `junctions_raw == sum(4 downstream buckets)`.

**Style nits also addressed in the same commit:**

- Removed a defensive `try/except Exception` around `pd.read_csv(stats_tsv)` in `_build_report_tsv` — `junction_filter_stats.tsv` is a *required* Snakemake input, so the error path can't fire under normal pipeline runs. Per CLAUDE.md "don't validate scenarios that can't happen."
- Removed a redundant outer `int()` cast on `sum(...)` in the same aggregation.
- Switched `getattr(snakemake.input, "junction_filter_stats", None)` to direct attribute access `snakemake.input.junction_filter_stats`. The `getattr(...,None)` pattern is reserved for genuinely optional inputs (`hla_qc`, `pdb`, `scores_tsv`); using it for a required input misleads the reader about optionality.

**Lab notebook process change captured as a memory:** user pushed back when I tried to *rewrite* the morning 10:48 entry to be a single self-contained "complete" account. Saved `shared/feedback_lab_notebook.md` "Entries are immutable" rule: each session writes its own entry, never edit/replace previously-committed entries — the lab notebook is a journal not a wiki, and the bug-catch narrative is exactly what's most useful for future-grep ("when did we add `mean_reads_filtered`?" → afternoon entry, not a retroactively-edited morning one).

**Verified:** 234/234 tests passing (5 in `TestClassifyJunctionsStats` now, up from 4); Snakemake dry-run clean; rule graph correctly triggers `filter_junctions` from the code change.

### 16:07 UTC — Editor: PM

#### Issue #235 skim — completed the Anthropic 2026 Agentic Coding Trends Report cross-checks

Pulled the PDF, extracted text + diagrams (text-only extraction missed the page-8 multi-agent architecture diagram, which turned out to be the densest content). Findings live in the [#235 comments](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235); two threads:

1. **Main writeup** — 5 cross-checks (4/5 convergent), 3-layer orchestration framing (human → agent-orchestrator → specialists), per-cross-check evidence table.
2. **Domain-bespoke follow-up** — sharpened cross-check #1: the report is *Coding Trends*, ours is *Bioinformatics Trends*; the orchestrator+specialists pattern is convergent but the specialist roster is domain-bespoke (we have Sci because we have scientific reasoning, not because we copied a template).

**Three framing insights captured as shared memories** (more durable than this entry):

- `feedback_cerebrum_vs_project.md` — Cerebrum is the meta-framework above all projects/roles; this project is one instance using it. Don't over-scope project-local issues to "Cerebrum".
- `feedback_multi_role_not_multi_agent.md` — Prefer "multi-role workflow" over "multi-agent" when describing this project. "Multi-agent" implies agent-to-agent autonomy we don't have (user is the message bus). Calling it multi-agent overclaims.
- `feedback_domain_bespoke_roles.md` — Our PM/Sci/Dev split is adapted to bioinformatics research. The orchestrator+specialists *pattern* is convergent across domains; the specific *roster* is domain-bespoke.

**Future direction (parking lot, no issue):** partial autonomy via cron jobs / scheduled agents (e.g. PM doing autonomous overnight triage, archive sweeps, scheduled lit reviews). Worth revisiting when a specific bottleneck warrants it.

**Process note:** this entry ships on its own time-suffixed branch (`docs/pm/lab-notebook-2026-05-02-1607`) per the convention introduced morning-of for journal-style entries that aren't tied to issue lifecycle. #235 is already closed; the lab notebook is a journal, not a deliverable, and shouldn't depend on issue state.

### 10:48 UTC — Editor: Scientist

#### Morning routine — AI-predicted TCR-pMHC structures paper

Surfaced [*"AI predicted TCR-pMHC structures differentiate immune interactions"*](https://www.biorxiv.org/content/10.64898/2026.02.24.707744v1) (bioRxiv 2026-02). Three findings relevant to us: AlphaFold2 most consistent among AI tools for TCR-pMHC multimer prediction; **structural features outperform sequence features** for binding discrimination; non-binders produce less stable conformations under MD simulation. Reinforces the structural step at the top of our funnel and gives independent evidence that energy/stability signals are meaningful — directly relevant to the [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES eval.

User discussion converged on framing the field as **hierarchical convergence**, not a sequence→structure paradigm shift: sequence models scale, structural models discriminate hard cases, and hybrids (structure-informed sequence models) are the emerging best-of-both. Our pipeline already runs MHCflurry (sequence-aware) → TCRdock (full structural) — the new finding doesn't argue against this, it argues *for* keeping a structural step.

#### Issue #232 — S7 manuscript lit review opened

PM's standalone S7 milestone (#16, `i3 - S7 - Publication - Splice Neoantigen Tooling Landscape (Lit Review)`) was greenlit yesterday for the 5 papers accumulated in the [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) lit review session. Opened [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232) (P2, S) with full paper inventory: Zotero keys + DOI links + per-section mapping (CNNeoPP→INTRODUCTION; splice2neo+AlphaGenome→METHODS; ENEO+SpliceMutr+AlphaGenome→DISCUSSION). PM's suggested title was INTRODUCTION+DISCUSSION only — corrected to INTRODUCTION, METHODS+DISCUSSION since the actual Zotero tags showed METHODS coverage too. AlphaGenome's DISCUSSION entry depends on #223/#224/#225 outcomes; other entries unblocked.

#### Issue #236 — Hybrid TCR-pMHC scoring eval (t2pmhc, TCRLens)

Today's lit search surfaced two strong hybrid candidates that fill a real gap in our TCR-scoring landscape (we had pure-sequence [#201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET, structure replacement [#188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2, structure confidence [#218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES, but no hybrid). Opened [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) (P2, S) with explicit framing of three inference-time modes (a/b/c: sequence-only / sequence+fast-structure / sequence+full-AF2) and three pipeline-fit roles (triage / replacement / confidence cross-check). The pipeline-fit decision is a *consequence* of the inference-mode characterisation, not a starting assumption — issue's first deliverable is mode classification per candidate.

t2pmhc abstract review (Zotero key `E3WRMAAH`) revealed it falls in **mode (c)** — uses predicted full TCR-pMHC complex structures as input. This narrows its pipeline-fit to confidence cross-check only (peer to HERMES); triage and replacement are off the table. TCRLens (Zotero key `272HU8FV`) inference mode still pending — paper read deferred to issue execution. Both papers added to Zotero with `manuscript-METHODS` + `manuscript-DISCUSSION` tags.

bioRxiv DOI crash in `zotero_add.py` reproduced for t2pmhc (Issue #229 — empty `container-title` from CrossRef). Worked around via direct Zotero API call.

#### News_log first Scientist entry under new convention

PM merged [PR #238](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/238) this morning establishing the `### HH:MM UTC — Editor: <Role>` sub-heading convention for `news_log.md`. First Scientist entry under that format went out as [PR #239](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/239) (merged 10:44 UTC) covering today's three TCR-pMHC scoring papers. PR hit a conflict with the parallel PM merge (both branches added a `## 2026-05-02` section); resolved by interleaving sub-headings newest-first under the shared date — Scientist 10:04 UTC above PM 09:41 UTC.

#### Memory adjustments — parking-lot leak + Zotero API reference

User pushed back on a single-line glossary entry triggering a full PR; saved `shared/feedback_batch_trivial_docs.md` codifying batch-then-flush for trivial docs. **Initial design leaked state**: I added a "Currently parked" section inside the memory file with the actual pending entries — PM read it cross-role via the `shared/` symlink and was confused. Corrected: memory files are for rules, working tree is for state. Pending entries (today's MD glossary entry) now live as uncommitted edits on a docs branch; the memory file just points there.

Also saved `reference_zotero_api.md` after an avoidable detour hunting for the `.env` location during Zotero tag fetches. Direct-query pattern (collection `Z38GTJNW`, `urllib.request` with `Zotero-API-Key` header) now documented for next time.

#### Standup status

Own message [`2026-05-01 13:57 UTC`] (S7 milestone request) marked Done — PM replied 14:02 UTC and created milestone #16 (now backing #232). No outstanding Pending messages addressed to Scientist.

### 10:48 UTC — Editor: Developer

#### Issue #214 — junction-funnel totals in report.tsv

Fast-ship slice of [Issue #104](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104), unblocks Scientist's RESULTS.md. Adds three patient-level rows to `report.tsv` under stage `junction_filtering`: `junctions_extracted_total`, `junctions_annotated_discarded`, `junctions_unannotated_total` (all sums across tumor samples; normals omitted by design — see brainstorming below).

**Architecture:** new artefact `junction_filter_stats.tsv` per patient, written by `filter_junctions.py` alongside the existing `novel_junctions.tsv`. Long-format schema (`sample_id, sample_type, category, count`) where `category ∈ {junctions_raw, annotated_discarded, normal_shared, tumor_exclusive}`. `generate_report.py` reads this file (declared as new input on both `generate_report` and `generate_report_with_structure` via the shared `_generate_report_input` helper), aggregates by category, and emits the 3 patient-level totals. The original per-sample rows are preserved untouched.

**Why a new artefact rather than re-reading BED:** `filter_junctions.py` already computes `n_annotated`, `n_normal_shared`, `n_tumor_exclusive` per tumor sample — they were just being logged and thrown away. Persisting them is a 5-line change. The alternative (re-reading BED in `generate_report` and re-classifying junctions) would duplicate the entire reference/normal-set classification logic.

**Why tumor-only:** the spec row `junctions_unannotated_total = normal_shared + tumor_exclusive` is a tumor-classification concept — only tumor samples get split this way (normal samples are used as a filter set, not classified). User confirmed (a) tumor-only over (b)/(c) mixed scopes; flagged that normal-sample `junctions_extracted_total` could be useful in future via a separate Issue if needed.

**Naming alignment with [Issue #215](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/215):** picked the long-format schema and category names (`junctions_raw`, `annotated_discarded`, etc.) to match #215's planned `filtering_stats.tsv`. When #215 ships, it can either supersede `junction_filter_stats.tsv` with the unified multi-step file, or aggregate it as one source among many. Migration cost stays small either way.

**Tests:** 4 new tests in `test_filter_junctions.py` (schema, normal-omission, multi-sample, back-compat for callers without stats path) + 3 in `test_generate_report.py` (totals correctness, back-compat, notes-field convention). 233/233 total. Snakemake dry-run validates the rule graph: `filter_junctions` correctly triggers with the new output, downstream rules unchanged.

### 10:06 UTC — Editor: PM

#### Morning routine — introduced PM news as Step 0

PM has been doing morning warmups (board recap, standup, triage) but never surfacing external context. Added Step 0 — PM news (web search before the board recap), scoped to GitHub Projects/Issues updates, PM tooling, methodology shifts. Each item logged to `research/news_log.md` to dedupe across roles. Updated `pm/feedback_morning_warmup.md` to codify the step + the new `## 📰 PM news` section header.

Test-ran on first day with three queries (GitHub Projects updates, PM tooling news, software engineering methodology trends). Filtered listicles; kept items with concrete signal. Three made it through → see `research/news_log.md` for the full log.

#### Issue #234 — GitHub MCP Server eval (XS, P2)

Today's news surfaced that the [official GitHub MCP server](https://github.blog/changelog/2026-01-28-github-mcp-server-new-projects-tools-oauth-scope-filtering-and-new-features/) (Jan 2026) exposes typed `project_v2` mutation tools (Status, Priority, Size, Target date) at lower token cost than raw `gh api graphql`. Today PM uses hand-rolled GraphQL with hardcoded field IDs (e.g. `PVTSSF_lAHOB17eGc4BSomPzhAHGh8`). Every triage / re-arrangement burns context on boilerplate. [Issue #234](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/234) opened to evaluate migration; assigned to `i2 - S1 - Tool Landscape Evaluations`, P2 (strategic DX win, not blocking).

#### Issue #235 — Anthropic 2026 Agentic Coding Trends Report skim (XS, P1)

Anthropic published industry data on multi-agent coding workflows. Our PM/Sci/Dev split, file-based memory, markdown-standup pattern, and scope-discipline rules were designed iteratively without surveying industry patterns. Low-cost opportunity to sanity-check before they ossify. [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235) — 30-min skim with 5 PM-flavoured cross-checks (role decomposition, coordination protocol, memory architecture, scope-discipline failure modes, handoff mechanics).

**Why P1 (vs P2):** This is rare for a research/eval issue. Justified because the cross-check informs every future PM decision (and indirectly Sci/Dev workflows) at low cost. Catching divergences early is much cheaper than after the patterns are entrenched. Not blocking active work, but high information value per minute spent.

#### Cerebrum vs project scope distinction

User flagged that I'd labelled #235 "Cerebrum cross-checks" when it should have been "PM practices cross-checks". **Cerebrum** = the meta multi-agent framework above all projects/roles; **this project** = one instance using Cerebrum. Conflating them inflates scope (project-local issues become "Cerebrum architecture" reviews). Saved as shared memory `feedback_cerebrum_vs_project.md` so all roles get the rule. Reframed #235's title and body accordingly.

#### News_log format extension — time + editor sub-heading

Originally the news_log had one date section with bullets. With PM joining Sci/Dev as a logger, editor attribution was previously implicit and going to break down. Mirrored the lab-notebook format: each session adds a `### HH:MM UTC — Editor: <Role>` sub-heading under the date. Historical entries (pre-2026-05-02) kept as-is. Documented in `shared/reference_news_log.md` and `pm/feedback_morning_warmup.md`. Shipping convention: `docs/<role>/news-log-YYYY-MM-DD-HHMM` time-suffixed branch (no issue link), mirroring the multi-session lab-notebook pattern.

#### Standup — Developer Pending re-raise cleared

[Re-raised the closing-run issues #193/#194/#195 message](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues) after 3 days of no response. Developer acknowledged at 08:56 UTC — same response as Scientist (will populate when scoping each iteration's work). Standup is now fully clear.

---

## 2026-05-01

### 18:30 UTC — Editor: Developer

#### CLAUDE.md — pin GitHub project board ID

Quick board scan surfaced that the project board ID isn't documented in-repo, and `gh api graphql` against `organization(login: "Jin-HoMLee")` silently returns null because it's a *user* project, not an org project. Added a one-liner to CLAUDE.md pinning project #9 + the correct GraphQL shape ([PR #231](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/231)).

### 17:52 UTC — Editor: Developer

#### Issue #90 / PR #179 — closed without merging; revival deferred

[PR #179](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/179) had Approve-with-nit on April 28; the nit was fixed in `51f6c43` this morning, but the PR then sat ~3 days during which [PR #210](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/210) and [PR #227](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/227) reshaped the report-rule outputs. By the time we returned, PR #179's consolidated rule was stale (only knew about the original single `report_html` output) and an attempted main-merge produced an incomplete `report.smk` (missing `report_3d_structure_tsv`) plus a real conflict in `structure.smk`. Closed without merging; left full revival context on [Issue #90](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/90#issuecomment-4360756527).

The Issue's premise is still valid: when TCRdock is enabled, both `generate_report` and `generate_report_with_structure` exist and emit overlapping outputs (`report_html`, `report_tsv`, `report_top_candidates_tsv` at identical paths) with no `ruleorder` to disambiguate. PR #210 added `report_3d_structure_tsv` only to the structure rule; that's the one output that distinguishes them. Three viable revival paths captured on the Issue: (1) consolidate with a sentinel/empty 3D-structure TSV when TCRdock is disabled, (2) function-driven `output:` unpack, (3) just add a `ruleorder` directive — option 3 is the smallest fix that addresses the actual ambiguity, options 1–2 are the cleaner consolidation.

The "log path nit" the reviewer flagged on PR #179 (`_LOGS/{patient_id}/analysis/report.log` → `report/report.log`) was a regression introduced by PR #179's own rename — main still uses `analysis.smk` so the path is internally consistent there. No follow-up fix needed on main.

### 14:51 UTC — Editor: Developer

#### Issue #221 — PR #210 review-fix follow-ups

Bundled the four lower-priority items deferred from PR #210 review into one PR. Reviewer flagged them as Low/Nit so the diff stays focused on the original refactor; this picks them up while context is fresh.

**Items addressed:**

1. **Test now actually verifies the artefact-driven path** (was Low-Medium). The original test claimed corrupted raw inputs shouldn't break HTML rendering but never corrupted anything — it ran `generate_report` once and asserted `<html` was in the output, a tautology that would pass whether the inversion was real or not. The reviewer's literal suggestion ("overwrite raw inputs with garbage, call generate_report with artefact paths") doesn't work because `generate_report` always reads raw inputs at the top of every call before writing artefacts — corrupting them just crashes the second call. Equivalent-strength approach: monkey-patch `_build_report_top_candidates_tsv` so that AFTER it writes the artefact we inject a sentinel peptide name ("SENTINELPEPTIDE") into the file. The HTML rendering happens after the writer and reads from the modified file. The sentinel would never appear in the HTML if the renderer fell back to raw `pred_df` — its presence proves the artefact path is genuinely active.

2. **`effective_patient_id` path-derived fallback now warns.** When `patient_id` isn't passed, the function falls back to `output_html.parts[-3]` assuming the layout `…/{patient_id}/reports/report.html`. If that layout ever changes, every artefact row records the wrong patient ID silently. Added a `log.warning` when the fallback fires — the Snakemake path always passes `wildcards.patient_id`, so the warning only surfaces on ad-hoc CLI runs (which is when it's most useful).

3. **`_build_report_3d_structure_tsv` gracefully handles column rename.** Was hardcoding `top.get("mhc", "")` for the allele. If TCRdock ever renames the column (or some local TCRdock version uses a different convention), the manifest silently records empty allele → 3D viewer shows "NA". One-liner change to `top.get("mhc") or top.get("allele", "")`.

4. **Truncation-notice nit kicked to a proper Issue.** `_build_strong_table_html_from_top_candidates` (artefact path) can't say "Showing N of M rows" the way the raw `_build_strong_table_html` can — when the artefact is written, the writer caps to `TOP_CANDIDATES_LIMIT=10` and the original count is discarded. A literal parity fix needs a schema change in `report.tsv` to surface the pre-cap total. That's not trivial and not appropriate for a "nit" PR — opened [Issue #226](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/226) for the proper fix and added a docstring note linking to it.

**Verified:** `pytest workflow/tests/test_generate_report.py` — 57/57 passing. The new sentinel test would fail if anyone broke the artefact-driven path, so it's a real regression guard now.

[Issue #221](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/221) closed by the PR carrying this entry.

### 14:00 UTC — Editor: Scientist

#### Issue #203 — rescope to AlphaGenome-only + validation strategy designed

Picked up [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (Backlog, P1) per "what's next" recommendation — only pure-Scientist research issue in the queue. Set Status to In progress and refreshed the issue body in two steps over the session.

**Rescope: population-panel axis closed.** #203 originally proposed two complementary axes for rethinking normal filtering — (1) GTEx population panel as always-on filter, (2) AlphaGenome predicted-normal fallback for unmatched-normal patients. Yesterday's [#191](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/191) decomposition session resolved (1) cleanly via [#126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126) → [#211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) + [#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212), with the key design choice of *pan-tissue* (vaccine systemic-CTL safety dominates) over *tissue-matched*. So #203 is now narrowed strictly to axis (2) — AlphaGenome predicted-normal fallback. Updated body explicitly strikes through axis (1) with pointers to where it lives, so any reader of #203 understands the boundary.

**Literature review (2024–2026) on unmatched-normal handling.** Searched three angles: clinical pipelines that handle unmatched-normal cases, AlphaGenome's published splice-junction validation, and prior art on predicted-normal approaches. Five papers added to the Zotero collection with `manuscript-INTRODUCTION` / `-METHODS` / `-DISCUSSION` tags so the manuscript-section relevance is captured up-front:

| Paper | Year | Key role for us |
|---|---|---|
| [splice2neo](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae080/7684965) | 2024 | Closest published methodological analog (sequence-based prediction → validate against healthy tissue) |
| [AlphaGenome](https://www.nature.com/articles/s41586-025-10014-0) | 2026 | Core technology candidate; predicts splice site presence + usage + junction connectivity from sequence |
| [ENEO](https://academic.oup.com/nargab/article/7/3/lqaf026/8196479) | 2025 | Closest unmatched-normal pipeline (Bayesian tumor-only); inspires probabilistic combination framing |
| [SpliceMutr](https://aacrjournals.org/cancerrescommun/article/4/12/3137/750561) | 2024 | Pan-cancer splice neoantigen burden reference |
| [CNNeoPP](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1722117/full) | 2026 | Newest end-to-end personalized pipeline (LLM-enhanced) |

**Convergent finding from lit review:** the published default for cancer-specific splicing is `positive sample rate < 1%` in GTEx (per splice2neo and the [Tumour-wide RNA splicing aberrations](https://www.nature.com/articles/s41586-024-08552-0) paper). Our [#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) config uses `min_read_count: 1` — strictly more aggressive than the 1%-positive-rate convention. Justified by the vaccine-safety framing (precision ≫ recall for 10–20 vaccine slots), but worth flagging so the deviation is read as deliberate, not arbitrary. Posted as a non-scope-changing comment on [#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212#issuecomment-4359624169) — natural relaxation point if future analysis shows we lose too many true tumor-exclusive junctions.

**Sharpening finding on AlphaGenome.** The lit review surfaced a meaningful gap in AlphaGenome's published benchmarking. The Nature 2026 paper covers (a) tissue-aggregate junction prediction across 7 GTEx tissues, (b) variant-effect on junction counts, (c) held-out genomic interval performance. But **per-individual normal junction prediction from WGS — our actual use case — is not directly benchmarked**. AlphaGenome is trained on healthy reference data; feeding it a patient's WGS may produce something closer to a "tissue prior" than a "patient-specific normal predictor". This isn't a deal-breaker, but it sharpens the validation question: we have to test whether patient-specificity is actually present, not assume it.

**Validation strategy (three experiments on patient_001).** Patient_001 has matched-normal RNA-seq, so it functions as the gold standard:

- **Experiment 1 — Predictive validity.** AlphaGenome on patient_001 WGS → precision/recall/F1 vs. observed matched-normal BED, stratified by tissue track.
- **Experiment 2 — Patient-specificity test (the discriminating experiment).** AlphaGenome on GRCh38 reference for the same regions → quantify how many predictions change when patient germline variants are included. **If delta is small, AlphaGenome is a tissue prior (redundant with #211/#212); if meaningful, it captures germline-driven splicing variation (genuine 3rd source).**
- **Experiment 3 — Comparative filter strength.** Apply matched-normal / GTEx pan-tissue / AlphaGenome separately on patient_001's tumor junction set → quantify overlap, unique contributions, total filtered.

**Decision rules locked into [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body:** F1 ≥ 0.8 AND patient-specific delta > 10% → adopt as 3rd always-on source; F1 ≥ 0.7 AND ≥ 5% unique vs GTEx → adopt as fallback only; F1 < 0.5 OR delta < 1% → no-go.

**Why this is publishable either way.** Both outcomes yield novel insights: a *positive* result is the first demonstration of foundation-model-driven per-individual normal filtering for splice neoantigens (methodological contribution); a *negative* result is a cautionary insight that genome foundation models operate as tissue priors not individual predictors on this task (informs the field's expectations + framing). The discriminating measurement is Exp 2's patient-specificity delta — that's the scientific novelty regardless of direction. Manuscript follow-up tracked under PM's standalone S7 [milestone #16](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/16).

**Sub-issues created and linked under [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203):**

- [Sub-Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223) — `feat(api): confirm AlphaGenome API access + cost feasibility` (Dev, XS) — gates B+C. Cohort upscale signal captured: if costs are low enough, expand beyond patient_001 for statistical power.
- [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) — `research(filter): patient_001 notebook for Experiments 1+2` (Sci, S, blocked by #223)
- [Sub-Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) — `research(filter): Experiment 3 — comparative filter strength` (Sci, XS, blocked by #224)

**Side products of the lit review session:**
- [Issue #222](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/222) — `research(splice): evaluate splice2neo as alternative or component` opened. splice2neo is the closest published methodological analog to our pipeline; the question is whether it's a better-engineered version of what we're building or whether our differentiation (GPS scoring, TCRdock structural validation, vaccine-safety pan-tissue filter) keeps it complementary. Same eval pattern as [#188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) (Boltz-2), [#218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) (HERMES), [#201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) (ImmSET) — but on the upstream junction-calling side rather than the structure side.
- [Threshold-stringency comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212#issuecomment-4359624169) on [#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) (see above).

[Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) stays In progress while sub-issues run; Sci closure happens after Exp 3 + decision write-up + S7 manuscript follow-up issue scoping.

### 11:52 UTC — Editor: Developer

#### Issue #219 — `research/news_log.md` shared news-deduplication file

Established a shared news log convention to prevent re-surfacing the same items across morning briefings (Developer + Scientist). Without it, repeated WebSearches over a single sprint hit the same tool releases / papers / deprecations and waste briefing time. Two cerebrum changes already landed this morning to wire the convention into the morning routine — this PR commits the actual file to the repo so both roles can read + append.

**Format (minimal — won't grow unwieldy):**

```
- **Item** (date) — keywords. → action. *Role. Signal.*
```

- *Role* = Developer / Scientist / Both (whoever surfaced it)
- *Signal* = pipeline-relevant / industry-standard / rising-trend / portfolio-differentiator (per `feedback_portfolio_lens.md`)
- *Action* = Issue number opened, or "no action" / "no action; revisit if X"

Initial entries cover items from this week's briefings: `uv 0.11.8 / ruff 0.15.7` (industry-standard, no action), HERMES (portfolio-differentiator, → [#218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218)), MHCflurry 2.2.0 (pipeline-relevant, no action — torch already compat), Snakemake 9.19.0 (pipeline-relevant, → [#200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200)), AlphaGenome (portfolio-differentiator, → [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)), ImmSET (portfolio-differentiator, → [#201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201)).

**Cerebrum side (already done, not in this PR):** `feedback_morning_routine.md` Phase 1 now starts with "read `research/news_log.md` before WebSearch"; new shared memory `reference_news_log.md` documents the convention; `developer/shared/MEMORY.md` lists the news log under a new "Always read before morning routine (Phase 1 — news)" section.

[Issue #219](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/219) closed by the PR carrying this entry.

### 11:37 UTC — Editor: Developer

#### Issue #79 / PR #210 — review-fix follow-up

[PR #210](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/210) reviewer (Claude code-review action) flagged 2 **Medium** items, 3 **Low** items, and 1 **Nit** on yesterday's Phase 2 refactor. Addressed both mediums in `98ed6fd` and deferred the rest to a single follow-up Issue to keep this PR's diff focused on the original refactor scope.

**Mediums fixed:**

- **Filter alignment.** `_build_report_top_candidates_tsv` was filtering `presentation_class == "strong"` while `_build_report_tsv` and `_resolve_top_candidate_for_structure` both used `isin(["strong", "weak"])`. The divergence meant a patient whose only top presenter was "weak" would have a populated `top_candidate` row in `report.tsv` but an empty `report_top_candidates.tsv` and a `("NA", "NA")` 3D viewer annotation. Aligned all three on `isin(["strong", "weak"])` — this matches the docstring's "all three surfaces stay aligned" claim and is consistent with the percentile gate already referencing the *weak* threshold. Updated the corresponding empty-output log message ("no strong presenters" → "no qualifying presenters") so it doesn't lie about the broader filter.
- **`_rank_presenters()` helper extracted.** The three-key sort (`genotype_presentation_score` desc → `n_strong_alleles` desc → `best_presentation_percentile` asc, with `presentation_percentile` fallback) was copy-pasted across three functions. Now defined once at the top of the report-generation section; all three call sites do `df = _rank_presenters(df)`.

**Verified:** `pytest workflow/tests/test_generate_report.py` — 57/57 passing.

**Deferred to follow-up Issue (PR #210 polish):**

- Low-Medium: `test_artefacts_drive_html_when_provided` doesn't actually corrupt raw inputs after artefact write — needs to overwrite raw inputs with garbage and assert correct patient ID is still in HTML to genuinely prove the artefact path is active.
- Low: silent `effective_patient_id` path-derived fallback — needs `log.warning` when the fallback fires.
- Low: `_build_report_3d_structure_tsv` hardcodes `"mhc"` column for allele — needs `top.get("mhc") or top.get("allele", "")` fallback.
- Nit: `_build_strong_table_html_from_top_candidates` no truncation notice — parity gap vs raw-path table; low value while `TOP_CANDIDATES_LIMIT=10` but worth a comment.

Reviewer LGTM at 10:56 UTC; deferral rationale accepted ("scope call is reasonable"). [PR #210](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/210) merging next.

### 11:23 UTC — Editor: Scientist

#### Issue #103 — folded into Issue #215 (filtering audit trail)

While running the "what's next" recommendation algorithm, [#103](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/103) (track per-file mean threshold + read-count distribution) surfaced as a top-3 candidate. On closer look it overlapped substantially with [#215](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/215) (full filtering audit trail, sub-issue of [#104](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104)).

**Column-by-column overlap:** 3/7 of #103's TSV columns (`n_junctions_raw`, `n_after_mean_filter`, `pct_retained`) duplicate #215's `junction-filter | junctions_raw` and `after_mean_filter` rows. The remaining 4/7 (`min_reads`, `mean_reads`, `median_reads`, `max_reads`) are uniquely additive — they answer the **threshold-tuning question** ("is the per-file mean appropriate for each sample?") rather than the **funnel question** ("where are junctions lost?") that #215 already covers. The mean filter is currently silent (the threshold value isn't recorded anywhere); on a deeply-sequenced sample it can sit at hundreds of reads and silently filter out rare-but-real junctions.

**Format compatibility ruling:** #215 is long-format (`step | category | count`). The 4 distribution stats are different *measurements* of the same `junction-filter` step on the same sample, so they fit naturally as 4 additional rows under `junction-filter`, not as a separate file or hybrid schema. No restructuring needed.

**Decision: fold #103 into #215.** Updated #215's body to add the 4 distribution categories under `junction-filter` and a new "Why include the 4 distribution stats" section that captures the threshold-tuning rationale (so the motivation isn't lost when #103 closes). Single coherent instrumentation pass for Developer, no sequencing dance, marginal scope add (~XS-worth on top of an S issue, well within the band). #103 closed with a comment explaining the overlap analysis and pointing to #215. FYI heads-up posted to Developer in standup — net effect is 4 extra stat rows per sample, no sequencing change.

### 09:30 UTC — Editor: Scientist

#### Issue #218 — HERMES (Visani et al., PNAS 2025) eval issue created

Developer flagged [HERMES (PNAS, Oct 2025)](https://doi.org/10.1073/pnas.2504783122) in standup as a structure-based ML scoring candidate for TCR-pMHC. Despite zero TCR-pMHC training, HERMES achieves 0.72 correlation with experimental binding affinities — it's a physics-guided 3D-equivariant model trained on the protein universe, predicting amino acid preferences from local structural environments. The [Feb 2026 bioRxiv follow-up](https://www.biorxiv.org/content/10.64898/2026.02.24.707744v1.full) applies AF2-predicted structures with the same scoring logic.

Framed for our pipeline as a **post-TCRdock confidence cross-check** rather than a replacement: HERMES could score the structural quality / interaction strength of the TCRdock-predicted complex, particularly appealing for neoantigens where domain-specific TCR-pMHC training data is scarce. A second potential use is **VDJdb panel pre-screening** before expensive co-folding ([#205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205)). Same eval pattern as [#188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) (Boltz-2) and [#201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) (ImmSET) — assess methodology, decide whether to integrate.

[Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) opened. Zotero key `MWZFINV6`. Follow-up reply to Developer in standup confirming.

---

## 2026-04-30

### 15:43 UTC — Editor: Scientist

#### Issue #191 — decomposition session for Issue #126 (GTEx pan-tissue junction filter)

Decomposed [parent Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126) into two implementation sub-issues. Most of the parent's scientific design was already locked (pan-tissue scope, GTEx V10 source, `filter_junctions.py` integration point), so the session was mostly mechanical — except for one important design pivot.

**Design pivot — adopted *always-on stacking*, rejected *tissue-matched scope*.**

The morning's Sci+Dev discussion (captured in [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203), Developer-created) had floated changing #126's design in two ways at once: (a) always-on filter that stacks with matched-normal (vs. opt-in fallback), and (b) tissue-matched scope (vs. pan-tissue). The morning framing motivated both via the GATK Panel-of-Normals analogue from somatic variant calling.

On reflection in this Sci session, **(a) was kept and (b) was rejected**:

- **Always-on stacking ✅ adopted.** Matched normal captures one observation of one individual in one tissue; GTEx pan-tissue captures population-level expression across ~54 tissues. They are complementary, not redundant. Stacking is the right default.
- **Tissue-matched scope ❌ rejected.** A vaccine-induced CTL response is *systemic* — patrols all tissues, not just the tumour's lineage. Restricting GTEx to e.g. bone tissue for an osteosarcoma patient leaves heart, liver, kidney, etc. unprotected against autoimmune off-tumour toxicity from cross-tissue junction expression. The PoN analogue from somatic variant calling is useful framing for *always-on stacking*, but the *tissue scope* choice is dominated by the vaccination-specific safety argument — same reasoning as the original [#126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126) spec, which I now re-confirm.

This is why the Scientist role exists separately from the Developer role: the GATK PoN intuition is sound for somatic SNV/indel calling, but breaks at the vaccination/autoimmunity level because the deliverable is no longer a list of mutations but a population of CTLs released into a patient's circulation. Captured this distinction explicitly in the resolution comment on [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203#issuecomment-4353911289).

**Locked-in defaults (`gtex_filter` config):**
- `enabled: true` (always-on; was `false` in #126)
- `min_read_count: 1` (most aggressive — precision over recall, vaccine 10–20 slot framing)
- `reference_bed: gs://splice-neoepitope-project/resources/gtex/v10/gtex_v10_pan_tissue_junctions.bed` (GCS-staged, consistent with rest of pipeline data flow)

**Output:**
- [Sub-Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) — `feat(filter): GTEx V10 pan-tissue junction reference set construction` (M) — pin GTEx V10 release, build pan-tissue union BED, stage to GCS
- [Sub-Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) — `feat(filter): integrate GTEx pan-tissue filter into filter_junctions.py — always-on, stacks with matched-normal` (M, blocked by #211) — wire into pipeline, validate on patient_001 (stacking) + patient_002 (sole filter beyond annotated)

Both linked under [parent #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126), milestone i2-S3. Cross-reference comment posted on [#126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126#issuecomment-4353911112) recording the default-flag change so the parent body's `enabled: false` snippet doesn't mislead future readers.

**Manuscript follow-up (not in this PR):** the DISCUSSIONS.md "Normal sample filtering → GTEx pan-tissue filter" section's *scientific reasoning* is unchanged, but if it currently describes the filter as "opt-in" or "fallback when no matched normal", that one-liner needs updating once [#212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) lands and we have validation numbers to cite.

[Issue #191](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/191) closed by the PR carrying this notebook entry.

### 14:15 UTC — Editor: Developer

#### Issue #79 Phase 2 — HTML now driven by report.tsv + two new artefacts

Completed the data/presentation decoupling for [Issue #79 (regen report.html from report.tsv)](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/79). Phase 1 yesterday wrote three machine-readable artefacts (`report.tsv`, `report_top_candidates.tsv`, `report_3d_structure.tsv`) but HTML still read the raw pipeline files. Phase 2 inverts that — `generate_report()` now writes the artefacts first, then reloads them and renders HTML from those projections.

**New helpers in `generate_report.py`:**

- `_load_report_tsv(path)` — long→wide pivot returning per-stage projections (junction_filtering DataFrame, mhc_prediction count dict, top_candidate dict, hla_typing dict, tcrdock dict). Non-lossless by design — exposes only what the HTML render needs.
- `_render_contig_peek(peek)` — parses the bracketed plain-text peek from `report_top_candidates.tsv` (e.g. `AAA[CC|GG]TTT`) into the same span-classed HTML `_render_contig` produces from raw sequence. Round-trip locked in by test.
- `_build_strong_table_html_from_top_candidates(df)` — consumes the wide TSV directly. No raw `pred_df` or contigs FASTA needed; the writer's quality gate, sort order, and `TOP_CANDIDATES_LIMIT=10` cap all flow through automatically.
- `_presenter_counts_html(mhc_prediction_dict)` — renders presentation-class counts from the loaded `report.tsv` projection.
- `_resolve_top_candidate_from_manifest()` / `_resolve_top_candidate_for_structure()` — split the two sourcing paths for the 3D viewer's annotation; `_build_structure_section()` signature simplified to `(pdb_path, peptide, allele)`.

**Vocabulary sweep** per `developer/shared/feedback_presenters_terminology.md`: "binder" → "presenter" / "presentation" everywhere user-facing and in internal naming. The IC50 column tooltip "Binding affinity in nM" stays — IC50 IS binding affinity per the rule's stated exception.

**Backward compat:** when `output_tsv` / `output_top_candidates_tsv` / `output_3d_structure_tsv` aren't passed (CLI-only path), the original raw-input renderers are used as fallbacks. No regression for ad-hoc `python generate_report.py` runs.

**Verified:**

- `pytest workflow/tests/`: **226/226 passing** (16 new tests added across `TestRenderContigPeek`, `TestBuildStrongTableHtmlFromTopCandidates`, `TestPresenterCountsHtml`, `TestGenerateReportEndToEnd`)
- `snakemake -n --configfile config/test_config.yaml`: 6 jobs, `generate_report` resolves
- `snakemake -n` with GPU overlay: 7 jobs, `generate_report_with_structure` resolves with the 3D manifest output

**Two commits on top of yesterday's three:**

- `3f55309` — feat(report): add `_load_report_tsv()` loader for HTML decoupling
- `0e1610f` — refactor(report): drive HTML from report.tsv + top-candidates + 3D manifest

**Follow-up still scoped under [Issue #198](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/198):** distributing artefact ownership across the producing scripts (move `report_top_candidates.tsv` and `report_3d_structure.tsv` writes to `run_mhcflurry.py` and `run_tcrdock.py` respectively). Today's PR keeps writers in `generate_report.py` for atomicity — `#198` will redistribute.

---

### 10:54 UTC — Editor: Scientist

#### Issue #190 — decomposition session for Issue #86 (patient-HLA-matched TCR panel)

Decomposed [parent Issue #86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86) into three implementation sub-issues. The big scope decision in this session was **collapsing the structural fan-out**: #86's original spec (top 10 TCRs/allele × top neoepitope) implied ~60 TCRdock jobs/patient (~30 h on P100 — Scientist-side estimate, not measured). After discussion we landed on a much narrower and more honest scope:

- **Panel itself stays at top 10/allele** as a *reference set* surfaced in the report (transparency about which TCRs were considered)
- **TCRdock structural prediction stays at 1 prediction/patient** (today's baseline) — single TCR × top neoepitope × top neoepitope's `best_allele`, but now HLA-matched instead of the hardcoded DMF5
- **Selection rule:** highest VDJdb confidence in the panel for the top neoepitope's `best_allele`, with VDJdb donor ID as deterministic tiebreaker

Why this works: today's pipeline runs DMF5 (HLA-A\*02:01-restricted) regardless of the patient's actual `best_allele`, which is structurally meaningless when the top neoepitope binds e.g. HLA-C\*07:01. Just swapping in an HLA-matched TCR — even a single one — is the core scientific gain of #86. Multi-TCR variance measurement is interesting but a separate experiment, not a production blocker.

**Supertype / 2-digit fallback decision.** The narrowed scope (1 TCR instead of 10) makes the "0 exact matches" case bite harder — we drop straight to DMF5. Considered moving supertype fallback in-scope here. Decided to defer (option a) but **instrument now**: [sub-issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) logs per-allele exact-match counts in a `vdjdb_panel_qc.tsv` sidecar. If patient_001/002's alleles show poor coverage in practice, we file a focused follow-up with the actual data. Avoids speculative scope creep + Sidney 2008 vs Lund 2004 supertype-table debates that don't yet have empirical motivation.

**Output:**
- [Sub-Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) — `feat(tcrdock): fetch_vdjdb_panel Snakemake rule + stitchr/VDJdb setup` (M)
- [Sub-Issue #205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205) — `feat(tcrdock): select HLA-matched top TCR from panel for structural prediction` (S, blocked by #204)
- [Sub-Issue #206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206) — `feat(report): surface VDJdb panel + matched TCR in report.html` (S, blocked by #205)

All three linked under [parent #86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86) and on the M5 - Modeling i2 milestone. [Issue #190](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/190) closed by the PR carrying this notebook entry.

---

### 08:35 UTC — Editor: Scientist

#### PR #199 — review nit: weak-presenter threshold definition tightened

Claude reviewer flagged that `(≤ 2.0%)` for weak presenters is ambiguous — it technically includes the strong-presenter pool (≤ 0.5%). Fixed both patient_001 KF#3 and patient_002 KF#7 in `CONCLUSIONS.md` to use the explicit range `(0.5% < percentile ≤ 2.0%)`, matching the definition already used in `RESULTS.md`. Counts unchanged; wording only.

---

## 2026-04-29

### 22:04 UTC — Editor: Scientist

#### Issue #197 — manuscript cross-doc cleanup (CONCLUSIONS + METHODS + DISCUSSIONS)

Bundled cross-doc cleanup PR addressing three threads surfaced during [PR #196](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/196) review:

**A. CONCLUSIONS.md patient_002 KF#5–7 refresh.** Numbers verified at source against `gs://splice-neoepitope-project/results/patient_002/reports/report.tsv`:

| Finding | Was (pre-#148, IC50) | Now (post-#148, GPS) |
|---|---|---|
| 5: tumor-exclusive | 55,912 from 347,046 (16.1%); 3 WES junctions removed | 58,914 from 364,168 (16.2%); 0 WES junctions removed (WES has no RNA splice reads) |
| 6: peptides | 781,424 9-mers (775,440 unique) | 2,330,687 across 8/9/10-mers (2,313,700 unique, 99.3%) |
| 7: predictions | 12,430 IC50-strong (0.32% of 3,907,120) | 67,935 percentile-strong (2.9% of 2,330,687) + top FADLRPLLL / HLA-C\*01:02 (GPS 0.9999, 4/5 alleles) |
| 8: HLA validation | OK (no change) | OK (no change) |

**B. METHODS.md §4 reconciliation.** Investigation in `workflow/scripts/aggregate_hla_alleles.py:5-22` showed the actual code cascade is **tumor → serology → normal → fallback**, with the docstring rationale: *"OptiType from Primary Tumor — reflects what the tumor cell actually presents, including any HLA loss-of-heterozygosity (LOH) events."* The "normal-first policy" wording was therefore fully out of date — pure doc-fix, not a config bug. METHODS.md §4 rewritten to describe the cascade with the LOH rationale; null-allele exclusion preserved. Output table caption `(normal-first)` → `(tumor-first cascade)`. Output table peptide row updated `9-mers` → `8/9/10-mers` (PR #99 alignment).

**C. CONCLUSIONS.md §Significance ProcessPoolExecutor stale claim.** Per CLAUDE.md, `ProcessPoolExecutor` was removed from `run_mhcflurry.py` (CUDA crash on GPU / 48 GB OOM on CPU). Replaced with the actual implementation: single genotype-level `Class1PresentationPredictor.predict()` call with peptides batched in-call. Also flipped "binding predictions" → "presentation predictions" per the [presenters terminology](feedback_presenters_terminology.md) rule.

**Cross-doc cascade fixes:**

- DISCUSSIONS.md "Impact of missing matched normal" paragraph rewritten — the pre-#148 narrative ("106,474 WES junctions, 3 overlapped") was an artifact of the chr-naming bug; post-#148 truth is simpler (zero RNA junctions extracted from WES alignment by design). Numbers refreshed.
- CONCLUSIONS.md Future Directions patient_002 longitudinal entry refreshed.
- patient_001_results.ipynb §7 item #4 marked fully resolved (was already partially resolved by Issue #178; METHODS reconciliation closed the loop).

Source-of-truth discipline reaffirmed throughout: `report.tsv` first, notebook post-calc only for what `report.tsv` doesn't expose, alignment intermediate (`alignment/{sample}/junctions.tsv` line count) only for the funnel-top rows pending [Issue #104](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104).

---

### 19:51 UTC — Editor: Scientist

#### PR #196 — review fix: `(one prediction per unique peptide)` parenthetical incorrect

Claude review on [PR #196](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/196) caught one substantive inconsistency: `RESULTS.md` Peptide Translation section said `(one prediction per unique peptide)` but `total_predictions = 1,286,492` matches the total peptide count, not the unique count (1,260,074). `Class1PresentationPredictor` predicts on all input peptides; duplicate sequences from different junctions or reading frames are each predicted separately.

Fixed parenthetical to `(one prediction per input peptide; duplicate sequences from different junctions or reading frames are each predicted separately)`.

Review also noted a pre-existing stale `ProcessPoolExecutor` claim in CONCLUSIONS.md §Significance — added to [Issue #197](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/197) scope alongside the patient_002 KF#5–7 staleness and METHODS.md reconciliation.

---

### 15:21 UTC — Editor: Scientist

#### Issue #178 — patient_001 results into RESULTS.md (consolidation of #171)

Mirroring the patient_002 RESULTS.md layout, replaced the legacy patient_001 section
(pre-#148, IC50-based Run 1 / Run 2 comparison) with the post-#148 valid-run numbers in
the 8-section structure: Dataset / HLA Typing / Junction Filtering / Peptide Translation
/ MHC Presentation / GPS / Top Candidate / TCRdock.

Source-of-truth hierarchy applied per [user clarification](
https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/178):

- **`report.tsv` first** — junction filtering (unannotated/normal_shared/tumor_exclusive),
  MHC class breakdown, top-candidate fields, HLA typing rows, TCRdock pdb flag.
- **Notebook post-calculations** only for what `report.tsv` doesn't expose — peptide
  length distribution, unique-peptide count, strong-presenters-per-allele, GPS
  distribution stats, per-allele dominance summary, read-support stats.
- **Alignment intermediate** for the two top-of-funnel rows (`junctions_extracted_total`
  and `junctions_annotated_discarded`) until [Issue #104](
  https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104) adds them to
  `report.tsv` directly. Computed from `gs://splice-neoepitope-project/results/patient_001/alignment/SRR9143066/junctions.tsv`
  line count = 146,647; annotated discarded = 146,647 − 30,029 = 116,618. Both rows
  are explicitly tagged with their derivation in the section.

Stale-number diff caught at source: current RESULTS.md said `146,644` total / `116,615`
annotated — both off by 3 from the post-#148 line count. Likely carried over from the
pre-#148 run. Patient_002's parallel rows verified separately (`alignment/BG003082_T0/junctions.tsv`
line count = 364,168, matches RESULTS.md exactly) — no fix needed there.

HLA typing section rewritten: the legacy "normal-first policy" wording was incorrect for
this run — `report.tsv` `hla_typing` rows show `source: tumor` for all three loci. The
MHCflurry predictions consumed the tumor calls (HLA-B = B\*18:01 / B\*15:63), not the
normal calls (B\*15:01 / B\*18:02) that the old RESULTS.md table listed. Documented the
B-locus tumor-vs-normal noise pattern factually without reasserting a policy.

Cross-patient framing kept tight: HLA-C dominance partial (57.2% vs ~69% in patient_002),
GPS top candidates landing at GPS ≈ 0.9999 / IC50 ≈ 33–34 nM in both patients despite
different best-allele HLA-C — suggests structural ceiling at the GPS-best slot.

Open follow-up: [Issue #104](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/104)
will collapse the alignment-intermediate rows into `report.tsv`; once it lands the
section can be re-pointed entirely at `report.tsv` and the lab-notebook footnote
removed.

---

### 14:09 UTC — Editor: Scientist

#### PR #189 — review round 2: stale cross-reference in patient_002

Second review on [PR #189](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/189) confirmed all round-1 fixes verified, with one trivial leftover: `patient_002_results.ipynb` §7 finding 6 read `"58,914 vs 27,347"` — the patient_001 cross-reference was a pre-fix value. After round 1 corrected patient_001's tumor_exclusive count to 27,348, this cross-reference stayed stale.

Updated to `27,348`. Functionally harmless (the "~2×" ratio was unchanged) but internally inconsistent across the two notebooks.

Pattern lesson recorded: when a number is corrected in one Jupyter notebook, every cross-reference to it from other notebooks/docs needs to be searched and updated — copy-paste of values between notebooks is a recurring vector for stale drift.

---

### 13:49 UTC — Editor: Scientist

#### PR #189 — review fixes (junction-count mislabel + minors)

[Claude review on PR #189](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/189) flagged a medium-severity factual error: §2.1 cell printed `"Total tumor_exclusive junctions: 30,029"` but `report.tsv` (loaded by the same notebook in §1) shows unannotated=30,029 / tumor_exclusive=27,348 / normal_shared=2,681. Root cause: silent hallucination from copy-paste — in patient_002 with no matched normal, unannotated == tumor_exclusive (both 58,914), so the label was numerically correct there. For patient_001 the difference surfaced.

Decision: filter at load to tumor_exclusive (clinically correct fix), and show all three counts in the print so the matched-normal value-add is visible:

```
Unannotated junctions:         30,029
  ├─ normal_shared (removed):   2,681
  └─ tumor_exclusive:          27,348   ← carried forward to MHC prediction
```

Applied the same template to patient_002 §2.1 for consistency (prints `0 ← no matched normal: nothing removed` — honest about absent filtering).

Cascading fixes in `patient_001_results.ipynb`:
- §7 finding 7 — count corrected to 27,348; removed contradictory comparison to "27,347 from earlier runs"; reframed as "after matched-normal filtering removed 2,681 of 30,029 (8.9%)"
- §7 finding 8 — read-support stats now computed on tumor_exclusive (mean 49 vs old 56; min/median/max effectively unchanged at 16 / 26 / 82,098 since the dropped normal_shared junctions skewed slightly higher than the kept tumor_exclusive ones)

Also addressed the two minor review notes:
- §6 GPS formula — added caveat that locus weighting is approximate; MHCflurry applies its own locus-aware calibration internally during scoring
- §5.2 — added inline source for the "HLA-C\*07:01 shared with patient_002" cross-patient claim

Re-ran both notebooks via `jupyter nbconvert --execute --inplace` (cache warm — no re-download).

---

### 10:32 UTC — Editor: Scientist

#### Sub-Issue #177 — patient_001 results analysis notebook

Created `research/notebooks/patient_001_results.ipynb` by copying the patient_002 template and adapting:
- BUCKET / CACHE_DIR → patient_001 paths
- Title block: SRR9143066 (gastric cancer) + SRR9143065 matched normal (vs. patient_002's BG003082 osteosarcoma + WES-only normal)
- Junction-analysis intro: rephrased to reflect matched RNA-seq normal applied (vs. patient_002's no-matched-normal caveat)
- HLA-allele dominance intro: listed the actual tumor-first HLA calls used by MHCflurry (HLA-A\*31:01/A\*26:01, HLA-B\*15:63/B\*18:01, HLA-C\*07:01/C\*03:03); reframed cross-patient question — does HLA-C dominance persist with patient_001's different HLA-C alleles?
- Top-candidate cell (5.2): rewrite around SQIPRTHSY / HLA-C\*07:01 (from `report.tsv`); cross-patient note that both top candidates land at GPS ≈ 0.9999, IC50 ≈ 33–34 nM despite different splice contexts and different HLA-C alleles (C\*01:02 vs C\*07:01)
- Inflation-check comment (6.3): generalised, with explicit cross-reference to patient_002's 174 cases
- Section 7 summary: converted to a fill-in-after-run template scaffolded for cross-patient comparison

Ran the notebook end-to-end via `jupyter nbconvert --execute --inplace` (1h timeout). All cells executed cleanly; the `parallel_process_count=1` fix held for the gsutil cp downloads (~633 MiB total: 97 MiB peptides + 536 MiB mhc_presentation), and the atomic `.tmp + os.rename` pattern left no `.gstmp` leftovers.

**Key first-look findings (full interpretation pending):**
- GPS inflation cases: **25 candidates** with `GPS > 0.9` and `n_strong_alleles = 0` (vs. 174 in patient_002 — large drop)
- HLA-C dominance partially recapitulated: HLA-C\*03:03 (33.7% as best) + HLA-C\*07:01 (23.5%) = **~57% of strong presenters** (less extreme than patient_002's ~69%)
- HLA-A\*31:01 is surprisingly active here (28.2% as best) — different role from patient_002's nearly-silent HLA-A\*01:01
- HLA-B\*18:01 is the quietest in patient_001 (median percentile 10.7 — analogous to HLA-A\*01:01 in patient_002)

**Doc note:** RESULTS.md still describes HLA typing under the legacy "normal-first policy" wording; the pipeline now uses tumor-first. This is a stale-doc issue to clean up in a separate PR — not blocking this one.

---

### 09:05 UTC — Editor: Scientist

#### Issue #186 — Discussion section on immune-pathway gene neoepitopes and presentation paradox

Morning routine surfaced [Coelho et al., *Cancer Cell* 2023](https://doi.org/10.1016/j.ccell.2022.12.009) (base-editing screens mapping IFN-γ signalling missense vulnerability). The discussion that followed produced a clinical-translation angle worth committing to the manuscript: splice-derived neoepitopes from immune-pathway genes (JAK1/2, STAT1, B2M, NLRC5, TAP1/2) sit at the intersection of "immunogenic" and "selected for during evasion" — vaccine-prioritisable but with a presentation paradox.

Added a new subsection to `research/manuscript/DISCUSSIONS.md` between *Allele breadth and immunodominance* and *Structural validation*. Section covers:

1. **Preventive vaccination** angle (public neoantigens; Kwok et al. 2024; constitutive HLA-I via NLRC5 + NK missing-self maintain presentation despite partial JAK LOF).
2. **Combination therapy** for established tumours: rescue must land *downstream* of the broken kinase. Upstream interventions (intratumoural IFN-γ, IFN-γR delivery) are futile because the signal terminates at the dead JAK. Effective bypasses: TLR agonists (NF-κB/IRF3), STING agonists (type I IFN; impaired in JAK1-LOF since JAK1+TYK2 is shared), HDAC/DNMT inhibitors (epigenetic derepression — addresses the transcriptional bottleneck directly), and NK-cell engagers.
3. **Mermaid diagram** (Figure 1) visualising the pathway with broken vs. rescue vs. constitutive nodes colour-coded.

Also added Coelho et al. 2023 to the Zotero collection (key M2HCXESA) following the new three-section note format. Initially miscited the venue as Nature 2024 — corrected to *Cancer Cell* 2023 (DOI 10.1016/j.ccell.2022.12.009) before committing to the manuscript.

---

## 2026-04-28

### 16:36 UTC — Editor: Scientist

#### PR #185 — review fixes + conflict resolution

[PR #185](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/185) (closes [Sub-Issue #161](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/161), local GCS cache for research notebooks). Reviewer flagged one medium-severity reliability bug:

- **Partial-download cache poisoning** — if `gsutil cp` fails midway (network timeout, credential expiry, disk full), the partial file sits at `local_path` and the next kernel restart skips the download via `os.path.exists()`, then feeds pandas a truncated TSV. Fixed in `1fc4a2c` with the standard atomic-rename pattern: download to `local_path + ".tmp"`, then `os.rename` to `local_path` on success; on any exception, remove the tmp file and re-raise.

Also resolved a merge conflict against main from PR #184 (Developer's Issue #180 lab notebook entry) — both PRs prepended entries to the 2026-04-28 section. Merged main into the branch, kept all entries ordered newest-first by timestamp. Merge commit `edd268d`.

CI re-triggered after the merge resolution and is now queued; PR back to MERGEABLE.

---

### 15:16 UTC — Editor: Scientist

#### Issue #161 — additional fixes from end-to-end notebook run

Ran the patient_002 notebook end-to-end to populate outputs before committing. Two issues surfaced and were fixed:

1. **macOS gsutil multiprocessing hang** — `gsutil cp` (unlike `gsutil cat`) hangs on macOS due to the Python fork bug ([Python issue #33725](https://bugs.python.org/issue33725)). The `peptides_novel.tsv` download stalled indefinitely after ~14 min. Fixed by passing `-o "GSUtil:parallel_process_count=1"` to the gsutil call in `gcs_read_tsv_cached()` — disables multiprocessing but keeps multithreading.
2. **Cell 6.4 (per-allele contribution) `KeyError: 'genotype_presentation_percentile'`** — the allele filter `not c.startswith("best")` excluded `best_presentation_score` but failed to exclude `genotype_presentation_score`, so `"genotype"` ended up in the alleles list. Fixed by inverting to a positive filter `c.startswith("HLA-")`. Pre-existing bug; surfaced now because the prior run's outputs predated the GPS column.

---

### 14:51 UTC — Editor: Developer

#### Issue #180 — verify TCRdock has no dependency on AFDB legacy API

Audit confirmed no dependency. The AlphaFold Database (AFDB) legacy API retires June 2026, so any pipeline component fetching from `alphafold.ebi.ac.uk` would break — this issue checked whether we do.

Findings:
- Zero AFDB / `alphafold.ebi.ac.uk` references repo-wide (`grep -rni` across `*.py`, `*.sh`, `Dockerfile*`, `*.yaml`, `*.smk`).
- `docker/Dockerfile.pipeline` downloads AlphaFold v2 params from Dropbox (TCRdock's curated hosting) at image build time, not from AFDB.
- `workflow/scripts/run_tcrdock.py` imports no HTTP clients (`requests`, `httpx`, `urllib`, `http.client`) and shells out to no `wget`/`curl`. It invokes Docker and reads/writes local files.
- Inference uses `--data_dir /opt/TCRdock/alphafold_params` (local container path, baked at build time).

Conclusion: the June 2026 AFDB retirement does not affect us. Closed [Issue #180](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/180) with the audit details.

---

### 13:49 UTC — Editor: Scientist

#### Issue #161 — added local GCS cache to `patient_002_results.ipynb`

Replaced `gcs_read_tsv()` (streamed from GCS via `gsutil cat` on every kernel restart) with `gcs_read_tsv_cached()` that downloads each TSV once to `/tmp/splice-neoepitope-cache/patient_002/` and reads locally on subsequent calls. Cache is `/tmp`-based (cleared on reboot) — no stale data risk, no git noise.

Eliminates the ~3 min wait on `mhc_presentation.tsv` (2.3M rows) that previously blocked every analysis session. Pattern carries forward to patient_001 notebook (#177).

Updated five cells: setup (function definition + `CACHE_DIR` constant) plus four call sites (`report.tsv`, `novel_junctions.tsv`, `peptides_novel.tsv`, `mhc_presentation.tsv`).

---

### 13:30 UTC — Editor: Developer

#### Issue #181 — add research/glossary.md

Added `research/glossary.md` as a living dictionary of project-relevant abbreviations. Seeded with three entries that came up in conversation today (AFDB, GKE, SLURM). Format: alphabetical sections, `**ABBREV** — Expansion. One-line context. *Domain: bio | ml | cloud | pipeline | stats.*` Going forward, new abbreviations are added when they come up in conversation (rule saved in memory). Merged via PR #182.

---

## 2026-04-27

### 17:17 UTC — Editor: Developer

#### Issue #158 — rename config/gpu.yaml → config/gpu_config.yaml

Renamed `config/gpu.yaml` to `config/gpu_config.yaml` for naming consistency with other `*_config.yaml` files and to avoid confusion with conda env files in `workflow/envs/*.yaml`. Updated all live references across 9 files: `config/config.yaml`, `scripts/run_cloud_gpu.sh`, `scripts/setup_vm.sh`, `.github/workflows/tests.yml`, `README.md`, `CLAUDE.md`, `docs/configuration.md`, `docs/google_cloud_guide.md`. Historical entries in this notebook (lines 173, 572, 576, 708) intentionally left unchanged — they are point-in-time records.

Evaluated Snakemake profile migration (Issue #158 step 2) and deferred: `gpu_config.yaml` holds pipeline-level config values (TCRdock settings, fallback HLA/TCR alleles), not executor/resource parameters — the `--configfile` overlay is the correct abstraction.

Merged via PR #174 (replaced PR #170 which closed during a branch rename operation).

---

### 16:58 UTC — Editor: Scientist

#### PR #173 — patient_002 Junction Filtering correction merged

Reviewed by Claude (automated, LGTM, math verified). Both CI checks passed. Squash merged, branch deleted. Closes #172.

---

### 15:43 UTC — Editor: Scientist

#### Issue #172 — RESULTS.md patient_002 Junction Filtering corrected to valid post-#148 run

The Junction Filtering table previously reflected the invalid pre-#148 run (total 347,046, annotated 291,131, unannotated 55,915, normal-shared 3, tumor-exclusive 55,912). Corrected to the valid post-#148 run numbers derived from `alignment/BG003082_T0/junctions.tsv` and `reports/report.tsv` on GCS:

| Stage | Old (invalid) | New (valid) |
|-------|--------------|-------------|
| Total extracted | 347,046 | 364,168 |
| Annotated | 291,131 (83.9%) | 305,254 (83.8%) |
| Unannotated | 55,915 (16.1%) | 58,914 (16.2%) |
| Normal-shared | 3 (~0.0%) | 0 (0.0%) |
| Tumor-exclusive | 55,912 (16.1%) | 58,914 (16.2%) |

`normal_shared = 0` is correct and expected: BG003082_N0_WES was used as the normal input and yielded valid HLA typing results, but WES data contains no RNA splice junctions so junction-level normal subtraction was not effective. Also updated Dataset section to reflect this accurately.

Developer standup (2026-04-26 21:31 UTC, Q4) confirmed the valid run used no RNA-seq normal; the WES-yielding-0-junctions behaviour was separately clarified by the user.

---

## 2026-04-26

### 22:40 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: acceptance criteria verified (both patients)

Verified all GCS outputs against the Issue #123 acceptance criteria. All criteria passed.

| Criterion | patient_001 | patient_002 |
|---|---|---|
| Pipeline completes end-to-end without errors | ✅ | ✅ |
| `mhc_presentation.tsv` has breadth columns (`best_allele`, `genotype_presentation_score`, `n_strong_alleles`, `best_presentation_percentile`) | ✅ | ✅ |
| `genotype_presentation_score` in [0, 1] and scientifically sensible | ✅ [0.0147, 0.9999] — 1,286,492 rows | ✅ [0.0122, 0.9999] — 2,330,687 rows |
| `report.tsv` and `report.html` generated | ✅ | ✅ |
| Results uploaded to GCS | ✅ | ✅ |

Note: two column names in the acceptance criteria are stale — `strong_alleles` (actual: `n_strong_alleles`) and `presentation_percentile_strong` (actual: `best_presentation_percentile`). The equivalent columns are present; these are pre-rename leftovers and do not block closure.

Issue #123 ready to close.

---

### 22:18 UTC — Editor: Scientist

#### PR #168 — review fixes (docs/scientist/issue-165-patient-002-documentation)

Three items addressed before merge:

1. **Weak threshold label clarified** — `Weak (percentile ≤ 2.0%)` → `Weak (0.5% < percentile ≤ 2.0%)` to make the exclusive band unambiguous.
2. **GPS phrasing fixed** — `exceeded GPS > 0.9` → `had GPS > 0.9` (redundant double-negative removed).
3. **9-mer discrepancy** — added WES proxy normal usage as a second possible cause alongside the #148 fix.

---

### 21:50 UTC — Editor: Scientist

#### PR #167 — review fixes (feat/scientist/issue-164-patient-002-notebook)

Three issues addressed in code review before merge:

1. **GPS subsection headings moved to markdown cells** — headings were embedded as comments inside code cells; split into proper markdown cells and renumbered 5.1–5.4 → 6.1–6.4 to match the parent Section 6.
2. **Inflation check comment corrected** — comment previously said candidates *are* caught by the quality gate; corrected to say they are *not* (best_percentile ~0.5–0.55%, below the 2% threshold — GPS inflates from allele breadth, not per-allele strength).
3. **Allele list derived programmatically** — removed hardcoded patient-specific allele list; now auto-detected from `_presentation_score` column names, making the notebook reusable across patients.

---

### 21:28 UTC — Editor: Scientist

#### PR #167 — patient_002 results analysis notebook (feat/scientist/issue-164-patient-002-notebook)

Added `research/notebooks/patient_002_results.ipynb` with full scientific analysis of the post-#148 valid run:

- **Junction analysis:** 55,912 tumor-exclusive junctions (no matched RNA-seq normal; WES proxy filtered 3 junctions)
- **Peptide translation:** 703,106 (8-mer) / 781,159 (9-mer) / 846,422 (10-mer) = 2,330,687 total, 2,313,700 unique
- **MHC presentation:** 67,935 strong presenters (2.9%); HLA-C\*01:02 + HLA-C\*07:01 account for ~69% of strong-presenting candidates
- **GPS validation:** mean 0.101, median 0.026 — discriminating; 174 inflation edge cases (GPS > 0.9 but n_strong_alleles = 0, best_percentile ~0.5–0.55%)
- **Top candidate:** FADLRPLLL / HLA-C\*01:02 (IC50 = 33.2 nM, percentile = 0.0045%, GPS = 0.9999, n_strong = 4)

Open questions posted to Developer in team standup: min read support filter, GPS n_strong_alleles gate, HLA-C calibration, WES normal confirmation.

---

### 21:24 UTC — Editor: Scientist

#### PR #166 — Jupyter notebook environment setup (feat/scientist/issue-163-notebook-env)

Merged notebook environment scaffold for `research/notebooks/`:

- `research/requirements.txt` — lower-bound pins: jupyterlab ≥ 4, pandas ≥ 2, matplotlib ≥ 3.7
- `.python-version` added to `.gitignore` (personal env config, not committed); `pyenv local 3.14.4` kept in README as the explicit setup step
- `research/README.md` updated with setup instructions: pyenv → venv → pip install → VSCode kernel selection; GCS auth note

---

### 18:15 UTC — Editor: Scientist

#### Issues #162–165 — patient_002 results analysis (valid run, post Issue #148)

First scientific analysis session for patient_002 (BG003082, osteosarcoma). Full analysis in `research/notebooks/patient_002_results.ipynb`.

**Key findings:**

- Top candidate: **FADLRPLLL / HLA-C\*01:02** — IC50 33.2 nM, presentation_percentile 0.0045%, GPS 0.9999, strong presenter across 4/5 alleles.
- **HLA-C dominance:** HLA-C\*01:02 + HLA-C\*07:01 account for ~69% of all 67,935 strong-presenting peptides despite 0.5× locus weight in GPS. Three hypotheses: (a) broader HLA-C groove/promiscuity, (b) osteosarcoma splice motifs enriched for HLA-C, (c) MHCflurry percentile calibration less precise for HLA-C (sparser training data). Sent to Developer for investigation.
- **HLA-A\*01:01 nearly silent:** median presentation percentile 8.5% among strong presenters — contributes almost nothing to GPS for this patient.
- **GPS is discriminating:** mean 0.101, median 0.026; only 1.1% of 2.3M predictions exceed GPS > 0.9. Formula works as designed.
- **GPS inflation edge case (minor):** 174 candidates (0.7% of GPS > 0.9) have `n_strong_alleles = 0` with best_percentile ~0.5–0.55%. Quality gate (2%) does not catch these. Flagged to Developer to consider adding `n_strong_alleles ≥ 1` as a hard filter.
- **Minimum junction read support is 174** — suspiciously high. Awaiting Developer confirmation of whether a min-reads filter is applied.
- **Peptide translation:** 2,330,687 total (703,106 × 8-mer, 781,159 × 9-mer, 846,422 × 10-mer), 2,313,700 unique. 9-mer count differs by 265 from a prior run (781,424) — untracked discrepancy, motivates run registry (Issue TBD).

**Actions taken:**

- Created `research/notebooks/patient_002_results.ipynb` (Issue #164).
- Set up `research/.venv` (Python 3.14.4 via pyenv) + `research/requirements.txt` + `research/README.md` notebook section (Issue #163).
- Updated `research/manuscript/RESULTS.md` patient_002 section: MHC presentation + GPS, peptide translation, top candidate FADLRPLLL, TCRdock. Junction Filtering table on hold pending Developer confirmation of WES normal usage (Issue #165).
- Opened Issue #161 (local GCS cache for notebooks); drafted run registry issue for Developer review in standup.
- Terminology convention established: "presenters/presenting" replaces "binders/binding" throughout.

**Open / pending:**
- Developer standup: min read filter, GPS gate, HLA-C calibration, WES normal confirmation, run registry engineering input.
- RESULTS.md Junction Filtering: on hold.

---

### 17:03 UTC — Editor: Developer

#### Issue #148 — PR #152 review fixes (assemble_contigs.py + alignment.smk)

Two issues raised in code review addressed before merge:

1. **Silent failure mode fixed (`assemble_contigs.py`)** — `bedtools getfasta` stderr is now always logged as WARNING (was only logged on non-zero exit, so chromosome-not-found warnings were silently discarded). Added a skip-rate check at the summary: if >10% of junctions are skipped due to short sequences, the log is upgraded to WARNING with an explicit message pointing to UCSC/ENSEMBL chr naming as the likely cause.

2. **Docstring clarified (`alignment.smk`)** — `hisat2_download_index` rule docstring and inline comment updated from "GRCh38" (ambiguous) to "hg38 (UCSC naming)" to remove the ambiguity that caused Issue #148 in the first place.

---

### 15:25 UTC — Editor: Developer

#### Issue #148 — patient_002 rerun results (after hg38_tran HISAT2 index fix)

Full pipeline rerun on patient_002 after merging `fix/issue-148-hisat2-chr-naming` (PR #152 — https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/152), which switched the HISAT2 prebuilt index from `grch38_tran` (ENSEMBL naming, no `chr` prefix) to `hg38_tran` (UCSC naming, `chr` prefix) to match the genome FASTA. Previous results were invalid due to the chromosome naming mismatch causing `bedtools getfasta` to return empty sequences for 56,735 of 56,774 junctions.

**Junction filtering**

| Metric | Count |
|---|---|
| Unannotated junctions | 58,914 |
| Tumor-exclusive | 58,914 |
| Normal-shared | 0 |

No matched normal sample for patient_002 → all unannotated junctions labeled tumor_exclusive (expected behaviour with a warning in the pipeline log).

**MHC predictions**

| Class | Count |
|---|---|
| Total predictions | 2,330,687 |
| Strong binders (≤ 0.5%) | 67,935 |
| Weak binders (≤ 2.0%) | 222,823 |
| Non-binders (> 2.0%) | 2,039,929 |

**HLA typing (OptiType, tumor)**

| Locus | Allele(s) |
|---|---|
| HLA-A | A\*01:01 |
| HLA-B | B\*08:01 / B\*27:05 |
| HLA-C | C\*01:02 / C\*07:01 |

**Top neoepitope candidate**

| Field | Value |
|---|---|
| Peptide | FADLRPLLL |
| Best allele | HLA-C\*01:02 |
| Presentation percentile | 0.0045% |
| Genotype presentation score | 0.9999 |
| IC50 | 33.2 nM |
| Presentation class | strong |

TCRdock PDB available: yes.

---

### ~14:54 UTC — Editor: Developer

#### Issue #118 — CI dry-run with GPU config overlay (PR #159)

Added a `pipeline-snakemake-dry-run` CI job to `.github/workflows/tests.yml` that runs `snakemake -n --configfile config/config.yaml config/gpu.yaml` on every push/PR to main. This catches config key errors in `structure.smk` TCRdock rule blocks — previously invisible because `gpu.yaml` was never loaded in CI. The root motivation was the `KeyError: 'ic50_strong'` on the first cloud prod run (fixed in #117), which this job would have caught. Placeholder input files (`touch`) are created before the dry-run so Snakemake can build the DAG without real data. Also renamed the existing `test` job to `pipeline-pytest` for clarity, and updated the branch protection required status checks accordingly.

---

### ~10:01 UTC — Editor: Scientist

#### Issue #153 — Add --update-note flag to zotero_add.py

Added `update_note()` helper and `--update-note ITEM_KEY` flag to `research/scripts/zotero_add.py`. Allows updating the note on an existing Zotero entry without re-adding the paper — needed for the morning reading routine when the note format evolves after initial add (e.g. adding a Limitations section). Uses optimistic locking via `If-Unmodified-Since-Version`. Falls back to creating a new note if none exists.

---

### ~09:30 UTC — Editor: Scientist

#### Issue #154 — Move zotero_add.py to research/scripts/ and add research/README.md

Moved `scripts/zotero_add.py` to `research/scripts/zotero_add.py` to clarify that the script is a research support tool, not part of the bioinformatics pipeline. Created `research/README.md` with a folder overview and full usage docs for `zotero_add.py` including the three-section note format convention (Results / Methods / Limitations).

---

### ~08:30 UTC — Editor: Developer

#### Issue #148 — patient_002 full rerun completed (first-pass success)

Full rerun of patient_002 completed successfully after merging the `hg38_tran` HISAT2 index fix (Issue #148, branch `fix/issue-148-hisat2-chr-naming`). All Snakemake steps finished without errors. At first view the results look good — junction counts and peptide counts appear in the expected range, unlike the previous invalid run (which produced only 783 peptides due to the ENSEMBL/UCSC chr naming mismatch). Detailed result review and top candidate recording deferred to next session.

---

## 2026-04-25

### ~16:56 UTC — Editor: Developer

#### Issue #136 — VM data cleanup + NVIDIA driver fixes for common-cu129 image

**Context:** patient_002 run restarted after morning disk resize. Left running while user went to lunch (~12:00 UTC).

**NVIDIA driver failures — three iterations (~13:00–15:00 UTC)**

The `ubuntu-accelerator-2204-amd64-with-nvidia-570` Deep Learning VM image was retired by Google; the replacement is `common-cu129-ubuntu-2204-nvidia-580`. The new image ships driver 580, which lacks GSP firmware support for P100 (Pascal SM 6.0). Three driver downgrade attempts were needed before a working configuration was found:

| Attempt | Package | Failure reason |
|---------|---------|----------------|
| 1 | `nvidia-driver-570-server` | Pulls display dependencies (`libgbm1`, `libxcb-*`, `libnvidia-egl-wayland1`) absent on headless VMs |
| 2 | `nvidia-headless-no-dkms-570-server` | Pre-compiled kernel modules don't match GCP kernel (`6.8.0-1053-gcp`) on `common-cu129` image — `nvidia-smi` fails with "driver not loaded" |
| 3 | `nvidia-headless-570-server` (DKMS) | **Works** — DKMS recompiles modules against the running kernel at install time |

`IMAGE_FAMILY` in `run_cloud_gpu.sh` updated to `common-cu129-ubuntu-2204-nvidia-580`. Driver downgrade block updated to install the DKMS variant. Three commits on `feat/issue-136-vm-data-cleanup`; PR #151 opened.

**patient_002 pipeline completed (~15:00–15:46 UTC)**

TCRdock rerun manually with `snakemake --forcerun run_tcrdock` after DKMS driver confirmed active. All 22 Snakemake steps finished. Results and logs uploaded to GCS at ~15:46 UTC.

Top candidate: **NSISRPSSL / HLA-C\*01:02, ic50\_nM=39.77, pLDDT=91.99**

**Issue #148 discovered — results invalid: chromosome naming mismatch (~16:00 UTC)**

Report showed only ~800 total neoepitopes (expected thousands). Root cause: `hisat2_prebuilt_url` in `config/config.yaml` used `grch38_tran.tar.gz` (ENSEMBL naming, no `chr` prefix) while the genome FASTA uses UCSC naming (`chr` prefix). `bedtools getfasta` returned empty sequences → 56,735 of 56,774 junctions skipped in `assemble_contigs.py` → only 18 contigs → 783 peptides. patient_001 was unaffected as it predates the config key and built the index locally from the UCSC FASTA. Fix: `hg38_tran.tar.gz` on `fix/issue-148-hisat2-chr-naming`. **Today's patient_002 results are invalid; full rerun required after Issue #148 is merged.**

---

### ~15:00 UTC — Editor: Scientist

#### TCR-pMHC binding prediction — field overview and structural improvement plan (Issue #86)

**Field overview:**

Two main paradigms:
- **Sequence/ML-based** (NetTCR-2.x, pMTnet, ERGO, MixTCRpred, TULIP): treat TCR+peptide as a sequence classification problem; generalisation to unseen epitopes remains a key limitation across all methods — to revisit in a future Scientist session
- **Structure-based** (TCRdock, AlphaFold3): model the 3D TCR-pMHC complex; TCRdock (Alam et al., *Science* 2023) is already integrated; AF3 (Abramson et al., *Nature* 2024) is a newer competitor

Key databases: VDJdb (curated TCR-pMHC specificity), IEDB, McPAS-TCR.

**Three axes for improving the structural approach — agreed execution order:**

1. **Better inputs (Issue #86):** patient-specific HLA-matched TCR panel from VDJdb — prerequisite for all downstream improvements
2. **Model upgrade:** benchmark AlphaFold3 vs. TCRdock on a known-binding panel once Axis 1 is in place
3. **Rescoring:** post-process PDB outputs with Rosetta `InterfaceAnalyzer` or FoldX `AnalyseComplex` for interface ΔΔG as a secondary ranking signal

**Issue #86 — VDJdb TCR panel design (first-pass, conservative):**

| Parameter | Value |
|---|---|
| HLA matching | Exact 4-digit |
| MHC class | Class I only |
| Paired α/β chains | Required |
| VDJdb confidence score | ≥ 2 |
| Panel size | Top 10 per allele |
| Redundancy reduction | None |
| Antigen category filter | None |

Key non-obvious dependency: **`stitchr`** (Peacock et al. 2022, *Bioinformatics*) — VDJdb stores CDR3 + V/J gene assignments only; TCRdock needs full α/β sequences; stitchr reconstructs them from IMGT germline references.

Pipeline integration: new Snakemake rule `fetch_vdjdb_panel` (input: patient HLA alleles; output: `results/{sample}/tcrdock/vdjdb_panel.tsv`) feeds into existing `run_tcrdock`. VDJdb data from `antigenomics/vdjdb-db` flat TSV (pinned release). Issue #86 body updated to reflect this design.

---

### ~12:00 UTC — Editor: Scientist

#### Zotero integration for morning science reading routine (Issue #137, PR #138)

Established daily science reading habit for Scientist sessions. Each morning warm-up ("good morning") now produces a Zotero entry rather than a markdown log.

**Setup:**
- Created Zotero collection "Splice Neoepitope Pipeline" (key `Z38GTJNW`, library `lee.jin-ho`, user ID 10082130)
- `scripts/zotero_add.py`: CLI tool to add a paper by DOI — fetches metadata from CrossRef (authors, journal, ISSN) with PubMed fallback for full date, structured abstract, and PMID; supports `--note` (child note), `--tags`, `--dry-run`
- Credentials in `.env` (gitignored)

**First entry:**
Weber et al. (2024) KEYNOTE-942, *The Lancet* 403:632–644. DOI: 10.1016/S0140-6736(23)02268-7. Phase 2b trial: personalised mRNA neoepitope vaccine (mRNA-4157/V940) + pembrolizumab → ~44% reduction in recurrence/death vs. pembrolizumab alone in resected melanoma. Targets SNV/indel neoantigens only; splice-junction neoepitopes (our focus) are absent — directly motivates this project.

**PR #138 merged.**

---

### ~10:30 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: sub-issue retrospective, VM disk fix, branch rebase

**patient_002 run failure — root cause investigation (~09:30 UTC)**

patient_002 run launched the previous evening (~22:12 UTC, 2026-04-24) got stuck overnight.
SSH inspection of the VM showed `samtools sort` hanging with no output: the HISAT2 index was
missing (`resources/hisat2_index/genome_tran` not found), causing HISAT2 to exit immediately.
`samtools sort` was reading from HISAT2's stdout pipe and hung instead of propagating the error.
`set -euo pipefail` was insufficient because the right-hand side of the pipe does not exit when
the left-hand side fails — it just blocks waiting for more input.

Immediate fix: killed the tmux session, deleted the stale `resources/hisat2_index/` directory
(contained `genome.*.ht2` from an old build, but `index.done` sentinel was present so Snakemake
skipped re-download). Pipeline restarted by user after disk issue resolved (see below).

Code fix (Issue #131, commit `e2b1a65`): added a pre-check at the top of `hisat2_align` shell
block — exits with a clear error and writes to `{log}` before `samtools` is ever launched:
```bash
if [[ ! -f "{params.index_prefix}.1.ht2" ]]; then
    echo "ERROR: HISAT2 index not found at {params.index_prefix}.*.ht2" | tee -a {log}
    exit 1
fi
```
This ensures the error propagates to `pipeline.log` and the orchestrator can detect it.

**VM disk full**

VM (100 GB SSD) was full — patient_001 data had not been cleaned up post-GCS upload. 300 GB
resize rejected (europe-west1 SSD quota is 250 GB; `neoepitope-pipeline` 200 GB +
`orchestrator` 10 GB = 210 GB in use). Resized to 200 GB (`gcloud compute disks resize`),
then `sudo resize2fs /dev/sda1` to extend the filesystem. patient_001 data deleted from VM
(already safe in GCS).

**Sub-issue retrospective — PRs #132–#135 merged**

10 commits accumulated on `feat/issue-123-prod-cloud-run` since patient_001 re-run were
retrospectively organised into 4 focused sub-issues, each cherry-picked to its own branch
and merged into main independently:

| Issue | PR | Commit(s) | Description |
|-------|----|-----------|-------------|
| #128 | #132 | `6dc0fdc` | fix(cloud): upload pipeline.log + orchestrator.log to GCS |
| #129 | #133 | `8d3ff15` | fix(tcrdock): match `_pdb_file` suffix to avoid pdbid false-match |
| #130 | #134 | `57c1f3f` | fix(cloud): remove `--conda-cleanup-envs` (deletes all envs on empty DAG) |
| #131 | #135 | `e2b1a65` | fix(alignment): fail fast on missing HISAT2 index |

`feat/issue-123-prod-cloud-run` now contains only 3 docs commits
(`bbd88c5`, `f2eb363`, `86e7e15`) on top of the latest main. Branch rebased and
force-pushed. PR for this branch will be opened after patient_002 run completes.

**patient_002 run status:** restarted after disk resize; in progress.

---

## 2026-04-24

### ~22:12 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: patient_002 started (WES normal, interim)

Launched production pipeline run for patient_002 (osteosarcoma IPISRC044, BG003082) via
`run_cloud_gpu.sh --detach` on `neoepitope-pipeline` VM. Branch: `feat/issue-123-prod-cloud-run`.

Normal sample: Blood Derived WES (`BG003082_N0_WES`) — interim proxy; yields near-zero junction
overlap (~3 junctions). All unannotated tumour junctions will be labeled `tumor_exclusive` with
a pipeline warning. Proper normal (Jan 2025 CD3+ PBMC Cell Ranger BAM) pending Issue #127
implementation.

Also includes two bug fixes committed since patient_001 re-run:
- `2625ec7` fix(tcrdock): match `_pdb_file` suffix in column search (avoids `pdbid` false-match)
- `117d174` fix(cloud): remove `--conda-cleanup-envs` (was deleting all 4 envs when DAG empty)

Issue #127 opened: support pre-aligned BAM as normal input for junction filtering.

---

### ~22:00 UTC

#### Issue #123 / #126 — patient_002 normal sample decision + GTEx filter design (Researcher session)

**Context:** patient_002 (osteosarcoma IPISRC044) has no matched RNA-seq normal. Current
production run uses Blood WES as proxy, yielding effectively no junction filtering (3 junctions
overlapping tumour set — consistent with WES spliced reads being alignment artefacts).

**Normal sample decision for current production run:**

Longitudinal single-cell RNA-seq (PBMC, sorted CD3+ T cells) is available from the Hudson Lab:

- Protocol: 10x Chromium (confirmed by Cell Ranger output + Seurat objects)
- Time points: Jan 2025 – Dec 2025 (growing dataset)
- Location: `hudson_lab/PBMC_scRNAseq/FASTQ` and `hudson_lab/PBMC_scRNAseq/cellranger`

Decision: use **Jan 2025 Cell Ranger BAM** as the normal input for the current patient_002 run.
Rationale:
- CD3+ T cells are the primary TIL population contaminating solid tumour RNA-seq; filtering their
  junctions directly targets the most relevant contamination source
- Jan 2025 is the earliest (most treatment-naive) time point — least treatment-modified T cell state
- 10x 3' capture gives sparse junction coverage, but still produces genuine biological splice junctions
  unlike WES (which gave near-zero overlap)
- This is an interim measure pending GTEx pan-tissue filter implementation (issue #126)

Developer action: run regtools on the Jan 2025 Cell Ranger BAM (already genome-aligned; no
re-alignment needed) and use resulting junction set as `normal_junctions` for patient_002.

**GTEx pan-tissue filter — scientific rationale documented (issue #126):**

A pan-tissue GTEx filter (not tissue-matched) is the correct long-term approach for a vaccination
application. Key argument: vaccine-induced CTLs are systemic — they patrol all tissues, not only
the tumour site. A junction present in any normal tissue could be presented on those cells, creating
off-tumour autoimmune toxicity risk. Pan-tissue filtering serves dual purpose:

1. Safety — excludes junctions expressed in any of ~54 GTEx tissue types
2. Candidate quality — pan-tissue absence is a stronger tumour-exclusivity claim; precision over
   recall is appropriate given 10–20 vaccine slot constraint

Full rationale added to `research/manuscript/DISCUSSIONS.md` (section "Normal sample filtering →
GTEx pan-tissue filter"). Issue #126 opened for Developer implementation.

---

### ~18:30 UTC

#### Issue #123 — M1 production cloud run: patient_001 started

Launched production pipeline run for patient_001 (gastric cancer, SRR9143066 T / SRR9143065 N) via `run_cloud_gpu.sh --detach` on `neoepitope-pipeline` VM (n1-highmem-8 + P100, europe-west1-b). Branch: `main`. patient_002 run to follow after patient_001 completes.

---

### ~17:30 UTC

#### Issue #124 — CLI flag polish (Developer session, PR #125, branch `feat/issue-124-cli-flag`)

**Goal:** Expose `--presentation-percentile-weak` via CLI parsers in both `generate_report.py` and `run_tcrdock.py`; fix docstring omission.

**Done:**

- Added `--presentation-percentile-weak` (type=float, default=2.0) to `_cli_main()` in both scripts and threaded through to the respective `generate_report()` / `run_structural_validation()` calls.
- Fixed `generate_report()` docstring: `presentation_percentile_weak` was present in the function signature but absent from the `Args:` section.
- 4 new `TestCLIParser` tests (default + custom value, one class per script). 200 tests passing.

**Commits (2, branch `feat/issue-124-cli-flag`):**
- `6cea0dd` feat(cli): add --presentation-percentile-weak flag to generate_report and run_tcrdock CLI parsers
- `7c6dd35` test(cli): add CLI parser tests for --presentation-percentile-weak flag

**PR #125 opened.** Awaiting merge.

---

### ~16:42 UTC

#### Issue #119 — PR #122 re-review fix + merge prep

**Re-review verdict:** Ready to merge. Three minor observations: (1) `presentation_notes` strings in `_build_report_tsv` still hardcoded `"2%"` for `weak`/`non` entries despite `presentation_percentile_weak` parameter being available; (2) CLI parsers missing `--presentation-percentile-weak` flag; (3) `generate_report()` docstring missing the parameter.

Fixed (1) now — `weak`/`non` entries converted to f-strings using `presentation_percentile_weak` (`18e1fac`). Items (2) and (3) deferred to follow-up issue.

**29 tests passing** (generate_report suite).

---

### ~16:15 UTC

#### Issue #119 — PR #122 review cycle (Developer session)

**Reviewer comments (6):** stale docstring, GPS out-of-range clamping, hardcoded quality-gate threshold, GPS test gap in `generate_report.py`, hardcoded 9-mer `end_nt_incl`, misleading log message on quality-gate rejection.

**Key design decisions made during review:**
- `hla_c_weight` outside [0,1] → hard `ValueError` (own config, must be correct).
- MHCflurry `presentation_score` outside [0,1] → warning only, no clamp (third-party output; out-of-range values cannot be reliably interpreted and should be fixed upstream).
- "binder" → "presenter" terminology throughout `run_tcrdock.py` and tests.

**Fix commits (8):**
- `4366af0` fix(mhc): fix stale docstring and add hla_c_weight validation + score warning
- `9f0b6c6` test(mhc): add tests for hla_c_weight validation and out-of-range score warning
- `b73b27b` fix(rules): thread presentation_percentile_weak into report and tcrdock params
- `94a480a` fix(report): replace hardcoded quality-gate 2.0 with configurable presentation_percentile_weak
- `34f7aa5` fix(tcrdock): replace hardcoded quality-gate 2.0 with configurable presentation_percentile_weak
- `1b08ebe` test(report): add GPS quality-gate and ranking tests for generate_report.py
- `fe6e188` fix(report): compute end_nt_incl from peptide length instead of hardcoded 9-mer
- `d0e7456` fix(tcrdock): distinguish no-presenters from quality-gate-empty in log messages

**196 tests passing.** Re-review requested.

---

### ~14:38 UTC

#### Issue #119 — Two-level MHC presentation prediction: implementation (Developer session, PR #122, branch `feat/issue-119-allele-breadth`)

**Goal:** Implement the allele breadth model designed in the Researcher session (~11:30 UTC) and open PR #122.

**Done:**

- Replaced `_compute_strong_alleles()` in `run_mhcflurry.py` with `_compute_per_allele_features()`: one `Class1PresentationPredictor.predict()` call per allele (single-element genotype list) to recover per-allele `presentation_score` and `presentation_percentile` discarded by the genotype API.
- New `genotype_presentation_score` formula: `1 − ∏(1 − wᵢ·pᵢ)`, with `w(HLA-A) = w(HLA-B) = 1.0`, `w(HLA-C) = hla_c_weight` (default 0.5). New supporting columns: `n_strong_alleles`, `best_presentation_percentile`.
- `allele` column renamed to `best_allele` throughout (`run_mhcflurry.py`, `generate_report.py`, `run_tcrdock.py`, tests).
- Updated ranking in `generate_report.py` and `run_tcrdock.py`: GPS ↓ → `n_strong_alleles` ↓ → `best_presentation_percentile` ↑. Quality gate: `best_presentation_percentile > 2%` excluded from top-candidates list.
- Added `mhcflurry.hla_c_weight: 0.5` to `config/config.yaml` and propagated through `mhc_affinity.smk` params.
- `METHODS.md` Section 6 rewritten: two-level architecture, GPS formula in LaTeX, quality gate, ranking rule, full output column table.
- `docs/configuration.md`: `hla_c_weight` documented with HLA-C surface density rationale.
- Test suite: 186 passed, 1 skipped (stale integration pipeline output — schema-version skip guard added).

**Key implementation detail:** The genotype API (`predict()` with all alleles at once) reports only the best-allele scores and silently discards all others. Per-allele scores must be recovered via separate single-allele calls. This dual-call design is why `_compute_per_allele_features` iterates over `resolved_alleles` individually.

**Column name change mid-session:** `breadth_score` → `genotype_presentation_score` after Researcher consultation. All code, tests, and docs use the new name.

**Local test run (chr22, patient_001_test):** Full pipeline completed. `mhc_presentation.tsv` has 24 columns (6 allele pairs + 3 breadth columns). Top hit: DVFGTPFSR / HLA-A\*33:03, GPS = 0.9996, best-allele percentile = 0.001%.

**Commits (9, all on `feat/issue-119-allele-breadth`):**
- `530ddfc` docs(manuscript): rename breadth_score → genotype_presentation_score
- `a13fd30` feat(config): add mhcflurry.hla_c_weight config key and Snakemake param
- `38b6ac8` feat(mhc): add per-allele genotype_presentation_score model (Issue #119)
- `49b8e96` feat(report): update ranking to genotype_presentation_score, add quality gate
- `b90c75d` feat(tcrdock): update candidate selection to genotype_presentation_score ranking
- `291fd23` test(mhc): update tests for new output schema and add breadth feature tests
- `3d8b903` test(report,tcrdock): update test data to best_allele; add breadth-aware tcrdock tests
- `260da3a` test(integration): update required columns for new schema with skip guards
- `f792ac1` docs(methods): update Section 6 for two-level MHC presentation prediction
- `ac89bd2` docs(config): document hla_c_weight and update mhcflurry section

---

### ~11:30 UTC

#### Issue #119 — Allele breadth model: scientific design (Researcher session, branch `feat/issue-119-allele-breadth`)

**Goal:** Develop and document the scientific rationale for a multi-allele breadth scoring
model for `mhc_presentation.tsv`. No code changes — manuscript only.

**Key scientific decisions:**

- **Two-level architecture:** MHCflurry is a molecular-level predictor (one peptide × one
  MHC allele → `presentation_score`). The breadth model is a separate genotype-level
  combiner: `breadth_score = 1 − ∏(1 − wᵢ·pᵢ)`. These are distinct modelling layers.
- **HLA-C weight:** enters at the genotype-combination level (surface density ~50% of
  HLA-A/B), not at the molecular prediction level. Configurable (`w_C`, default 0.5).
- **`presentation_score` in formula (not `presentation_percentile`):** `presentation_score`
  is a calibrated absolute probability — the correct input for a complementary probability
  formula. `presentation_percentile` is a rank statistic; using it would require an
  arbitrary mapping. The two metrics can disagree (high breadth_score with all alleles in
  non-binder percentile territory is possible for promiscuous alleles) — this is the
  intended behaviour of the two-tier system.
- **Committed ranking for cancer vaccine application:**
  1. `breadth_score` — primary (multi-allele coverage, LOH robustness, vaccine slot efficiency)
  2. `n_strong_alleles` — secondary (count of alleles at ≤ 0.5% percentile)
  3. `best_presentation_percentile` — minimum quality gate only (not a ranking dimension);
     filters candidates where no allele reaches weak-binder territory
- **Immunodominance acknowledged but not dominant in vaccine context:** in natural
  anti-tumour immunity, one very strong allele can dominate via intramolecular MHC groove
  competition and immunodomination (Yewdell & Bennink, *Annu Rev Immunol* 1999; Chen &
  McCluskey, *Adv Cancer Res* 2006). In therapeutic vaccination, immunodomination is largely
  bypassed; HLA LOH and vaccine slot efficiency make breadth the primary criterion.

**Manuscript changes (committed `c9d5eb3`, branch `feat/issue-119-allele-breadth`):**
- `INTRODUCTION.md`: new section "HLA Genotype, Surface Expression, and Allele Breadth"
- `DISCUSSIONS.md`: new section "Allele breadth and immunodominance: two complementary
  ranking signals" (subsections: two-level architecture, immunodominance mechanisms,
  vaccine application with committed table, calibration note)

**METHODS.md:** not updated — Developer session responsibility after implementation.

**Design spec for Developer session:** `memory/project_allele_breadth_design.md`

---

### ~10:02 UTC

#### Issue #115 — Rename `mhc_affinity.tsv` → `mhc_presentation.tsv` (PR #121, branch `feat/issue-115-rename-mhc-affinity-tsv`)

**Motivation:** The output file was named after `Class1AffinityPredictor`, but the pipeline now uses `Class1PresentationPredictor` (since Issue #85). The name `mhc_affinity.tsv` was misleading — the primary ranking metric is `presentation_percentile`, not IC50.

**Changes:** Pure mechanical rename across 9 files — output path in `mhc_affinity.smk`, Snakemake output key `mhc_affinity_tsv` → `mhc_presentation_tsv` in `mhc_affinity.smk`, `run_mhcflurry.py`, `structure.smk`; path references in `analysis.smk`, `generate_report.py`, `run_tcrdock.py`, `test_integration.py`; docs in `METHODS.md`, `configuration.md`, `README.md`. The `.smk` rule module filename (`mhc_affinity.smk`) is unchanged — it is the rule module for the MHC prediction step, not named after the output.

**Testing:** 25/25 integration tests pass. Local test result file renamed from `mhc_affinity.tsv` → `mhc_presentation.tsv`.

---

### ~08:49 UTC

#### Hotfix — rename `ic50_strong` → `presentation_percentile_strong` in `structure.smk` (PR #117, branch `fix/structure-smk-ic50-key`)

**Problem:** First prod cloud run after Issue #85 failed immediately with `KeyError: 'ic50_strong'` in `workflow/rules/structure.smk` line 77. The `generate_report_with_structure` rule was still referencing the removed config key.

**Root cause:** The TCRdock rule block in `structure.smk` is guarded by `if _TCRDOCK_ENABLED:`, which is only true when `gpu.yaml` is merged in. Local tests and CI never enable TCRdock, so the stale key was never evaluated during development.

**Fix:** One-line rename of `ic50_strong` → `presentation_percentile_strong` in the `params:` block of `generate_report_with_structure`, matching the already-corrected `analysis.smk`.

**Follow-up filed:** Issue #118 — add a Snakemake dry-run with `gpu.yaml` overlay to CI to catch this class of bug automatically.

**Prod cloud run** (patient_001, full genome, GPU): completed successfully after fix. UTC: 2026-04-24T00:12:39Z. Key results vs. baseline (pre-Issue #85):
- Junction counts unchanged: 30,029 unannotated → 27,348 tumor-exclusive
- Total predictions: 7,718,952 → 1,286,492 (6× reduction — one row per peptide, not per peptide×allele)
- Strong presenters: 15,880 (IC50 < 50 nM) → 44,916 (presentation_percentile ≤ 0.5%)
- Top candidate changed: EVAETLSLF/HLA-A\*26:01 (IC50 16.6 nM) → FAFPFAQTL/HLA-C\*03:03 (percentile 0.0007%)

---

### ~00:30 UTC

#### Issue #85 — Switch to Class1PresentationPredictor (PR #114, branch `feat/issue-85-prediction-mode`)

**Goal:** Replace `Class1AffinityPredictor` with `Class1PresentationPredictor`, which integrates binding affinity with antigen processing (proteasomal cleavage, TAP transport) trained on mass-spec MHC ligand data.

**API surprises encountered:**
1. `predict_to_dataframe()` does not exist on `Class1PresentationPredictor` (affinity predictor only) — use `predict()` instead.
2. `Class1PresentationPredictor.predict()` is a *genotype-level* API: takes all patient HLA alleles at once (≤6) and returns **one best-allele prediction per peptide**. Passing the same allele repeated N times (the `AffinityPredictor` convention) raises `ValueError: alleles list must have at most 6 elements`.

**Design decisions:**
- Dropped `affinity_percentile` / `binder_class`: `predict()` does not return `affinity_percentile`, and a second inference pass would double runtime with no clear biological gain.
- Single `presentation_class` label from `presentation_percentile`: strong ≤ 0.5%, weak ≤ 2%, non > 2% — thresholds from Jiang et al. 2024 (*Communications Biology*).
- Merged `_load_predictor_for_gpu()` / `_load_predictor_for_cpu()` into `_load_predictor()`: both were identical after parallel-worker removal.

**Final output schema (`mhc_affinity.tsv`):** `contig_key, start_nt, peptide, allele, ic50_nM, processing_score, presentation_score, presentation_percentile, presentation_class`. One row per peptide; `allele` = best presenter in patient genotype.

**Local test run** (patient_001_test, chr22, 6-allele genotype): 4045 unique peptides → 147 strong / 401 weak / 3571 non. UTC: 2026-04-23T21:43:51Z.

182 tests pass (unit + integration).

---

## 2026-04-23

### ~20:00 UTC

#### Issue #110 — Fix stale integration test assertions (PR #113, branch `fix/issue-110-contig-length-test`)

**Goal:** Update two integration test assertions that became stale after PR #99 changed flank size to `3 × (max_length − 1)` and introduced multi-length peptides.

**Changes:**
- `test_all_contigs_are_50nt` → `test_all_contigs_are_54nt`: contig length is now 54 nt (= 27 + 27 flanks with `peptide_lengths=[8, 9, 10]`, `max_length=10`)
- `test_all_peptides_are_9mers` → `test_all_peptides_are_8_to_10mers`: pipeline produces 8-, 9-, and 10-mer peptides per config

No functional code changes; tests updated to match current pipeline behaviour.

**Note on branch base:** this branch was originally created off `feat/issue-97-cds-reading-frame` before #112 merged, making the PR appear to contain ~300 lines of #97 reading frame code. After #112 merged to main, the branch was rebased onto main — the diff is now exactly the 4-line test fix.

---

### ~19:30 UTC

#### PR #112 review fixes — issue #97

- Removed redundant `gene_cds` intermediate `defaultdict` from `_build_cds_donor_lookup`; donor frames now computed in a single GTF parsing pass (halves peak memory on large GTFs)
- Fixed misleading comment on GTF end coordinate: `# GTF end is 1-based inclusive = 0-based exclusive (no adjustment needed)`
- Added `test_minus_strand_reading_frame_annotated` to `TestClassifyJunctionsReadingFrame` — covers the `donor_coord = junction.end` path through the full `classify_junctions` pipeline for a − strand junction
- 157 tests pass

### ~18:30 UTC

#### Issue #97 — CDS reading frame annotation for novel junctions (branch `feat/issue-97-cds-reading-frame`)

**Goal:** Annotate each tumor-exclusive junction in `novel_junctions.tsv` with its canonical CDS reading frame, derived from the GENCODE GTF. This supports downstream biological interpretation without restricting translation.

**Key design decisions:**

1. **Annotation only, not restriction.** The `reading_frame` column is informational metadata. All three sense-strand frames are still translated downstream. Restricting translation to the CDS-derived frame would introduce false negatives in hypermutated tumors (MSI-high, POLE-mutant) where upstream frameshift indels, SVs, or chained novel junctions shift the active frame — precisely the tumors most likely to harbour actionable junction neoepitopes. These frame-altering events cannot be inferred from RNA-seq junction data alone.

2. **Protein-coding transcripts only.** CDS records are filtered to `transcript_type "protein_coding"` in the GENCODE GTF. NMD and other non-translated isoforms are excluded to avoid spurious frame assignments from transcripts that are not translated by the ribosome.

3. **Union of frames across transcripts.** When a donor site is annotated in multiple protein-coding transcripts with different reading frames, all attested frame offsets are retained (e.g. `"0,2"`). This is strictly more informative than falling back to all three frames when ambiguous.

4. **No chained frame propagation.** If a second upstream novel junction shifts the reading frame before reaching the junction of interest, the propagated frame could in principle be computed as `(phase_B + L) % 3`. However, phasing two junctions to the same transcript is impossible from short-read RNA-seq without transcript assembly — left for a future improvement.

5. **Sense-strand translation only.** For strand-specific RNA-seq (dUTP / first-strand), `bedtools getfasta -s` already reverse-complements minus-strand contigs so they represent the correct transcript sequence. Antisense translation of a strand-corrected contig has no established biological basis. Confirmed strand-specific for patient_001 (KAPA RNA HyperPrep Kit with RiboErase HMR, dUTP second-strand marking). **patient_002 strandedness not yet verified.**

**Frame offset formula:**
```
phase_at_donor = (exon_length − gtf_frame) % 3
frame_offset   = (−phase_at_donor) % 3          → 0, 1, or 2
```
`upstream_nt` is always divisible by 3 by pipeline design, so `frame_offset` equals the start position within the upstream codon context (0 = in-frame, 1 = +1 shift, 2 = +2 shift). Donor coordinate: `junction.start` for + strand, `junction.end` for − strand.

**Changes:**
- `workflow/scripts/filter_junctions.py`: new `_build_cds_donor_lookup(gtf_path)` + `reading_frame` column in output TSV; `classify_junctions()` accepts optional `gencode_gtf` arg; `--gencode-gtf` CLI flag added
- `workflow/rules/filter_junctions.smk`: `gencode_gtf` added as explicit input to `filter_junctions` rule
- `research/manuscript/METHODS.md` §5: updated to describe reading frame annotation subsection
- `research/manuscript/DISCUSSIONS.md`: new section "Reading frame annotation: why translation is not restricted to the canonical frame"
- `workflow/tests/test_filter_junctions.py`: 13 new tests (9 × `_build_cds_donor_lookup`, 4 × `reading_frame` column integration); all 156 tests pass

### ~13:30 UTC

#### Issue #16 — Pre-built HISAT2 index (PR #111)

**Goal:** Replace the 60–90 min `hisat2-build` step with a download of the official pre-built `grch38_tran` index (~10–15 min).

**Changes:**
- Added `alignment.hisat2_prebuilt_url` config key; production config points to `grch38_tran.tar.gz` on the HISAT2 S3 mirror; test config sets it empty to keep building from the chr22 FASTA
- Replaced `hisat2_index` rule with conditional `hisat2_download_index` (URL set) / `hisat2_index` (URL empty) branching; introduced `_HISAT2_INDEX_PREFIX` (`genome_tran` vs `genome`) used by `hisat2_align`
- `grch38_tran` includes GENCODE splice sites baked in — strictly better than the plain genome index previously built from scratch

**Validation:**
- Tarball URL confirmed reachable; structure verified locally (`tar -tz`): top-level `grch38_tran/` with `genome_tran.{1..8}.ht2`
- Snakemake dry runs confirmed correct rule selection for both empty and non-empty URL
- Cloud extraction test on `hisat2-index-test` (e2-standard-4, europe-west1-b): all 8 `.ht2` files extracted correctly via `--strip-components=1`
- 143 unit tests pass

**Post-review fixes:** wrapped `curl | tar` in subshell for correct log capture; added `--retry 3`; added `resources: mem_mb=1000`; added re-download comment in config.

### ~12:52 UTC

#### Known limitation — contig assembly uses reference genome, not patient reads

`assemble_contigs.smk` uses `bedtools getfasta` on GRCh38 to extract flanking sequences around novel junction coordinates. It does not assemble from patient RNA-seq reads.

**Implication:** somatic SNVs or indels in the exon flanks are ignored. Predicted peptide sequences assume wild-type exonic context around the junction.

**Why acceptable for now:** the neoepitope signal is primarily from the novel exon–exon junction itself. Flank mutations are a second-order effect; reference extraction is standard in comparable published pipelines.

**Future improvement:** in hypermutated tumors (MSI-high, POLE-mutant), flank mutations could meaningfully alter predicted peptides. A future direction would be to extract junction-spanning reads from the BAM and assemble the local haplotype directly.

### ~12:41 UTC

#### Issue #107 — Single-VM consolidation merged (PR #109)

**Goal:** Eliminate the two-VM / three-phase architecture (CPU alignment VM + GPU TCRdock VM with GCS handoff) in favour of a single consolidated VM.

**Changes:**
- Deleted `scripts/setup_cloud.sh` and `scripts/setup_tcrdock_vm.sh`; replaced with unified `scripts/setup_vm.sh` (8-step idempotent setup: deps → GPU check → Docker + NVIDIA Container Toolkit → Miniforge3 → snakemake env → repo clone → TCRdock image build → reference data)
- Rewrote `scripts/run_cloud_gpu.sh` to manage a single `neoepitope-pipeline` VM (n1-highmem-8 + P100, 200 GB); Snakemake runs alignment → MHCflurry (GPU) → TCRdock (Docker/GPU) → report in one session with no GCS handoff
- Renamed `config/tcrdock_gpu.yaml` → `config/gpu.yaml` (overlay now enables MHCflurry GPU acceleration in addition to TCRdock); introduced `GPU_CONFIG_FILE` variable in `run_cloud_gpu.sh`; factored out mode-independent settings (`RESULTS_DIR`, `GCS_PATH`) from `case` block
- Fixed `--conda-cleanup-envs` call missing the GPU overlay; added `nvidia-container-toolkit` apt install fallback; forwarded `--keep-vm` through detach re-invocation

**Validated:** end-to-end cloud run on `neoepitope-pipeline` VM (2026-04-23, alignment → MHCflurry 7.72M predictions → TCRdock → HTML report).

### 09:40 UTC

#### Issue #105 — GPU MHCflurry validation (patient_001)

**Goal:** Validate that MHCflurry runs correctly on the P100 GPU and measure the speedup vs CPU baseline.

**Run:** patient_001 (SRR37781424, Luminal A breast cancer) on `neoepitope-pipeline` VM (n1-highmem-8 + P100, europe-west1-b).

| Timestamp (UTC) | Event |
|---|---|
| 09:16:25 | Snakemake started |
| 09:18:42 | Conda env activated, alleles loaded |
| 09:18:46 | GPU detected — sequential execution path selected |
| 09:18:54 | HLA-A\*31:01 started |
| 09:21:40 | HLA-A\*26:01 started (2m46s) |
| 09:24:24 | HLA-B\*18:01 started (2m44s) |
| 09:27:09 | HLA-B\*15:63 started (2m45s) |
| 09:29:53 | HLA-C\*07:01 started (2m44s) |
| 09:32:37 | HLA-C\*03:03 started (2m44s) |
| 09:36:15 | All 6 alleles complete (3m38s for last) |
| 09:36:33 | Report generated, pipeline done |

**Results:**

- Total predictions: 7,718,952 (1,260,074 unique peptides × 6 alleles)
- Strong binders (IC50 ≤ 50 nM): 15,880
- Weak binders (IC50 ≤ 500 nM): 149,662
- **Total MHCflurry time: ~17m21s vs CPU baseline ~1.5 hours → ~5.2× speedup on P100**

**Why not 10–50× as estimated?** Sequential allele execution (one CUDA context, one model in VRAM) — can't parallelise alleles because each would require a separate model copy in the 16 GB P100 VRAM. GPU parallelism applies within each allele's `predict_to_dataframe()` call (all 1.26M peptides at once), not across alleles.

**Infrastructure issues resolved during development (branch `feat/issue-105-gpu-mhcflurry`):**

1. **`hisat2.yaml` samtools/libdeflate conflict:** `regtools >=1.0.0` requires `libdeflate >=1.26`; no bioconda linux-64 samtools/htslib build is compiled against it. Removed samtools from conda env; pipeline uses system apt samtools 1.13 instead.
2. **`ProcessPoolExecutor` OOM on CPU:** 6 workers × full model (~8 GB RAM each) = ~48 GB on a 52 GB VM → BrokenProcessPool. Fixed by removing the process pool entirely — sequential execution on both CPU and GPU paths.
3. **NVIDIA driver 580 incompatible with P100:** Driver 580 requires GSP firmware; P100 (Pascal, SM 6.0) lacks it. Fixed by in-place downgrade to driver 570 (`apt-get install nvidia-driver-570-server && purge *580* && reboot`) — no VM deletion needed.
4. **PyTorch SM 6.0 mismatch:** PyTorch 2.5+ dropped SM 6.0 (P100/Pascal) support. Pinned `torch>=2.0,<2.5` in `python.yaml` (installs 2.4.1). Also rewrote `_has_gpu()` to use a PyTorch smoke-test kernel instead of TensorFlow — TF reported GPU available even when PyTorch kernels would fail on SM mismatch.
5. **Orphan sentinel file:** `.snakemake/conda/080a2daa..._.env_setup_done` existed without the actual env directory → Snakemake skipped rebuild → `ModuleNotFoundError: No module named 'mhcflurry'`. Fixed by deleting the orphan sentinel.

---

## 2026-04-22

### Issue #105 — GPU MHCflurry: CUDA architecture notes

**Host driver vs CUDA toolkit — how they relate:**

The **host NVIDIA driver** lives on the VM OS and is the only software that directly controls the GPU hardware. Its version determines the highest CUDA toolkit version it can support (e.g. driver 570 → CUDA 12.8).

The **CUDA toolkit** (cuDNN, cuBLAS, etc.) is what a specific application was compiled against. It lives wherever the application lives — inside a Docker container, inside a pip package — completely independent of the host.

**NVIDIA's backward-compatibility rule:** a host driver supports all CUDA toolkit versions ≤ its own. Driver 570 supports CUDA 12.8, 12.0, 11.8, 10.x, etc.

In our pipeline, two applications run on the same VM with different CUDA requirements:

```
Host (Deep Learning VM image)
  └─ NVIDIA driver 570 (CUDA 12.8 capable)
       ├─ conda env: python.yaml
       │    └─ tensorflow[and-cuda] pip package → bundles CUDA 12.x libs → MHCflurry
       └─ Docker container: tcrdock:latest → bundles CUDA 11.8 libs → TCRdock
```

Neither application depends on the host CUDA toolkit — each carries its own. The host driver just needs to be ≥ the highest toolkit version in use. 12.8 covers both.

**Why `ProcessPoolExecutor` crashes on GPU:**

MHCflurry's CPU path spawns multiple worker processes (one per allele), each importing TensorFlow and initialising its own CUDA context. Multiple CUDA contexts on the same GPU conflict → `BrokenProcessPool` crash. Fix: detect GPU at startup, load the predictor once in the main process, and iterate over alleles sequentially. The GPU handles massive internal parallelism within each `predict_to_dataframe()` call — sequential allele iteration adds negligible overhead.

---

### Issue #98 — Proteome k-mer filter (PR #106, merged)

Rewrote `blastp_filter.py` as `proteome_filter.py`. Instead of running blastp as a subprocess, the script parses the Swiss-Prot FASTA once to build a `set[str]` of all k-mers from every canonical human protein, then checks each query peptide via O(1) set lookup. Drops the BLAST conda env entirely — runtime goes from ~2 hours (424K peptides × 4 threads) to seconds.

The `matched_accessions` column in `peptides_excluded.tsv` carries all source proteins for each excluded peptide (semicolon-separated), replacing the separate `blastp_hits.tsv` audit file. Rule, script, test, and config key all renamed `blastp_filter` → `proteome_filter`. 19 unit tests passing.

**Local test results (chr22, patient_001_test, 8/9/10-mers):**

| | Count | % |
|---|---|---|
| Total peptides into filter | 4,163 | 100% |
| Novel (passed) | 4,119 | 98.9% |
| Excluded (self-peptides) | 44 | 1.1% |

---

### PR #101 — Move research artefacts to `research/` (Issue #100)

Separated operational documentation from research outputs. `docs/lab_notebook.md` → `research/lab_notebook.md`; `docs/manuscript/` → `research/manuscript/`. `docs/` now contains only software/infrastructure documentation (installation, configuration, cloud guide, etc.).

---

### PR #102 — Remove automatic Claude PR review CI workflow

Removed `.github/workflows/claude-code-review.yml`, which triggered a full Claude review on every push. In a solo project this was noisy and consumed usage unnecessarily. On-demand review is still available via `@claude` mentions through `claude.yml`.

---

### Issue #82 — Multi-length peptide extraction (8-mer, 9-mer, 10-mer)

**Motivation:** MHC-I presents 8–10-mer peptides. Extracting only 9-mers misses a significant fraction of true binders.

**Design:**
- Replaced hardcoded 9-mer extraction in `translate_peptides.py` with a loop over a configurable `translation.peptide_lengths: [8, 9, 10]` list.
- Contig flank is now auto-derived: `flank_nt = 3 * (max_length - 1)`. For max length 10 this gives symmetric 27/27 nt flanks (54 nt contigs), up from 26/24 (50 nt). The old 50 nt contigs were 3 nt too short to cover all 10-mer junction-spanning windows on the downstream side.
- MHCflurry natively accepts 8–15-mers; blastp operates on whatever peptides it receives. No changes needed in either step.
- 16 unit tests updated/added covering 8/9/10-mer start ranges and known-start sequences.

**Local test results (chr22, patient_001_test):**

| | 9-mer only | 8/9/10-mer |
|---|---|---|
| Unique peptides into MHCflurry | 1,339 | 4,045 (~3×) |
| Total predictions (6 alleles) | 8,172 | 24,714 (~3×) |
| Strong binders (≤50 nM) | 39 | 37 |
| Weak binders (≤500 nM) | 290 | 502 (+73%) |

**Follow-on issues filed:**
- **Issue #97** — Use GENCODE CDS reading frame to restrict translation to the biologically relevant frame. Currently all 3 reading frames are used per junction; the correct frame can be determined from the GENCODE GTF for junctions where (a) the upstream donor maps unambiguously to a single CDS exon end, and (b) no other tumor-exclusive junction is detected upstream in the same gene. Both conditions must hold; overlap between genes and upstream frame-shifting junctions are edge cases requiring fallback to all 3 frames.
- **Issue #98** — Replace blastp exact-match filter with a Python set-based k-mer lookup. On the production patient (424K unique peptides), blastp ran for >2 hours. A set of all k-mers from the Swiss-Prot FASTA built once at startup reduces this to seconds via O(1) lookup per peptide. Also supports near-self filtering in future (generate all 1-mismatch variants of query, check against exact set).

---

### Issue #89 — blastp filter against human reference proteome (PR #96, merged)

**Overview:** Added `blastp_filter_peptides` rule between `translate_peptides` and `run_mhcflurry`. Excludes peptides with a 100% identity / 100% query coverage match in UniProt Swiss-Prot (self-peptides the immune system is tolerized to).

**Cloud run debugging:**
- *Run 1 — blastp silently skipped:* `mhc_affinity.tsv` existed from a pre-blastp run. With `--rerun-triggers mtime`, Snakemake never propagated upward to check for the missing `peptides_novel.tsv`. Fixed by deleting downstream outputs before re-running.
- *Run 2 — `Error: Unknown argument: "perc_identity"`:* `-perc_identity` is a `blastn`-only flag; blastp 2.16.0 rejects it. Fixed: removed flag, now uses `-outfmt "6 qseqid sseqid pident"` + `-qcov_hsp_perc 100` (full coverage at search time) + Python `int(float(pident)) == 100` filtering. Added 15 unit tests.

**Local test (9-mers):** 1,371 peptides → 9 self-peptides excluded (0.7%) → 1,362 novel passed to MHCflurry.

---

### PR #93 — Remove `gcs:` block and BAM upload fix

Removed the `gcs:` block from `config/config.yaml` and `config/test_config.yaml` (BAM upload had been moved to `scripts/run_cloud_cpu.sh` and no longer read config). Also removed `temp()` + `upload_bam` from `alignment.smk` which was preventing BAM/BAI/BED files from reaching GCS.

---

### PR #94 — Snakefile refactor (Issue #92)

Reduced the Snakefile to a minimal entry point: patient ID derivation and shared config moved into `common.smk`. Added a zero-rows guard in `_read_samples_tsv` that raises a clear `ValueError` on empty or missing TSV rather than a bare `IndexError`.

---

### PR #95 — Claude GitHub Actions workflows

Added Claude Code Review and Claude PR Assistant GitHub Actions workflows (`.github/workflows/`).

---

## 2026-04-21

### Patient_001 (gastric cancer) — GPU run completed

Patient_001 TCRdock run completed successfully overnight. Top candidate EVAEYNASF / HLA-A\*26:01 (IC50 = 16.5 nM) was run through TCRdock on the P100 GPU VM. Outputs archived to `gs://splice-neoepitope-project/results/patient_001/`: `docking_scores.tsv`, `top_candidate.pdb`, `report.html`. Both CPU and GPU VMs stopped cleanly (TERMINATED).

This is the second confirmed end-to-end run for patient_001 (first was on the old `tcrdock-handoff` bucket; this run used the new `splice-neoepitope-project` bucket with the updated `run_cloud_gpu.sh`).

### Infrastructure fix — `.snakemake/metadata` sync via GCS (commit `3ddfe96`)

**Problem:** The GPU VM was re-running all CPU pipeline steps (alignment, filtering, MHCflurry) on every invocation, even when results were already present. Root cause: without `.snakemake/metadata/` on the GPU VM, Snakemake has no baseline for `--rerun-triggers code params` and triggers re-runs of every rule.

**Fix:** Extended `run_cloud_gpu.sh` to sync `.snakemake/metadata/` from the CPU VM to GCS (Phase 2 upload) and then from GCS to the GPU VM (Phase 3 download), alongside `results/` and `logs/`. Only `.snakemake/metadata/` is synced — the `conda/` subdir is gigabytes and is not transferred.

**Related investigation — "Incomplete files" warning:**
The GPU VM also showed an `IncompleteFilesException` on startup despite all result files being present. Root cause: stale `.snakemake/incomplete/` markers from a previous Spot VM preemption. These are separate from `.snakemake/metadata/` — incomplete markers are only cleared by `--rerun-incomplete`, `--cleanup-metadata`, or manual removal. Resolved by manually running `rm -rf .snakemake/incomplete/` on the GPU VM as a one-time fix. Decision not to add this to the script: (a) the VM is now STANDARD (not Spot), so future preemptions are unlikely; (b) clearing incomplete markers blindly would mask legitimately partial TCRdock outputs from future interrupted runs, since `gcloud storage cp` does not delete local files absent from GCS.

### Issue #73 — patient_002 manuscript results (PR #74, merged)

Added patient_002 osteosarcoma results to `research/manuscript/`:
- **RESULTS.md:** full patient_002 section — dataset, HLA typing + serology validation, junction funnel (347,046 raw → 55,912 tumor-exclusive), peptide translation (781,424 9-mers), MHC predictions (12,430 strong binders), top candidate TTDPVQALY / HLA-A\*01:01 (IC50 = 23.9 nM), TCRdock caveat (fallback allele bug).
- **CONCLUSIONS.md:** added patient_002 Key Findings (points 5–8); updated Limitations HLA section (confirmed match, no longer future tense); updated Future Directions (T0 complete, focus on T1/T2).
- **DISCUSSIONS.md:** updated "Impact of missing matched normal: patient_002" section with real WES proxy counts (106,474 apparent junctions, only 3 overlap tumor).

PR #74 (`docs/issue-73-manuscript-patient002` → `main`) opened and merged.

### Issue #59 — Known HLA serology input via samples TSV (PR #69, merged)

**Motivation:** Clinical HLA serotyping (Red Cross / WGS) is more reliable than OptiType from RNA-seq for germline alleles. Patient_002 has confirmed Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01) which should take priority over the blood WES OptiType call when no tumor call is available.

**Design decision:** Following nf-core convention (samplesheet as single source of truth), serology alleles are stored as six inline columns in each patient's samples TSV (`serology_A1/A2`, `serology_B1/B2`, `serology_C1/C2`). The tumor row is left empty; only the normal/blood row carries the germline values.

**Priority order** implemented in `aggregate_hla_alleles.py`:
1. **Tumor OptiType** — preferred because HLA LOH can occur in the tumor, altering what is actually presented
2. **Serology** — germline ground truth from clinical typing
3. **Normal/blood OptiType**
4. **Config fallback**

**Null allele handling:** A\*01:11N is a null allele (not expressed at the cell surface). It is excluded from the prediction allele list but retained in `hla_qc.tsv` for transparency. The serology validation result is `match (null allele excluded)` when OptiType finds only the expressed allele.

**QC output extended:** `hla_qc.tsv` now includes `serology_allele1`, `serology_allele2`, and `serology_validation` columns. The HTML report's HLA section shows a "Known alleles / validation" column when serology data is present.

---

### Issue #78 / PR #77 — `gcloud storage rsync` + `mtime` rerun trigger (merged)

**Problem:** `gcloud storage cp` fresh-timestamps all downloaded files. With `--rerun-triggers code params` on the GPU VM, this is fine as long as only code/param changes matter. But if `mhc_affinity.tsv` changes content (e.g. HLA alleles change), TCRdock would not re-run under `code params` alone because TCRdock's script and params are unchanged.

**Fix:** Switched Phase 3 GPU download from `gcloud storage cp` to `gcloud storage rsync --recursive` (CRC32C checksum comparison, only downloads changed files, preserves mtime of unchanged files). Also added `mtime` to GPU Snakemake's `--rerun-triggers` so file content changes cascade correctly. Metadata download remains `cp` (always overwrite with CPU VM's authoritative `.snakemake/metadata/`).

---

### Issue #57 Phase 1 — `report.tsv` structured summary artifact (PR #80 merged, PR #81 open)

Added `_build_report_tsv()` to `generate_report.py`. Emitted alongside `report.html` whenever `output_tsv` is provided. Schema: `patient_id | stage | metric | value | notes`. Stages:

| Stage | Metrics written |
|---|---|
| `junction_filtering` | `unannotated`, `tumor_exclusive`, `normal_shared` per sample |
| `mhc_prediction` | `total_predictions`, `strong`, `weak`, `non` binder counts |
| `top_candidate` | `peptide`, `allele`, `ic50_nM`, `binder_class` |
| `hla_typing` | `HLA-A/B/C` alleles with source and read count |
| `tcrdock` | `pdb_available` true/false |

`report_tsv` added as a declared Snakemake output in both `analysis.smk` (non-GPU path) and `structure.smk` (GPU path with TCRdock). The `structure.smk` output declaration was missing from PR #80, causing an `AttributeError` on the GPU VM. Fixed in PR #81.

Phase 2 (refactor `report.html` to read from `report.tsv` instead of recomputing) deferred to Issue #79.

---

### Patient_001 (gastric cancer) — re-run with tumor-first HLA (successful)

Re-ran patient_001 with all merged changes (tumor-first HLA priority, serology columns, rsync, report.tsv). TCRdock re-ran despite existing results because `aggregate_hla_alleles.py` code changed → HLA re-typed → `mhc_affinity.tsv` regenerated with fresh mtime → `mtime` trigger correctly fired TCRdock (~4 min GPU time).

**HLA-B alleles updated** (tumor-first policy now active): B\*15:01/B\*18:02 → **B\*15:63/B\*18:01**. The A alleles were unchanged so the top candidate is the same.

**`report.tsv` read directly from GCS** (`gcloud storage cat gs://splice-neoepitope-project/results/patient_001/reports/report.tsv`):

| Stage | Metric | Value |
|---|---|---|
| junction_filtering | tumor_exclusive | 27,348 |
| junction_filtering | normal_shared | 2,681 |
| mhc_prediction | total_predictions | 2,598,882 |
| mhc_prediction | strong | 13,139 |
| top_candidate | peptide | EVAEYNASF |
| top_candidate | allele | HLA-A\*26:01 |
| top_candidate | ic50_nM | 16.51 |
| hla_typing | HLA-B | HLA-B\*18:01 / HLA-B\*15:63 |

patient_002 re-run planned next (with all latest changes including PR #81).

---

### Documentation update (PR #87)

Updated three docs to reflect recent pipeline changes:
- `docs/data_preparation.md` — serology columns added to sample manifest format and column table; HLA typing roles table updated with tumor-first priority order
- `docs/configuration.md` — HLA section: added allele priority order note and serology column cross-reference
- `docs/google_cloud_guide.md` — fixed stale bucket name (`<PROJECT_ID>-tcrdock-handoff` → `splice-neoepitope-project`); added `{patient_id}` to result paths; added `report.tsv` to retrieval example; noted `gcloud storage rsync` in "How it works"

---

## 2026-04-20

### Patient_002 (osteosarcoma BG003082) — first full production run

**Patient:** BG003082 T0 tumor (Boston Gene, Nov 2022, paired-end RNA-seq ~10 GB) + BG003082 N0 WES normal (blood-derived, used for HLA typing only).

**HLA typing:** A\*01:01/A\*01:01, B\*08:01/B\*27:05, C\*07:01/C\*01:02 — confirmed match to Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01). First patient with ground-truth HLA validation.

**Results:** Run completed end-to-end: alignment → HLA typing → MHCflurry → TCRdock → HTML report with Mol\* 3D viewer. Final outputs archived to `gs://splice-neoepitope-project/results/patient_002/`.

**Infrastructure bugs discovered and fixed (PR #69, branch `feat/issue-65-patient002-cloud-run`):**

- **mtime re-run cascade:** Re-downloaded temp FASTQs had newer mtime than existing `junctions.tsv`, causing unnecessary re-alignment and OptiType re-runs. Fixed with `ancient()` on FASTQ inputs in both `hisat2_align` and `run_optitype`.
- **OptiType OOM:** razers3 peaks at ~36 GB on full RNA-seq FASTQs. CPU VM upgraded from n1-standard-8 → n1-highmem-8 (52 GB) → n2-highmem-8 (64 GB; n1 unavailable in zone). OptiType threads capped at 5 to force sequential sample execution and prevent concurrent OOM.
- **samtools sort OOM:** `-m 3G` caused OOM on the WES normal sample (8 threads × 3 GB). Reduced to `-m 1G`.
- **Boot disk:** 50 GB → 100 GB pd-ssd to handle reference index + paired-end FASTQ staging.
- **`--rerun-incomplete`:** Added to orchestrator snakemake invocation so killed runs resume cleanly instead of raising `IncompleteFilesException`.
- **MHCflurry re-run on GPU VM:** `resources/mhcflurry_models.done` sentinel absent on GPU VM → cascade triggered `run_mhcflurry` re-run. Fixed by running `snakemake resources/mhcflurry_models.done --use-conda` before the main TCRdock run (uses the correct `python.yaml` env; direct `mhcflurry-downloads fetch` in the `snakemake` bootstrap env would fail on a fresh VM).
- **GPU VM GCS upload permission:** Existing TERMINATED GPU VM lacked `--scopes=cloud-platform`. `gcloud storage cp` failed with permission denied; TCRdock results were not uploaded. Fixed by adding `gcloud compute instances set-service-account ... --scopes=cloud-platform` before `instances start` in orchestrator.
- **tmux not installed on GPU VM:** Deep Learning VM image does not include tmux by default. Added idempotent `apt-get install -y -q tmux` to GPU provisioning block.
- **VM auto-stop:** CPU and GPU VMs now unconditionally stop on pipeline exit (success or failure).

**TCRdock result:** pLDDT 92.25 for top candidate `FMSGFLYFV` on `HLA-A*02:01` (fallback allele — note: for patient_002 the actual alleles are A\*01:01, not A\*02:01; fallback was used due to a sentinel issue, now fixed for future runs).

---

### Patient_001 (gastric cancer) — cloud run started

Updated `config/samples/patient_001.tsv` to use ENA HTTPS URLs instead of local `data/` paths, enabling cloud runs without pre-staging FASTQs. Run started; `n2-highmem-8` required after `n1-highmem-8` hit `ZONE_RESOURCE_POOL_EXHAUSTED` in `europe-west1-b`.

---

### Documentation update

README slimmed from ~600 to 337 lines. Detailed content moved to new dedicated docs:
- `docs/installation.md` — full conda/Snakemake setup
- `docs/data_preparation.md` — aligner selection, FASTQ sources, manifest format
- `docs/configuration.md` — full `config.yaml` parameter reference
- `docs/google_cloud_guide.md` — added manual TCRdock run section

---

## 2026-04-17

### Patient_001 (gastric cancer) — first full production run with HLA typing + TCRdock

**Patient:** SRR9143065 (Solid Tissue Normal) / SRR9143066 (Primary Tumor) — gastric cancer surgical section, single-end Illumina HiSeq 3000.

**Junction filtering results:**
- Unannotated junctions: 30,029
- Normal-shared (filtered out): 2,682 (8.9%)
- Tumor-exclusive candidates: 27,347

**HLA typing (OptiType):** HLA-A\*26:01, HLA-A\*31:01, HLA-B\*15:05, HLA-B\*18:20, HLA-C\*03:21, HLA-C\*07:01 — no ground-truth alleles available for this patient, so typing cannot be validated.

**MHCflurry predictions:** 8,226 peptide × allele pairs (1,371 9-mers across 6 alleles), 54 strong binders (IC50 ≤ 50 nM), 317 weak binders (IC50 ≤ 500 nM).

**Top neoepitope candidate:** EVAEYNASF / HLA-A\*26:01, IC50 16.5 nM (strong binder). TCRdock structure predicted on P100 GPU VM; Mol\* viewer renders ternary complex with correct chain labels (A=MHC, B=peptide, C=TCR-α, D=TCR-β). Results archived to `gs://splice-neoepitope-project-tcrdock-handoff/results/`.

**Infrastructure issues resolved during this run:**
- P100 GPU VM required proprietary nvidia kernel modules (`linux-modules-nvidia-570-server-6.8.0-1053-gcp`); open-source modules (`linux-modules-nvidia-570-server-open-*`) do not support Pascal (P100) GPUs.
- GCS download of pipeline outputs onto the GPU VM fresh-stamped all files, causing Snakemake's `--rerun-triggers mtime` to re-run MHCflurry unnecessarily. Fixed by switching the GPU phase to `--rerun-triggers code params`.
- Mol\* viewer silently broken by unpinned CDN URL (`unpkg.com/molstar` → v5.8.0, breaking API). Pinned to `molstar@4.9.0`.

---

### Issue #50 — parallel MHCflurry allele predictions (ProcessPoolExecutor)

**Problem:** Serial per-allele loop was the bottleneck for patients with 6+ HLA alleles.

**First attempt (ThreadPoolExecutor) — failed:** `predict_to_dataframe()` is not thread-safe due to shared TensorFlow state. Testing revealed identical IC50 values (267.12 nM) across different alleles for the same peptide — confirmed state corruption. Threads release the GIL during I/O but TensorFlow's internal state is mutated during inference.

**Fix:** Switched to `ProcessPoolExecutor` with `initializer=_worker_init`. Each worker process loads its own predictor copy, sets `TF_NUM_INTRAOP_THREADS=1` / `OMP_NUM_THREADS=1` before TensorFlow imports (prevents CPU oversubscription when running multiple workers), and returns a lean per-allele DataFrame (no `peptides_df` pickling across process boundaries).

**Test result:** 8,226 rows, 54 strong, 317 weak — max IC50 diff vs. serial baseline = 0.00e+00. Local run: 6 alleles, 4 workers, ~12 s total.

**Also renamed:** `predict.smk` → `mhcflurry.smk`, `predictions.tsv` → `mhc_affinity.tsv` (tool-agnostic naming).

---

### Patient_002 planning — osteosarcoma IPISRC044

**Dataset:** Publicly available osteosarcoma dataset (https://osteosarc.com/data/). Patient IPISRC044, multi-institutional (UCLA / UCSF / Boston Gene / Tempus). GCS bucket: `gs://osteosarc-genomics`.

**Plan:** Start with T0 tumor (Boston Gene, Nov 2022) as the baseline timepoint.
- FASTQs: `rna-seq/fastq/bostongene_2022/202211_bostongene_tumor_rna_BG003082_R1.fastq.gz` + `_R2.fastq.gz` (paired-end, ~10 GB)
- Run OptiType for HLA typing — ground-truth Class I alleles are known from Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01), giving us a validation opportunity we didn't have for patient_001.

**No matched RNA-seq normal available.** Blood WGS DNA cannot be used as a substitute — `regtools junctions extract` requires spliced RNA-seq reads (`N` CIGAR operations); DNA-seq reads map continuously without splicing and produce no junctions. Pipeline runs in warning mode, labelling all unannotated junctions `tumor_exclusive`. Based on patient_001 statistics, approximately 9% of unannotated junctions may be normal-shared and would be misclassified. The downstream MHCflurry + TCRdock funnel is expected to absorb most of this noise.

---

## 2026-04-16

**Goal:** Implement HLA typing with OptiType and stabilise the cloud pipeline for a production run.

**Done:**
- Implemented HLA typing step (issue #42): OptiType runs on each sample's FASTQ, typing results aggregated per patient into `alleles.tsv`, passed to MHCflurry in place of fallback alleles when `hla.enabled: true`.
- Added configurable CBC ILP solver for OptiType (issue #49) — reduces OptiType runtime significantly on VMs with multiple cores.
- Fixed `PYTHONUNBUFFERED=1` for OptiType log flushing (issue #52): without it, log output was buffered and appeared to stall.
- Fixed HLA concordance display in report (issue #53): loci with no normal/tumor discrepancy now show ✓ concordant instead of blank.
- Upgraded CPU VM to n1-standard-8 and made thread counts config-driven (issue #40).
- Fixed Snakemake unlock on interrupted restarts (issue #47).
- Suppressed interactive "Next steps" prompts when `run_cloud_gpu.sh` calls sub-scripts (issue #48).
- Dropped all GDC/TCGA code (issue #44) — pipeline is now fully open-access, no registration-gated data sources.
- Renamed junction origin labels to `tumor_exclusive` / `normal_shared` (issue #38).

---

## 2026-04-14

**Goal:** Refactor pipeline to be patient-centric and unify alignment rules.

**Done:**
- Replaced `{cancer_type}` wildcard with `{patient_id}` throughout (issue #26). Config key `cancer_types: [local]` → `patient_id` string. Routing in `filter.smk` now switches on `config["data_source"]` rather than hardcoded `"local"` path segments.
- Unified STAR and HISAT2 alignment into a single `alignment.smk` module (issue #35). Both aligners produce junction TSVs in the same format; downstream rules are aligner-agnostic.
- `PATIENT_IDS` now derived from `samples_tsv` at DAG construction time; sample IDs renamed to SRR accessions throughout.

---

## 2026-04-13

**Goal:** Merge TCRdock structural validation and get a clean production baseline.

**Done:**
- Merged PR #27 (issue #25): TCRdock step, Mol\* 3D viewer, `run_cloud_gpu.sh` three-phase lifecycle (CPU → GCS → GPU Spot VM). Full end-to-end test passed.

---

## 2026-04-10

**Goal:** Get TCRdock structural validation running end-to-end via Docker on the GCP GPU Spot VM.

**Done:**
- Fixed Docker build failure: the official TCRdock Dockerfile fails at `conda install openmm=7.7.0` (no longer available). Created `docker/Dockerfile.pipeline` which uses pip only and omits openmm/pdbfixer — these are only needed for optional Amber relaxation, not for `run_prediction.py`.
- Successfully built `tcrdock:latest` on GPU VM and ran TCRdock end-to-end. Pipeline completed 3/3 steps, report generated.
- Identified visual issue in report: AlphaFold outputs all residues as a single chain A, so Mol* rendered the ternary complex without chain-colour distinction. Added `relabel_pdb_chains()` to `run_tcrdock.py`, which reassigns chain IDs (A=MHC, B=peptide, C=TCR-alpha, D=TCR-beta) using the chain lengths from TCRdock's `alphafold_setup/targets.tsv`.
- Updated `setup_tcrdock_vm.sh` to build from `docker/Dockerfile.pipeline` instead of cloning TCRdock separately for the (broken) official Dockerfile.

### 2026-04-10 15:30 — Mol* COMPND fix

Fixed Mol* sequence panel showing generic "Polymer 1/2/3/4" instead of meaningful chain names. Root cause: PDB COMPND records had off-by-one column positions (continuation number in cols 9–11 instead of PDB-standard 8–10) and the first line incorrectly included a continuation number. Also padded lines to 80 chars and transliterated Unicode α/β to ASCII for PDB compliance. Confirmed fix locally — Mol* now shows "MHC heavy chain", "Peptide", "TCR alpha", "TCR beta".

### 2026-04-10 17:00 — Pre-PR refactoring

Code cleanup before creating PR for #25:
- `generate_report.py`: fixed `html` module shadowing (local variable named `html` overrode the stdlib import → renamed to `report_html`, import aliased to `html_mod`). Moved `import json` from function body to top-level. Extracted COMPND record building into `_build_compnd_records()` helper.
- `run_tcrdock.py`: replaced `assert` with `raise ValueError` for input validation.
- All 57 tests passing.

### 2026-04-10 18:00 — Documentation updates

Updated all project documentation for the branch:
- README.md: pipeline diagram (now 7 steps), TOC, config table with `tcrdock.*` parameters, output tree with `tcrdock/` directory, project structure with `docker/` and new files, citations for TCRdock and Mol*.
- `docs/google_cloud_guide.md`: new "Automated GPU Pipeline (TCRdock)" section with quick start, retrieval, how-it-works, detached mode, cost estimate. Reordered TOC so Troubleshooting is last.
- CLAUDE.md: added "TCRdock via Docker" and "PDB chain relabelling" decision notes.

### 2026-04-10 19:00 — Final cloud test

Kicked off `run_cloud_gpu.sh` for end-to-end validation on GCP. Pending result before creating PR for #25.

---

## 2026-04-09

**Goal:** Implement TCRdock structural validation (issue #25) and automate the full CPU→GPU cloud pipeline.

**Done:**
- Merged issues #20 (filter junction-spanning 9-mers) and #22 (move junction-spanning filter to translation step; output peptides as TSV). Closed #21.
- Implemented TCRdock step (`workflow/rules/tcrdock.smk`, `workflow/scripts/run_tcrdock.py`) and Mol* 3D viewer in the HTML report. Verified TCRdock input column format and two-step workflow (`setup_for_alphafold.py` → `run_prediction.py`) against the real API before writing code.
- Wrote `scripts/run_cloud_gpu.sh`: three-phase lifecycle (CPU VM steps 1–5 → GCS handoff → GPU Spot VM TCRdock). Used GCS bucket (`tcrdock-handoff`) for VM-to-VM transfer.
- Hit a series of dependency issues with the conda-based TCRdock env: wrong script name (`predict.py` vs `run_prediction.py`), missing tensorflow, Python 3.8 incompatibility with dm-haiku, cuDNN 8.9 vs 9.10 mismatch on the Deep Learning VM image.
- Decided to switch to Docker to sidestep the dependency issues entirely. Rewrote `run_tcrdock.py` to call TCRdock via `docker run --gpus all` instead of a conda env; updated `tcrdock_gpu.yaml` and `tcrdock.smk` accordingly.
- Official TCRdock Dockerfile fails at `conda install openmm=7.7.0`. Confirmed openmm/pdbfixer are only used in `alphafold/relax/` (Amber relaxation), not in `run_prediction.py`. Session ended with this as the open problem.

**Key decisions:**
- Use Docker for TCRdock rather than a conda env — eliminates host-side cuDNN/JAX version management.
- Use `--new_docking` flag (1 AlphaFold run per target instead of 3) to reduce GPU time.
- GCS bucket for VM-to-VM result transfer rather than direct SCP, so the two VMs don't need to be up simultaneously.

---

## 2026-04-07

**Goal:** Run full pipeline on real cancer data; add local test dataset for macOS development.

**Done:**
- Full cloud pipeline run succeeded on SRR37781424 (Luminal A breast cancer, GEO GSE119889). All steps 1–5 completed on `splice-prod-test` VM.
- Added `scripts/prepare_test_data.sh` and chr22 test config for local macOS runs (M1, 8 GB RAM). Downloads chr22 reference + 500K-read subsets of a matched gastric cancer pair (SRR9143066 tumor / SRR9143065 normal) via ENA HTTPS to avoid sra-tools arm64 issues.
- Merged PR #1 (Copilot-assisted modernisation baseline).

---

## 2026-04-05

**Goal:** Fix pipeline bugs found during first cloud run.

**Done:**
- Fixed `auto_stop.sh`: `gcloud compute instances stop` fails with `ACCESS_TOKEN_SCOPE_INSUFFICIENT` from inside the VM. Switched to `sudo shutdown -h now`.
- Fixed conda env dependency conflicts (samtools/libdeflate, regtools version pinning).
- Fixed mhcflurry 2.2.0 API change: `predict()` now returns a numpy array, not a DataFrame. Switched to `predict_to_dataframe()`.
- Fixed `statistical_analysis.py`: extract `sample_type` from `source_header` correctly.
- Auto-download MHCflurry models as a pipeline step rather than requiring manual setup.

---

## 2026-04-04

**Goal:** Document GCP deployment so the pipeline can be handed off or reproduced.

**Done:**
- Wrote `docs/google_cloud_guide.md`: full step-by-step from project creation to pipeline run, including VM setup, conda, sra-tools version pinning (3.1.1 — newer versions segfault), and regtools argument order gotcha.

---

## 2026-04-03

**Goal:** Make alignment work on the GCP VM (8 GB RAM) without STAR.

**Done:**
- Added HISAT2 as a low-memory alternative to STAR. STAR requires ~30 GB RAM for hg38; HISAT2 indexes fit in ~8 GB.
- Added DAG diagram (PDF) for workflow visualisation.

---

## 2026-04-02

**Goal:** Replace NetMHCPan with an open-source MHC binding predictor.

**Done:**
- Replaced NetMHCPan (registration-gated, no programmatic access) with MHCflurry 2.x. MHCflurry is fully open-source, pip-installable, and produces IC50 predictions compatible with the 500 nM strong-binder threshold used in the original 2015 paper.

---

## 2026-03-25

**Goal:** Implement the modernised pipeline from scratch.

**Done:**
- Initial Snakemake pipeline implementing all steps: STAR/HISAT2 alignment → regtools junction extraction → GENCODE annotation filtering → normal-sample filtering → peptide translation → MHCflurry binding prediction → HTML report.
- Junction origin classification hierarchy: annotated → discard; unannotated + in normal → patient-specific (discard); unannotated + absent in normal → tumor-specific (predict). This replaces the original Fisher's exact test with an upstream biological filter.

---
