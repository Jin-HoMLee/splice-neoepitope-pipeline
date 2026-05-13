# Lab Notebook — Scientist

Per-role lab notebook for Scientist sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

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
