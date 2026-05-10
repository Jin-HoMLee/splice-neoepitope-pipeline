# Lab Notebook — Scientist

Per-role lab notebook for Scientist sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-10

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
