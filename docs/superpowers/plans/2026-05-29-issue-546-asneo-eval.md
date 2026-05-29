# ASNEO Eval Implementation Plan (Issue #546)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a desk-eval verdict on ASNEO (Zhang et al., *Aging* 2020) as a splice-neoepitope candidate generator peer to our pipeline — a verified Zotero note, a 5-mode pipeline-fit decision on #546, and a tool-primer Quarto deck.

**Architecture:** A multi-agent research Workflow (gather → adversarial-verify → return structured findings) feeds a single-author synthesis pass that writes the Zotero note, the deck, and the decision comment. No code execution / no ASNEO run (desk eval).

**Tech Stack:** Claude Workflow tool (research engine); Quarto reveal.js (deck); `research/scripts/zotero_add.py` (Zotero); `gh` CLI + GraphQL (board/PR); `quarto` + headless Chrome (render + visual verify).

**Domain note — there is no unit-test loop here.** This plan produces research + documentation artifacts, not testable code. The TDD "red/green" gates are replaced by **verification gates**: (1) the Workflow's adversarial-verify phase rejects unsourced claims, (2) a headless per-slide overflow check, (3) the AC closure-ritual audit. Each is called out explicitly in its task.

**Spec:** [docs/superpowers/specs/2026-05-29-issue-546-asneo-eval-design.md](../specs/2026-05-29-issue-546-asneo-eval-design.md)

---

## File Structure

| Path | Responsibility | Action |
|---|---|---|
| `research/evals/issue_546_asneo/slides.qmd` | Tool-primer deck source (title + 9 content slides) | Create |
| `research/evals/issue_546_asneo/refs.bib` | Deck bibliography (ASNEO + any cited comparators) | Create |
| `research/evals/issue_546_asneo/slides.html` + `slides_files/` | Rendered deck (kept in tree, gitignored from *commit*? no — committed per HERMES precedent) | Create (rendered) |
| `research/lab_notebook/scientist.md` | Dated lab-notebook entry (after review) | Modify (append) |
| Zotero collection `Z38GTJNW` (external) | 3-section HTML note on DOI `10.18632/aging.103581` | Create (via API) |
| Issue #546 (external) | 5-mode decision comment + AC ticks | Modify |
| Issue #258 (external) | Verdict-comment update *iff* calculus changes | Conditionally modify |

Already done (pre-plan): #546 → In progress (`47fc9ee4`); branch `research/scientist/issue-546-asneo-eval` checked out; spec committed (`6c9ad06`).

---

## Task 1: Research Workflow → verified structured findings

**Files:** none committed (findings live in-context for the synthesis tasks).

This is the **gather + adversarial-verify** engine (spec approach C). It returns ONE validated object the rest of the plan consumes. The schema field names defined here are referenced verbatim in Task 3 — keep them consistent.

- [ ] **Step 1: Author the Workflow script**

Save inline via the Workflow tool (it persists to the session dir). Script:

```javascript
export const meta = {
  name: 'asneo-eval-research',
  description: 'Gather + adversarially verify ASNEO facts and our-pipeline backbone for the Issue #546 desk eval',
  phases: [
    { title: 'Gather' },
    { title: 'Verify' },
  ],
}

const FINDINGS = {
  type: 'object',
  required: ['asneo', 'ours', 'modes', 'open_questions', 'caveats'],
  properties: {
    asneo: {
      type: 'object',
      required: ['citation','license','repo_url','repo_last_activity','language',
                 'code_available','weights_available','junction_backbone',
                 'normal_filtering','peptide_assembly','frame_handling',
                 'mhc_predictor','validation_cohort','validation_n'],
      properties: {
        citation: {type:'string'}, license: {type:'string'},
        repo_url: {type:'string'}, repo_last_activity: {type:'string'},
        language: {type:'string'}, code_available: {type:'string'},
        weights_available: {type:'string'}, junction_backbone: {type:'string'},
        normal_filtering: {type:'string'}, peptide_assembly: {type:'string'},
        frame_handling: {type:'string'}, mhc_predictor: {type:'string'},
        validation_cohort: {type:'string'}, validation_n: {type:'string'},
      },
    },
    ours: {
      type: 'object',
      required: ['junction_backbone','normal_filtering','peptide_assembly','frame_handling','mhc_predictor'],
      properties: {
        junction_backbone: {type:'string'}, normal_filtering: {type:'string'},
        peptide_assembly: {type:'string'}, frame_handling: {type:'string'},
        mhc_predictor: {type:'string'},
      },
    },
    modes: {
      type: 'object',
      required: ['triage','replacement','cross_check','component_reuse','reject'],
      properties: {
        triage: {type:'object', required:['verdict','rationale'], properties:{verdict:{type:'string'},rationale:{type:'string'}}},
        replacement: {type:'object', required:['verdict','rationale'], properties:{verdict:{type:'string'},rationale:{type:'string'}}},
        cross_check: {type:'object', required:['verdict','rationale'], properties:{verdict:{type:'string'},rationale:{type:'string'}}},
        component_reuse: {type:'object', required:['verdict','rationale'], properties:{verdict:{type:'string'},rationale:{type:'string'}}},
        reject: {type:'object', required:['verdict','rationale'], properties:{verdict:{type:'string'},rationale:{type:'string'}}},
      },
    },
    verified_claims: {
      type: 'array',
      items: {type:'object', required:['claim','sources','status'],
        properties:{claim:{type:'string'}, sources:{type:'array',items:{type:'string'}}, status:{type:'string'}}},
    },
    open_questions: {type:'array', items:{type:'string'}},
    caveats: {type:'array', items:{type:'string'}},
  },
}

const CLAIM = {
  type:'object', required:['claim','status','sources','confidence'],
  properties:{ claim:{type:'string'}, status:{type:'string'},
    sources:{type:'array',items:{type:'string'}}, confidence:{type:'string'} },
}

phase('Gather')
const [asneoPaper, asneoRepo, ours] = await parallel([
  () => agent(`Read the ASNEO paper (Zhang et al., Aging 2020, DOI 10.18632/aging.103581) and any accessible mirrors. Extract, with source URLs: junction-calling backbone (STAR-only? extra filtering/motif selection?), normal-junction filtering (matched-normal? GTEx panel? both?), peptide assembly + flank handling + frameshift detection, MHC presentation predictor (netMHCpan / MHCflurry / in-house?), validation cohort + N. Mark any value you cannot source as "UNVERIFIED".`,
        {label:'asneo-paper', phase:'Gather', schema: FINDINGS.properties.asneo}),
  () => agent(`Survey the ASNEO GitHub/code repository. Report with URLs: repo URL, last commit/activity date, language (Python 2 vs 3), whether code is available, whether trained weights/models are available, and the license (verify against the actual LICENSE file, not just the README). Mark unsourced values "UNVERIFIED".`,
        {label:'asneo-repo', phase:'Gather', schema: {type:'object', required:['repo_url','repo_last_activity','language','code_available','weights_available','license'], properties:{repo_url:{type:'string'},repo_last_activity:{type:'string'},language:{type:'string'},code_available:{type:'string'},weights_available:{type:'string'},license:{type:'string'}}}}),
  () => agent(`Read these files in the current repo and summarize OUR splice-neoepitope backbone per axis: workflow/scripts/bed12_to_junctions.py, workflow/scripts/star_sj_to_junctions.py, workflow/scripts/filter_junctions.py, workflow/scripts/assemble_contigs.py, workflow/scripts/translate_peptides.py, workflow/scripts/run_mhcflurry.py. Axes: junction_backbone, normal_filtering, peptide_assembly, frame_handling, mhc_predictor. Quote file:line for each claim.`,
        {label:'our-backbone', phase:'Gather', agentType:'Explore', schema: FINDINGS.properties.ours}),
])

phase('Verify')
const loadBearing = [
  `ASNEO MHC predictor identity: claimed "${asneoPaper.mhc_predictor}"`,
  `ASNEO validation cohort: claimed "${asneoPaper.validation_cohort}" (N=${asneoPaper.validation_n})`,
  `ASNEO license: claimed "${asneoRepo.license}"`,
  `ASNEO normal-junction filtering: claimed "${asneoPaper.normal_filtering}"`,
]
const verdicts = await parallel(loadBearing.map(c => () =>
  agent(`Adversarially verify this claim against >=2 INDEPENDENT sources (paper text, repo, third-party). Default to status="unverified" if you cannot find >=2 concordant sources. Claim: ${c}`,
        {label:`verify`, phase:'Verify', schema: CLAIM})))

return {
  asneo: { citation: asneoPaper.citation || 'Zhang et al., Aging 2020, 10.18632/aging.103581', ...asneoPaper, ...asneoRepo },
  ours,
  modes: {
    triage: {verdict:'', rationale:''},
    replacement: {verdict:'', rationale:''},
    cross_check: {verdict:'', rationale:''},
    component_reuse: {verdict:'', rationale:''},
    reject: {verdict:'', rationale:''},
  },
  verified_claims: verdicts.filter(Boolean),
  open_questions: [],
  caveats: [],
}
```

- [ ] **Step 2: Run the Workflow**

Invoke the Workflow tool with the script above. It runs in the background and notifies on completion.

- [ ] **Step 3 (VERIFICATION GATE): Audit the returned findings**

Read the returned object. For EVERY field used in the deck (Task 3), confirm it is NOT the literal string `UNVERIFIED` and that any numeric/cohort claim appears in `verified_claims` with `status:"verified"` and ≥2 `sources`. The `modes` verdicts come back empty by design — the Scientist fills them in synthesis (Step 4), they are not delegated.

Expected: a populated `asneo` + `ours` object; `verified_claims` with the four load-bearing claims each `verified` or `unverified`.

- [ ] **Step 4: Synthesize the 5-mode verdicts + decision (Scientist, in main loop)**

Using ONLY verified findings, fill `modes.{triage,replacement,cross_check,component_reuse,reject}` (verdict ∈ {✅,⚠️,❌} + one-line rationale), pick the headline decision mode, and draft `open_questions` + `caveats`. Hedge or drop any `unverified` claim (memory: never quote unreachable sources). No commit — these feed Tasks 2/3/5.

---

## Task 2: Zotero entry (dedup-check → 3-section note)

**Files:** Zotero collection `Z38GTJNW` (external; no git change). `.env` at repo root holds creds.

- [ ] **Step 1 (GATE): Dedup pre-check the DOI**

`zotero_add.py` does not dedup. Query existing items first:

```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-scientist
python - <<'PY'
import os, requests
from dotenv import load_dotenv
load_dotenv()
uid=os.environ["ZOTERO_USER_ID"]; key=os.environ["ZOTERO_API_KEY"]
r=requests.get(f"https://api.zotero.org/users/{uid}/items/top",
    headers={"Zotero-API-Key":key}, params={"q":"10.18632/aging.103581","qmode":"everything","limit":50})
hits=[(i["key"], i["data"].get("DOI",""), i["data"].get("title","")) for i in r.json()
      if "aging.103581" in (i["data"].get("DOI","") or "")]
print("EXISTING:", hits or "none")
PY
```

Expected: `EXISTING: none` (then add) or a `(KEY, doi, title)` tuple (then `--update-note KEY`).

- [ ] **Step 2: Surface the entry to the user before adding**

Per `feedback_zotero_surface_extras.md`: show the user the title + the drafted 3-section note and confirm before writing to Zotero. Do not silently add.

- [ ] **Step 3: Add (or update) with the 3-section HTML note**

Note format = **Findings / Methods / vs. our pipeline**, bold-keyword leads, telegraphic bullets (`feedback_zotero_note_format.md`). Tags are **space-separated**:

```bash
python research/scripts/zotero_add.py 10.18632/aging.103581 \
  --tags splice-neoepitope alternative-splicing rna-seq tool-eval \
  --note '<p><b>Findings.</b></p><ul><li>...</li></ul><p><b>Methods.</b></p><ul><li>...</li></ul><p><b>vs. our pipeline.</b></p><ul><li>...</li></ul>'
# If Step 1 found an existing key K:  ... --update-note K --note '<...>'
```

- [ ] **Step 4 (VERIFICATION GATE): Confirm the entry + capture its Zotero key**

Re-query `items/top` for the DOI; record the returned item key for `refs.bib` (`note = {Zotero: KEY}`).

---

## Task 3: Scaffold the deck (`slides.qmd` + `refs.bib`)

**Files:**
- Create: `research/evals/issue_546_asneo/refs.bib`
- Create: `research/evals/issue_546_asneo/slides.qmd`

- [ ] **Step 1: Write `refs.bib`**

```bibtex
@article{zhang2020asneo,
  title = {{ASNEO}: Identification of personalized alternative splicing based neoantigens with {RNA}-seq},
  author = {Zhang, Zhanbo and others},
  journal = {Aging},
  year = {2020},
  doi = {10.18632/aging.103581},
  note = {Zotero: <KEY from Task 2 Step 4>},
}
```

Add a second `@article`/`@misc` entry only for any comparator actually cited on a slide (e.g. the #258 NeoGuider paper, Zotero `Z8FJSDVT`).

- [ ] **Step 2: Write `slides.qmd` — front matter (mirror HERMES verbatim except title/date)**

```yaml
---
title: "ASNEO as a splice-neoepitope candidate generator — peer or upgrade?"
subtitle: "Tool evaluation · Issue #546 · 2026-05-29"
author: "Jin-Ho Lee"
institute: "Splice Neoepitope Pipeline · JH M Lee Lab"
date: 2026-05-29
date-format: "YYYY-MM-DD"
format:
  revealjs:
    theme: [simple, ../../slides/custom.scss]
    footer: "Splice Neoepitope Pipeline · JH M Lee Lab"
    incremental: false
    slide-number: c/t
    progress: true
    chalkboard: false
    code-line-numbers: false
    fig-align: center
    fig-cap-location: bottom
    transition: fade
    background-transition: fade
    width: 1280
    height: 720
# PDF (beamer) render disabled — Quarto/pandoc longtable issue (per CLAUDE.md).
bibliography: refs.bib
csl: ../../slides/nature.csl
execute:
  echo: false
  warning: false
---
```

- [ ] **Step 3: Write the 9 content slides**

Each `## Heading` + `---` separator. Fill bracketed `<from findings.*>` slots from the Task 1 object (these are data-dependent values, NOT deferred decisions). Slide map:

1. **The question** — centered: *"Should ASNEO change how we generate splice neoepitopes?"* + a `callout-note` "Why now": surfaced from the #258 NeoGuider eval as the actual splice tool it delegates to.
2. **ASNEO at a glance** `{.smaller}` — `<findings.asneo.citation>`, `<license>`, `<repo_url>` + `<repo_last_activity>` + `<language>`, code/weights status. `callout-tip` one-liner of its CLI if known.
3. **How it works** — STAR → `SJ.out.tab` → `<junction_backbone>` filtering → `<peptide_assembly>` / `<frame_handling>` → `<mhc_predictor>` → ranking. Use an `.incremental` list.
4. **Head-to-head** `{.smaller}` — a 5-row table, columns *Axis | ASNEO | Ours*, rows = junction calling / normal filtering / peptide+frame / MHC predictor / validation. Values from `findings.asneo.*` and `findings.ours.*`.
5. **Where it would plug in** — a `{mermaid}` flowchart showing our DAG (junction extract → filter → assemble → translate → MHCflurry) with the candidate ASNEO insertion point(s) highlighted per the surviving mode(s).
6. **Reasons by mode** `{.smaller}` — table: *Mode | Verdict | Rationale* for Triage/Replacement/Cross-check/Component-reuse/Reject, from `findings.modes.*`.
7. **Decision** — centered big text = the headline mode + one-line framing; `. . .` then the sub-issue title (if integrate) or decline-note (if not); `callout-note` on whether the #258 calculus changed.
8. **Open questions & caveats** `{.smaller}` — `.incremental` list from `findings.open_questions` + `findings.caveats`.
9. **References** — `::: {#refs}\n:::` and cite entries inline elsewhere with `[@zhang2020asneo]`.

- [ ] **Step 4: Commit the source**

```bash
git add research/evals/issue_546_asneo/slides.qmd research/evals/issue_546_asneo/refs.bib
git commit -m "research(eval): ASNEO tool-primer deck source (Issue #546)"
```

(Commit only — surface SHA, do not push yet.)

---

## Task 4: Render + visual-verify the deck

**Files:** Create `research/evals/issue_546_asneo/slides.html` + `slides_files/`.

- [ ] **Step 1: Render**

```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-scientist/research/evals/issue_546_asneo
quarto render slides.qmd
```

Expected: `slides.html` written, exit 0, no unresolved-citation warnings.

- [ ] **Step 2 (VERIFICATION GATE): Headless per-slide overflow check**

reveal.js silently lets content overflow the 1280×720 canvas (`feedback_slide_visual_check.md`). Screenshot each slide index in a throwaway `.check_shots/` dir and page through:

```bash
for i in $(seq 0 9); do
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
    --headless --screenshot=".check_shots/slide_$i.png" --window-size=1280,720 \
    "file://$PWD/slides.html#/$i" 2>/dev/null
done
```

Read each PNG; confirm no clipped text/tables. `.check_shots/` is a throwaway debug dir — safe to delete after (the *exception* in `feedback_keep_render_artifacts.md`).

- [ ] **Step 3: Fix any overflow**

If a slide clips, add `{.smaller}` / split content / trim bullets in `slides.qmd`, re-render (Step 1), re-check (Step 2). Repeat until clean.

- [ ] **Step 4: Commit the rendered deck**

```bash
git add research/evals/issue_546_asneo/slides.html research/evals/issue_546_asneo/slides_files
git commit -m "research(eval): render ASNEO deck (Issue #546)"
```

Keep `slides.html` + `slides_files/` in the tree (`feedback_keep_render_artifacts.md`). Delete only `.check_shots/`.

---

## Task 5: Record the 5-mode decision + board-back the actionables

**Files:** Issue #546 (comment + AC ticks); Issue #258 (conditional); board #9 (any new follow-up Issues).

- [ ] **Step 1: Post the decision comment on #546**

Render-for-humans: a `Mode | Verdict | Rationale` table + the headline decision + the verified-claims provenance. End with **Created by: Scientist**. Use the link+prefix+keyword ref style for every Issue/PR mention.

```bash
gh issue comment 546 --body-file <(cat <<'MD'
## Verdict — <headline mode>
... 5-mode table + rationale + verified-claims list ...
**Created by:** Scientist
MD
)
```

- [ ] **Step 2: Conditionally update the #258 verdict comment**

If the ASNEO eval changes the NeoGuider component-reuse calculus (#546 AC), append a note to the #258 verdict; otherwise state explicitly in the #546 comment that #258 is unaffected (so it is not re-litigated).

- [ ] **Step 3 (GATE): Board-back every slide-surfaced actionable**

Per `feedback_slide_findings_on_board.md`: each item on the Decision/Open-questions/Caveats slides that is actionable (a promised follow-up, an open scientific question, a forward-looking caveat) must exist as an Issue or an Issue comment on board #9 BEFORE the deck PR merges. File them now; capture their #numbers for the PR Test plan.

---

## Task 6: Open the PR + request review

**Files:** none (git + GitHub).

- [ ] **Step 1: Push the branch (after surfacing)**

Surface the full commit list + diffstat, get the user's OK, THEN:

```bash
git push -u origin research/scientist/issue-546-asneo-eval
```

- [ ] **Step 2: Create the PR**

Body includes a **Test plan** with checkboxes: deck renders clean / per-slide overflow verified / Zotero note added / 5-mode decision posted on #546 / slide actionables board-backed (list the #s from Task 5 Step 3). Reference `closes #546` is wrong for a desk eval that stays open until merge — use `advances #546` or link without auto-close, matching the #258 pattern.

```bash
gh pr create --title "research(eval): ASNEO splice-neoepitope candidate generator — desk eval (#546)" --body-file <path>
```

- [ ] **Step 3: Flip PR Status → Ready for review**

A PostToolUse hook (#558) may auto-board the new PR and set Status. VERIFY the PR's board Status; only set `8bf9192f` (`Ready for review`) via GraphQL if the hook did not.

- [ ] **Step 4: Offer `@claude review`**

Per `feedback_github_workflow.md`, offer the review pass. On approval only: `gh pr comment <N> --body "@claude review"` (the exact canonical trigger the hook allows).

---

## Task 7: Lab notebook + close ritual + merge

**Files:** Modify `research/lab_notebook/scientist.md`; Issue #546 + PR (AC/Test-plan ticks).

- [ ] **Step 1: Write the lab-notebook entry (AFTER review, before merge)**

Append a dated entry (`## 2026-05-29` / `### HH:MM UTC — Editor: Scientist` / `#### [Issue #546](url) ASNEO desk eval`) capturing the decision + rationale + deck path + Zotero key + that it reflects the post-review final state. Reference the PR # (closure-audit bot checks for it). Commit:

```bash
git add research/lab_notebook/scientist.md
git commit -m "research(notebook): ASNEO eval decision entry (Issue #546)"
git push
```

- [ ] **Step 2 (GATE): Tick the closure checklists**

Verify each #546 AC and each PR Test-plan box, then tick `- [x]` (or comment-defer with a link). Never leave aspirational unchecked boxes.

- [ ] **Step 3: Merge via the closure-ritual gate**

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER> --squash --delete-branch
```

The script blocks the merge if any `- [ ]` remains on the PR Test plan or any linked Issue's ACs. On clean audit it forwards to `gh pr merge`.

- [ ] **Step 4: Confirm board auto-flip**

#546 + the PR should auto-flip to **Done** on merge (never set manually). If #546 used `advances` (no auto-close), tick-with-defer-link or close manually per the #258 closure path.

---

## Self-Review (against the spec)

**Spec coverage:**
- Zotero verified 3-section note → Task 2. ✅
- Repo state (maintenance/lang/weights/license) → Task 1 Step 1 (asneo-repo agent). ✅
- Head-to-head across 5 axes → Task 1 (gather) + Task 3 slide 4. ✅
- 5-mode decision recorded → Task 1 Step 4 + Task 5 Step 1. ✅
- Tool-primer deck, visually verified → Tasks 3 + 4. ✅
- #258 conditional update → Task 5 Step 2. ✅
- Lab-notebook entry → Task 7 Step 1. ✅
- Board-backed actionables → Task 5 Step 3. ✅
- Lifecycle (status/branch/commit-push-separate/PR/@claude/close) → Tasks 6 + 7. ✅

**Placeholder scan:** The `<from findings.*>` slots in Task 3 are data-dependent values resolved by Task 1 (legitimate, not deferred decisions); the schema field names are concrete and defined in Task 1. The `<KEY>`, `<N>`, `<path>` slots are values captured by an explicit earlier step. No "TBD"/"add error handling"-class deferrals.

**Type consistency:** Deck slots in Task 3 reference `findings.asneo.{citation,license,repo_url,repo_last_activity,language,junction_backbone,normal_filtering,peptide_assembly,frame_handling,mhc_predictor,validation_cohort,validation_n}`, `findings.ours.{junction_backbone,normal_filtering,peptide_assembly,frame_handling,mhc_predictor}`, and `findings.modes.{triage,replacement,cross_check,component_reuse,reject}` — all defined in the Task 1 `FINDINGS` schema. Consistent.
