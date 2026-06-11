# Arc Work-Structuring — Foundation (Plan 1 of 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the arc dimension as GitHub labels, apply the validated v1 taxonomy to the 75 open issues, set the opening active slate, and make arcs queryable from the board tool.

**Architecture:** Arcs are orthogonal labels (`arc:<slug>`) plus a focus marker (`arc-phase:active|next|later`), per the design spec. A committed TSV manifest (`scripts/pm/arc_taxonomy.tsv`) is the source of truth; two idempotent bash scripts create the labels and apply them; `board_open_items.py` gains `--arc`/`--arc-phase` filters so the taxonomy is queryable (and so Plan 2's pull rule has a query layer). All work here is in the **project repo** (`splice-neoepitope-pipeline`) + GitHub-side label/issue edits — no personas-repo or memory changes (those are Plan 2).

**Tech Stack:** `gh` CLI (labels + issue edits + ProjectV2 GraphQL), Bash, Python 3 (stdlib only), pytest (via `workflow/tests/.venv`).

**Spec:** [`docs/superpowers/specs/2026-06-05-arc-work-structuring-design.md`](../specs/2026-06-05-arc-work-structuring-design.md) (anchored to Issue #633).

---

## File Structure

| Path | Repo | Responsibility |
|---|---|---|
| `scripts/pm/arc_labels.sh` | project | Create/update the 8 `arc:*` + 3 `arc-phase:*` labels (idempotent, `--force`) |
| `scripts/pm/arc_taxonomy.tsv` | project | v1 taxonomy — source of truth: `arc_slug  phase  issue…` per row |
| `scripts/pm/apply_arc_labels.sh` | project | Apply arc + arc-phase labels to issues per the manifest (idempotent) |
| `scripts/board_open_items.py` | project | **Modify** — derive `arc`/`arc_phase` fields + add `--arc`/`--arc-phase` filters |
| `scripts/tests/test_board_open_items_arc.py` | project | **Create** — unit tests for the arc derivation + filters |

**Out of scope for Plan 1** (see "Follow-on Plans"): the shared pull rule, board-hygiene un-arced check, milestone cleanup, the arc-review process.

---

## Task 1: Arc + arc-phase label definitions

**Files:**
- Create: `scripts/pm/arc_labels.sh`

- [ ] **Step 1: Write the label-creation script**

```bash
#!/usr/bin/env bash
# scripts/pm/arc_labels.sh
# Idempotently create the arc:* (narrative throughline) and arc-phase:* (focus
# marker) labels. Re-runnable: --force updates colour/description if the label
# already exists. Spec: docs/superpowers/specs/2026-06-05-arc-work-structuring-design.md
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"

# name | hex colour (no '#') | description
ARC_LABELS=(
  "arc:aligner-junctions|0F766E|Arc: aligner & junction extraction (STAR/HISAT2 verification, SJ/annotation correctness)"
  "arc:junction-filtering|047857|Arc: junction filtering & tumor-specificity (GTEx + matched-normal + AlphaGenome, integrity)"
  "arc:variant-cohort|15803D|Arc: variant calling & cohort expansion (somatic calling, SpliceAI/MMSplice, new patients)"
  "arc:scoring-tcr-pmhc|B45309|Arc: neoepitope scoring & TCR-pMHC modeling (MHCflurry calibration, scorers, structures)"
  "arc:results-reporting|A16207|Arc: results, reporting & manuscript (HLA panel runs, report layer, decks, writeup)"
  "arc:cloud-reproducibility|1D4ED8|Arc: cloud execution & reproducibility (Batch, run registry, GPU, run_cloud_gpu.sh)"
  "arc:memory-methodology|7C3AED|Arc: memory & methodology (MEMORY.md slimming/audit, MM role, persona framing)"
  "arc:board-governance|9D174D|Arc: board governance & enforcement (Kanban governance, hooks/guards, PM-tooling evals)"
)
PHASE_LABELS=(
  "arc-phase:active|2DA44E|Arc focus: actively pulled now (cap 3)"
  "arc-phase:next|D4A72C|Arc focus: queued next"
  "arc-phase:later|8C959F|Arc focus: parked"
)

create() {
  local spec="$1" name color desc
  IFS='|' read -r name color desc <<<"$spec"
  gh label create "$name" --repo "$REPO" --color "$color" --description "$desc" --force
}
for spec in "${ARC_LABELS[@]}" "${PHASE_LABELS[@]}"; do create "$spec"; done
echo "Created/updated ${#ARC_LABELS[@]} arc + ${#PHASE_LABELS[@]} arc-phase labels on $REPO."
```

- [ ] **Step 2: Make it executable and run it**

Run: `chmod +x scripts/pm/arc_labels.sh && bash scripts/pm/arc_labels.sh`
Expected: 11 lines of `gh` output (created/updated) + the final summary line.

- [ ] **Step 3: Verify all 11 labels exist**

Run:
```bash
gh label list --repo Jin-HoMLee/splice-neoepitope-pipeline --limit 200 \
  --json name -q '[.[]|select(.name|startswith("arc"))]|length'
```
Expected: `11`

- [ ] **Step 4: Commit**

```bash
git add scripts/pm/arc_labels.sh
git commit -m "feat(pm): arc:* + arc-phase:* label definitions (#633)"
```

---

## Task 2: The v1 arc taxonomy manifest

**Files:**
- Create: `scripts/pm/arc_taxonomy.tsv`

- [ ] **Step 1: Write the manifest** (whitespace-delimited; `#` lines are comments; column 1 = arc slug, column 2 = phase, remaining tokens = issue numbers)

```
# scripts/pm/arc_taxonomy.tsv — v1 arc taxonomy (2026-06-05)
# Source of truth for apply_arc_labels.sh. Spec §5:
# docs/superpowers/specs/2026-06-05-arc-work-structuring-design.md
# Columns (whitespace-delimited): arc_slug  arc-phase  issue...
arc:aligner-junctions     active  297 375 377 378 411 636
arc:junction-filtering    later   126 212 304 381 594 663
arc:variant-cohort        later   413 416 436 437 438 440
arc:scoring-tcr-pmhc      active  433 492 547 566 585 601 659
arc:results-reporting     later   193 198 233 435 455
arc:cloud-reproducibility next    66 183 195 310 630 651 658 664 669 673 674
arc:memory-methodology    later   248 265 324 326 346 353 527 538 539 540 541 542 672
arc:board-governance      active  234 250 294 295 406 445 498 499 533 553 569 578 617 626 633 655 665
# unfiled (intentionally NO arc label): 446 451 570 641
```

- [ ] **Step 2: Verify the manifest covers exactly 71 issues across 8 arcs (75 open − 4 unfiled)**

Run:
```bash
awk '!/^#/ && NF {for(i=3;i<=NF;i++)c++} END{print c}' scripts/pm/arc_taxonomy.tsv
```
Expected: `71`

- [ ] **Step 3: Verify no issue is double-assigned**

Run:
```bash
awk '!/^#/ && NF {for(i=3;i<=NF;i++)print $i}' scripts/pm/arc_taxonomy.tsv | sort | uniq -d
```
Expected: (empty output — no duplicates)

- [ ] **Step 4: Commit**

```bash
git add scripts/pm/arc_taxonomy.tsv
git commit -m "feat(pm): v1 arc taxonomy manifest — 71 issues across 8 arcs (#633)"
```

---

## Task 3: Apply the taxonomy to the issues

**Files:**
- Create: `scripts/pm/apply_arc_labels.sh`

- [ ] **Step 1: Write the apply script**

```bash
#!/usr/bin/env bash
# scripts/pm/apply_arc_labels.sh
# Apply arc:* + arc-phase:* labels to issues per scripts/pm/arc_taxonomy.tsv.
# Idempotent: gh add-label of an already-present label is a no-op.
# Prereq: run scripts/pm/arc_labels.sh first (labels must exist).
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:-$SCRIPT_DIR/arc_taxonomy.tsv}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

while read -r arc phase rest; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  for n in $rest; do
    echo "  #$n -> $arc + arc-phase:$phase"
    gh issue edit "$n" --repo "$REPO" \
      --add-label "$arc" --add-label "arc-phase:$phase"
  done
done < "$MANIFEST"
echo "Done applying arc taxonomy from $MANIFEST."
```

- [ ] **Step 2: Dry-run-read the manifest to confirm parsing (no GitHub writes yet)**

Run:
```bash
while read -r arc phase rest; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  echo "$arc [$phase]: $(echo $rest | wc -w) issues"
done < scripts/pm/arc_taxonomy.tsv
```
Expected:
```
arc:aligner-junctions [active]: 6 issues
arc:junction-filtering [later]: 6 issues
arc:variant-cohort [later]: 6 issues
arc:scoring-tcr-pmhc [active]: 7 issues
arc:results-reporting [later]: 5 issues
arc:cloud-reproducibility [next]: 11 issues
arc:memory-methodology [later]: 13 issues
arc:board-governance [active]: 17 issues
```

- [ ] **Step 3: Run the apply (writes labels to 71 issues)**

Run: `chmod +x scripts/pm/apply_arc_labels.sh && bash scripts/pm/apply_arc_labels.sh`
Expected: 71 `#N -> arc:… + arc-phase:…` lines + final "Done" line.

- [ ] **Step 4: Verify per-arc issue counts match the manifest**

Run:
```bash
for a in aligner-junctions:6 junction-filtering:6 variant-cohort:6 scoring-tcr-pmhc:7 \
         results-reporting:5 cloud-reproducibility:11 memory-methodology:13 board-governance:17; do
  slug=${a%:*}; want=${a#*:}
  got=$(gh issue list --repo Jin-HoMLee/splice-neoepitope-pipeline --state open \
        --label "arc:$slug" --json number -q 'length')
  printf 'arc:%-22s got=%s want=%s %s\n' "$slug" "$got" "$want" \
    "$([ "$got" = "$want" ] && echo OK || echo MISMATCH)"
done
```
Expected: all 8 lines end in `OK`.

- [ ] **Step 5: Verify the active slate (3 arcs = 30 issues) and that unfiled stay unlabeled**

Run:
```bash
echo -n "active: "; gh issue list --repo Jin-HoMLee/splice-neoepitope-pipeline --state open \
  --label arc-phase:active --json number -q 'length'    # expect 30
echo -n "next:   "; gh issue list --repo Jin-HoMLee/splice-neoepitope-pipeline --state open \
  --label arc-phase:next --json number -q 'length'      # expect 11
for n in 446 451 570 641; do
  c=$(gh issue view "$n" --repo Jin-HoMLee/splice-neoepitope-pipeline \
      --json labels -q '[.labels[].name|select(startswith("arc"))]|length')
  echo "unfiled #$n arc-labels=$c"                       # expect 0 each
done
```
Expected: `active: 30`, `next: 11`, and `arc-labels=0` for all four unfiled issues.

- [ ] **Step 6: Commit**

```bash
git add scripts/pm/apply_arc_labels.sh
git commit -m "feat(pm): apply_arc_labels.sh — sync arc taxonomy to issues (#633)"
```

---

## Task 4: Make arcs queryable in `board_open_items.py`

**Files:**
- Modify: `scripts/board_open_items.py` (`normalize()`, `matches_filter()`, argparse in `main()`)
- Create: `scripts/tests/test_board_open_items_arc.py`

- [ ] **Step 1: Write the failing unit test**

```python
# scripts/tests/test_board_open_items_arc.py
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402


def _item(labels):
    return {
        "content": {
            "__typename": "Issue", "number": 1, "title": "t", "url": "u", "state": "OPEN",
            "createdAt": None, "updatedAt": None, "closedAt": None,
            "labels": {"nodes": [{"name": n} for n in labels]},
        },
        "fieldValues": {"nodes": [{"name": "Ready", "field": {"name": "Status"}}]},
    }


def _args(**kw):
    base = dict(role=None, status=None, priority=None, size=None, arc=None, arc_phase=None)
    base.update(kw)
    return argparse.Namespace(**base)


def test_normalize_derives_arc_and_phase():
    n = b.normalize(_item(["role:pm", "arc:board-governance", "arc-phase:active"]))
    assert n["arc"] == "arc:board-governance"
    assert n["arc_phase"] == "active"


def test_normalize_arc_none_when_absent():
    n = b.normalize(_item(["role:pm"]))
    assert n["arc"] is None and n["arc_phase"] is None


def test_arc_phase_label_not_mistaken_for_arc():
    n = b.normalize(_item(["arc-phase:later"]))
    assert n["arc"] is None and n["arc_phase"] == "later"


def test_filter_by_arc_phase():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc_phase="active"))
    assert not b.matches_filter(it, _args(arc_phase="later"))


def test_filter_by_arc_slug_accepts_short_and_full_form():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc="scoring-tcr-pmhc"))
    assert b.matches_filter(it, _args(arc="arc:scoring-tcr-pmhc"))
    assert not b.matches_filter(it, _args(arc="cloud-reproducibility"))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest scripts/tests/test_board_open_items_arc.py -v`
Expected: FAIL — `KeyError: 'arc'` in `test_normalize_*` and `AttributeError`/`KeyError` on `arc`/`arc_phase` (the fields and filters don't exist yet).

- [ ] **Step 3: Add arc derivation in `normalize()`**

In `scripts/board_open_items.py`, find (inside `normalize()`):
```python
    role = next((l for l in labels if l.startswith("role:")), None)
    return {
```
Replace with:
```python
    role = next((l for l in labels if l.startswith("role:")), None)
    arc = next((l for l in labels if l.startswith("arc:")), None)
    arc_phase = next(
        (l.removeprefix("arc-phase:") for l in labels if l.startswith("arc-phase:")),
        None,
    )
    return {
```

Then in the same returned dict, find:
```python
        "role": role,
        "labels": labels,
```
Replace with:
```python
        "role": role,
        "arc": arc,
        "arc_phase": arc_phase,
        "labels": labels,
```

- [ ] **Step 4: Add the filters in `matches_filter()`**

Find:
```python
    if args.size and it["size"] != args.size:
        return False
    return True
```
Replace with:
```python
    if args.size and it["size"] != args.size:
        return False
    if args.arc:
        want = args.arc if args.arc.startswith("arc:") else f"arc:{args.arc}"
        if it["arc"] != want:
            return False
    if args.arc_phase and it["arc_phase"] != args.arc_phase:
        return False
    return True
```

- [ ] **Step 5: Add the argparse flags in `main()`**

Find:
```python
    p.add_argument("--size", choices=list(SIZE_ORDER), help="Filter by Size")
```
Insert immediately after it:
```python
    p.add_argument("--arc", help='Filter by arc label (e.g. "scoring-tcr-pmhc" or "arc:scoring-tcr-pmhc")')
    p.add_argument("--arc-phase", dest="arc_phase", choices=["active", "next", "later"],
                   help="Filter by arc focus phase")
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest scripts/tests/test_board_open_items_arc.py -v`
Expected: PASS (5 passed).

- [ ] **Step 7: Smoke-test against the live board**

Run: `python3 scripts/board_open_items.py --arc-phase active --json 2>/dev/null | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))'`
Expected: a number close to `30` (open board items in the active slate; equals the active issue count unless an issue isn't on the board).

- [ ] **Step 8: Commit**

```bash
git add scripts/board_open_items.py scripts/tests/test_board_open_items_arc.py
git commit -m "feat(pm): board_open_items --arc/--arc-phase filters + arc derivation (#633)"
```

---

## Task 5: Open the PR

- [ ] **Step 1: Push the branch**

Run: `git push -u origin design/pm/issue-633-arc-work-structuring-spec`

- [ ] **Step 2: Open the PR** (note: `@-claude` hyphen form is intentional — never the literal trigger)

```bash
gh pr create --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --title "feat(pm): arc work-structuring — foundation (labels + taxonomy + queryability)" \
  --body "$(cat <<'BODY'
Implements Plan 1 of the arc work-structuring rollout (design spec + plan under docs/superpowers/).

- arc:* + arc-phase:* label definitions (scripts/pm/arc_labels.sh)
- v1 taxonomy manifest, 71 issues / 8 arcs (scripts/pm/arc_taxonomy.tsv)
- apply script (scripts/pm/apply_arc_labels.sh) — taxonomy synced to live issues
- board_open_items.py --arc/--arc-phase filters + unit tests

Opening active slate: arc:aligner-junctions, arc:scoring-tcr-pmhc, arc:board-governance.
Plan 2 (shared pull rule + board hygiene, personas repo) and Plan 3 (milestone cleanup) follow.

Related to the Issue #633 board-governance review.

## Test plan
- [ ] `gh label list` shows 11 arc/arc-phase labels
- [ ] per-arc issue counts match the manifest (8x OK)
- [ ] active=30, next=11, unfiled (446/451/570/641) carry no arc label
- [ ] `pytest scripts/tests/test_board_open_items_arc.py` passes
BODY
)"
```

- [ ] **Step 3: Auto-board the PR** if not auto-added (the post-create hook normally handles this; verify it landed on project #9).

- [ ] **Step 4: Merge via the closure-ritual gate** (after review per the routine-PR bot-review rule if applicable)

Run: `bash scripts/audit_and_merge.sh <PR_NUMBER> --squash --delete-branch`

---

## Follow-on Plans (not implemented here)

**Plan 2 — Arc-aware pull + board hygiene (personas repo, MM-landed).** Depends on Plan 1's labels existing.
- `shared/feedback_best_next_issue.md` Step 1: restrict the candidate pool to `arc-phase:active` (one edit → all four roles get arc-aware pull). This is the spec's *core deliverable* (§8).
- `shared/feedback_board_hygiene.md`: add the "un-arced `Ready` issue = triage gap" sweep check + surface the active slate in the PM morning routine.
- New shared rule: the cadenced **arc review** (split/merge/retire + re-pick ≤3 active; §6).
- Update `CLAUDE.md` board-governance section + `feedback_milestones.md` for the three-axis model.
- Landed via the personas repo PR flow (MM / role-self-commit per Issue #672).

**Plan 3 — Milestone de-overloading & cleanup (project repo + GitHub).** Independent; can follow Plan 1.
- Delete dead `M1/M2/M7`; close/populate empty placeholder milestones.
- Drop the arc-suffix from live milestone names; resolve Milestones-vs-Iteration-field (Issue #633 open question).

---

## Self-Review

**Spec coverage (Plan 1 scope only):** §4.3 carriers (arc + arc-phase labels) → Tasks 1,3. §5 taxonomy → Task 2 (manifest) + Task 3 (application). §5.1 active slate → Task 2/3 (phase column). Queryability for §8 pull → Task 4. Milestone de-overloading (§7), pull rule (§8), arc review (§6 governance ops), board hygiene (§9 Phase 1 sweep) → **deferred to Plans 2/3** (flagged, not dropped).

**Placeholder scan:** none — every script and edit shows full content; `<PR_NUMBER>` in Task 5 is a runtime value, not a content placeholder.

**Type/name consistency:** arc slugs identical across `arc_labels.sh`, `arc_taxonomy.tsv`, `apply_arc_labels.sh`, and the verification loops. `arc`/`arc_phase` field names consistent between `normalize()`, the returned dict, `matches_filter()`, and the test. `--arc-phase` → `dest="arc_phase"` matches `args.arc_phase`.

**Count reconciliation:** 6+6+6+7+5+11+13+17 = 71 arc-labeled + 4 unfiled = 75 open. Active = 6+7+17 = 30; next = 11; later = 6+6+5+13 = 30.
