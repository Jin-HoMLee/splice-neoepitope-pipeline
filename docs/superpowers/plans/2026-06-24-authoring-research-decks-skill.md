# authoring-research-decks Skill Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract the four research-artifact convention sections out of always-on `AGENTS.md` into a first project skill (`.claude/skills/authoring-research-decks/SKILL.md`), leaving a one-line pointer-stub behind, proving the slimming pattern from epic #859.

**Architecture:** A pure content relocation - no code paths touched. The four sections (lines 52-106 of `AGENTS.md`: experiment-notebook layout + the three Quarto slide-deck tiers) move **verbatim** into a new `SKILL.md` body under a description-triggered frontmatter. `AGENTS.md` keeps a single pointer-stub so the topic stays discoverable in resident context. Fidelity is verified by a git-diff cross-check: the lines removed from `AGENTS.md` must equal the lines added to the skill body.

**Tech Stack:** Markdown, YAML frontmatter, git, Claude Code project-skill loading (`.claude/skills/`).

## Global Constraints

- The moved content is **byte-for-byte verbatim** - no paraphrasing, no reformatting, no heading-level changes. The `##`/`###` headings travel as-is into the skill body.
- One sentence per physical line is **not** retrofitted onto the moved content - it moves exactly as it currently reads (verbatim rule wins).
- No em dash `-` is introduced in any new prose (stub, frontmatter); use a plain dash. (The moved content keeps whatever it already contains - verbatim.)
- Commit messages carry **no** `Co-Authored-By` agent trailer.
- Branch is `docs/pm/issue-860-authoring-decks-skill-pilot` (already checked out). Do not create a new branch.
- Reference other Issues/PRs as genuine `#N` links only; never use a bare `#N` for a positional/step number.
- The four target sections are, in file order: `## Experiment notebooks live under \`research/experiments/\``, `## Slide decks for experiment Issues`, `## Slide decks for eval Issues`, `## Slide decks for research-decision Issues`. They are contiguous (lines 52-106), bounded above by the `### PDB chain relabelling` block (ends line 50) and below by `## MHC Presentation Vocabulary` (line 108).

---

## File Structure

- **Create:** `.claude/skills/authoring-research-decks/SKILL.md` - the extracted skill. Frontmatter (`name`, `description`) + the four sections verbatim. This is the first entry under `.claude/skills/` (the directory does not exist yet).
- **Modify:** `AGENTS.md` - delete the four sections (lines 52-106), insert the one-line pointer-stub in their place.
- **Modify (final task only):** Issue #860 acceptance criteria (via `gh`, not a file) + open the PR.

No tests directory - this is a documentation move. Verification is structural (file existence, content equality, section absence, skill load) rather than a pytest cycle.

---

## Task 1: Create the skill file with the four sections verbatim

**Files:**
- Create: `.claude/skills/authoring-research-decks/SKILL.md`

**Interfaces:**
- Consumes: nothing (first task).
- Produces: a skill named `authoring-research-decks` whose body is the verbatim text currently at `AGENTS.md` lines 52-106. Task 2 relies on that body matching the removed `AGENTS.md` block exactly.

- [ ] **Step 1: Capture the exact source block to a scratch file**

This is the source of truth for the verbatim move. Extract lines 52-106 of `AGENTS.md` (the four sections, ending at the blank line 107 before `## MHC Presentation Vocabulary`).

Run:
```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
sed -n '52,106p' AGENTS.md > /tmp/decks_block.md
head -1 /tmp/decks_block.md; echo "..."; tail -1 /tmp/decks_block.md; wc -l /tmp/decks_block.md
```
Expected: first line is `## Experiment notebooks live under \`research/experiments/\``, last line is the `## Slide decks for research-decision Issues` Format line ending in `New decks use \`issue_NNN_<short>\`.`, and the count is 55 lines.

- [ ] **Step 2: Create the skill directory and write the frontmatter + body**

Create `.claude/skills/authoring-research-decks/SKILL.md`. The frontmatter is exactly as the spec specifies; the body is the captured block, unchanged.

```bash
mkdir -p .claude/skills/authoring-research-decks
{
cat <<'EOF'
---
name: authoring-research-decks
description: >-
  Conventions for authoring research artifacts in this repo - experiment
  notebooks (research/experiments/issue_NNN_<short>/ layout, cross-experiment
  data sharing, size bands) and the three Quarto slide-deck tiers (experiment /
  eval / research-decision). Use when creating or editing a slide deck, an
  experiment notebook, or their outputs/ folder.
---

EOF
cat /tmp/decks_block.md
} > .claude/skills/authoring-research-decks/SKILL.md
```

- [ ] **Step 3: Verify the skill body equals the source block (fidelity check)**

The body (everything after the frontmatter) must be byte-identical to the captured block.

Run:
```bash
# Strip frontmatter (lines through the second '---' and the blank line after) and compare
awk 'f{print} /^---$/{c++} c==2 && !f {f=1}' .claude/skills/authoring-research-decks/SKILL.md | sed '1{/^$/d}' > /tmp/skill_body.md
diff /tmp/decks_block.md /tmp/skill_body.md && echo "VERBATIM OK"
```
Expected: `VERBATIM OK` (no diff output). If diff shows anything, the body was altered - fix until the diff is empty.

- [ ] **Step 4: Verify the frontmatter is valid YAML with the required keys**

Run:
```bash
conda activate snakemake
python -c "
import yaml, pathlib, re
t = pathlib.Path('.claude/skills/authoring-research-decks/SKILL.md').read_text()
fm = re.match(r'^---\n(.*?)\n---\n', t, re.S).group(1)
d = yaml.safe_load(fm)
assert d['name'] == 'authoring-research-decks', d.get('name')
assert 'Use when' in d['description'], 'description must state WHEN to use'
print('FRONTMATTER OK:', d['name'])
"
```
Expected: `FRONTMATTER OK: authoring-research-decks`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/authoring-research-decks/SKILL.md
git commit -m "feat(skills): add authoring-research-decks skill (verbatim deck/notebook conventions) (#859, #860)"
```

---

## Task 2: Replace the four sections in AGENTS.md with the pointer-stub

**Files:**
- Modify: `AGENTS.md` (delete lines 52-106, insert the stub)

**Interfaces:**
- Consumes: the skill body from Task 1 (used as the equality reference for the removed block).
- Produces: an `AGENTS.md` with the four `##` sections gone and one stub line present; downstream tasks rely on the four headings being absent.

- [ ] **Step 1: Confirm the boundary lines before editing**

Line numbers can shift if Task 1 somehow touched `AGENTS.md` (it should not have). Re-confirm.

Run:
```bash
grep -nE '^## (Experiment notebooks live under|Slide decks for|MHC Presentation Vocabulary)' AGENTS.md
```
Expected: `52:## Experiment notebooks live under ...`, `82:## Slide decks for experiment Issues`, `88:## Slide decks for eval Issues`, `98:## Slide decks for research-decision Issues`, `108:## MHC Presentation Vocabulary`. If 52 and 108 differ, recompute the delete range as `<first> .. <108-line minus 2>` (the blank line before MHC).

- [ ] **Step 2: Delete the four sections and insert the stub**

Replace lines 52-106 (the four sections + their trailing content, leaving the blank line 107 intact as the separator before MHC) with the single stub line. Use an explicit, auditable edit:

```bash
# Delete lines 52-106, then insert the stub at line 52
sed -i '' '52,106d' AGENTS.md
sed -i '' '52i\
**Authoring an experiment notebook or Quarto slide deck?** Layout, cross-experiment sharing, size bands, and the 3 deck tiers (experiment / eval / research-decision) live in the `authoring-research-decks` skill.
' AGENTS.md
```
(Note: this is macOS `sed` - the `-i ''` and the `i\` newline form are BSD-specific, matching the dev platform.)

- [ ] **Step 3: Verify the four headings are gone and the stub is present**

Run:
```bash
echo "--- headings that must be ABSENT ---"
grep -nE '^## (Experiment notebooks live under|Slide decks for (experiment|eval|research-decision))' AGENTS.md || echo "ALL FOUR ABSENT (good)"
echo "--- stub that must be PRESENT ---"
grep -n 'authoring-research-decks` skill' AGENTS.md
echo "--- context around the stub ---"
sed -n '49,55p' AGENTS.md
```
Expected: `ALL FOUR ABSENT (good)`; the stub line is found at line 52; the context shows `### PDB chain relabelling` content, a blank line, the stub, a blank line, then `## MHC Presentation Vocabulary`.

- [ ] **Step 4: Verify the removed block exactly equals the skill body (no content lost)**

The strongest fidelity gate: what left `AGENTS.md` must be exactly what landed in the skill.

Run:
```bash
# Reconstruct the removed lines from the diff against the parent commit's AGENTS.md
git show HEAD~1:AGENTS.md | sed -n '52,106p' > /tmp/removed_block.md
diff /tmp/removed_block.md /tmp/decks_block.md && echo "REMOVED == SKILL SOURCE (good)"
```
Expected: `REMOVED == SKILL SOURCE (good)`. (`HEAD~1` is the Task 1 commit, whose `AGENTS.md` still had the sections.)

- [ ] **Step 5: Verify the line-count drop**

Run:
```bash
echo "AGENTS.md now: $(wc -l < AGENTS.md) lines (was 468; expect ~414 = 468 - 55 removed + 1 stub - 0... allow 412-415)"
```
Expected: roughly 414 lines (468 minus 55 section lines plus 1 stub line; small variance from blank-line handling is fine).

- [ ] **Step 6: Confirm CLAUDE.md still resolves through the symlink**

`CLAUDE.md` is a symlink to `AGENTS.md`; the edit must be visible through it.

Run:
```bash
test -L CLAUDE.md && grep -q 'authoring-research-decks` skill' CLAUDE.md && echo "CLAUDE.md symlink reflects the stub (good)"
```
Expected: `CLAUDE.md symlink reflects the stub (good)`.

- [ ] **Step 7: Commit**

```bash
git add AGENTS.md
git commit -m "docs(spec): replace 4 deck/notebook sections with skill pointer-stub (#859, #860)"
```

---

## Task 3: Verify the skill loads, update Issue #860, open the PR

**Files:**
- Modify: Issue #860 acceptance criteria (via `gh`)
- No file changes in this task beyond the PR.

**Interfaces:**
- Consumes: the committed skill (Task 1) and the stubbed `AGENTS.md` (Task 2).
- Produces: a verified-loading skill and an open PR off the pilot branch.

- [ ] **Step 1: Reload plugins and confirm the skill is discovered**

This step is interactive - it runs Claude Code slash commands in the live session.

Run in the Claude Code session:
```
/reload-plugins
```
Then verify the skill appears. In a fresh check, ask the session to list skills or invoke the Skill tool's discovery; the `authoring-research-decks` skill must be listed with its description.
Expected: `authoring-research-decks` shows up as an available project skill. If it does not appear, stop - the project-skill-loading mechanism is the pilot's core risk; do not proceed to PR until it loads.

- [ ] **Step 2: Confirm config health**

Run in the Claude Code session:
```
/doctor
```
Expected: no errors introduced by the new `.claude/skills/` directory; hooks still listed under `/hooks` (the skill add must not have disturbed `settings.json`).

- [ ] **Step 3: Tick the Issue #860 acceptance criteria that are now met**

Fetch the current body, tick the boxes the pilot satisfies (skill created verbatim; sections replaced by stub; diff accounts for moved lines; skill loads after `/reload-plugins`), and leave any genuinely-deferred box with a linked carrier. Use `gh issue edit 860 --body-file <edited>` after editing locally.

Run:
```bash
gh issue view 860 --json body --jq .body > /tmp/issue860.md
# edit /tmp/issue860.md: change the satisfied "- [ ]" lines to "- [x]"
gh issue edit 860 --body-file /tmp/issue860.md
```
Expected: the closure-ritual-gated boxes on #860 are ticked, so `scripts/audit_and_merge.sh` will not block on them later.

- [ ] **Step 4: Write the PM lab-notebook entry**

Per the closure ritual, a role-tagged Issue needs a `## <merge-date>` entry in `research/lab_notebook/pm.md` referencing the PR/Issue. Add it now (it can reference the PR number after Step 5; write a minimal entry, finalize the PR ref post-create).

Add to `research/lab_notebook/pm.md` a dated entry summarizing: first project skill extracted from `AGENTS.md`; pattern proven for epic #859; verbatim move + stub; resident `AGENTS.md` dropped ~55 lines.

- [ ] **Step 5: Open the PR**

```bash
git push -u origin docs/pm/issue-860-authoring-decks-skill-pilot
gh pr create --fill --base main \
  --title "feat(skills): pilot authoring-research-decks skill - first AGENTS.md slimming extraction (#860)" \
  --body "$(cat <<'EOF'
Pilots the AGENTS.md slimming pattern (epic #859) by extracting the four research-artifact convention sections into the first project skill.

## What changed
- New `.claude/skills/authoring-research-decks/SKILL.md` - the experiment-notebook layout + three Quarto deck tiers, moved **verbatim**.
- `AGENTS.md`: those four sections replaced by a one-line pointer-stub (~55 resident lines removed).

## Verification
- Removed block == skill body (git-diff cross-check, no content lost).
- Skill loads after `/reload-plugins`; `/doctor` clean; hooks intact.

## Test plan
- [x] Skill body byte-identical to the removed AGENTS.md block
- [x] Frontmatter valid YAML with name + when-to-use description
- [x] Four `##` sections absent from AGENTS.md; stub present
- [x] CLAUDE.md symlink reflects the stub
- [x] Skill appears after /reload-plugins

Closes #860.
EOF
)"
```
Expected: PR opens against `main`, links #860. Do not merge here - merge runs later via `scripts/audit_and_merge.sh 860` after the bot-review offer.

- [ ] **Step 6: Offer the bot review**

After the PR is open, offer a `@-claude review` (the non-trivial-PR convention). Post the canonical trigger comment if the user wants the bot pass before merge.

---

## Self-Review

**Spec coverage** (against `docs/superpowers/specs/2026-06-24-agents-md-skill-extraction-design.md`):
- Pilot design (`authoring-research-decks`, the four sections, the frontmatter) -> Task 1. ✓
- Pointer-stub replacing the four sections -> Task 2 (stub text matches spec line). ✓
- "Verify the skill loads after `/reload-plugins`" (rollout step 1 + acceptance criterion) -> Task 3 Step 1. ✓
- "git diff accounts for every moved line; the move is verbatim" (Verification section) -> Task 1 Step 3 + Task 2 Step 4. ✓
- "Resident-line count drops by ~56 lines" (Verification) -> Task 2 Step 5 (55 lines moved + 1 stub). ✓
- Out of scope here (later epic #859 follow-ups): the skill *candidates* needing evals, the reference-file extractions, the trim-to-pointer pass. Correctly excluded - this plan is the pilot only.

**Placeholder scan:** No TBD/TODO. The one judgment step (Task 3 Step 3, "tick the satisfied boxes") is inherent to a live Issue edit and is bounded by naming the specific boxes. Acceptable.

**Type/name consistency:** Skill name `authoring-research-decks` and path `.claude/skills/authoring-research-decks/SKILL.md` are identical across Tasks 1, 2, 3 and the stub. Line range 52-106 is consistent and re-confirmed before the destructive edit (Task 2 Step 1). ✓

**Note on verbatim vs the one-sentence-per-line house rule:** the global "one sentence per line" rule is intentionally *not* applied to the moved content (it already lives as-is in `AGENTS.md`); enforcing it would violate the stronger verbatim-move constraint and break the diff-equality gate. Flagged in Global Constraints so the implementer does not "helpfully" reflow.
