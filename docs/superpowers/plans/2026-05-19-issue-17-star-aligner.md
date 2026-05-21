# Issue #17 STAR Aligner Switch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Switch production aligner default from HISAT2 to STAR with the minimum-viable parameter changes that make the switch correctness-equivalent to the Nature paper command. Defer sensitivity tuning to a follow-up Issue.

**Architecture:** Four-file change — flip `aligner` default in `config/config.yaml`, pin `aligner: hisat2` in `config/test_config.yaml` to prevent the M1 8 GB local dev environment from inheriting the new default (STAR's 32 GB index won't fit), modify the `star_align` shell block in `workflow/rules/alignment.smk` to add `--twopassMode Basic` + `--limitSjdbInsertNsj 2000000` and remove three silent-throttle filters, and add a snapshot pytest that asserts the rendered shell command contains/lacks the expected flags.

**Tech Stack:** Snakemake 8, STAR aligner, pytest, conda envs (`workflow/envs/star.yaml`, already present).

**Reference:** Design spec at [docs/superpowers/specs/2026-05-19-issue-17-star-aligner-design.md](../specs/2026-05-19-issue-17-star-aligner-design.md).

---

## File Structure

**Files to modify:**
- `config/config.yaml` — flip `alignment.aligner` from `"hisat2"` to `"star"`
- `config/test_config.yaml` — add explicit `alignment.aligner: "hisat2"` override
- `workflow/rules/alignment.smk` — modify `star_align` shell block (~5-line diff)

**Files to create:**
- `workflow/tests/test_alignment_star_command.py` — pytest snapshot test of the rendered STAR shell command

**Files indirectly affected:**
- `dag.pdf` — regenerated via `bash scripts/visualize_dag.sh` to resolve a stray unstaged deletion and satisfy the DAG-visualization rule for Snakemake rule changes

---

## Task 0: Branch setup + spec commit

**Files:**
- Modify: project board (set Issue #17 Status → In Progress)
- Create: feature branch `feat/developer/issue-17-star-aligner`
- Commit: `docs/superpowers/specs/2026-05-19-issue-17-star-aligner-design.md` (already written, uncommitted)

- [ ] **Step 1: Sync local main with origin**

Run:
```bash
git fetch origin && git checkout main && git pull origin main
```
Expected: working tree on `main`, in sync with `origin/main`. If `D dag.pdf` is unstaged, leave it for Task 5 (regenerated then).

- [ ] **Step 2: Set Issue #17 Status to In Progress on the board**

Run:
```bash
ISSUE_ITEM_ID=$(gh api graphql -f query='
{ user(login: "Jin-HoMLee") { projectV2(number: 9) { items(first: 100) {
  nodes { id content { ... on Issue { number } } }
} } } }' --jq '.data.user.projectV2.items.nodes[] | select(.content.number==17) | .id')

gh api graphql -f query='
mutation($item: ID!) {
  updateProjectV2ItemFieldValue(input: {
    projectId: "PVT_kwHOB17eGc4BSomP",
    itemId: $item,
    fieldId: "PVTSSF_lAHOB17eGc4BSomPzhAHFf8",
    value: { singleSelectOptionId: "47fc9ee4" }
  }) { projectV2Item { id } }
}' -f item="$ISSUE_ITEM_ID"
```
Expected: mutation returns the item ID without an `errors` block.

If pagination misses #17 (board has 371+ items), increase the search by re-running with later pages or use the GraphQL `node(id: ...)` form to fetch directly. Field ID + option ID are from the shared MEMORY.md Always-in-effect block.

- [ ] **Step 3: Create the feature branch**

Run:
```bash
gh issue develop 17 --name feat/developer/issue-17-star-aligner --checkout
```
Expected: branch created, linked to Issue #17, currently checked out. `git status` shows clean tree on the new branch (still with `D dag.pdf` unstaged — handled in Task 5).

- [ ] **Step 4: Commit the spec doc (docs-first)**

Run:
```bash
git add docs/superpowers/specs/2026-05-19-issue-17-star-aligner-design.md
git commit -m "$(cat <<'EOF'
docs(spec): Issue #17 STAR aligner switch design

Spec for the minimum-viable STAR aligner switch: flips config default,
adds --twopassMode Basic + --limitSjdbInsertNsj 2000000, removes silent
multi-mapper throttle. Defers paper-fidelity sensitivity tuning to a
follow-up Issue with per-flag benchmarks.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands on `feat/developer/issue-17-star-aligner`.

---

## Task 1: Write the failing snapshot test

**Files:**
- Create: `workflow/tests/test_alignment_star_command.py`

The test renders the `star_align` shell command via `snakemake --dry-run --printshellcmds` and asserts the expected flag set. It must FAIL on current `main` (because `--twopassMode Basic` is absent and `--outSJfilterReads Unique` is present) so Task 2's rule edit can make it pass.

- [ ] **Step 1: Write the test file**

Create `workflow/tests/test_alignment_star_command.py` with this content:

```python
"""Snapshot test of the rendered star_align shell command.

Asserts the production STAR alignment command contains the flags required
for novel-junction sensitivity (--twopassMode Basic, --limitSjdbInsertNsj)
and does not silently throttle multi-mappers
(--outSJfilterReads Unique, --outSJfilterCountUniqueMin, --outSJfilterCountTotalMin).

The test runs `snakemake --dry-run --printshellcmds` programmatically with an
aligner=star config override and stub inputs in a tmp_path. The captured
shell output is asserted against substring presence/absence.

This is the first shell-level rule snapshot test in the codebase. Existing
tests cover pure Python helpers (test_strandness.py, test_bed12_to_junctions.py,
test_star_sj_to_junctions.py, etc.); this complements them at the rule layer.
"""
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def star_dry_run_output(tmp_path):
    """Render the star_align shell command and return captured stdout."""
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not on PATH — activate the snakemake conda env")

    # Stub FASTQs (must exist for `ancient(_get_fastq*)` to skip producer-rule lookup)
    fq1 = tmp_path / "test_R1.fq.gz"
    fq2 = tmp_path / "test_R2.fq.gz"
    fq1.touch()
    fq2.touch()

    # Stub STAR index (sentinel + dir must exist for star_align inputs to resolve)
    star_index_dir = tmp_path / "star_index"
    star_index_dir.mkdir()
    (star_index_dir / "index.done").touch()

    # Stub samples.tsv
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "patient_id\tsample_id\tsample_type\tfastq1\tfastq2\n"
        f"testpatient\ttestsample\tPrimary Tumor\t{fq1}\t{fq2}\n"
    )

    # Config override: aligner=star, point at stub paths
    override = tmp_path / "override.yaml"
    override.write_text(textwrap.dedent(f"""\
        samples_tsv: "{samples}"
        alignment:
          aligner: "star"
          star_index_dir: "{star_index_dir}"
    """))

    target = "results/testpatient/alignment/testsample/junctions.tsv"

    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "--printshellcmds",
            "--configfile", "config/config.yaml", str(override),
            "--",
            target,
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            f"snakemake dry-run failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )
    return result.stdout


def test_star_align_uses_twopass_mode(star_dry_run_output):
    """2-pass mode is the flag that gives STAR its novel-junction edge over HISAT2."""
    assert "--twopassMode Basic" in star_dry_run_output


def test_star_align_raises_sjdb_insert_limit(star_dry_run_output):
    """Default 1M cap can abort runs when 2-pass + multi-mappers discover more."""
    assert "--limitSjdbInsertNsj 2000000" in star_dry_run_output


def test_star_align_does_not_filter_to_unique_reads_only(star_dry_run_output):
    """Paper keeps multi-mappers; --outSJfilterReads Unique silently drops them."""
    assert "--outSJfilterReads Unique" not in star_dry_run_output


def test_star_align_does_not_throttle_min_unique_count(star_dry_run_output):
    """--outSJfilterCountUniqueMin throttles output below paper baseline."""
    assert "--outSJfilterCountUniqueMin" not in star_dry_run_output


def test_star_align_does_not_throttle_min_total_count(star_dry_run_output):
    """--outSJfilterCountTotalMin throttles output below paper baseline."""
    assert "--outSJfilterCountTotalMin" not in star_dry_run_output
```

- [ ] **Step 2: Run the new test and verify it FAILS as expected**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_alignment_star_command.py -v
```
Expected: tests collected. The fixture should produce captured stdout successfully. Then the 5 test functions split as follows on current `main`:
- `test_star_align_uses_twopass_mode` → **FAIL** (current rule lacks `--twopassMode Basic`)
- `test_star_align_raises_sjdb_insert_limit` → **FAIL** (current rule lacks `--limitSjdbInsertNsj 2000000`)
- `test_star_align_does_not_filter_to_unique_reads_only` → **FAIL** (current rule has `--outSJfilterReads Unique`)
- `test_star_align_does_not_throttle_min_unique_count` → **FAIL** (current rule has `--outSJfilterCountUniqueMin 1 1 1 1`)
- `test_star_align_does_not_throttle_min_total_count` → **FAIL** (current rule has `--outSJfilterCountTotalMin 1 1 1 1`)

If the fixture itself fails (e.g. `snakemake dry-run` errored), inspect the captured `STDOUT` / `STDERR` in the failure message and fix the stub setup before proceeding. Most likely cause: a different rule's input is unresolved → add another stub file or extend the config override.

- [ ] **Step 3: Commit the failing test**

Run:
```bash
git add workflow/tests/test_alignment_star_command.py
git commit -m "$(cat <<'EOF'
test(alignment): snapshot test for star_align shell command flags

Asserts --twopassMode Basic + --limitSjdbInsertNsj 2000000 present, and
the three --outSJfilter* throttles absent. First shell-level rule snapshot
test in the codebase — complements existing pure-Python helper tests.

Currently fails on 5 assertions; fixed in the next commit (Issue #17).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands. Failing tests are a deliberate intermediate state (TDD red).

---

## Task 2: Update `star_align` rule to make the test pass

**Files:**
- Modify: `workflow/rules/alignment.smk` lines 338-349 (`star_align` shell block)

- [ ] **Step 1: Edit the STAR command in alignment.smk**

The current shell block at [workflow/rules/alignment.smk:338-349](../../../workflow/rules/alignment.smk#L338-L349):

```bash
            STAR \\
                --runMode alignReads \\
                --runThreadN {threads} \\
                --genomeDir {input.index_dir} \\
                --readFilesIn $FASTQ_FILES \\
                $READCMD \\
                --outFileNamePrefix {params.output_prefix} \\
                --outSAMtype None \\
                --outSJfilterReads Unique \\
                --outSJfilterCountUniqueMin 1 1 1 1 \\
                --outSJfilterCountTotalMin 1 1 1 1 \\
                2>&1 | tee {log}
```

Replace with:

```bash
            STAR \\
                --runMode alignReads \\
                --runThreadN {threads} \\
                --genomeDir {input.index_dir} \\
                --readFilesIn $FASTQ_FILES \\
                $READCMD \\
                --outFileNamePrefix {params.output_prefix} \\
                --outSAMtype None \\
                --twopassMode Basic \\
                --limitSjdbInsertNsj 2000000 \\
                2>&1 | tee {log}
```

Net diff:
- ADD: `--twopassMode Basic \` (line after `--outSAMtype None \`)
- ADD: `--limitSjdbInsertNsj 2000000 \` (line after `--twopassMode Basic \`)
- REMOVE: `--outSJfilterReads Unique \`
- REMOVE: `--outSJfilterCountUniqueMin 1 1 1 1 \`
- REMOVE: `--outSJfilterCountTotalMin 1 1 1 1 \`

- [ ] **Step 2: Run the snapshot test and verify all 5 assertions PASS**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_alignment_star_command.py -v
```
Expected: 5 passed, 0 failed.

- [ ] **Step 3: Run the snakemake dry-run check for the star branch**

Run:
```bash
snakemake -n --use-conda --configfile config/test_config.yaml \
  --config samples_tsv=config/samples/patient_001_test.tsv alignment='{"aligner": "star", "star_index_dir": "/tmp/stub_star_idx"}' \
  -- results/patient_001_test/alignment/SRR9143066/junctions.tsv 2>&1 | tail -30
```
Expected: snakemake parses successfully, prints a dry-run plan including the star_align shell command. No `MissingInputException` is OK to ignore for the index_dir input — we're just verifying parse-time correctness of the Snakefile, not running the full DAG. Any *parse* errors must be fixed.

Alternative if the `--config` JSON form is finicky in the test shell: use the same fixture pattern as the pytest (write an override.yaml to /tmp and pass two configfiles).

- [ ] **Step 4: Commit the rule change**

Run:
```bash
git add workflow/rules/alignment.smk
git commit -m "$(cat <<'EOF'
feat(alignment): enable STAR 2-pass mode + remove silent multi-mapper throttle

Issue #17 — minimum-viable STAR config changes for production switch:
- Add --twopassMode Basic (novel junction sensitivity, the actual reason
  to switch from HISAT2)
- Add --limitSjdbInsertNsj 2000000 (defensive raise from default 1M;
  2-pass + multi-mappers can exceed and cause a fatal abort mid-run)
- Remove --outSJfilterReads Unique (silently drops multi-mapper-derived
  junctions vs Nature paper baseline)
- Remove --outSJfilterCountUniqueMin / --outSJfilterCountTotalMin
  (throttles output below paper baseline)

Paper-fidelity sensitivity flags (--outFilterMatchNminOverLread,
--alignIntronMax, --alignSJDBoverhangMin, etc.) deferred to a follow-up
Issue with per-flag benchmarks — avoids cargo-culting the recipe wholesale.

Test plan covered by workflow/tests/test_alignment_star_command.py (5 passing
assertions). End-to-end cloud verification deferred to Issue #378's full
patient_002 re-run.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands; tests pass.

---

## Task 3: Flip `config/config.yaml` aligner default

**Files:**
- Modify: `config/config.yaml` line ~47 (`alignment.aligner`)

- [ ] **Step 1: Edit config.yaml**

In [config/config.yaml](../../../config/config.yaml) under the `alignment:` block, change:

```yaml
alignment:
  aligner: "hisat2"                           # "star" (32 GB RAM) or "hisat2" (8 GB RAM)
```

to:

```yaml
alignment:
  aligner: "star"                             # "star" (32 GB RAM) or "hisat2" (8 GB RAM)
```

(Only the value changes; the inline comment lists both options and stays accurate.)

- [ ] **Step 2: Verify the production config parses with the new default**

Run:
```bash
snakemake -n --use-conda --configfile config/config.yaml \
  --config samples_tsv=config/samples/patient_001_test.tsv \
  -- results/patient_001_test/alignment/SRR9143066/junctions.tsv 2>&1 | tail -20
```
Expected: dry-run plan includes the `star_align` rule, not `hisat2_align`. Any parse errors must be fixed before commit. MissingInputException for genome FASTA / star_index is OK in this environment.

- [ ] **Step 3: Commit the config flip**

Run:
```bash
git add config/config.yaml
git commit -m "$(cat <<'EOF'
config(alignment): flip production default from hisat2 to star

Issue #17 — production aligner is now STAR by default. HISAT2 path
remains available via aligner: "hisat2" override (used in test_config.yaml
for the chr22 / M1 local dev environment, where STAR's 32 GB index won't
fit on 8 GB RAM).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands.

---

## Task 4: Pin `aligner: hisat2` in `config/test_config.yaml`

**Files:**
- Modify: `config/test_config.yaml` `alignment:` block

Without this pin, `snakemake --configfile config/test_config.yaml` on M1 would inherit the new `star` default and try to build a 32 GB STAR index on 8 GB RAM.

- [ ] **Step 1: Edit test_config.yaml**

In [config/test_config.yaml](../../../config/test_config.yaml), the `alignment:` block currently reads:

```yaml
alignment:
  hisat2_index_dir: "resources/test/hisat2_index"   # separate from production index
  hisat2_prebuilt_url: ""                            # build from chr22 FASTA (no pre-built available)
  star_index_dir: "resources/test/star_index"        # separate from production index
  threads: 4      # M1 MacBook Air (8 GB RAM) — leave headroom for other processes
```

Add an `aligner` key at the top of this block:

```yaml
alignment:
  aligner: "hisat2"                                  # explicit pin: M1 8 GB can't fit STAR's 32 GB index
  hisat2_index_dir: "resources/test/hisat2_index"   # separate from production index
  hisat2_prebuilt_url: ""                            # build from chr22 FASTA (no pre-built available)
  star_index_dir: "resources/test/star_index"        # separate from production index
  threads: 4      # M1 MacBook Air (8 GB RAM) — leave headroom for other processes
```

- [ ] **Step 2: Verify the chr22 test config parses + uses hisat2**

Run:
```bash
snakemake -n --use-conda --configfile config/test_config.yaml \
  -- results/patient_001_test/alignment/SRR9143066/junctions.tsv 2>&1 | tail -20
```
Expected: dry-run includes `hisat2_align`, NOT `star_align`. Confirms the override defeats inheritance.

- [ ] **Step 3: Commit the pin**

Run:
```bash
git add config/test_config.yaml
git commit -m "$(cat <<'EOF'
config(test): pin aligner=hisat2 to prevent star inheritance on M1

Issue #17 — test_config.yaml inherits from config.yaml via Snakemake's
recursive merge. After flipping production to STAR, the chr22 local dev
environment would inherit the new default and attempt to build a 32 GB
STAR index on the M1's 8 GB RAM. Explicit pin to hisat2 preserves the
local dev experience.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands.

---

## Task 5: Regenerate `dag.pdf`

**Files:**
- Modify: `dag.pdf` (currently unstaged-deleted on main; regenerate to satisfy the DAG-visualization rule for Snakemake rule changes)

The `dag.pdf` was unstaged-deleted on `main` at session start — likely a leftover from a previous session. The DAG-visualization memory rule says to regenerate after any Snakemake rule change before opening a PR. This task does both at once.

- [ ] **Step 1: Run the DAG visualization script**

Run:
```bash
bash scripts/visualize_dag.sh
```
Expected: a fresh `dag.pdf` is written to the repo root. The script may also produce intermediate `.dot` artifacts — leave those gitignored if they are.

If `scripts/visualize_dag.sh` fails (e.g. missing `graphviz` on macOS), install via `brew install graphviz` and retry. The script is mentioned in `memory/feedback_dag_visualization.md` as the canonical DAG render path.

- [ ] **Step 2: Diff vs the previous version**

Run:
```bash
git diff --stat dag.pdf
```
Expected: file present (resolves the unstaged deletion). The PDF binary diff may show byte-level changes even if structure is unchanged — that's fine. Structural change is not expected since we only modified shell-block content inside an existing rule, not added/removed rules or edges.

- [ ] **Step 3: Commit the regenerated DAG**

Run:
```bash
git add dag.pdf
git commit -m "$(cat <<'EOF'
chore(dag): regenerate dag.pdf after star_align rule update

Issue #17 — required by the DAG-visualization rule for any Snakemake
rule change before PR. Also resolves the stray unstaged deletion on main
from a previous session.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```
Expected: commit lands; `git status` clean.

---

## Task 6: Run the full pipeline pytest suite

**Files:** (no edits; verification only)

- [ ] **Step 1: Run the full pytest suite**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/ -v 2>&1 | tail -30
```
Expected: all tests pass. The new `test_alignment_star_command.py` reports 5 passed. No regression in existing tests (`test_strandness`, `test_star_sj_to_junctions`, `test_bed12_to_junctions`, etc.).

If any unrelated tests fail, investigate before opening the PR (per the "test before PR" memory).

- [ ] **Step 2: Run the snakemake dry-run CI check explicitly for both configs**

Run:
```bash
echo "=== Production (star) ===" && \
  snakemake -n --use-conda --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001_test.tsv \
    -- results/patient_001_test/alignment/SRR9143066/junctions.tsv 2>&1 | grep -E "(error|Error|ERROR)" | head -5 || echo "  OK"
echo "=== Test (hisat2) ===" && \
  snakemake -n --use-conda --configfile config/test_config.yaml \
    -- results/patient_001_test/alignment/SRR9143066/junctions.tsv 2>&1 | grep -E "(error|Error|ERROR)" | head -5 || echo "  OK"
```
Expected: both `OK` lines printed; no errors. Confirms both configs parse cleanly.

---

## Task 7: Push, open PR, update Issue #378 ACs, set PR Status

**Files:** (no file edits; GitHub state mutations)

- [ ] **Step 1: Push the feature branch**

Run:
```bash
git push -u origin feat/developer/issue-17-star-aligner
```
Expected: branch pushed, tracking origin.

- [ ] **Step 2: Open the PR**

Run:
```bash
gh pr create --base main --title "feat(alignment): switch production default to STAR with 2-pass + multi-mapper inclusion (Issue #17)" --project "JH M Lee Lab" --body "$(cat <<'EOF'
**Created by:** Developer

## Summary

Switches production aligner default from HISAT2 to STAR. Adds the two STAR
parameters that materially affect novel-junction sensitivity (`--twopassMode
Basic`, `--limitSjdbInsertNsj 2000000`) and removes three filters that
silently throttle multi-mapper-derived junctions (`--outSJfilterReads Unique`,
`--outSJfilterCountUniqueMin`, `--outSJfilterCountTotalMin`).

Paper-fidelity sensitivity tuning (`--outFilterMatchNminOverLread`,
`--alignIntronMax`, `--alignSJDBoverhangMin`, etc.) is **deferred** to a
follow-up Issue with per-flag benchmarks — avoids cargo-culting the Nature
paper command wholesale.

Closes [Issue #17](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17).

Design spec: [docs/superpowers/specs/2026-05-19-issue-17-star-aligner-design.md](docs/superpowers/specs/2026-05-19-issue-17-star-aligner-design.md).

## Test plan

- [x] New pytest `workflow/tests/test_alignment_star_command.py` (5 assertions) passes locally
- [x] Full pytest suite (`workflow/tests/.venv/bin/python -m pytest workflow/tests/`) green — no regressions
- [x] `snakemake -n` with `aligner: star` config parses cleanly
- [x] `snakemake -n` with `config/test_config.yaml` confirms hisat2 path still resolves (M1 8 GB stays on hisat2)
- [x] `dag.pdf` regenerated via `bash scripts/visualize_dag.sh`
- [x] Three new ACs added to [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) for the cloud verification (see Task 7 Step 3)

## Out of scope — deferred follow-ups

Paper-fidelity tuning Issue to be filed post-merge (per Task 7 Step 5).
EOF
)"
```
Expected: PR created, linked to Issue #17. Note the URL.

- [ ] **Step 3: Add three new ACs to Issue #378**

Run:
```bash
ORIG_BODY=$(gh issue view 378 --repo Jin-HoMLee/splice-neoepitope-pipeline --json body --jq .body)

NEW_BODY="$ORIG_BODY

---

## Additional ACs from [PR #](TBD) (Issue #17, STAR aligner switch)

The full patient_002 re-run inherits the new STAR defaults from PR #17. Verify in the run output:

- [ ] \`SJ.out.tab\` contains at least one row with \`multi_reads > 0 AND unique_reads == 0\` — confirms multi-mapper inclusion is wired through end-to-end
- [ ] 2-pass log line present in \`star_align.log\` (look for \`STAR: 2nd pass:\` or equivalent)
- [ ] No \`limitSjdbInsertNsj exceeded\` fatal error in stderr — confirms the 2M raise was sufficient

**Created by:** Developer (PR for Issue #17 wrap-up)
"

gh issue edit 378 --repo Jin-HoMLee/splice-neoepitope-pipeline --body "$NEW_BODY"
```

Then **replace `(TBD)` with the actual PR number** from Step 2:
```bash
PR_NUMBER=<from Step 2>
gh issue view 378 --repo Jin-HoMLee/splice-neoepitope-pipeline --json body --jq .body \
  | sed "s/PR #](TBD)/PR #${PR_NUMBER}](https:\\/\\/github.com\\/Jin-HoMLee\\/splice-neoepitope-pipeline\\/pull\\/${PR_NUMBER})/" \
  | gh issue edit 378 --repo Jin-HoMLee/splice-neoepitope-pipeline --body-file -
```
Expected: Issue #378 body now has three new unchecked ACs linked to this PR.

- [ ] **Step 4: Set PR Status to "Ready for review" on the board**

Run:
```bash
PR_NUMBER=<from Step 2>
PR_ITEM_ID=$(gh api graphql -f query='
{ user(login: "Jin-HoMLee") { projectV2(number: 9) { items(first: 100) {
  nodes { id content { ... on PullRequest { number } } }
} } } }' --jq ".data.user.projectV2.items.nodes[] | select(.content.number==$PR_NUMBER) | .id")

gh api graphql -f query='
mutation($item: ID!) {
  updateProjectV2ItemFieldValue(input: {
    projectId: "PVT_kwHOB17eGc4BSomP",
    itemId: $item,
    fieldId: "PVTSSF_lAHOB17eGc4BSomPzhAHFf8",
    value: { singleSelectOptionId: "8bf9192f" }
  }) { projectV2Item { id } }
}' -f item="$PR_ITEM_ID"
```
Expected: mutation returns the PR item ID. PR Status on the board is now "Ready for review".

- [ ] **Step 5: File the follow-up paper-fidelity tuning Issue**

Run:
```bash
gh issue create --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --project "JH M Lee Lab" \
  --label role:developer --label enhancement \
  --title "STAR sensitivity tuning — benchmark Nature paper command flags individually" \
  --body "$(cat <<'EOF'
**Created by:** Developer (follow-up to Issue #17)

## Motivation

PR for [Issue #17](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17) (STAR aligner switch) ships the minimum-viable parameter changes — only those that are correctness issues (silent multi-mapper throttle removed, 2-pass mode added). Several sensitivity flags from the Nature paper command were deliberately deferred to avoid cargo-culting the recipe wholesale without benchmarks.

This Issue tracks per-flag benchmarks for the deferred flags. Each flag should be evaluated individually on a patient_002 STAR run, measured against:
- Total `SJ.out.tab` junction count
- % annotated (vs novel)
- Runtime delta
- Multi-mapper utilization

Adopt a flag if its delta on novel-junction recall justifies its complexity; reject otherwise.

## Flags to evaluate

- [ ] `--outFilterMatchNminOverLread 0.33` (vs default 0.66) — more permissive match-length fraction
- [ ] `--outFilterScoreMinOverLread 0.33` (vs default 0.66) — more permissive score fraction
- [ ] `--outFilterMultimapNmax 20` (vs default 10) — raise multi-mapper count cap
- [ ] `--alignIntronMax 500000` (vs STAR default) — intron length cap
- [ ] `--alignMatesGapMax 1000000` (vs STAR default)
- [ ] `--alignSJDBoverhangMin 1` (vs default 3) — permissive overhang for novel junctions
- [ ] `--sjdbGTFfile` removal from `star_index` (vs current GTF-in-index)
- [ ] BAM emission via `--outSAMtype BAM SortedByCoordinate` + SAM-formatting flags (`--outSAMattributes NH HI NM MD AS XS`, `--outSAMstrandField intronMotif`, `--outSAMunmapped None`, `--outSAMmultNmax 1`) — required only if a downstream rule needs the BAM

## Acceptance criteria

- [ ] Each flag benchmarked on patient_002 STAR run; numbers recorded in `research/lab_notebook/developer.md`
- [ ] PR(s) adopt flags whose delta is meaningful; reject the rest with one-line justification
EOF
)"
```
Expected: new issue created, added to project board #9.

---

## Self-review (do these before reporting plan complete)

- [ ] **Spec coverage:** Each change in [the spec's "Changes — four files" section](../specs/2026-05-19-issue-17-star-aligner-design.md) is implemented by a task here. ✓ (config.yaml = Task 3, test_config.yaml = Task 4, alignment.smk = Task 2, new pytest = Task 1)
- [ ] **Placeholder scan:** No "TODO", "TBD" in step bodies. (`(TBD)` in the PR # for Issue #378's AC update is a deliberate sentinel that Step 3's follow-up sed replaces.)
- [ ] **Type consistency:** Test function names match across the test file definition and the assertion descriptions. ✓
- [ ] **Verification ordering:** TDD red→green sequence (Task 1 writes failing tests; Task 2 makes them pass). ✓
- [ ] **Branch + board hygiene:** Task 0 sets Issue Status In Progress; Task 7 Step 4 sets PR Status Ready for review. ✓
- [ ] **Closure ritual:** PR body Test plan + Issue #378 ACs both use `- [x]`/`- [ ]` boxes that `scripts/audit_and_merge.sh` will check at merge time. ✓
- [ ] **Lab notebook coda:** Lab notebook entry happens between review and merge (post-review-feedback, pre-`audit_and_merge.sh`) per the shared rule — not pre-commit. Not a task in this plan; will surface at merge time.

---

## Post-merge follow-ups (not in this plan's scope)

1. Write `research/lab_notebook/developer.md` entry **after** any review feedback is incorporated, **before** `bash scripts/audit_and_merge.sh <PR#>` is run.
2. Merge via `bash scripts/audit_and_merge.sh <PR#>` — the closure-ritual gate audits PR body + Issue #17 ACs (note Issue #17 currently has none; the gate only enforces what exists).
3. Wait for [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) full patient_002 cloud re-run to verify the three new ACs end-to-end.
