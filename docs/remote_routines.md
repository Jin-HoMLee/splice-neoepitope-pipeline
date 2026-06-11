# Remote Routines (scheduled cloud agents)

> **Status: established 2026-06-03.** Validated by the overnight one-shot batch — #632 pilot ([PR #647](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/647)) + env-probe smoke test; #641 ([PR #648](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/648)) confirmed the corrected `pip install`, the stray-closer guard, and the env-ferry, and surfaced the **network egress allowlist**; #375 ([PR #649](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/649)) validated the heavier `requirements-test.txt` install (pandas/biopython/numpy) **and** the conda-absent → **`snakemake -n` deferred-to-CI** pattern (#647/#648/#649 all CI-green incl. `pipeline-snakemake-dry-run`); #435 ([PR #650](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/650)) was the 4th. Re-probe the sandbox table if the platform image changes.

Remote routines are scheduled or one-time Claude Code agents that run in an **isolated Anthropic-cloud sandbox** (CCR) — *not* on a local machine and *not* an OS cron. They are created via the `/schedule` skill (the RemoteTrigger API) or at <https://claude.ai/code/routines>. Two uses:

- **Recurring hygiene sweeps** — board/milestone/PR/dependency drift the morning routine only catches when a human runs a session (see the routine-candidate analysis).
- **Overnight one-shot work-dispatch** — a bounded, machine-verifiable Issue → a **draft PR you review at breakfast**.

## CCR sandbox environment — reference facts

Empirically probed in the `github-pat` environment (`env_01FffKZvNhx6EVwk7nEqMSkM`) on 2026-06-03. Re-probe if the platform image changes.

| Fact | Value | Implication for a routine |
|---|---|---|
| OS / CPU | Linux 6.18.5 x86_64, 4 cores | fine for tests/dry-runs; not for heavy compute |
| Python | 3.11.15 (`/usr/local/bin`) | |
| `gh` | 2.65.0, authed via `$GH_TOKEN` — fine-grained PAT, **repo + project** scope | can push, open PRs, comment, **and** move project-board fields |
| PyPI egress | ✅ works | `pip install` reaches pypi.org |
| **network egress** | **allowlisted** — PyPI + GitHub reachable; other hosts (e.g. `api.datacite.org`) return `403 Host not in allowlist` | an in-sandbox live API-shape probe works **only** for GitHub + PyPI; a mapper/fixture that depends on another external API can't be confirmed in-sandbox → build from the documented schema and **defer live confirmation to human/CI** (leave that AC unchecked). Surfaced by #641/[PR #648](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/648). |
| **pip pyyaml trap** | `pip install --upgrade … pyyaml` **fails** — the debian-managed PyYAML 6.0.1 can't be uninstalled, which **aborts the whole transaction → pytest never installs** | install with **`python3 -m pip install pytest`** (no `--upgrade`, don't touch pyyaml — it's already present). Capture `pip`'s rc directly (`rc=$?` with no pipe; `\| tail` swallows it). |
| conda / mamba / snakemake | ❌ absent (`command not found`) | a `snakemake -n` gate **cannot run in-sandbox** → defer the dry-run to **CI on the PR**; pytest is the in-sandbox gate |
| git clone | **shallow**, detached HEAD, served through a `local_proxy` git proxy | the proxy **403s any push to `main`** → push only to **`claude/*`** branches |
| git identity | preset to `Claude <noreply@anthropic.com>` | commit author is already Claude (so the `Co-Authored-By: Claude …` trailer is technically redundant — **omit the trailer**) |
| local files | a fresh **single-repo** checkout — **no** per-clone venvs (`workflow/tests/.venv`, `research/.venv`), **no** gitignored data (`indices/`, `references/`, `*_outputs/`), **no** GCP VM, **no** secrets (`ZOTERO_*`, GCS) | VM-bound / data-bound / venv-bound work is out of scope; rebuild only light deps via pip |
| human-in-loop | **none** — any interactive permission prompt **hangs the run forever** | never write under `.claude/` (the Claude Code ≥2.1.78 protected-dir gate prompts → hang); broaden/enumerate `allowed_tools` so no command prompts |
| observability | the run log / stdout lives **only** in the gated claude.ai session UI — unreachable by `gh`/API/WebFetch | **ferry** any facts you need out to a durable channel (a `## Environment` PR-body section or a `gh pr comment`) *during* the run |

## Hardened one-shot dispatch — prompt checklist

The sandbox boots with **zero memory and no skills**, so the prompt must encode every convention explicitly:

1. **Env self-report first** → `tee /tmp/run.log`. Install the gate deps with `python3 -m pip install pytest` (capture `pip-rc` with no pipe). If egress fails, **STOP and post the log to the Issue** — never ship an unverified PR.
2. **Read the Issue + target files in full** (no assumed context). Restate conventions (commit style, draft-only, no bare `#N`).
3. **One live data-shape probe** before fixturing if the change depends on an external API/endpoint shape — don't mock a fictional schema. **But only GitHub + PyPI hosts are reachable** (the egress allowlist 403s others). If the host is blocked, build from the documented schema and **defer live confirmation to human/CI** (leave that AC unchecked) — never claim a probe you couldn't run.
4. **TDD red→green**; paste the **real** green gate output into the PR body (never claim green without it).
5. **Branch `claude/issue-NNN-…` off `origin/main`** (plain `git switch -c`, *not* `gh issue develop` — untested through the sandbox's proxied shallow checkout, and the local parent-guard hook that protects `gh issue develop` doesn't load in-sandbox); open a **draft** PR; **never merge, never force-push**.
6. **`Closes #NNN` in the PR body** (the only link to the Issue for a plain-branch PR's merge gate) — then a **stray-closer final scan**: neutralize any *other* `close/fix/resolve`+`#N`, including incidental prose like `closed #192`, and verify `closingIssuesReferences == {NNN}`.
7. **Ferry the env/run log out** into a `## Environment` PR-body section (and/or a PR comment) — the only durable record.
8. **Leave creds/network-gated ACs unchecked** for human ratify; **don't write the lab-notebook entry** (the human writes it at merge); never the literal `@`+`claude` mention (use `@-claude`).
9. **`allowed_tools`:** for a first/uncertain task, list the **`Bash` tool with no argument-level filter** + file tools (e.g. `["Bash", "Read", "Edit", "Write", "Glob", "Grep"]`) — an un-allow-listed *command* otherwise prompts → hangs the unattended run.

## When to dispatch a one-shot

**Good fit** = *specifiable* (a complete brief) ∧ *machine-verifiable* (a gate the agent can run: pytest, `snakemake -n` deferred to CI, lint) ∧ *reviewable-as-artifact* (a draft PR you ratify).
**Bad fit** = needs the GCP VM / live GCS / local data / a human judgment call on direction / ambiguous scope.

## Handoff — routine-drafted PRs go to the owning role

A routine-drafted draft PR is reviewed, lab-notebooked, and merged by the **role persona named on its Issue's `role:` label** — *not* by the PM who dispatched it, and *not* by another routine. Dispatch deliberately leaves three human-in-loop steps the unattended run must never do:

- the **lab-notebook entry** (`research/lab_notebook/<role>.md`) — reflective judgment; enforced at merge by `scripts/audit_and_merge.sh` (lab-notebook gate);
- any **creds/network-gated ratification** the sandbox allowlist blocked — e.g. a live external-API `--dry-run` (the routine builds from the documented schema and leaves that AC unchecked);
- the **merge** itself (the closure-ritual gate).

Routing is mostly automatic: each PR carries its `role:` label, sits at *Ready for review* on board #9, and its body documents its own residuals. But because a `Claude`-authored draft at Ready-for-review can read as *"did I make this?"* to the picking-up persona, the dispatcher posts a brief **handoff comment** on each PR naming the owner role + the pre-merge residuals. Coordinate via the board/PR, **not** a standup message (per the [Issue #569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569) work-vs-message split).

**In short: PM dispatches + routes; the owning role finishes.** Always delegate to the role *persona* (a human-in-loop session), never to another unattended routine — those three residuals are exactly the steps kept human by design. (Optional future step: a routine could *draft* the lab-notebook entry for the persona to ratify, mirroring how the code PR is itself a draft — but the merge gate stays human.)

## Provenance

- #632 pilot (the first overnight code dispatch) → [PR #647](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/647); env-probe smoke test → the `## Environment` comment on that PR.
- In-session-vs-one-shot difference analysis and the routine-candidate board scan: PM session 2026-06-03.
