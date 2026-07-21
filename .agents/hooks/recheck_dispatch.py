#!/usr/bin/env python3
"""PostToolUse hook: trigger milestone-capacity AND parent-status rechecks on capacity-change events.

Watches Bash commands for trigger shapes:
  1. gh issue close N                              -> recheck issue's milestone + parent-status chain
  2. gh issue edit N ... --milestone X             -> recheck old + new milestones (via /events history)
  3. project board Size mutation (Size field ID)   -> recheck affected issue's milestone
  4. gh api .../milestones/N PATCH (due_on edit)   -> recheck milestone N
  5. project board Status mutation (Status field ID) -> recheck affected issue's parent-status chain

For each match, invokes the relevant script and emits its output as
`additionalContext` so it appears in the next prompt.

Scope filter (Issue #454): each check declares a "scope" in HOOK_CONFIG and
runs only when that scope is active for this invocation (`--scope shared|pm|all`).
Committed `.agents/settings.json` invokes `--scope shared`; the PM-local
`.agents/settings.local.json` invokes `--scope pm`. So Scientist/Developer
sessions get only the shared, broadly-actionable checks, while PM-only checks
(e.g. capacity rechecks) stay local. A flagless invocation runs everything.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import graphql_meter  # noqa: E402

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
SIZE_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHGiA"
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
PROJECT_ID = "PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER = 9
TARGET_DATE_FIELD_ID = "PVTF_lAHOB17eGc4BSomPzhAHGiM"
TARGET_DATE_FIELD_NAME = "Target date"
SCRIPT = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "pm" / "recheck_milestone.py")
PARENT_STATUS_SCRIPT = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "pm" / "recheck_parent_status.py")

# ---------------------------------------------------------------------------
# Fire-log infrastructure (Issue #453)
# ---------------------------------------------------------------------------

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"


def _log_fire(hook_name: str, issue: int | None = None, **metadata) -> None:
    """Append one JSONL line to LOG_PATH recording a hook fire.

    POSIX guarantees each write() syscall on an O_APPEND regular file is
    atomic (no concurrent-writer interleave) regardless of length. The 4 KB
    guard is sanity hygiene — if a future schema change pushes a line above
    ~4 KB, the write is skipped rather than persisted, so the log structure
    stays parseable.
    """
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": hook_name,
        "issue": issue,
        **metadata,
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    if len(line.encode("utf-8")) >= 4096:
        print(f"hook_fires: oversized line (>=4KB) for {hook_name}, skipped", file=sys.stderr)
        return
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line)


def _count_fires(hook_name: str) -> int:
    """Count lines in LOG_PATH where hook_name appears in the JSONL record.

    Substring match against the serialized form (`"hook":"<name>"`); cheap
    at K=3 scale, no JSON parsing needed.
    """
    if not LOG_PATH.exists():
        return 0
    needle = f'"hook":"{hook_name}"'
    return sum(1 for line in LOG_PATH.read_text(encoding="utf-8").splitlines() if needle in line)


HOOK_CONFIG: dict[str, dict] = {
    # scope: "shared" → runs in all role sessions (committed settings.json);
    #        "pm"     → PM-local only (settings.local.json).
    "target_sync_check":     {"threshold": 3, "dock": 454, "scope": "shared"},   # re-sync actionable by whoever moved the milestone (any role)
    "recheck_milestone":     {"threshold": 3, "dock": 618, "scope": "pm"},      # PM-local by design (dock #618); capacity rebalance is PM-coordinated + fires on every close → noise for Sci/Dev
    "recheck_parent_status": {"threshold": 3, "dock": 617, "scope": "pm"},      # keep PM-local: #617 proves-out concluded PM-local (A2 mooted the shared-promotion); its fires measured the child→parent mirror that #794 retired, so they can't justify promoting the narrowed residual — reassess scope against the residual once fresh fires accrue
}


# ---------------------------------------------------------------------------
# Scope filter (Issue #454)
# ---------------------------------------------------------------------------
ACTIVE_SCOPES: set[str] = {"shared", "pm"}  # default: run everything (flagless / --scope all)


def _parse_scope(argv: list[str]) -> set[str]:
    """Map a --scope {shared,pm,all} CLI arg to the active-scope set. Unknown/absent → all (fail-open)."""
    scope = "all"
    for i, arg in enumerate(argv):
        if arg == "--scope" and i + 1 < len(argv):
            scope = argv[i + 1]
        elif arg.startswith("--scope="):
            scope = arg.split("=", 1)[1]
    if scope == "shared":
        return {"shared"}
    if scope == "pm":
        return {"pm"}
    return {"shared", "pm"}


def _in_scope(hook_name: str) -> bool:
    """True if hook_name's scope is active this invocation. No scope key → always (fail-open)."""
    cfg = HOOK_CONFIG.get(hook_name)
    scope = cfg.get("scope") if cfg else None
    return scope is None or scope in ACTIVE_SCOPES


def _threshold_prompt(hook_name: str) -> str | None:
    """Return the 🎯 promotion-review prompt when count == K (exactly); else None.

    Equals-once (not >=) so the prompt fires on exactly the K-th fire and
    never again. No nagging.
    """
    cfg = HOOK_CONFIG.get(hook_name)
    if cfg is None:
        return None
    if _count_fires(hook_name) != cfg["threshold"]:
        return None
    dock = f"Issue #{cfg['dock']}" if cfg["dock"] else "(no dock Issue filed yet)"
    return (
        f"🎯 {hook_name} has now fired {cfg['threshold']} times — review for promotion: "
        f"bash scripts/check_hook_health.sh ; {dock}"
    )


# Per-hook no-fire sentinels. Each subprocess script may emit one of several
# strings to signal "nothing actionable here"; if any sentinel is present in
# output, the predicate returns False (no fire). Sentinels are verified
# against the live script source in scripts/pm/.
_NO_FIRE_SENTINELS = {
    "recheck_milestone": (
        "Status: [No change]",       # no drift OR no remaining capacity
        "has no milestone",          # issue lost its milestone (stderr)
        # Post-move listing not yet converged (Issue #406). Surfaced to the PM as
        # context, but it is NOT an actionable capacity drift (a "come back
        # later" note), so it must not inflate the promotion fire-log.
        "[stale state, verification pending]",
    ),
    "recheck_parent_status": (
        "Status: [No change]",       # no drift detected on parent chain
        "has no parent",             # leaf issue (the common close case)
    ),
}


def _is_fire(hook_name: str, output: str) -> bool:
    """Return True if this output represents a real warning (not 'no change' / not 'no action needed')."""
    if hook_name == "target_sync_check":
        # Only the manual-param fallback ("[target re-sync needed …]") is a real fire — it
        # means the auto-mutation failed and a human must finish the sync. A successful
        # auto-sync confirmation ("[target auto-synced …]", Route A #782) is the mechanism
        # working as intended, NOT a promotion signal, so it must not inflate the fire-log.
        return bool(output) and not output.startswith("[target auto-synced")
    sentinels = _NO_FIRE_SENTINELS.get(hook_name)
    if sentinels is None:
        return False
    return not any(s in output for s in sentinels)


def _wrap_warning(hook_name: str, issue: int | None, warning: str, outputs: list[str]) -> None:
    """Append warning to outputs; if it's a real fire, log it + append threshold prompt."""
    if not _in_scope(hook_name):
        return  # out of scope — suppress output (the dispatch() guard prevents the subprocess before reaching here; this is a defensive backstop)
    outputs.append(warning)  # always emit so user sees recheck context
    if not _is_fire(hook_name, warning):
        return
    _log_fire(hook_name, issue=issue)
    prompt = _threshold_prompt(hook_name)
    if prompt:
        outputs.append(prompt)


# ---------------------------------------------------------------------------
PATTERN_CLOSE = re.compile(r"\bgh\s+issue\s+close\s+(\d+)")
PATTERN_MOVE = re.compile(r"\bgh\s+issue\s+edit\s+(\d+)\b[^|;&]*--milestone\b")
PATTERN_GH_API = re.compile(r"\bgh\s+api\b")
PATTERN_MILESTONE_PATH = re.compile(r"/milestones/(\d+)\b")
PATTERN_PATCH_METHOD = re.compile(r"-X\s+PATCH|--method\s+PATCH|method=PATCH")
PATTERN_ITEMID = re.compile(r'itemId:\s*"(PVTI_[A-Za-z0-9_-]+)"')


def _run_script(script: str, label: str, *args: str) -> str:
    if not Path(script).is_file():
        return f"({label} error: script not found at {script})"
    try:
        result = subprocess.run(
            ["python3", script, *args],
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        return f"({label} error: {exc})"
    out = result.stdout
    if result.stderr:
        out += result.stderr
    return out


def run_recheck(*args: str) -> str:
    return _run_script(SCRIPT, "recheck", *args)


def run_parent_status_recheck(*args: str) -> str:
    return _run_script(PARENT_STATUS_SCRIPT, "parent-status recheck", *args)


def lookup_milestone_number_by_title(title: str) -> int | None:
    # per_page=100 covers the current ~17 milestones with headroom. If this repo
    # ever exceeds 100 (closed + open), switch to --paginate.
    result = subprocess.run(
        ["gh", "api", f"repos/{REPO}/milestones?state=all&per_page=100"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    for m in data:
        if m.get("title") == title:
            return m.get("number")
    return None


def lookup_issue_for_item(item_id: str) -> int | None:
    query = f'query {{ rateLimit {{ cost remaining }} node(id: "{item_id}") {{ ... on ProjectV2Item {{ content {{ ... on Issue {{ number }} }} }} }} }}'
    result = subprocess.run(
        ["gh", "api", "graphql", "-f", f"query={query}"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    graphql_meter.log_graphql_spend("recheck_dispatch", data, query_name="lookup_issue_for_item")
    node = (data.get("data") or {}).get("node") or {}
    num = (node.get("content") or {}).get("number")
    return int(num) if isinstance(num, int) else None


def get_issue_milestone(issue: int) -> tuple[str | None, str | None]:
    """Return (milestone_title, due_on YYYY-MM-DD) for the issue's current milestone."""
    result = subprocess.run(
        ["gh", "api", f"repos/{REPO}/issues/{issue}",
         "--jq", "{title: .milestone.title, due_on: .milestone.due_on}"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return (None, None)
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return (None, None)
    title = data.get("title")
    due_on = data.get("due_on")
    if due_on:
        due_on = due_on[:10]
    return (title, due_on)


def get_issue_target_date(issue: int) -> tuple[str | None, str | None]:
    """Return (target_date YYYY-MM-DD, project_item_id) for the issue's entry on project #9.

    Returns (None, None) if the issue is not on the project board. Returns
    (None, item_id) if the issue is on the board but Target date is unset.
    """
    query = (
        'query { rateLimit { cost remaining } repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") '
        '{ issue(number: ' + str(issue) + ') { projectItems(first: 10) { nodes { '
        'id project { number } fieldValues(first: 20) { nodes { '
        '... on ProjectV2ItemFieldDateValue { date field { ... on ProjectV2FieldCommon { name } } } '
        '} } } } } } }'
    )
    result = subprocess.run(
        ["gh", "api", "graphql", "-f", f"query={query}"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return (None, None)
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return (None, None)
    graphql_meter.log_graphql_spend("recheck_dispatch", data, query_name="get_issue_target_date")
    issue_node = ((data.get("data") or {}).get("repository") or {}).get("issue") or {}
    items = (issue_node.get("projectItems") or {}).get("nodes") or []
    for item in items:
        if ((item.get("project") or {}).get("number")) != PROJECT_NUMBER:
            continue
        item_id = item.get("id")
        for fv in (item.get("fieldValues") or {}).get("nodes") or []:
            if not fv:
                continue
            field_name = (fv.get("field") or {}).get("name")
            if field_name == TARGET_DATE_FIELD_NAME:
                return (fv.get("date"), item_id)
        return (None, item_id)
    return (None, None)


def target_sync_check(issue: int) -> str | None:
    """Compare issue's Target date against current milestone due_on. Returns warning or None."""
    ms_title, ms_due_on = get_issue_milestone(issue)
    target_date, item_id = get_issue_target_date(issue)

    if item_id is None:
        return None

    if ms_title:
        if target_date == ms_due_on:
            return None
        return (
            f"[target re-sync needed — move on #{issue}, new milestone \"{ms_title}\"]\n"
            f"Project board Target date: {target_date or '(unset)'}\n"
            f"New milestone due_on:      {ms_due_on or '(unset)'}\n"
            f"Re-sync via GraphQL `updateProjectV2ItemFieldValue`:\n"
            f"  projectId: {PROJECT_ID}\n"
            f"  itemId:    {item_id}\n"
            f"  fieldId:   {TARGET_DATE_FIELD_ID}  ({TARGET_DATE_FIELD_NAME})\n"
            f"  date:      {ms_due_on}"
        )

    if target_date is not None:
        return (
            f"[target re-sync needed — demilestone on #{issue}]\n"
            f"Issue is un-milestoned but Target date is still: {target_date}\n"
            f"Clear via GraphQL `clearProjectV2ItemFieldValue`:\n"
            f"  projectId: {PROJECT_ID}\n"
            f"  itemId:    {item_id}\n"
            f"  fieldId:   {TARGET_DATE_FIELD_ID}  ({TARGET_DATE_FIELD_NAME})"
        )

    return None


def _mutate_target_date(item_id: str, due_on: str | None) -> bool:
    """Execute the Target-date mutation on project #9. `due_on=None` clears it.

    Returns True on success, False on ANY failure (non-zero gh exit, GraphQL
    `errors` payload, timeout, or `gh` missing). The False path is the fail-open
    trigger: the caller (apply_target_sync) falls back to surfacing the manual
    mutation params, so a token/scope/network problem never loses the sync.
    """
    if not re.match(r"PVTI_[A-Za-z0-9_-]+$", item_id):
        return False  # malformed id → fail-open rather than build a malformed mutation string
    if due_on:
        query = (
            'mutation { updateProjectV2ItemFieldValue(input:{'
            f'projectId:"{PROJECT_ID}", itemId:"{item_id}", fieldId:"{TARGET_DATE_FIELD_ID}", '
            f'value:{{date:"{due_on}"}}'
            '}){ projectV2Item{ id } } }'
        )
    else:
        query = (
            'mutation { clearProjectV2ItemFieldValue(input:{'
            f'projectId:"{PROJECT_ID}", itemId:"{item_id}", fieldId:"{TARGET_DATE_FIELD_ID}"'
            '}){ projectV2Item{ id } } }'
        )
    try:
        result = subprocess.run(
            ["gh", "api", "graphql", "-f", f"query={query}"],
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False
    if result.returncode != 0:
        return False
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return False
    try:
        probe = subprocess.run(
            ["gh", "api", "graphql", "-f", f"query={graphql_meter.RATE_LIMIT_PROBE_QUERY}"],
            capture_output=True, text=True, check=False,
        )
        graphql_meter.log_graphql_probe("recheck_dispatch", json.loads(probe.stdout), query_name="mutate_target_date")
    except Exception:
        pass
    return "errors" not in data


def apply_target_sync(issue: int) -> str | None:
    """Auto-apply the Target-date sync for `issue` (Route A, Issue #782).

    Executes the deterministic derivation Target = current-milestone `due_on`
    (or CLEARS Target when the issue is un-milestoned), so the previously-manual
    follow-up step can't slip. Returns:
      - a confirmation line on a successful mutation (surfaced so the move is visible);
      - the manual-param warning from target_sync_check() on mutation failure (fail-open);
      - None when no sync is needed (already in sync / not on the board).

    Only this deterministic Target derivation is auto-applied. The capacity
    recheck stays advisory (a PM judgment surface) and is never auto-mutated.
    """
    ms_title, ms_due_on = get_issue_milestone(issue)
    target_date, item_id = get_issue_target_date(issue)

    if item_id is None:
        return None  # not on the project board — nothing to sync

    if ms_title:
        if target_date == ms_due_on:
            return None  # already in sync — idempotent no-op
        if _mutate_target_date(item_id, ms_due_on):
            shown = ms_due_on or "(cleared — milestone has no due date)"
            return f'[target auto-synced — #{issue} → milestone "{ms_title}", Target {shown}]'
        return target_sync_check(issue)  # fail-open: surface manual params, unchanged

    if target_date is not None:
        # un-milestoned but Target still set → clear it
        if _mutate_target_date(item_id, None):
            return f"[target auto-synced — #{issue} demilestoned, Target cleared]"
        return target_sync_check(issue)  # fail-open

    return None


def prior_milestones_for_issue(issue: int) -> list[int]:
    """Return milestone numbers from the 2 most recent milestoned/demilestoned events."""
    result = subprocess.run(
        ["gh", "api", f"repos/{REPO}/issues/{issue}/events?per_page=100",
         "--jq",
         'map(select(.event == "milestoned" or .event == "demilestoned")) '
         '| sort_by(.created_at) | reverse | .[0:2] | .[].milestone.title'],
        capture_output=True, text=True, check=False,
    )
    titles = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    numbers: list[int] = []
    for title in titles:
        ms = lookup_milestone_number_by_title(title)
        if ms is not None and ms not in numbers:
            numbers.append(ms)
    return numbers


def dispatch(cmd: str) -> list[str]:
    outputs: list[str] = []

    m = PATTERN_CLOSE.search(cmd)
    if m:
        issue_n = int(m.group(1))
        if _in_scope("recheck_milestone"):
            _wrap_warning(
                "recheck_milestone",
                issue_n,
                f"[milestone recheck — close #{issue_n}]\n{run_recheck('--issue', str(issue_n))}",
                outputs,
            )
        if _in_scope("recheck_parent_status"):
            _wrap_warning(
                "recheck_parent_status",
                issue_n,
                f"[parent-status recheck — close #{issue_n}]\n"
                f"{run_parent_status_recheck('--issue', str(issue_n))}",
                outputs,
            )

    m = PATTERN_MOVE.search(cmd)
    if m:
        issue = int(m.group(1))
        if _in_scope("recheck_milestone"):
            ms_numbers = prior_milestones_for_issue(issue)
            if not ms_numbers:
                _wrap_warning(
                    "recheck_milestone",
                    issue,
                    f"[milestone recheck — move on #{issue} "
                    f"(milestone history empty; rechecking current only)]\n"
                    f"{run_recheck('--issue', str(issue), '--moved-issue', str(issue))}",
                    outputs,
                )
            else:
                # --moved-issue enables post-move eventual-consistency
                # reconciliation of the laggy listing endpoint (Issue #406):
                # each source/dest recheck confirms the move propagated before
                # recomputing capacity, instead of a false [UPDATE NEEDED].
                for ms in ms_numbers:
                    _wrap_warning(
                        "recheck_milestone",
                        issue,
                        f"[milestone recheck — move on #{issue}, milestone {ms}]\n"
                        f"{run_recheck('--milestone', str(ms), '--moved-issue', str(issue))}",
                        outputs,
                    )
    # Target-date auto-sync (Route A, Issue #782): a dedicated finditer loop so
    # EVERY issue in a batched `gh issue edit N --milestone … && …` is synced,
    # not just the first match (the batch-miss that triggered #782). The capacity
    # recheck above intentionally stays first-match + advisory; only this
    # deterministic Target derivation auto-applies.
    if _in_scope("target_sync_check"):
        seen_for_sync: set[int] = set()
        for mv in PATTERN_MOVE.finditer(cmd):
            issue = int(mv.group(1))
            if issue in seen_for_sync:
                continue
            seen_for_sync.add(issue)
            sync_warn = apply_target_sync(issue)
            if sync_warn:
                _wrap_warning("target_sync_check", issue, sync_warn, outputs)

    if SIZE_FIELD_ID in cmd and _in_scope("recheck_milestone"):
        item_match = PATTERN_ITEMID.search(cmd)
        if item_match:
            issue = lookup_issue_for_item(item_match.group(1))
            if issue:
                _wrap_warning(
                    "recheck_milestone",
                    issue,
                    f"[milestone recheck — size change on #{issue}]\n"
                    f"{run_recheck('--issue', str(issue))}",
                    outputs,
                )

    if STATUS_FIELD_ID in cmd and _in_scope("recheck_parent_status"):
        item_match = PATTERN_ITEMID.search(cmd)
        if item_match:
            issue = lookup_issue_for_item(item_match.group(1))
            if issue:
                _wrap_warning(
                    "recheck_parent_status",
                    issue,
                    f"[parent-status recheck — Status change on #{issue}]\n"
                    f"{run_parent_status_recheck('--issue', str(issue))}",
                    outputs,
                )

    if PATTERN_GH_API.search(cmd) and PATTERN_PATCH_METHOD.search(cmd) and _in_scope("recheck_milestone"):
        m = PATTERN_MILESTONE_PATH.search(cmd)
        if m:
            _wrap_warning(
                "recheck_milestone",
                None,
                f"[milestone recheck — due_on PATCH on milestone {m.group(1)}]\n"
                f"{run_recheck('--milestone', m.group(1))}",
                outputs,
            )

    return outputs


def main() -> int:
    global ACTIVE_SCOPES
    ACTIVE_SCOPES = _parse_scope(sys.argv[1:])
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (payload.get("tool_input") or {}).get("command", "")
    if not cmd:
        return 0

    outputs = dispatch(cmd)
    if not outputs:
        return 0

    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": "\n\n".join(outputs),
        }
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
