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
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
SIZE_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHGiA"
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
PROJECT_ID = "PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER = 9
TARGET_DATE_FIELD_ID = "PVTF_lAHOB17eGc4BSomPzhAHGiM"
TARGET_DATE_FIELD_NAME = "Target date"
SCRIPT = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "pm" / "recheck_milestone.py")
PARENT_STATUS_SCRIPT = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "pm" / "recheck_parent_status.py")

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
    query = f'query {{ node(id: "{item_id}") {{ ... on ProjectV2Item {{ content {{ ... on Issue {{ number }} }} }} }} }}'
    result = subprocess.run(
        ["gh", "api", "graphql", "-f", f"query={query}",
         "--jq", ".data.node.content.number"],
        capture_output=True, text=True, check=False,
    )
    out = result.stdout.strip()
    return int(out) if out.isdigit() else None


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
        'query { repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") '
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
        outputs.append(
            f"[milestone recheck — close #{m.group(1)}]\n{run_recheck('--issue', m.group(1))}"
        )
        outputs.append(
            f"[parent-status recheck — close #{m.group(1)}]\n"
            f"{run_parent_status_recheck('--issue', m.group(1))}"
        )

    m = PATTERN_MOVE.search(cmd)
    if m:
        issue = int(m.group(1))
        ms_numbers = prior_milestones_for_issue(issue)
        if not ms_numbers:
            outputs.append(
                f"[milestone recheck — move on #{issue} "
                f"(milestone history empty; rechecking current only)]\n"
                f"{run_recheck('--issue', str(issue))}"
            )
        else:
            for ms in ms_numbers:
                outputs.append(
                    f"[milestone recheck — move on #{issue}, milestone {ms}]\n"
                    f"{run_recheck('--milestone', str(ms))}"
                )
        sync_warn = target_sync_check(issue)
        if sync_warn:
            outputs.append(sync_warn)

    if SIZE_FIELD_ID in cmd:
        item_match = PATTERN_ITEMID.search(cmd)
        if item_match:
            issue = lookup_issue_for_item(item_match.group(1))
            if issue:
                outputs.append(
                    f"[milestone recheck — size change on #{issue}]\n"
                    f"{run_recheck('--issue', str(issue))}"
                )

    if STATUS_FIELD_ID in cmd:
        item_match = PATTERN_ITEMID.search(cmd)
        if item_match:
            issue = lookup_issue_for_item(item_match.group(1))
            if issue:
                outputs.append(
                    f"[parent-status recheck — Status change on #{issue}]\n"
                    f"{run_parent_status_recheck('--issue', str(issue))}"
                )

    if PATTERN_GH_API.search(cmd) and PATTERN_PATCH_METHOD.search(cmd):
        m = PATTERN_MILESTONE_PATH.search(cmd)
        if m:
            outputs.append(
                f"[milestone recheck — due_on PATCH on milestone {m.group(1)}]\n"
                f"{run_recheck('--milestone', m.group(1))}"
            )

    return outputs


def main() -> int:
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
