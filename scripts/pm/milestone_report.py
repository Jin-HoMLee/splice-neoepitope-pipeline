#!/usr/bin/env python3
"""Generate a self-contained HTML milestone closure report.

Per-milestone artifact produced *before* a milestone closes, serving three
layered audiences (PM retrospective · lab seminar · portfolio showcase) from one
file. Authored via the author-editor-critic triad — see the design spec at
``docs/superpowers/specs/2026-06-16-milestone-closure-report-design.md`` and
Issue #752.

Four layers (this module):
  1. Data        — pull every issue in the milestone (open + closed) from board #9
                   + ``gh issue list --milestone`` + per-issue labels/timestamps.
  2. Metrics     — pure, unit-tested functions over the data structs (no I/O).
  3. Aggregation — seed a first-draft narrative from lab-notebook + closing
                   comments into an author-owned ``<slug>.narrative.md`` sidecar.
  4. Render      — Jinja2 -> one self-contained HTML file (inline CSS).

Layers 3 and 4 lazy-import ``markdown`` / ``jinja2`` so the metrics functions
stay importable in the bare ``ci-tools-pytest`` env (pytest + pyyaml only).

Usage:
  scripts/pm/milestone_report.py "pm-i6 - PM Tooling & Memory Methodology II"
  scripts/pm/milestone_report.py "<milestone>" --out-dir docs/pm/milestone_reports
  scripts/pm/milestone_report.py "<milestone>" --dry-run   # print metrics, no write
"""

import argparse
import json
import re
import statistics
import subprocess
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional

OWNER = "Jin-HoMLee"
REPO = "splice-neoepitope-pipeline"
PROJECT_NUMBER = 9
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent  # scripts/pm -> scripts -> repo root
DEFAULT_OUT_DIR = REPO_ROOT / "docs" / "pm" / "milestone_reports"
DEFAULT_SDR_OUT_DIR = REPO_ROOT / "docs" / "pm" / "sdr_reports"
DEFAULT_TREND_WEEKS = 4
TEMPLATE_PATH = HERE / "templates" / "milestone_report.html.j2"
# Arrival marker (Issue #811): the item shipped without crossing Backlog -> Ready.
UNPLANNED_LABEL = "unplanned"


# --- slug -------------------------------------------------------------------

def slugify(name: str) -> str:
    """Sanitize a milestone name into a filename slug.

    Lowercase, collapse every non-alphanumeric run to a single hyphen, strip
    leading/trailing hyphens. ``"pm-i6 - PM Tooling & Memory Methodology II"``
    -> ``"pm-i6-pm-tooling-memory-methodology-ii"``.
    """
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


# --- time helpers (pure) ----------------------------------------------------

def parse_iso(ts: Optional[str]) -> Optional[datetime]:
    """Parse an ISO-8601 timestamp (``...Z`` accepted); None-safe."""
    if not ts:
        return None
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))


# --- metrics layer (pure, unit-tested) --------------------------------------

def closed_issues(issues: list[dict]) -> list[dict]:
    return [i for i in issues if i.get("state") == "CLOSED"]


# GitHub's IssueStateReason enum closes that are NOT deliveries: NOT_PLANNED
# (won't-do / superseded) and DUPLICATE (tracked elsewhere). Defined as an
# explicit set rather than a "!= COMPLETED" complement so a closed issue with
# no recorded reason still falls through to *delivered* (legacy data is never
# dropped), and so a future enum addition doesn't silently land in either bucket
# without a deliberate edit here. Issue #851 (DUPLICATE gap caught in review).
DESCOPED_REASONS = frozenset({"NOT_PLANNED", "DUPLICATE"})


def is_descoped(issue: dict) -> bool:
    """True if the issue closed as a descope (``NOT_PLANNED`` / ``DUPLICATE``)
    rather than a delivery. Single source for both the metrics split and the
    inventory badge."""
    return issue.get("state") == "CLOSED" and issue.get("state_reason") in DESCOPED_REASONS


def descoped_issues(issues: list[dict]) -> list[dict]:
    """Closed issues that did NOT ship (descoped / superseded / wontfix / dup).

    These are *not* deliverables — counting them as such inflates the delivered
    headline and masks dropped scope (the machine analogue of the parent-rollup
    close-reason rule). Issue #851.
    """
    return [i for i in issues if is_descoped(i)]


def delivered_issues(issues: list[dict]) -> list[dict]:
    """Closed issues that actually shipped (closed, and not a descope reason).

    A closed issue with no recorded reason is treated as delivered, so legacy
    data without ``stateReason`` is never silently dropped from the count.
    """
    return [i for i in issues if i.get("state") == "CLOSED" and not is_descoped(i)]


def is_unplanned(issue: dict) -> bool:
    """True when the issue shipped without crossing the Backlog -> Ready commitment.

    This is the **arrival** axis, not urgency (Issue #811). Priority (P0-P3) says
    how costly delay is; it cannot say whether the item was committed or barged
    in mid-session, so a P2 same-day fix and a P1 committed item are
    indistinguishable by band. The marker is the only thing that separates them,
    which is why the throughput split cannot be derived from priority.
    """
    return UNPLANNED_LABEL in issue.get("labels", [])


def unplanned_issues(issues: list[dict]) -> list[dict]:
    return [i for i in issues if is_unplanned(i)]


def committed_issues(issues: list[dict]) -> list[dict]:
    return [i for i in issues if not is_unplanned(i)]


def is_marker_in_use(repo_label_count: int) -> bool:
    """Is the arrival marker actually in use, repo-wide? Pure given the count.

    The whole point of the arrival axis is to be falsifiable, and this is the fact
    that makes it so. Split out as a pure function over an injected count (rather than
    doing the ``gh`` call inline) so the fabricating case is unit-testable without a
    network, which is exactly the property the metric was missing.

    A count of 0 means the label has never been applied to anything - so the arrival
    axis is **not in use**, and no delivered item can be classified by it. See
    ``throughput_breakdown`` for why that must not render as 0%.
    """
    return repo_label_count > 0


def throughput_breakdown(issues: list[dict], *, marker_in_use: bool) -> dict[str, Any]:
    """Split *delivered* work by arrival: committed vs unplanned.

    The Kanban Throughput Breakdown Chart, reduced to the one split we actually
    need. Keyed off delivered (not raw closed) for the same reason throughput and
    cycle time are: a descoped issue is not shipped work.

    **Two degenerate cases, both of which must NOT print a confident zero.**

    1. *Empty window* - no delivered issues at all. ``pct_unplanned`` is None (not
       0.0), so a zero-ship week reads as "no data". This guard already existed.

    2. *The marker is not in use* (``marker_in_use=False``) - i.e. the ``unplanned``
       label has never been applied to anything, repo-wide. This returned **0.0**,
       and 0.0 is indistinguishable from a real, hard-won, measured zero:

           | we genuinely absorbed no unplanned work | -> 0% |
           | nobody ever applied the label           | -> 0% |
           | the label-reading code is broken        | -> 0% |

       Three different worlds, one confident output. The number could only ever come
       back one way. That is not a measurement, it is a fabrication wearing a
       measurement's clothes - the most dangerous state for a metric (Issue #1180).

       When the marker is not in use, **nothing delivered can be classified by
       arrival at all**, so ``n_committed`` is None too. Printing "25 committed"
       would be exactly the same false assertion as "0% unplanned", one column over.
       The delivered count is reported as ``n_unclassifiable`` instead, so a degraded
       input is *visible* rather than silently absorbed.

    ``marker_in_use`` is a **required keyword**, deliberately: it cannot be derived
    from ``issues`` (from that list alone, "no delivered issue carries the label" is
    information-theoretically identical in both worlds), and a default would let a
    caller silently fall back into the fabricating behavior. Making it required means
    a new call site must *decide*.
    """
    delivered = delivered_issues(issues)
    n_delivered = len(delivered)

    if not marker_in_use:
        return {
            "n_committed": None,
            "n_unplanned": None,
            "n_unclassifiable": n_delivered,
            "pct_unplanned": None,
            "marker_in_use": False,
        }

    n_unplanned = len(unplanned_issues(delivered))
    return {
        "n_committed": n_delivered - n_unplanned,
        "n_unplanned": n_unplanned,
        "n_unclassifiable": 0,
        "pct_unplanned": (100.0 * n_unplanned / n_delivered) if n_delivered else None,
        "marker_in_use": True,
    }


def cycle_time_days(issue: dict) -> Optional[float]:
    """(closed_at - created_at) in days for a closed issue; None otherwise."""
    if issue.get("state") != "CLOSED":
        return None
    created = parse_iso(issue.get("created_at"))
    closed = parse_iso(issue.get("closed_at"))
    if not created or not closed:
        return None
    return (closed - created).total_seconds() / 86400.0


def cycle_times(issues: list[dict]) -> list[float]:
    return [t for t in (cycle_time_days(i) for i in issues) if t is not None]


def avg_cycle_time(issues: list[dict]) -> Optional[float]:
    ts = cycle_times(issues)
    return statistics.fmean(ts) if ts else None


def median_cycle_time(issues: list[dict]) -> Optional[float]:
    ts = cycle_times(issues)
    return statistics.median(ts) if ts else None


def per_role_counts(issues: list[dict]) -> dict[str, int]:
    """Delivered-issue counts grouped by ``role:*`` (multi-role issues count
    once per role). Descoped (``NOT_PLANNED``) closes are excluded — a dropped
    issue is not a per-role deliverable (Issue #851)."""
    counts: dict[str, int] = {}
    for issue in delivered_issues(issues):
        for role in issue.get("roles", []):
            if role.startswith("role:"):
                counts[role] = counts.get(role, 0) + 1
    return counts


def milestone_duration_days(
    issues: list[dict],
    milestone_closed_at: Optional[str],
    milestone_created_at: Optional[str] = None,
) -> Optional[float]:
    """Earliest issue ``created_at`` (or milestone creation, if no issues) to the
    milestone ``closed_at``, in days. None when either endpoint is unknowable."""
    end = parse_iso(milestone_closed_at)
    if not end:
        return None
    created = [parse_iso(i.get("created_at")) for i in issues if i.get("created_at")]
    created = [c for c in created if c]
    start = min(created) if created else parse_iso(milestone_created_at)
    if not start:
        return None
    return (end - start).total_seconds() / 86400.0


def throughput_per_week(n_delivered: int, duration_days: Optional[float]) -> Optional[float]:
    """Delivered issues per week; None when duration is zero/unknown (guarded)."""
    if not duration_days or duration_days <= 0:
        return None
    return n_delivered / (duration_days / 7.0)


# --- window-mode metrics (weekly SDR; pure, unit-tested) --------------------
# The meta-work Service Delivery Review (Issue #915) is a cadence-based
# retrospective decoupled from any milestone (per Issue #902 facet 2, meta-work
# flows milestone-free, so the per-milestone report never fires for it). These
# pure functions bucket delivered issues into trailing weeks so the report can
# show a throughput + cycle-time TREND, not just a single-window snapshot.

def normalize_until(until: datetime) -> datetime:
    """Anchor a window end on the last instant of its UTC day.

    The trend is bucketed with datetime precision but *sourced* from GitHub's
    ``closed:A..B`` search operator, which is day-granular. A mid-day anchor (what
    ``datetime.now()`` yields in a live run) therefore puts every bucket edge at a
    time-of-day the source data cannot resolve, which both displaces issues into
    adjacent buckets and strands those past the final edge. Snapping the anchor to
    the day boundary makes the bucket grid agree with the granularity of the data
    it buckets. Idempotent. Issue #1099.

    ``until`` must be timezone-aware. A naive datetime would be read as *local*
    time by ``astimezone``, silently shifting every bucket edge by the host's UTC
    offset - the exact class of bug this function exists to remove.
    """
    if until.tzinfo is None:
        raise ValueError(
            "normalize_until requires a timezone-aware datetime; a naive one would "
            "be interpreted as local time and shift every bucket edge."
        )
    return until.astimezone(timezone.utc).replace(
        hour=23, minute=59, second=59, microsecond=999999
    )


def _first_covered_day(span_start: datetime) -> date:
    """The first whole day a ``(span_start, ...]`` bucket actually holds.

    ``span_start`` is an **exclusive** edge sitting at the last instant of its own
    day, so that day belongs to no bucket. Both the fetch range and the rendered
    row label must key off the day *after* it, and they must agree - hence one
    helper rather than the same ``+ 1 day`` written twice. Issue #1099.
    """
    return (span_start + timedelta(days=1)).date()


def week_windows(until: datetime, n_weeks: int) -> list[tuple[datetime, datetime]]:
    """``n_weeks`` trailing 7-day ``(start, end]`` windows ending at ``until``.

    Chronological (oldest first) so a rendered series reads left-to-right;
    ``windows[-1]`` is the reporting week ``(until - 7d, until]``, earlier
    entries are the trend history.
    """
    windows = [
        (until - timedelta(days=7 * (i + 1)), until - timedelta(days=7 * i))
        for i in range(max(n_weeks, 1))
    ]
    return list(reversed(windows))


def fetch_span(until: datetime, n_weeks: int) -> tuple[str, str]:
    """ISO date range whose whole-day coverage exactly matches ``week_windows``.

    Buckets are **half-open at the start** (``(span_start, end]``), and with a
    normalized ``until`` the ``span_start`` edge sits at the *last instant* of its
    own day - so that day belongs to no bucket. Fetching from ``span_start.date()``
    would pull a full day of issues that nothing can count, which is the bottom-edge
    leak behind Issue #1099. Start the fetch the day after instead, so every issue
    the fetch returns lands in exactly one bucket (conservation).
    """
    span_start = week_windows(until, n_weeks)[0][0]
    return _first_covered_day(span_start).isoformat(), until.date().isoformat()


def closed_in_window(issues: list[dict], start: datetime, end: datetime) -> list[dict]:
    """Closed issues whose ``closed_at`` falls in ``(start, end]``. Pure."""
    out = []
    for i in closed_issues(issues):
        c = parse_iso(i.get("closed_at"))
        if c and start < c <= end:
            out.append(i)
    return out


def weekly_series(issues: list[dict], until: datetime, n_weeks: int,
                  *, marker_in_use: bool) -> list[dict]:
    """Per-week throughput + median cycle-time trend over the trailing weeks.

    Each entry: ``{week_start, week_end, n_delivered, median_cycle_time_days}``
    plus the arrival split ``throughput_breakdown`` contributes
    (``{n_committed, n_unplanned, pct_unplanned}``) - dates as ISO ``YYYY-MM-DD``,
    chronological. Throughput + cycle time key off *delivered* issues (a descoped
    close is not shipped work), consistent with the milestone-mode metrics.

    ``week_start`` names the **first day the bucket actually covers**, not the
    half-open ``(start`` edge - which, with a normalized ``until``, is the last
    instant of the *preceding* day. Rendered as ``week_start..week_end`` the pair
    reads as an inclusive day range, so it must be one: a reader who cross-checks a
    row with ``closed:<week_start>..<week_end>`` has to get that row's number back.
    Naming the raw edge instead put a day the bucket does not hold into the label,
    and that mismatch is precisely what produced the bogus reference counts in
    Issue #1099.
    """
    series = []
    for start, end in week_windows(until, n_weeks):
        in_window = closed_in_window(issues, start, end)
        delivered = delivered_issues(in_window)
        series.append({
            "week_start": _first_covered_day(start).isoformat(),
            "week_end": end.date().isoformat(),
            **throughput_breakdown(in_window, marker_in_use=marker_in_use),
            "n_delivered": len(delivered),
            "median_cycle_time_days": median_cycle_time(delivered),
        })
    return series


def compute_window_metrics(
    week_issues: list[dict], all_issues: list[dict], until: datetime, n_weeks: int,
    *, marker_in_use: bool
) -> dict[str, Any]:
    """Headline (reporting-week) metrics + the weekly trend for SDR window mode.

    ``week_issues`` are the closed meta-work issues in the reporting week (the
    headline); ``all_issues`` span the full trend lookback (the series). Mirrors
    ``compute_metrics``' keys so one template renders both modes, plus the extra
    ``weekly_series``. There is no carried-forward concept in window mode (only
    closed issues are fetched), so it is fixed at 0.
    """
    delivered = delivered_issues(week_issues)
    return {
        "n_total": len(week_issues),
        "n_closed": len(closed_issues(week_issues)),
        "n_delivered": len(delivered),
        "n_descoped": len(descoped_issues(week_issues)),
        "n_carried_forward": 0,
        "duration_days": 7.0,
        "throughput_per_week": throughput_per_week(len(delivered), 7.0),
        "avg_cycle_time_days": avg_cycle_time(delivered),
        "median_cycle_time_days": median_cycle_time(delivered),
        "per_role_counts": per_role_counts(week_issues),
        **throughput_breakdown(week_issues, marker_in_use=marker_in_use),
        "weekly_series": weekly_series(all_issues, until, n_weeks,
                                       marker_in_use=marker_in_use),
    }


def compute_metrics(issues: list[dict], milestone: dict,
                    *, marker_in_use: bool) -> dict[str, Any]:
    """Assemble the headline metrics block from the data layer."""
    closed = closed_issues(issues)
    delivered = delivered_issues(issues)
    descoped = descoped_issues(issues)
    duration = milestone_duration_days(
        issues, milestone.get("closed_at"), milestone.get("created_at")
    )
    return {
        "n_total": len(issues),
        "n_closed": len(closed),
        "n_delivered": len(delivered),
        "n_descoped": len(descoped),
        "n_carried_forward": len(issues) - len(closed),
        "duration_days": duration,
        # Throughput + cycle time key off *delivered*, not raw closed — a
        # descoped issue is not shipped work, so neither the rate nor the
        # created→closed cycle should count it (Issue #851).
        "throughput_per_week": throughput_per_week(len(delivered), duration),
        "avg_cycle_time_days": avg_cycle_time(delivered),
        "median_cycle_time_days": median_cycle_time(delivered),
        "per_role_counts": per_role_counts(issues),
        **throughput_breakdown(issues, marker_in_use=marker_in_use),
    }


# --- data layer (gh / board #9) ---------------------------------------------

def _gh(args: list[str]) -> str:
    try:
        return subprocess.run(
            ["gh", *args], check=True, capture_output=True, text=True
        ).stdout
    except subprocess.CalledProcessError as exc:
        # Surface gh's own error text (auth/network/unknown-resource) instead of
        # the opaque CalledProcessError the caller would otherwise see.
        if exc.stderr:
            print(exc.stderr.rstrip(), file=sys.stderr)
        raise


def fetch_milestone(name: str) -> dict:
    """Resolve a milestone by full name -> {title, created_at, closed_at, state}.

    Uses ``--jq '.[]'`` so ``--paginate`` emits NDJSON (one object per line)
    across pages — concatenated REST array pages are NOT valid JSON, so a bare
    ``json.loads`` over a multi-page ``--paginate`` body crashes once the repo
    exceeds one page of milestones.
    """
    raw = _gh([
        "api", f"repos/{OWNER}/{REPO}/milestones",
        "--paginate", "-X", "GET", "-f", "state=all", "-f", "per_page=100",
        "--jq", ".[]",
    ])
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        ms = json.loads(line)
        if ms.get("title") == name:
            return {
                "title": ms["title"],
                "number": ms["number"],
                "state": ms.get("state"),
                "created_at": ms.get("created_at"),
                "closed_at": ms.get("closed_at"),
                "due_on": ms.get("due_on"),
                "description": ms.get("description"),
            }
    raise SystemExit(f"milestone not found: {name!r}")


def _board_fields_by_number() -> dict[int, dict]:
    """Best-effort map issue number -> board Status/Priority/Size/arc.

    Reuses scripts/board_open_items.py (the paginated helper). Returns {} if the
    helper is unavailable or errors — the report degrades to label-only data.
    """
    try:
        scripts_dir = str(HERE.parent)
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        import board_open_items as boi  # type: ignore
        items = boi.fetch_all_items()
        out: dict[int, dict] = {}
        for raw in items:
            norm = boi.normalize(raw)
            if norm and norm.get("number"):
                out[int(norm["number"])] = norm
        return out
    except Exception as exc:  # pragma: no cover - best-effort enrichment
        print(f"  (board enrichment unavailable: {exc})", file=sys.stderr)
        return {}


def _normalize_issue(it: dict, board: dict[int, dict]) -> dict:
    """Normalize one ``gh issue list --json`` record into the report's issue dict.

    Shared by the milestone and window fetchers so both carry an identical shape.
    """
    labels = [lbl["name"] for lbl in it.get("labels", [])]
    num = it["number"]
    b = board.get(num, {})
    # stateReason: COMPLETED | NOT_PLANNED | None (open). Upper-cased so the
    # metrics layer can compare against the canonical NOT_PLANNED token.
    reason = it.get("stateReason")
    issue = {
        "number": num,
        "title": it["title"],
        "url": it.get("url"),
        "state": it["state"].upper(),
        "state_reason": reason.upper() if reason else None,
        "created_at": it.get("createdAt"),
        "closed_at": it.get("closedAt"),
        "labels": labels,
        "roles": [l for l in labels if l.startswith("role:")],
        "arcs": [l for l in labels if l.startswith("arc:")],
        "status": b.get("status") or ("Done" if it["state"].upper() == "CLOSED" else "n/a"),
        # "n/a", not an em-dash: this string is emitted into the generated HTML,
        # where the house no-em-dash rule applies and no guard can see it (the
        # PreToolUse guard scans Claude's edits, not a script's output). Issue #1180.
        "priority": b.get("priority") or "n/a",
        "size": b.get("size") or "n/a",
    }
    # Single-source the descoped flag for the inventory badge (covers
    # NOT_PLANNED + DUPLICATE without the template hardcoding the set).
    issue["is_descoped"] = is_descoped(issue)
    return issue


def fetch_marker_in_use() -> bool:
    """Has the ``unplanned`` label EVER been applied, repo-wide?

    The falsifier the arrival metric was missing. Deliberately repo-wide and
    all-states: the question is not "did this window have unplanned work" (which the
    issue list already answers, ambiguously) but "is the marker being used at all" -
    a fact that no window of issues can supply about itself.

    ``--limit 1`` because the count is irrelevant; only zero-vs-nonzero matters.
    """
    raw = _gh([
        "issue", "list", "--repo", f"{OWNER}/{REPO}",
        "--label", UNPLANNED_LABEL, "--state", "all", "--limit", "1",
        "--json", "number",
    ])
    return is_marker_in_use(len(json.loads(raw)))


def fetch_milestone_issues(name: str) -> list[dict]:
    """All issues in the milestone (open + closed), normalized for the report."""
    raw = _gh([
        "issue", "list", "--repo", f"{OWNER}/{REPO}",
        "--milestone", name, "--state", "all", "--limit", "1000",
        "--json", "number,title,state,stateReason,labels,createdAt,closedAt,url",
    ])
    board = _board_fields_by_number()
    return [_normalize_issue(it, board) for it in json.loads(raw)]


def fetch_window_issues(since: str, until: str) -> list[dict]:
    """Closed *meta-work* issues in the ``[since, until]`` date window (ISO dates).

    Meta-work = closed issues carrying **no milestone** (per Issue #902 facet 2:
    meta-work flows milestone-free; lifecycle work is milestoned and covered by
    the per-milestone report). Milestoned closes are filtered out here so the SDR
    covers only flow work.
    """
    raw = _gh([
        "issue", "list", "--repo", f"{OWNER}/{REPO}",
        "--state", "closed", "--search", f"closed:{since}..{until}", "--limit", "1000",
        "--json", "number,title,state,stateReason,labels,createdAt,closedAt,url,milestone",
    ])
    board = _board_fields_by_number()
    return [
        _normalize_issue(it, board)
        for it in json.loads(raw)
        if not it.get("milestone")  # skip lifecycle (milestoned) work
    ]


# --- aggregation layer (narrative auto-seed) --------------------------------

_SEED_HEADLINE_MAX = 180
# Byline/timestamp sub-headers to skip when digesting an entry — e.g. the PM
# notebook's "### HH:MM UTC — Editor: PM" line, which is metadata, not content.
_BYLINE_RE = re.compile(r"^\d{1,2}:\d{2}\b|\b(editor|author|role)\s*:", re.IGNORECASE)


def _first_prose_line(body: str) -> str:
    """The first human-prose line of an entry body — skips blank lines, markdown
    bullet/heading markers, HTML comments, and byline/timestamp sub-headers
    (so the digest picks the descriptive title, not "Editor: PM"). One line."""
    for line in body.strip().splitlines():
        stripped = line.strip().lstrip("#*->–— ").strip()
        if not stripped or stripped.startswith("<!--"):
            continue
        if _BYLINE_RE.search(stripped):
            continue
        return stripped
    return ""


def _lab_notebook_seed(roles: set[str], window: tuple[Optional[datetime], Optional[datetime]]) -> str:
    """Best-effort: digest lab-notebook entries dated within the milestone window
    for the involved roles into one pointer bullet each (date + headline), NOT
    verbatim bodies — the sidecar is meant for a human to expand. '' if none."""
    start, end = window
    bullets: list[str] = []
    for role in sorted(roles):
        short = role.split(":", 1)[-1]
        nb = REPO_ROOT / "research" / "lab_notebook" / f"{short}.md"
        if not nb.exists():
            continue
        text = nb.read_text(encoding="utf-8", errors="replace")
        # Entries are "## <YYYY-MM-DD>" headed; digest those inside the window.
        for m in re.finditer(r"^## (\d{4}-\d{2}-\d{2})\b(.*?)(?=^## \d{4}-\d{2}-\d{2}\b|\Z)",
                             text, flags=re.DOTALL | re.MULTILINE):
            day = parse_iso(m.group(1) + "T00:00:00Z")
            if (start and day and day < start) or (end and day and day > end):
                continue
            headline = _first_prose_line(m.group(2))
            if len(headline) > _SEED_HEADLINE_MAX:
                headline = headline[: _SEED_HEADLINE_MAX - 1].rstrip() + "…"
            bullets.append(f"- **{m.group(1)}** ({short}): {headline}")
    return "\n".join(bullets)


def seed_narrative(milestone: dict, issues: list[dict], metrics: dict) -> str:
    """First-draft narrative markdown for the author-owned sidecar."""
    roles = {r for i in issues for r in i.get("roles", [])}
    window = (parse_iso(milestone.get("created_at")), parse_iso(milestone.get("closed_at")))
    seed = _lab_notebook_seed(roles, window)
    delivered = delivered_issues(issues)
    descoped = descoped_issues(issues)
    carried = [i for i in issues if i.get("state") != "CLOSED"]

    lines = [
        f"<!-- Author-owned narrative for {milestone['title']}. Sections 3/4/5 only.",
        "     The script regenerates the HTML from this file + fresh board data;",
        "     it never overwrites this sidecar once it exists. -->",
        "",
        "## Deliverables (Review layer)",
        "",
        "<!-- Lead role: what shipped, grouped by deliverable, with PR + slide links.",
        f"     The {len(delivered)} delivered issues are listed in the Inventory appendix",
        "     below - narrate the highlights here, don't re-list them. -->",
        "",
        "_Seeded from the lead role's lab-notebook entries in the milestone window - "
        "replace with the deliverables narrative._",
        "",
    ]
    lines += [seed or "<!-- no auto-seed found; author from scratch -->"]
    if descoped:
        lines += [
            "",
            "## Descoped (closed NOT_PLANNED)",
            "",
            "<!-- These closed WITHOUT shipping (superseded / YAGNI / wontfix). They are",
            "     excluded from the delivered count - record why each was dropped + where",
            "     the need (if any) was routed, so the descope isn't silently masked. -->",
            "",
        ]
        lines += [f"- [#{i['number']}]({i['url']}) {i['title']} - _why descoped: TBD_" for i in descoped]
    lines += [
        "",
        "## Carried-forward & routing",
        "",
        "<!-- PM: issues that didn't close + where they went (carve / arc) + the",
        "     closure-routing decision (a/b/c/d). -->",
        "",
    ]
    lines += [f"- [#{i['number']}]({i['url']}) {i['title']} - _route: TBD_" for i in carried] or ["- _(none carried forward)_"]
    lines += [
        "",
        "## Retrospective (process/health)",
        "",
        "<!-- PM: was it healthy? what to improve? WIP/aging observations. -->",
        "",
    ]
    return "\n".join(lines) + "\n"


# --- render layer (jinja2; lazy import) -------------------------------------

def _md_to_html(md_text: str) -> str:
    import markdown  # lazy
    return markdown.markdown(md_text, extensions=["tables", "fenced_code"])


def render_html(
    milestone: dict, issues: list[dict], metrics: dict, narrative_md: str,
    mode: str = "milestone",
) -> str:
    from jinja2 import Environment, FileSystemLoader, select_autoescape  # lazy

    env = Environment(
        loader=FileSystemLoader(str(TEMPLATE_PATH.parent)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    env.filters["pct"] = lambda v: f"{v:.0%}" if v is not None else "n/a"
    template = env.get_template(TEMPLATE_PATH.name)
    # Milestone-level arc(s) = the distinct arc labels across its issues (arc is
    # a per-issue label under the three-axis model, not a milestone property).
    arcs = sorted({a for i in issues for a in i.get("arcs", [])})
    return template.render(
        milestone=milestone,
        issues=issues,
        closed=closed_issues(issues),
        carried=[i for i in issues if i.get("state") != "CLOSED"],
        metrics=metrics,
        arcs=arcs,
        mode=mode,
        narrative_html=_md_to_html(narrative_md),
        generated_at=milestone.get("closed_at") or milestone.get("due_on") or "",
    )


def _fmt(v: Optional[float], unit: str = "") -> str:
    return "n/a" if v is None else f"{v:.1f}{unit}"


def print_metrics(milestone: dict, metrics: dict) -> None:
    state = milestone.get("state")
    if state and state != "n/a":  # milestone mode
        print(f"Milestone: {milestone['title']}  [{state}]")
    else:  # window (SDR) mode: the title already says what it is
        print(milestone["title"])
    print(f"  total / closed / carried-forward : "
          f"{metrics['n_total']} / {metrics['n_closed']} / {metrics['n_carried_forward']}")
    print(f"  delivered / descoped (closed)    : "
          f"{metrics['n_delivered']} / {metrics['n_descoped']}")
    pct = metrics.get("pct_unplanned")
    if not metrics.get("marker_in_use", True):
        # The marker has never been applied repo-wide, so NOTHING delivered can be
        # classified by arrival. Printing "N / 0  (0% unplanned)" here would be a
        # fabrication wearing a measurement's clothes (Issue #1180).
        print(f"  arrival: UNCLASSIFIABLE          : "
              f"{metrics['n_unclassifiable']} delivered, none classifiable "
              f"(the `{UNPLANNED_LABEL}` marker is not in use repo-wide)")
    else:
        print(f"  arrival: committed / unplanned   : "
              f"{metrics['n_committed']} / {metrics['n_unplanned']}"
              f"{f'  ({pct:.0f}% unplanned)' if pct is not None else '  (no delivered work)'}")
    print(f"  duration (days)                  : {_fmt(metrics['duration_days'])}")
    print(f"  throughput (delivered/week)      : {_fmt(metrics['throughput_per_week'])}")
    print(f"  cycle time avg / median (days)   : "
          f"{_fmt(metrics['avg_cycle_time_days'])} / {_fmt(metrics['median_cycle_time_days'])}")
    print(f"  per-role (delivered)             : {metrics['per_role_counts']}")
    series = metrics.get("weekly_series")
    if series:
        print("  weekly trend (delivered [committed+unplanned] / median cycle days):")
        for w in series:
            print(f"    {w['week_start']}..{w['week_end']} : "
                  f"{w['n_delivered']} [{_fmt_arrival(w)}] / "
                  f"{_fmt(w['median_cycle_time_days'])}")


def _fmt_arrival(w: dict) -> str:
    """Render one weekly row's arrival split, or say it cannot be rendered.

    Mirrors the headline `arrival:` line. Without this the row printed the literal
    ``[None+None]`` whenever the marker was not in use - not a fabrication (None is at
    least honest), but a leak of an internal sentinel into a human-facing report, and
    inconsistent with every other surface this fix touched.

    Caught in review, and worth recording HOW it survived: the live check that was
    supposed to verify this fix grepped the output for `arrival` - i.e. for the line
    already known to be fixed. A check aimed only at what you changed cannot show you
    what you missed.
    """
    if w["n_committed"] is None:
        return "unclassifiable"
    return f"{w['n_committed']}+{w['n_unplanned']}"


# --- orchestration ----------------------------------------------------------

def _generate(milestone: dict, issues: list[dict], metrics: dict,
              out_dir: Path, dry_run: bool, mode: str) -> int:
    """Shared metrics-print + seed + render path for both modes."""
    print_metrics(milestone, metrics)
    if dry_run:
        return 0
    slug = slugify(milestone["title"])
    out_dir.mkdir(parents=True, exist_ok=True)
    sidecar = out_dir / f"{slug}.narrative.md"
    if not sidecar.exists():  # never overwrite author edits
        sidecar.write_text(seed_narrative(milestone, issues, metrics), encoding="utf-8")
        print(f"  seeded narrative -> {sidecar}")
    else:
        print(f"  narrative exists (preserved) -> {sidecar}")
    html = render_html(milestone, issues, metrics, sidecar.read_text(encoding="utf-8"), mode=mode)
    html_path = out_dir / f"{slug}.html"
    html_path.write_text(html, encoding="utf-8")
    print(f"  rendered report  -> {html_path}")
    return 0


def _run_milestone_mode(name: str, out_dir: Optional[Path], dry_run: bool) -> int:
    milestone = fetch_milestone(name)
    issues = fetch_milestone_issues(name)
    metrics = compute_metrics(issues, milestone, marker_in_use=fetch_marker_in_use())
    return _generate(milestone, issues, metrics, out_dir or DEFAULT_OUT_DIR, dry_run, "milestone")


def _parse_date_arg(ap: argparse.ArgumentParser, flag: str, value: str) -> datetime:
    """Validate a ``YYYY-MM-DD`` CLI date -> tz-aware UTC datetime; ap.error else."""
    try:
        return datetime.strptime(value, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    except ValueError:
        ap.error(f"{flag} must be YYYY-MM-DD, got: {value!r}")


def _run_window_mode(ap: argparse.ArgumentParser, since: Optional[str], until: Optional[str],
                     trend_weeks: int, out_dir: Optional[Path], dry_run: bool) -> int:
    """Weekly meta-work SDR. Reporting week = the most recent 7-day window; the
    trend spans ``n_weeks`` back. Skips a reporting week with nothing closed."""
    # Normalize BOTH paths to the end of the UTC day. The explicit --until already
    # did; the default (live `now()`) path did not, so a real run bucketed on
    # mid-day edges the day-granular source data cannot resolve (Issue #1099).
    until_dt = normalize_until(
        _parse_date_arg(ap, "--until", until) if until else datetime.now(timezone.utc)
    )
    # --since sets the TREND SPAN: derive the week count from since..until so the
    # trend actually reaches back to --since (overriding --trend-weeks). Deriving
    # a single span-start without this would only clip the fetch and could
    # silently under-count the headline (PR #960 review).
    if since:
        since_dt = _parse_date_arg(ap, "--since", since)
        if since_dt >= until_dt:
            ap.error("--since must be before --until")
        span_days = (until_dt.date() - since_dt.date()).days
        n_weeks = max(1, -(-span_days // 7))  # ceil-divide into whole weeks
    else:
        n_weeks = trend_weeks
    windows = week_windows(until_dt, n_weeks)
    report_start, report_end = windows[-1]

    all_issues = fetch_window_issues(*fetch_span(until_dt, n_weeks))
    week_issues = closed_in_window(all_issues, report_start, report_end)
    metrics = compute_window_metrics(week_issues, all_issues, until_dt, n_weeks,
                                     marker_in_use=fetch_marker_in_use())

    if metrics["n_total"] == 0:  # nothing closed in the reporting week -> no artifact
        print(f"Meta-work SDR - week ending {report_end.date().isoformat()}: "
              "no meta-work issues closed in the reporting week; skipping (no artifact).")
        return 0

    pseudo = {
        "title": f"Meta-work SDR - week ending {report_end.date().isoformat()}",
        "number": None,
        "state": "n/a",
        "created_at": report_start.isoformat(),
        "closed_at": report_end.isoformat(),
        "due_on": None,
        "description": None,
    }
    return _generate(pseudo, week_issues, metrics, out_dir or DEFAULT_SDR_OUT_DIR, dry_run, "window")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Generate a milestone closure report, or a weekly meta-work "
                    "Service Delivery Review (SDR) when no milestone is given.")
    ap.add_argument("milestone", nargs="?",
                    help="full milestone name (e.g. 'pm-i6 - PM Tooling & Memory Methodology II'); "
                         "omit (or pass --weekly) for the weekly meta-work SDR window mode")
    ap.add_argument("--weekly", action="store_true",
                    help="force weekly SDR window mode even if a milestone arg is present")
    ap.add_argument("--since", help="SDR trend-span start (YYYY-MM-DD); sets the trend length "
                                    "from since..until, overriding --trend-weeks. Reaches back AT "
                                    "LEAST this far, rounded outward to whole 7-day weeks, so the "
                                    "first row may begin a few days before this date")
    ap.add_argument("--until", help="SDR window end (YYYY-MM-DD); default = today. Anchored to the "
                                    "end of that UTC day, matching the day-granular `closed:` search "
                                    "the trend is sourced from")
    ap.add_argument("--trend-weeks", type=int, default=DEFAULT_TREND_WEEKS,
                    help=f"SDR trailing weeks shown in the trend (default {DEFAULT_TREND_WEEKS})")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output dir (default: milestone_reports for milestone mode, "
                         "sdr_reports for weekly SDR mode)")
    ap.add_argument("--dry-run", action="store_true", help="print metrics only; write nothing")
    args = ap.parse_args()

    if args.milestone and not args.weekly:
        return _run_milestone_mode(args.milestone, args.out_dir, args.dry_run)
    if args.trend_weeks < 1:
        ap.error("--trend-weeks must be >= 1")
    return _run_window_mode(ap, args.since, args.until, args.trend_weeks, args.out_dir, args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())
