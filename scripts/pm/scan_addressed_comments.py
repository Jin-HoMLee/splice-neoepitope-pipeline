#!/usr/bin/env python3
"""Board-wide scan for comments addressed `**To:** <role>` (Daily Stand-up Beat 2).

A team-coordination comment opens `**From:** <Role> -> **To:** <Role(s)>`. A ping
`To:` your role can land on ANY Issue/PR regardless of its `role:*` label (e.g. a
Developer replies `To: PM` on a `role:developer` Issue), so a role-scoped scan of
only your own `role:*` Issues structurally misses every cross-role ping - which is
exactly how PM missed the Developer's `To: PM` reply on `role:developer`
[Issue #887] (2026-06-29). There is no @mention/notification for cross-role pings
(all roles are one GitHub user), so detection depends entirely on scanning the
right surface board-wide.

This replaces the error-prone daily hand-roll (the inline jq/zsh form broke 3x in
five minutes on 2026-06-29). Usage:

  scripts/pm/scan_addressed_comments.py --role pm
  scripts/pm/scan_addressed_comments.py --role developer --days 3
  scripts/pm/scan_addressed_comments.py --role memory_manager --since 2026-07-01T00:00:00Z

Window: the last-session watermark (`.agents/last_session_marker.json`, written by
the Stop hook) minus a 1-day backward overlap; a conservative 7-day floor when the
marker is absent (fresh clone / first run). `--since` / `--days` override it.

Exit 0 on a completed scan - finding nothing is not a failure. A hard `gh` /
network failure on the board-listing call surfaces loudly (non-zero traceback)
rather than as a false empty result. Issue #901.

[Issue #887]: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/887
"""
import argparse
import json
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
MARKER_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "last_session_marker.json"

# Backward overlap on the watermark + the floor when no marker exists.
OVERLAP_DAYS = 1
FLOOR_DAYS = 7

# Role -> display aliases as they appear in a `**To:** <Role>` field. Full names
# plus the short forms the team uses on the board (Sci / Dev / MM). A `To: all` /
# `To: team` broadcast is deliberately NOT matched: only a *named* addressee owes
# an ack (per feedback_team_coordination.md), and surfacing every broadcast to
# every role would be noise.
ROLE_ALIASES = {
    "pm": ["PM"],
    "scientist": ["Scientist", "Sci"],
    "developer": ["Developer", "Dev"],
    "memory_manager": ["Memory Manager", "MM"],
}

# The literal marker after which the addressee list appears. Mirrors the proven
# jq predicate `contains("To:** <Role>")`: the addressing convention bolds the
# label (`**To:** PM`), so the tail after `To:**` is the recipient list. We match
# on this literal (no `test()` regex) because the regex form's `\*` / `\b`
# escaping silently matched nothing when hand-rolled (2026-06-29).
_TO_MARKER = "To:**"


# --- pure helpers (unit-tested, no I/O) ---


def display_names(role):
    """Aliases to look for in a `To:` field for `role`. Raises on unknown role."""
    if role not in ROLE_ALIASES:
        raise ValueError(f"unknown role: {role!r} (choices: {', '.join(ROLE_ALIASES)})")
    return ROLE_ALIASES[role]


def _to_field_tails(body):
    """Yield the text after each `To:**` marker, one per line that carries it.

    Only the recipient list (the tail of the `**To:**` line) is searched, so a
    role name mentioned elsewhere in the comment body never false-matches.
    """
    for line in (body or "").splitlines():
        idx = line.find(_TO_MARKER)
        if idx != -1:
            yield line[idx + len(_TO_MARKER):]


def _word_bounded(alias_lower, tail):
    """True if `alias_lower` appears as a whole word/phrase in `tail`.

    Punctuation (commas between multiple recipients, arrows, emphasis) is
    flattened to spaces so `**To:** Scientist, PM` matches `pm`, while `PMx` or a
    substring hit inside another word does not. Pure string ops - no regex
    predicate, so the jq `test()` escaping footgun cannot recur.
    """
    flat = "".join(c if c.isalnum() else " " for c in tail.lower())
    padded = f" {' '.join(flat.split())} "
    return f" {alias_lower} " in padded


def body_addresses_role(body, names):
    """True if any `To:` field in `body` names one of `names` (the role aliases)."""
    aliases = [n.lower() for n in names]
    return any(
        _word_bounded(a, tail) for tail in _to_field_tails(body) for a in aliases
    )


def compute_since(marker, now):
    """Window floor as an aware UTC datetime.

    `marker` is the parsed `.agents/last_session_marker.json` dict (or None). Uses
    `last_session_end_utc` minus OVERLAP_DAYS when present + parseable, else a
    FLOOR_DAYS floor before `now`.
    """
    ts = (marker or {}).get("last_session_end_utc")
    if ts:
        try:
            end = datetime.fromisoformat(ts.replace("Z", "+00:00"))
            if end.tzinfo is None:
                end = end.replace(tzinfo=timezone.utc)
            return end - timedelta(days=OVERLAP_DAYS)
        except (ValueError, TypeError):
            pass
    return now - timedelta(days=FLOOR_DAYS)


def snippet(body, maxlen=140):
    """One-line preview of a comment body, whitespace-collapsed and truncated."""
    flat = " ".join((body or "").split())
    return flat if len(flat) <= maxlen else flat[: maxlen - 3] + "..."


def parse_ts(ts):
    """Parse a GitHub ISO-8601 timestamp to an aware UTC datetime (or None)."""
    if not ts:
        return None
    try:
        dt = datetime.fromisoformat(ts.replace("Z", "+00:00"))
        return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    except (ValueError, TypeError):
        return None


def select_pings(comments, names, since):
    """Filter raw REST comments to those created >= `since` AND addressing the role.

    Returns [(login, created_at, snippet)], oldest first. A comment with an
    unparseable/absent timestamp is kept (fail-open: better a stray line than a
    dropped ping).
    """
    out = []
    for c in comments or []:
        created = parse_ts(c.get("created_at"))
        if created is not None and created < since:
            continue
        if not body_addresses_role(c.get("body"), names):
            continue
        login = ((c.get("user") or {}).get("login")) or "?"
        out.append((login, c.get("created_at") or "?", snippet(c.get("body"))))
    return out


# --- gh I/O ---


def gh(*args, parse_json=True):
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def read_marker(path=MARKER_PATH):
    """Parse the session watermark file, or None if absent/malformed."""
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def recently_updated(since_date):
    """Board-wide Issues + PRs updated on/after `since_date` (YYYY-MM-DD).

    Two calls (issues + PRs) rather than the board's Done-first project query, so
    the pagination trap does not apply and PRs are covered.
    """
    search = f"updated:>={since_date}"
    items = {}
    for kind, sub in (("Issue", "issue"), ("PR", "pr")):
        rows = gh(
            sub, "list", "--repo", REPO, "--state", "all", "--limit", "1000",
            "--search", search, "--json", "number,title",
        )
        for r in rows:
            items[r["number"]] = {**r, "kind": kind}
    return [items[n] for n in sorted(items)]


def fetch_comments(number):
    """All REST comments on Issue/PR `number` (REST serves both). Empty on error."""
    try:
        return gh(
            "api", f"repos/{REPO}/issues/{number}/comments", "--paginate",
            # `.user` is null for a fully-deleted account; guard it so a single
            # such comment doesn't error the whole `--jq` page (which would exit
            # gh non-zero and silently drop EVERY comment on this issue).
            "--jq", '.[] | {body, created_at, user: {login: (.user.login? // "?")}}',
            parse_json=False,
        ).strip()
    except subprocess.CalledProcessError:
        return ""


def _parse_jsonl(text):
    """Parse newline-delimited JSON objects (gh --jq streams one per line)."""
    rows = []
    for line in (text or "").splitlines():
        line = line.strip()
        if line:
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return rows


# --- orchestration ---


def render(role, groups, since):
    lines = [
        f"Board pings **To:** {role} since {since.strftime('%Y-%m-%d %H:%MZ')} "
        f"(board-wide, all recently-updated Issues/PRs):"
    ]
    if not groups:
        lines.append("  (none)")
        return "\n".join(lines) + "\n"
    for item, pings in groups:
        lines.append(f"\n#{item['number']} [{item['kind']}] {item['title']}")
        for login, ts, snip in pings:
            lines.append(f"  - @{login} {ts}: {snip}")
    return "\n".join(lines) + "\n"


def scan(role, since):
    names = display_names(role)
    since_date = since.strftime("%Y-%m-%d")
    groups = []
    for item in recently_updated(since_date):
        comments = _parse_jsonl(fetch_comments(item["number"]))
        pings = select_pings(comments, names, since)
        if pings:
            groups.append((item, pings))
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--role", required=True, choices=sorted(ROLE_ALIASES),
                        help="scan for comments addressed To: this role")
    parser.add_argument("--since", help="ISO-8601 window floor (overrides the watermark)")
    parser.add_argument("--days", type=int, help="window = last N days (overrides the watermark)")
    args = parser.parse_args()

    now = datetime.now(timezone.utc)
    if args.since:
        since = parse_ts(args.since)
        if since is None:
            parser.error(f"unparseable --since: {args.since!r}")
    elif args.days is not None:
        if args.days < 0:
            parser.error(f"--days must be non-negative, got {args.days}")
        since = now - timedelta(days=args.days)
    else:
        since = compute_since(read_marker(), now)

    print(render(args.role, scan(args.role, since), since), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
