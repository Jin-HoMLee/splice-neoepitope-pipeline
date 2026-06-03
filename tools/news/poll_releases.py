#!/usr/bin/env python3
"""Morning-news release-feed poller (Issue #639).

Polls upstream feeds for the project's dependency basket (``tools.yaml``) and
emits only version DELTAS relative to a per-role watermark
(``watermark_<role>.yaml``). This replaces evergreen free-text web search as the
primary engine for the morning-routine News phase, ending the cross-day
repetition diagnosed in Issue #639.

Design:
  * PURE delta. A tool surfaces only when its upstream version differs from the
    last version we reported. The horizon is "since I last told you" (unbounded),
    NOT a fixed time window -- so a quiet tool stays quiet forever, and a multi-day
    gap yields one clean delta instead of a stale repeat.
  * First sighting seeds the watermark silently (no first-run flood).
  * Guards (frozen / watch_only / max_version) prevent false "you're behind" noise
    for intentionally-pinned or regression-risk deps.

Run:
    conda activate snakemake
    python tools/news/poll_releases.py --role developer

Network (urllib) and ``gh`` calls are dependency-injected so the unit tests run
fully offline. Only the IO helpers import PyYAML (lazily), so the pure-logic
functions import without it.
"""

import argparse
import json
import re
import subprocess
import sys
import urllib.request
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_BASKET = HERE / "tools.yaml"
REPO = "Jin-HoMLee/splice-neoepitope-pipeline"

# feed_types with no machine-readable feed (human watch only)
MANUAL_FEEDS = {"manual", "docker_image"}


# --------------------------------------------------------------------------- #
# Pure logic (no network, no yaml) -- this is the unit-tested core.
# --------------------------------------------------------------------------- #
def version_tuple(v):
    """Parse a version string into a comparable tuple of ints.

    Strips a leading ``v`` and any local-version segment (``+cu126``), then takes
    the leading integer of each dot-chunk. Good enough for cap checks and
    equality on the version strings we track (it is NOT a full PEP 440 parser).
    """
    main = str(v).split("+", 1)[0].lstrip("vV").strip()
    parts = []
    for chunk in main.split("."):
        digits = ""
        for ch in chunk:
            if ch.isdigit():
                digits += ch
            else:
                break
        parts.append(int(digits) if digits else 0)
    return tuple(parts)


def version_gt(a, b):
    """True if version ``a`` is strictly greater than ``b``."""
    return version_tuple(a) > version_tuple(b)


def apply_guards(tool, latest):
    """Decide how a fetched version should be treated.

    Returns ``(surfaced_version, kind)`` where kind is ``"bump"`` or ``"watch"``,
    or ``(None, reason)`` when the version must be suppressed
    (reason in: ``fetch-failed``, ``frozen``, ``above-cap``).
    """
    if latest is None:
        return None, "fetch-failed"
    if tool.get("frozen"):
        return None, "frozen"
    cap = tool.get("max_version")
    if cap and version_gt(latest, cap):
        return None, "above-cap"
    if tool.get("watch_only"):
        return latest, "watch"
    return latest, "bump"


def compute_delta(tool, latest, watermark_version):
    """Compute the briefing record for one tool.

    Returns a dict describing the outcome, or ``None`` when there is nothing to
    record at all (suppressed by a guard, or fetch failed). The returned dict
    carries ``surface`` (whether to show it in the briefing) and ``watermark``
    (the version to persist, or ``None`` to leave the watermark untouched).
    """
    name = tool["tool"]
    surfaced, kind = apply_guards(tool, latest)

    if surfaced is None:
        # Guard-suppressed or fetch failure: never surface, never touch watermark.
        return {
            "tool": name, "kind": kind, "surface": False, "watermark": None,
            "from": watermark_version, "to": latest, "suppressed": True,
        }

    if watermark_version is None:
        # First sighting: seed the baseline silently.
        return {
            "tool": name, "kind": "baseline", "surface": False, "watermark": surfaced,
            "from": None, "to": surfaced, "suppressed": False,
        }

    if surfaced != watermark_version:
        # Genuine delta: surface once. Pure-delta -> surfaced even if tracked.
        return {
            "tool": name, "kind": kind, "surface": True, "watermark": surfaced,
            "from": watermark_version, "to": surfaced, "suppressed": False,
        }

    # Unchanged since last report: muted (this is where the repetition dies).
    return {
        "tool": name, "kind": kind, "surface": False, "watermark": surfaced,
        "from": watermark_version, "to": surfaced, "suppressed": False,
    }


def parse_cu126_index(html):
    """Extract the highest ``torch==X.Y.Z+cu126`` version from the wheel index HTML.

    Returns the version string WITH the ``+cu126`` tag preserved, or None.
    """
    versions = re.findall(r"torch-(\d+\.\d+\.\d+)\+cu126", html or "")
    if not versions:
        return None
    best = max(versions, key=version_tuple)
    return f"{best}+cu126"


# --------------------------------------------------------------------------- #
# Network / gh IO (injectable; not exercised by the default offline test run).
# --------------------------------------------------------------------------- #
def _http_get_json(url, timeout=10):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def _http_get_text(url, timeout=10):
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def _gh_runner(args):
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, check=True
    ).stdout


def fetch_pypi(pkg, get_json=_http_get_json):
    return get_json(f"https://pypi.org/pypi/{pkg}/json")["info"]["version"]


def fetch_github(repo, gh=_gh_runner):
    tag = gh(["api", f"repos/{repo}/releases/latest", "--jq", ".tag_name"]).strip()
    return tag.lstrip("vV") or None


def fetch_bioconda(pkg, get_json=_http_get_json):
    return get_json(f"https://api.anaconda.org/package/bioconda/{pkg}").get("latest_version")


def fetch_pytorch_cu126(_feed, get_text=_http_get_text):
    return parse_cu126_index(get_text("https://download.pytorch.org/whl/cu126/torch/"))


def fetch_latest(tool, *, get_json=_http_get_json, get_text=_http_get_text, gh=_gh_runner):
    """Dispatch a single tool to its feed fetcher. Returns version str or None.

    Never raises: a fetch failure becomes ``None`` (-> guard ``fetch-failed``),
    so one flaky feed never aborts the whole briefing.
    """
    ft = tool.get("feed_type")
    if ft in MANUAL_FEEDS:
        return None
    try:
        if ft == "pypi":
            return fetch_pypi(tool["feed"], get_json=get_json)
        if ft == "github":
            return fetch_github(tool["feed"], gh=gh)
        if ft == "bioconda":
            return fetch_bioconda(tool["feed"], get_json=get_json)
        if ft == "pytorch_cu126":
            return fetch_pytorch_cu126(tool["feed"], get_text=get_text)
    except Exception as exc:  # noqa: BLE001 - intentional fail-soft per feed
        print(f"  ! fetch failed for {tool['tool']} ({ft}): {exc}", file=sys.stderr)
        return None
    return None


def find_tracking_issue(topic, gh=_gh_runner):
    """Return the number of an open Issue matching ``topic``, or None. Fails soft."""
    try:
        out = gh([
            "issue", "list", "--repo", REPO, "--state", "open",
            "--search", topic, "--limit", "1", "--json", "number",
        ])
        items = json.loads(out)
        return items[0]["number"] if items else None
    except Exception:  # noqa: BLE001 - annotation is best-effort
        return None


# --------------------------------------------------------------------------- #
# YAML IO (lazy import so pure-logic callers don't need PyYAML).
# --------------------------------------------------------------------------- #
def _yaml():
    import yaml  # noqa: PLC0415 - lazy on purpose
    return yaml


def load_basket(path=DEFAULT_BASKET):
    with open(path) as f:
        return _yaml().safe_load(f)


def watermark_path(role, here=HERE):
    return Path(here) / f"watermark_{role}.yaml"


def load_watermark(path):
    p = Path(path)
    if not p.exists():
        return {"tools": {}, "last_briefing": None}
    data = _yaml().safe_load(p.read_text()) or {}
    data.setdefault("tools", {})
    data.setdefault("last_briefing", None)
    return data


def save_watermark(path, data):
    with open(path, "w") as f:
        _yaml().safe_dump(data, f, sort_keys=True, default_flow_style=False)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(basket, watermark, *, check_tracked=True, fetcher=fetch_latest, tracker=find_tracking_issue):
    """Process the basket against the watermark. Returns (surfaced, updated_watermark).

    ``surfaced`` is the list of records to show the user (genuine deltas + watch
    notes). The watermark dict is updated in place and also returned.
    """
    wm_tools = watermark.setdefault("tools", {})
    surfaced = []

    for tool in basket.get("software", []) + basket.get("reference_data", []):
        name = tool["tool"]
        latest = fetcher(tool)
        rec = compute_delta(tool, latest, wm_tools.get(name))
        if rec is None:
            continue
        if rec["watermark"] is not None:
            wm_tools[name] = rec["watermark"]
        if rec["surface"]:
            if check_tracked and rec["kind"] != "watch":
                rec["tracking_issue"] = tracker(name)
            surfaced.append(rec)

    return surfaced, watermark


def format_briefing(surfaced, basket):
    """Render the surfaced deltas + the static watch-list as a text briefing."""
    lines = []
    if surfaced:
        lines.append("## Dependency deltas (since last briefing)")
        for rec in surfaced:
            tag = "watch" if rec["kind"] == "watch" else "→"
            arrow = f"{rec['from']} {tag} {rec['to']}" if rec["from"] else rec["to"]
            issue = rec.get("tracking_issue")
            suffix = f"  (tracked: Issue #{issue})" if issue else ""
            lines.append(f"- **{rec['tool']}** — {arrow}{suffix}")
    else:
        lines.append("## Dependency deltas: none — every tracked dep is unchanged since last briefing.")

    watch = basket.get("watch", [])
    if watch:
        names = ", ".join(w["tool"] for w in watch)
        lines.append(f"\n_Watch-list (candidate-to-adopt, muted): {names}_")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Morning-news dependency release-feed poller (Issue #639).")
    ap.add_argument("--role", default="developer", help="role whose watermark to use (default: developer)")
    ap.add_argument("--basket", default=str(DEFAULT_BASKET), help="path to tools.yaml")
    ap.add_argument("--watermark", default=None, help="override watermark path (default: watermark_<role>.yaml)")
    ap.add_argument("--no-write", action="store_true", help="do not persist the updated watermark")
    ap.add_argument("--no-check-tracked", action="store_true", help="skip the open-Issue annotation lookup")
    ap.add_argument("--json", action="store_true", help="emit surfaced deltas as JSON instead of text")
    args = ap.parse_args(argv)

    basket = load_basket(args.basket)
    wm_path = args.watermark or watermark_path(args.role)
    watermark = load_watermark(wm_path)

    surfaced, watermark = run(basket, watermark, check_tracked=not args.no_check_tracked)
    watermark["last_briefing"] = date.today().isoformat()

    if args.json:
        print(json.dumps(surfaced, indent=2))
    else:
        print(format_briefing(surfaced, basket))

    if not args.no_write:
        save_watermark(wm_path, watermark)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
