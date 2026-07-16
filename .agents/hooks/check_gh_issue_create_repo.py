#!/usr/bin/env python3
"""Pre-flight hook: warn-and-confirm when `gh issue create` targets a repo that
disagrees with the Issue's shape (Issue #1083).

Board 9 aggregates two repos with independent Issue numbering:
`splice-neoepitope-pipeline` (project) and
`claude-personas-splice-neoepitope-pipeline` (personas/memory). That collision
keeps producing a misfile-on-create: a memory/MM-shaped Issue lands in the
project repo, or a pipeline-shaped Issue lands in personas.

Mechanism class (Issue #1150, ratified 2026-07-15): this is a **best-effort
convenience control, not a guarantee** - a fuzzy, pattern-based guard, not an
enforcement boundary. So it emits `permissionDecision: "ask"` (confirm-or-cancel),
NEVER `deny`: a false positive costs one keystroke, not a blocked create. It is
deliberately **narrow and precision-first** (web best practice: an over-eager
guard is worse than none - false positives destroy trust and get disabled), and
it **fails OPEN** on every uncertain path. The fire-log is the "violations
dashboard" that would justify widening the surface later.

Shape signal is deterministic-first: the `role:memory_manager` label is the spine
(-> personas); title/body keywords only corroborate, and only a clearly one-sided
signal that contradicts the target repo fires. Anything ambiguous, unlabelled, or
mixed allows.

Sibling of `check_gh_issue_develop_parent.py` (same command-start parsing via
`_shell_parse.normalize_command`, same fail-open gh I/O, same fire-log infra).

Escape hatch: set `CLAUDE_ALLOW_REPO_MISFILE=1` to bypass (verbatim cross-repo
filing is occasionally legitimate).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1083.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_PUNCT = _shell_parse.PUNCT
_CREATE_PREFIX = ("gh", "issue", "create")

# Repo NAMES (owner is always Jin-HoMLee for both). The guard compares by name.
PROJECT_REPO = "splice-neoepitope-pipeline"
PERSONAS_REPO = "claude-personas-splice-neoepitope-pipeline"
_REPO_FOR_SHAPE = {"project": PROJECT_REPO, "personas": PERSONAS_REPO}

# Shape keywords (lowercased substring match on title+body). Kept tight on
# purpose: precision over recall. Widen only from fire-log evidence.
PERSONAS_KW = (
    "memory manager", "role:memory_manager", "memory.md", "episode", "episodic",
    "post-it", "lab-notebook", "lab notebook", "feedback_", "shared/", "persona",
    "semantic memory", "drain",
)
PROJECT_KW = (
    "snakemake", "pipeline", "workflow", "alignment", "junction", "neoepitope",
    "regtools", "bedtools", "mhcflurry", "tcrdock", "hisat2", "star", "conda",
    "fastq", "chr22", "gencode",
)

# `gh issue create` flags whose value is the following token.
_VALUE_FLAGS = {
    "-t", "--title", "-b", "--body", "-F", "--body-file", "-l", "--label",
    "-a", "--assignee", "-m", "--milestone", "-p", "--project", "-R", "--repo",
    "-T", "--template",
}


# --- pure helpers (unit-tested, no I/O) ---


def create_args(cmd: str) -> list[str] | None:
    """Return the token list AFTER a real `gh issue create`, or None.

    Tokenizes with shlex over the NORMALIZED command (heredoc bodies stripped,
    unquoted newlines -> separators, Issue #1142) so a `gh issue create` buried
    in a quoted body - e.g. a comment documenting this hook - does NOT match.
    Only a `gh issue create` at a command start counts; args are collected up to
    the next shell separator. Untokenizable input fails safe -> None.
    """
    import shlex

    try:
        lex = shlex.shlex(_shell_parse.normalize_command(cmd),
                          posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None
    at_command_start = True
    for i, tok in enumerate(tokens):
        if tok and all(ch in _PUNCT for ch in tok):
            at_command_start = True
            continue
        if at_command_start and tuple(tokens[i:i + 3]) == _CREATE_PREFIX:
            rest = tokens[i + 3:]
            args: list[str] = []
            for t in rest:
                if t and all(ch in _PUNCT for ch in t):
                    break
                args.append(t)
            return args
        at_command_start = False
    return None


def parse_fields(args: list[str]) -> tuple[str, str, list[str]]:
    """Extract (title, body, labels) from `gh issue create` args.

    Handles `--flag value` and `--flag=value` forms; `--label` may repeat and/or
    carry a comma list. `--body-file` content is NOT read (fail-open: title +
    labels are signal enough); only an inline `--body` is inspected.
    """
    title: list[str] = []
    body: list[str] = []
    labels: list[str] = []
    i = 0
    while i < len(args):
        tok = args[i]
        nxt = args[i + 1] if i + 1 < len(args) else None
        if tok in ("-t", "--title") and nxt is not None:
            title.append(nxt); i += 2; continue
        if tok.startswith("--title="):
            title.append(tok.split("=", 1)[1]); i += 1; continue
        if tok in ("-b", "--body") and nxt is not None:
            body.append(nxt); i += 2; continue
        if tok.startswith("--body="):
            body.append(tok.split("=", 1)[1]); i += 1; continue
        if tok in ("-l", "--label") and nxt is not None:
            labels.extend(x.strip() for x in nxt.split(",") if x.strip()); i += 2; continue
        if tok.startswith("--label="):
            labels.extend(x.strip() for x in tok.split("=", 1)[1].split(",") if x.strip())
            i += 1; continue
        i += 1
    return " ".join(title), " ".join(body), labels


def repo_from_args(args: list[str]) -> tuple[str, str] | None:
    """Parse an explicit `-R`/`--repo [HOST/]OWNER/REPO`, else None (cwd fallback)."""
    val = None
    for i, tok in enumerate(args):
        if tok in ("-R", "--repo"):
            val = args[i + 1] if i + 1 < len(args) else None
            break
        if tok.startswith("-R=") or tok.startswith("--repo="):
            val = tok.split("=", 1)[1]
            break
    if not val:
        return None
    parts = val.split("/")
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None


def repo_signal(title: str, body: str, labels: list[str]) -> tuple[str | None, str]:
    """Infer the expected repo from Issue shape: 'personas' | 'project' | None.

    Deterministic-first + precision-tuned:
      - `role:memory_manager` label -> personas (the deterministic spine).
      - else a clearly one-sided keyword signal (>=2 hits on one side, 0 on the
        other) -> that side.
      - else None (ambiguous / mixed / unlabelled -> caller allows).
    """
    label_set = {ln.lower() for ln in labels}
    if "role:memory_manager" in label_set:
        return "personas", "role:memory_manager label"
    text = f"{title} {body}".lower()
    p_hits = [k for k in PERSONAS_KW if k in text]
    j_hits = [k for k in PROJECT_KW if k in text]
    if len(p_hits) >= 2 and not j_hits:
        return "personas", f"personas keywords {p_hits}"
    if len(j_hits) >= 2 and not p_hits:
        return "project", f"project keywords {j_hits}"
    return None, ""


def ask_payload(target_repo: str, expected_repo: str, reason: str) -> dict:
    """The PreToolUse ask (confirm-or-cancel) decision for a shape/target mismatch."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "ask",
            "permissionDecisionReason": (
                f"This `gh issue create` targets `{target_repo}`, but the Issue's "
                f"shape looks like it belongs in `{expected_repo}` ({reason}). "
                "Board 9 spans two repos with independent numbering, so a misfile "
                "lands the Issue in the wrong tracker. Confirm to create here, or "
                "cancel and re-target with `--repo Jin-HoMLee/"
                f"{expected_repo}`. This is a best-effort convenience check "
                "(Issue #1150); set CLAUDE_ALLOW_REPO_MISFILE=1 to silence it."
            ),
        }
    }


# --- gh I/O (fail-open) ---


def resolve_target_repo(args: list[str]) -> tuple[str, str] | None:
    """(owner, name) for the create target: explicit `-R` wins, else cwd repo."""
    explicit = repo_from_args(args)
    if explicit:
        return explicit
    try:
        res = subprocess.run(
            ["gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner"],
            capture_output=True, text=True, timeout=6, check=True,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return None
    owner, _, name = res.stdout.strip().partition("/")
    return (owner, name) if owner and name else None


def _log_fire(target: str, expected: str, reason: str) -> None:
    """Append one fire-log line on an ask (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_gh_issue_create_repo",
        "target": target,
        "expected": expected,
        "action": f"ask-repo-misfile:{reason}",
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    try:
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(line)
    except OSError:
        pass


# --- orchestration ---


def main() -> int:
    if os.environ.get("CLAUDE_ALLOW_REPO_MISFILE") == "1":
        return 0  # explicit escape hatch

    try:
        data = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (data.get("tool_input") or {}).get("command", "")
    args = create_args(cmd)
    if args is None:
        return 0  # not a `gh issue create`

    title, body, labels = parse_fields(args)
    expected_shape, reason = repo_signal(title, body, labels)
    if expected_shape is None:
        return 0  # ambiguous shape -> fail open

    target = resolve_target_repo(args)
    if target is None:
        return 0  # can't resolve target repo -> fail open
    _, target_name = target

    expected_name = _REPO_FOR_SHAPE[expected_shape]
    if target_name not in (PROJECT_REPO, PERSONAS_REPO):
        return 0  # target is a third repo -> out of scope, allow
    if target_name == expected_name:
        return 0  # shape and target agree -> allow

    _log_fire(target_name, expected_name, reason)
    print(json.dumps(ask_payload(target_name, expected_name, reason)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
