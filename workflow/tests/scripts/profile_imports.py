"""Profile cold-subprocess import time for each test file and each project script.

Diagnostic for [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364):
local pytest collection hangs under memory pressure because module-level imports
get amplified by OS paging on the M1 8 GB dev box. This script establishes the
static cost surface so future refactors (e.g. lazy-pandas) can be measured.

Usage:
    workflow/tests/.venv/bin/python workflow/tests/scripts/profile_imports.py

Each entry is timed by exec-ing the file's module-level body (imports +
defs/classes only — `if __name__ == "__main__":` and any trailing statements
are stripped) in a fresh subprocess. This isolates *import* cost from test
execution, and matches what pytest pays during the collection phase.
"""
from __future__ import annotations

import ast
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
TESTS_DIR = REPO_ROOT / "workflow" / "tests"
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
PY = REPO_ROOT / "workflow" / "tests" / ".venv" / "bin" / "python"


def module_level_imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text())
    out: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            out.extend(n.name.split(".")[0] for n in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            out.append(node.module.split(".")[0])
    return sorted(set(out))


def truncate_to_imports_and_defs(src: str) -> str:
    """Return the module's top-level imports + defs/classes/assigns, stopping
    at the first `if __name__ == "__main__":` block or any other statement."""
    tree = ast.parse(src)
    body = []
    for node in tree.body:
        if isinstance(node, ast.If) and ast.unparse(node.test).startswith("__name__"):
            break
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                             ast.ClassDef, ast.AsyncFunctionDef, ast.Expr,
                             ast.Assign, ast.AnnAssign)):
            body.append(node)
            continue
        break
    return ast.unparse(ast.Module(body=body, type_ignores=[]))


def cold_import_time(path: Path, extra_path: Path) -> tuple[float, bool]:
    truncated = truncate_to_imports_and_defs(path.read_text())
    runner = (
        "import sys; "
        f"sys.path.insert(0, {str(extra_path)!r}); "
        f"exec(compile({truncated!r}, {str(path)!r}, 'exec'))"
    )
    t0 = time.perf_counter()
    r = subprocess.run([str(PY), "-c", runner], capture_output=True, text=True, timeout=180)
    return time.perf_counter() - t0, r.returncode == 0


def profile_dir(label: str, files: list[Path], extra_path: Path) -> None:
    print(f"\n## {label}")
    rows: list[tuple[str, float, bool, list[str]]] = []
    for f in files:
        dt, ok = cold_import_time(f, extra_path)
        imps = module_level_imports(f)
        rows.append((f.name, dt, ok, imps))
        print(f"  {'OK ' if ok else 'ERR'} {dt:6.2f}s  {f.name}", file=sys.stderr)
    rows.sort(key=lambda r: -r[1])
    print(f"\n| {'file':<42s} | {'sec':>6s} | module-level imports |")
    print(f"| {'-'*42} | {'-'*6} | --- |")
    for name, dt, ok, imps in rows:
        marker = "" if ok else " ⚠️"
        print(f"| `{name}` | {dt:6.2f}{marker} | {', '.join(imps)} |")


def main() -> int:
    if not PY.exists():
        print(f"ERROR: {PY} not found. Set up the test venv first (see workflow/tests/README.md).",
              file=sys.stderr)
        return 1
    test_files = sorted(TESTS_DIR.glob("test_*.py"))
    script_files = sorted(SCRIPTS_DIR.glob("*.py"))
    profile_dir("Per-test-file cold module-load (seconds)", test_files, SCRIPTS_DIR)
    profile_dir("Per-project-script cold module-load (seconds)", script_files, SCRIPTS_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
