"""Unit + invariant tests for the project-map graph extractor.

Run with the project's pytest venv:
    workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py -v

CI: collected by the `ci-tools-pytest` job in `.github/workflows/tests.yml`
(Issue #713). `build_graph()` discovers files via `git ls-files` (Issue #780),
so the atlas is reproducible across working-tree states: a local clone with
populated gitignored artifacts (`references/`, `logs/`, `data/`, `results/`,
`indices/`) produces the identical graph as a clean checkout, and the suite
passes regardless of whether the pipeline has been run locally. (Running it
requires a git working tree.)
"""
import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location(
    "extract_graph", pathlib.Path(__file__).parent / "extract_graph.py")
eg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(eg)


def test_classify_path_relative_paths():
    # These are the relative paths os.walk actually passes (no leading slash).
    assert eg.classify_path("research/manuscript/METHODS.md") == "manuscript"
    assert eg.classify_path("docs/google_cloud_guide.md") == "doc"
    assert eg.classify_path("research/evals/issue_188_boltz2/slides.qmd") == "slide"
    assert eg.classify_path("research/lab_notebook/pm.md") == "labnotebook"
    assert eg.classify_path("workflow/scripts/run_mhcflurry.py") == "script"
    assert eg.classify_path("workflow/envs/python.yaml") == "env"
    assert eg.classify_path("config/config.yaml") == "config"
    assert eg.classify_path("tools/ci/closure_audit.py") == "ci"
    assert eg.classify_path("research/notebooks/patient_001_results.ipynb") == "notebook"


def test_resolve_link():
    assert eg.resolve_link("docs/a.md", "b.md") == "docs/b.md"
    assert eg.resolve_link("docs/a.md", "sub/c.md#section") == "docs/sub/c.md"
    assert eg.resolve_link("research/evals/x/README.md",
                           "../../slides/nature.csl") == "research/slides/nature.csl"
    assert eg.resolve_link("docs/a.md", "https://x.com/y.md") is None  # external
    assert eg.resolve_link("docs/a.md", "#anchor-only") is None


def test_research_half_is_connected():
    g = eg.build_graph()
    ids = {n["id"] for n in g["nodes"]}
    research_ids = {n["id"] for n in g["nodes"] if n["group"] == "research"}
    real = [e for e in g["edges"]
            if e["type"] != "contained_in"
            and (e["source"] in research_ids or e["target"] in research_ids)]
    assert len(real) > 10, f"research still nearly disconnected: {len(real)} real edges"
    # every edge endpoint is a real node (D3 forceLink safety)
    assert all(e["source"] in ids and e["target"] in ids for e in g["edges"])


def test_issue_tagging():
    g = eg.build_graph()
    by_path = {n.get("path"): n for n in g["nodes"] if n.get("path")}
    hit = [n for p, n in by_path.items() if "issue_218" in p and n.get("issue") == 218]
    assert hit, "expected nodes under issue_218_* to carry issue==218"
    # nodes outside any issue folder must not be tagged
    assert by_path.get("workflow/scripts/run_mhcflurry.py", {}).get("issue") is None


def test_resource_blob_is_gone():
    g = eg.build_graph()
    by_group = {}
    for n in g["nodes"]:
        by_group.setdefault(n["group"], []).append(n["type"])
    for grp, types in by_group.items():
        frac = types.count("resource") / len(types)
        assert frac < 0.40, f"group '{grp}' is {frac:.0%} unclassified 'resource'"


def test_gitignored_runtime_dirs_contribute_no_nodes():
    """Regression for #780: discovery must exclude gitignored runtime/data dirs
    (`references/`, `results/`, `logs/`, `data/`, `indices/`) so the atlas is
    identical between a clean checkout and a clone that has run the pipeline.
    All five are fully gitignored (zero tracked files), so a git-tracked walk
    must emit no node — dir or file — rooted under any of them."""
    g = eg.build_graph()
    gitignored_roots = {"references", "results", "logs", "data", "indices"}
    offenders = []
    for n in g["nodes"]:
        path = (n.get("path") or "").replace("\\", "/")
        if path and path.split("/", 1)[0] in gitignored_roots:
            offenders.append(path)
    assert not offenders, (
        f"gitignored runtime dirs leaked {len(offenders)} node(s) into the atlas: "
        f"{sorted(offenders)[:10]}")
