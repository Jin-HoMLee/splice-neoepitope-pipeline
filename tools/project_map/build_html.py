#!/usr/bin/env python3
"""Assemble the self-contained project_map.html — idempotent.

project_map.html is both the template (hand-authored D3 / CSS / interaction code)
and the shipped artifact. This script inlines two things so the output needs **no
network at all** — open it straight off disk and it renders:

  * the graph data from graph.json (one <script id="graph-data"> block), and
  * the D3 v7 library from vendor/d3.v7.min.js (one <script id="d3-inlined">
    block, replacing the d3js.org CDN <script src>).

It is safe to run any number of times: existing graph-data / d3-inlined blocks and
the CDN tag are stripped first, and duplicate loadGraph() definitions left by an
older non-idempotent build are collapsed to one.

Usage:
    python tools/project_map/build_html.py
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
HTML_PATH = HERE / "project_map.html"
JSON_PATH = HERE / "graph.json"
D3_PATH = HERE / "vendor" / "d3.v7.min.js"

GRAPH_DATA_RE = re.compile(
    r'[ \t]*<script id="graph-data"[^>]*>.*?</script>\n?', re.DOTALL)
D3_INLINED_RE = re.compile(
    r'[ \t]*<script id="d3-inlined">.*?</script>\n?', re.DOTALL)
# The CDN tag we replace with the inlined library.
D3_CDN_RE = re.compile(
    r'[ \t]*<script src="https://d3js\.org/d3\.v7[^"]*"></script>\n?')


def collapse_duplicate_function(html, name):
    """Keep only the first of N byte-identical `function <name>() {...}` blocks."""
    needle = f"function {name}("
    spans, bodies = [], []
    idx = html.find(needle)
    while idx != -1:
        brace = html.find("{", idx)
        if brace == -1:
            break
        depth, end = 0, brace
        for i in range(brace, len(html)):
            if html[i] == "{":
                depth += 1
            elif html[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        spans.append((idx, end))
        bodies.append(html[idx:end])
        idx = html.find(needle, end)

    if len(spans) <= 1:
        return html, 0
    first_body, removed = bodies[0], 0
    for (start, end), body in reversed(list(zip(spans[1:], bodies[1:]))):
        if body == first_body:
            tail = end + 1 if html[end:end + 1] == "\n" else end
            html = html[:start] + html[tail:]
            removed += 1
    return html, removed


def build():
    for label, path in (("project_map.html", HTML_PATH), ("graph.json", JSON_PATH),
                        ("vendor/d3.v7.min.js", D3_PATH)):
        if not path.exists():
            hint = " — run extract_graph.py first" if label == "graph.json" else ""
            print(f"Error: {label} not found at {path}{hint}")
            return 1

    graph_data = json.loads(JSON_PATH.read_text())
    d3_src = D3_PATH.read_text()
    html = HTML_PATH.read_text()

    if "</head>" not in html:
        print("Error: no </head> in HTML — cannot embed")
        return 1

    # 1) Strip prior injected blocks + the CDN tag (idempotent; collapses dups).
    html, n_data = GRAPH_DATA_RE.subn("", html)
    html, _ = D3_INLINED_RE.subn("", html)
    html, n_cdn = D3_CDN_RE.subn("", html)

    # 2) Collapse duplicate loadGraph() definitions from older builds.
    html, n_fns = collapse_duplicate_function(html, "loadGraph")

    # 3) Inline D3 right after <head> so it's defined before any script runs.
    d3_block = f'<script id="d3-inlined">{d3_src}</script>\n'
    html = html.replace("<head>", "<head>\n" + d3_block, 1)

    # 4) Inline one fresh graph-data block before </head>.
    data_block = ('<script id="graph-data" type="application/json">'
                  + json.dumps(graph_data) + "</script>\n")
    html = html.replace("</head>", data_block + "</head>", 1)

    HTML_PATH.write_text(html)

    size_kb = len(html.encode("utf-8")) / 1024
    print(f"✅ Built {HTML_PATH.name} ({size_kb:.1f} KB) — fully self-contained, no network needed")
    print(f"   {len(graph_data['nodes'])} nodes, {len(graph_data['edges'])} edges + inlined D3 v7")
    if n_data > 1:
        print(f"   collapsed {n_data} duplicate data block(s) → 1")
    if n_cdn:
        print("   replaced d3js.org CDN <script> with vendored D3")
    if n_fns:
        print(f"   collapsed {n_fns} duplicate loadGraph() definition(s)")
    return 0


if __name__ == "__main__":
    sys.exit(build())
