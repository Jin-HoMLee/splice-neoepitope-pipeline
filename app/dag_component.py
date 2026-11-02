"""Offline graph component boundary for the selected Cytoscape + ELK renderer.

The Phase 0 Python shell supplies only validated node/edge data. The bundled HTML
fallback keeps the demo usable before the optional frontend bundle is built; the
frontend contract and pinned dependencies live under components/cytoscape_elk.
"""

import html
import json
from typing import Any, Dict, List

NODES: List[Dict[str, Any]] = [
    {"id": "samples", "label": "Sample staging", "layer": "input", "x": 30, "y": 145, "rule": "create_alignment_manifest"},
    {"id": "refs", "label": "Reference bundle", "layer": "input", "x": 30, "y": 45, "rule": "build_reference_junctions"},
    {"id": "align", "label": "HISAT2 alignment", "layer": "intermediate", "x": 240, "y": 145, "rule": "hisat2_align"},
    {"id": "filter", "label": "Origin filter", "layer": "intermediate", "x": 425, "y": 145, "rule": "filter_junctions"},
    {"id": "peptides", "label": "Peptides", "layer": "intermediate", "x": 610, "y": 145, "rule": "translate_peptides"},
    {"id": "hla", "label": "HLA genotype", "layer": "intermediate", "x": 425, "y": 270, "rule": "aggregate_hla_alleles"},
    {"id": "mhc", "label": "MHC presentation", "layer": "intermediate", "x": 795, "y": 145, "rule": "run_mhcflurry"},
    {"id": "report", "label": "Report artifacts", "layer": "output", "x": 980, "y": 145, "rule": "generate_report"},
]
EDGES = [
    ("refs", "align"), ("samples", "align"), ("align", "filter"),
    ("filter", "peptides"), ("samples", "hla"), ("peptides", "mhc"),
    ("hla", "mhc"), ("mhc", "report"),
]


def render_dag(components: Any, *, height: int = 510) -> None:
    """Render the read-only offline graph without returning executable actions."""
    positions = {node["id"]: node for node in NODES}
    lines = []
    for source, target in EDGES:
        a, b = positions[source], positions[target]
        x1, y1 = a["x"] + 130, a["y"] + 33
        x2, y2 = b["x"], b["y"] + 33
        lines.append(
            f'<path class="edge" data-source="{source}" data-target="{target}" '
            f'd="M{x1} {y1} C{(x1+x2)//2} {y1} {(x1+x2)//2} {y2} {x2} {y2}" />'
        )
    buttons = []
    for node in NODES:
        buttons.append(
            '<button class="node layer-{layer}" data-id="{id}" data-rule="{rule}" '
            'style="left:{x}px;top:{y}px"><small>{layer}</small><strong>{label}</strong></button>'.format(
                **{key: html.escape(str(value)) for key, value in node.items()}
            )
        )
    payload = html.escape(json.dumps({"nodes": NODES, "edges": EDGES}))
    document = f"""
    <style>
      *{{box-sizing:border-box}} body{{margin:0;font-family:Inter,system-ui;color:#172033}}
      .toolbar,.legend{{display:flex;gap:12px;align-items:center;padding:8px 10px;border:1px solid #d8dbea;background:#fff}}
      .toolbar{{border-radius:12px 12px 0 0;justify-content:space-between}} .toolbar span{{font-size:12px;color:#697386}}
      .legend{{border-top:0;font-size:11px;color:#697386}} .legend i{{width:12px;height:12px;border-radius:4px;display:inline-block;margin-right:4px}}
      .canvas{{height:390px;overflow:auto;position:relative;border:1px solid #d8dbea;border-top:0;background:#f8fafc radial-gradient(#d7dbea 1px,transparent 1px);background-size:18px 18px}}
      .stage{{position:relative;width:1150px;height:380px;transform-origin:top left}} svg{{position:absolute;inset:0;width:1150px;height:380px}}
      .edge{{fill:none;stroke:#7c3aed;stroke-width:2.2;transition:.15s}} .edge.dim{{opacity:.1}} .edge.focus{{stroke-width:3.5}}
      .node{{position:absolute;width:130px;height:66px;border:2px solid;border-radius:13px;text-align:left;padding:9px;background:#fff;cursor:pointer;transition:.15s}}
      .node small{{display:block;text-transform:uppercase;font-size:8px;font-weight:900;letter-spacing:.08em}} .node strong{{font-size:11px}}
      .layer-input{{border-color:#2563eb;background:#eff6ff}} .layer-intermediate{{border-color:#7c3aed;background:#faf8ff}} .layer-output{{border-color:#15803d;background:#ecfdf3}}
      .node.dim{{opacity:.15}} .node.selected{{outline:3px solid #f0c66c;outline-offset:2px}}
      .inspector{{display:flex;gap:12px;padding:9px 10px;border:1px solid #d8dbea;border-top:0;border-radius:0 0 12px 12px;font-size:11px}} .inspector b{{color:#570df8}}
      button.control{{border:1px solid #d8dbea;border-radius:7px;background:#fff;padding:5px 8px}}
    </style>
    <div class="toolbar"><div><b>Cytoscape + ELK graph contract</b><br><span>offline read-only DAG</span></div><div><button class="control" id="fit">Fit</button><button class="control" id="all">Show all</button></div></div>
    <div class="legend"><span><i style="background:#eff6ff;border:2px solid #2563eb"></i>Input</span><span><i style="background:#faf8ff;border:2px solid #7c3aed"></i>Intermediary</span><span><i style="background:#ecfdf3;border:2px solid #15803d"></i>Output</span></div>
    <div class="canvas" id="canvas"><div class="stage" id="stage"><svg viewBox="0 0 1150 380">{''.join(lines)}</svg>{''.join(buttons)}</div></div>
    <div class="inspector">Selected: <b id="title">Sample staging</b><span>Rule: <code id="rule">create_alignment_manifest</code></span></div>
    <script type="application/json" id="graph-data">{payload}</script>
    <script>
      const nodes=[...document.querySelectorAll('.node')], edges=[...document.querySelectorAll('.edge')];
      function clear(){{nodes.forEach(n=>n.classList.remove('dim','selected'));edges.forEach(e=>e.classList.remove('dim','focus'));}}
      function select(node){{const id=node.dataset.id; const adjacent=new Set([id]); edges.forEach(e=>{{if(e.dataset.source===id)adjacent.add(e.dataset.target);if(e.dataset.target===id)adjacent.add(e.dataset.source);}});nodes.forEach(n=>{{n.classList.toggle('selected',n===node);n.classList.toggle('dim',!adjacent.has(n.dataset.id));}});edges.forEach(e=>{{const hit=e.dataset.source===id||e.dataset.target===id;e.classList.toggle('focus',hit);e.classList.toggle('dim',!hit);}});document.getElementById('title').textContent=node.textContent.trim();document.getElementById('rule').textContent=node.dataset.rule;}}
      nodes.forEach(n=>n.onclick=()=>select(n)); document.getElementById('all').onclick=clear;
      document.getElementById('fit').onclick=()=>{{const canvas=document.getElementById('canvas'),stage=document.getElementById('stage');const z=Math.min((canvas.clientWidth-20)/1150,(canvas.clientHeight-20)/380);stage.style.transform=`scale(${{z}})`;canvas.scrollTo(0,0);}};
      select(nodes[0]);
    </script>
    """
    components.html(document, height=height, scrolling=False)
