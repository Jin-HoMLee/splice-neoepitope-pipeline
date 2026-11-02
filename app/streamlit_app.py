"""Local-only Streamlit shell for pipeline artifact exploration and run preview."""

from pathlib import Path
import sys
from typing import Any, Dict, List

import streamlit as st
import streamlit.components.v1 as components

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.command_builder import build_snakemake_argv, display_command
from app.config_contract import ConfigValidationError, LocalRunConfig
from app.dag_component import render_dag
from app.demo_data import load_demo

st.set_page_config(page_title="Splice Workbench", page_icon="🧬", layout="wide")

st.markdown(
    """
<style>
[data-testid="stHeader"] { background: rgba(244,247,251,.94); }
.block-container { padding-top: 1.2rem; max-width: 1440px; }
.demo { border:1px solid #d9c676; background:#fff8dc; padding:.75rem 1rem; border-radius:14px; }
.card { border:1px solid #d9dceb; background:#fff; padding:1rem; border-radius:16px; }
.safe { border-left:4px solid #0d9488; background:#ecfdf5; padding:.8rem; }
.small { color:#697386; font-size:.85rem; }
</style>
""",
    unsafe_allow_html=True,
)

data = load_demo()
metadata = data["metadata"]

head_left, head_right = st.columns([3, 1])
with head_left:
    st.title("Splice Workbench")
    st.caption("Local report-aware adapter over the existing Snakemake workflow")
with head_right:
    mode = st.segmented_control("Mode", ["Explore", "Configure"], default="Explore")

st.markdown(
    f'<div class="demo"><b>Offline synthetic demo</b> · {metadata["prediction_status"]} · '
    "no upload · no network · no clinical use</div>",
    unsafe_allow_html=True,
)

if mode == "Configure":
    st.header("Configure local run")
    st.info("Phase 0 validates and previews only. It does not launch Snakemake.")
    with st.form("run-config"):
        manifest_text = st.text_input("Local sample manifest", "config/samples/patient_001_test.tsv")
        output_text = st.text_input("Run output root", ".streamlit/runs/preview")
        depth = st.radio("Junction-read preset", [10, 5], horizontal=True, format_func=lambda x: "10 · recommended" if x == 10 else "5 · sensitivity")
        acknowledge = st.checkbox("I acknowledge that 5 reads is sensitivity-oriented and can increase false positives", disabled=depth != 5)
        threads = st.slider("Threads", 1, 8, 4)
        st.text_input("Aligner", "HISAT2 · native local smoke", disabled=True)
        submitted = st.form_submit_button("Validate and preview")
    if submitted:
        config = LocalRunConfig(
            manifest=REPO_ROOT / manifest_text,
            output_root=REPO_ROOT / output_text,
            depth=depth,
            threads=threads,
            sensitivity_acknowledged=acknowledge,
        )
        try:
            config.validate(require_paths=True)
        except ConfigValidationError as exc:
            st.error(str(exc))
        else:
            overlay = config.overlay()
            overlay_path = config.output_root / "config.overlay.yaml"
            target = config.output_root / "results" / "demo" / "reports" / "report.html"
            argv = build_snakemake_argv(
                configfiles=[REPO_ROOT / "config/config.yaml", overlay_path],
                target=target,
                cores=threads,
                dry_run=True,
            )
            st.success("Configuration is within the Phase 0 allowlist")
            st.subheader("Effective overlay preview")
            st.json(overlay)
            st.subheader("Exact dry-run command preview")
            st.code(display_command(argv), language="shell")
    st.stop()

section = st.radio(
    "Workspace",
    ["Overview", "Candidates", "3D structure", "Provenance", "Run history"],
    horizontal=True,
    label_visibility="collapsed",
)

if section == "Overview":
    st.header("Pipeline overview and filtering funnel")
    cols = st.columns(4)
    for column, (label, value) in zip(cols, [("Raw junctions", 124), ("Tumor-exclusive", 38), ("Translated peptides", 16), ("Illustrative presenters", 6)]):
        column.metric(label, value)
    scope, depth_view, display = st.columns(3)
    with scope:
        st.selectbox("Evidence scope", ["All evidence", "Tumor arm", "Normal arm"])
    with depth_view:
        st.selectbox("Depth view", ["10 reads · recommended", "5 reads · sensitivity"])
    with display:
        st.selectbox("Display", ["Counts", "Percent of raw"])
    st.bar_chart(data["funnel"], x="stage", y="count", horizontal=True, color="#570df8")
    selected_stage = st.selectbox("Inspect stage", [row["stage"] for row in data["funnel"]])
    stage = next(row for row in data["funnel"] if row["stage"] == selected_stage)
    st.markdown(f'<div class="safe"><b>{selected_stage}</b><br>Rule: <code>{stage["rule"]}</code><br>Artifact: <code>{stage["artifact"]}</code></div>', unsafe_allow_html=True)
    st.subheader("Junction origin policy")
    st.dataframe(data["origins"], width="stretch", hide_index=True)

elif section == "Candidates":
    st.header("Presentation results and candidate context")
    query = st.text_input("Search candidate, sequence, or allele")
    classes = st.multiselect("Presentation class", ["strong", "weak", "non"], default=["strong", "weak", "non"])
    candidates: List[Dict[str, Any]] = [row for row in data["candidates"] if row["presentation_class"] in classes]
    if query:
        needle = query.lower()
        candidates = [row for row in candidates if needle in " ".join(str(value) for value in row.values()).lower()]
    event = st.dataframe(
        candidates,
        width="stretch",
        hide_index=True,
        selection_mode="single-row",
        on_select="rerun",
        column_order=("candidate", "presentation_class", "best_allele", "presentation_percentile", "genotype_presentation_score", "ic50_nM"),
    )
    selected_index = event.selection.rows[0] if event.selection.rows else 0
    if candidates:
        selected = candidates[min(selected_index, len(candidates) - 1)]
        a, b, c = st.columns(3)
        a.metric("Selected", selected["candidate"])
        b.metric("Presentation percentile", selected["presentation_percentile"])
        c.metric("Best allele", selected["best_allele"])
        with st.expander("Evidence, sequence, and model context", expanded=True):
            st.code(selected["sequence"])
            st.write("Presentation percentile is primary. GPS and IC50 are secondary context.")

elif section == "3D structure":
    st.header("TCR-peptide-MHC structure")
    st.warning("Structure unavailable in the offline demo. No PDB is bundled and no TCRdock/GPU execution is offered.")
    candidate = st.selectbox("Candidate context", [row["candidate"] for row in data["candidates"][:3]])
    representation = st.segmented_control("Representation", ["Cartoon", "Lines", "Surface"], default="Cartoon")
    st.markdown('<div class="card" style="height:300px;display:grid;place-items:center;background:#17243a;color:#dbe8f4"><div><b>Mol* local-PDB slot</b><br><span class="small">Awaiting top_candidate.pdb</span></div></div>', unsafe_allow_html=True)
    st.caption(f"Selected {candidate} · {representation} · chain focus and pLDDT/pTM/PAE appear only when emitted")

elif section == "Provenance":
    st.header("Provenance explorer")
    st.json({"metadata": metadata, "provenance": data["provenance"]})
    tabs = st.tabs(["Output manifest", "Effective config", "Exact command", "Environment"])
    with tabs[0]:
        st.dataframe([
            {"artifact": "junctions", "path": "results/demo-v2/junctions/novel_junctions.tsv", "state": "synthetic"},
            {"artifact": "presentation", "path": "results/demo-v2/mhcflurry/presentation_calibrated.tsv", "state": "synthetic"},
            {"artifact": "structure", "path": "results/demo-v2/reports/report_3d_structure.tsv", "state": "not available"},
        ], width="stretch", hide_index=True)
    with tabs[1]:
        st.code("filtering:\n  min_junction_reads: 10\nalignment:\n  aligner: hisat2\ntcrdock:\n  enabled: false", language="yaml")
    with tabs[2]:
        st.code("snakemake --cores 4 --use-conda --configfile config/config.yaml overlay.yaml -- results/demo-v2/reports/report.html", language="shell")
    with tabs[3]:
        st.write("Local CPU-only · network disabled in demo · synthetic fixture")

else:
    st.header("Run history · Cytoscape + ELK DAG")
    run, status = st.columns(2)
    with run:
        st.selectbox("Run", ["demo-v2 · complete", "chr22 sensitivity · complete", "local draft · interrupted"])
    with status:
        st.selectbox("Path status", ["All paths", "Active only", "Optional only"])
    render_dag(components)

st.divider()
st.caption("Prediction-only artifact review. This local adapter is not a clinical, diagnostic, or treatment system.")
