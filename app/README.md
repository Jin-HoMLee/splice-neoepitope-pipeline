# Splice Workbench Phase 0

A local, read-only Streamlit adapter over the existing Snakemake pipeline artifacts.

## Scope

Phase 0 provides:

- an offline synthetic demo with no network or subprocess;
- Explore and Configure modes;
- drillable funnel, candidate, structure-unavailable, provenance, and DAG surfaces;
- a typed allowlist for the 10-read default and acknowledged 5-read sensitivity preset;
- exact Snakemake dry-run argv preview using one `--configfile` list and `--` before the target;
- native HISAT2 as the only executable-profile preview;
- no browser uploads, arbitrary YAML/shell, GPU, cloud, remote storage, or TCRdock execution.

The synthetic presenter rows are illustrative UI fixtures. They are not MHCflurry predictions and make no clinical claim.

## Run

Create a dedicated UI environment, separate from the Snakemake and per-rule environments:

```bash
python -m venv app/.venv
app/.venv/bin/pip install -r app/requirements.txt
app/.venv/bin/streamlit run app/streamlit_app.py --server.address 127.0.0.1
```

The demo is read-only. Configure mode validates and previews; it does not launch Snakemake in Phase 0.

## Test

```bash
workflow/tests/.venv/bin/python -m pytest app/tests
```

The tests cover pure contracts and do not require Streamlit.

## Cytoscape and ELK boundary

`app/dag_component.py` is the offline component boundary and Phase 0 fallback. It receives validated nodes and edges only and implements the approved read-only graph behavior. The production custom-component bundle must vendor and version-pin Cytoscape.js and ELK.js locally before replacing the fallback; runtime CDN loading is not allowed.
