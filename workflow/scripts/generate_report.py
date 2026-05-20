#!/usr/bin/env python3
"""generate_report.py — Generate a summary HTML report for pipeline results.

Produces a self-contained HTML page showing:
  1. Junction origin summary — counts of tumor_exclusive vs normal_shared
     junctions per sample, so results can be quickly verified.
  2. Top strong presenters — the highest-presentation predicted neoepitopes
     (presentation_percentile <= strong threshold).

Usage (standalone):
  python generate_report.py \\
      --novel-junctions results/junctions/local/novel_junctions.tsv \\
      --predictions results/predictions/local/mhc_presentation.tsv \\
      --output results/reports/local/report.html

Usage (Snakemake):
  Called automatically by the ``generate_report`` rule.
"""

from __future__ import annotations

import argparse
import html as html_mod
import json
import logging
from pathlib import Path
from typing import Any

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

_PIPELINE_DIAGRAM = """\
<div class="mermaid">
flowchart TD
    reads(["RNA-seq reads<br/>tumor + normal"])
    align["HISAT2 alignment<br/>junction extraction with regtools"]
    hla["OptiType HLA typing<br/>normal sample - optional"]
    ref["GENCODE reference filter<br/>remove annotated junctions"]
    cmp["Normal sample comparison"]
    ps["patient_specific<br/>also in normal tissue"]
    ts["tumor_specific<br/>absent from normal tissue"]
    asm["Contig assembly<br/>50 nt per junction"]
    trans["In-silico translation<br/>3 reading frames - 9-mers"]
    pred["MHCflurry prediction<br/>IC50 + presentation_score per peptide (best allele)"]
    out(["Neoepitope candidates<br/>strong / weak presenters"])

    reads --> align
    align --> hla
    align --> ref
    ref -->|unannotated junctions| cmp
    cmp --> ps
    cmp --> ts
    ts --> asm
    asm --> trans
    trans --> pred
    hla -->|patient-specific alleles| pred
    pred --> out

    style reads fill:#eafaf1,stroke:#27ae60
    style out   fill:#eafaf1,stroke:#27ae60
    style hla   fill:#f5eef8,stroke:#8e44ad
    style ps    fill:#fdecea,stroke:#e74c3c
    style ts    fill:#eafaf1,stroke:#27ae60
    style ref   fill:#fef9e7,stroke:#f39c12
</div>
"""

_MOLSTAR_VIEWER = """\
<div id="molstar-container" style="width:100%;height:500px;position:relative;border:1px solid #ddd;border-radius:6px;overflow:hidden;margin:1em 0;"></div>
<script type="text/javascript">
  // Mol* viewer — loads the inlined PDB and renders the TCR-peptide-MHC complex.
  // Mol* is MIT-licensed and loaded from the official CDN (no registration required).
  // Version is pinned to avoid breaking API changes on unpkg "latest" redirects.
  document.addEventListener("DOMContentLoaded", function() {{
    molstar.Viewer.create("molstar-container", {{
      layoutIsExpanded: false,
      layoutShowControls: true,
      layoutShowRemoteState: false,
      layoutShowSequence: true,
      layoutShowLog: false,
      layoutShowLeftPanel: true,
    }}).then(function(viewer) {{
      var pdbData = {pdb_data};
      viewer.loadStructureFromData(pdbData, "pdb", {{ dataLabel: "TCR-pMHC complex" }});
    }});
  }});
</script>
"""

# Chain ID → biological component name (standard TCRdock output order)
_VIEWER_CHAIN_NAMES = {
    "A": "MHC heavy chain",
    "B": "Peptide",
    "C": "TCR \u03b1-chain",
    "D": "TCR \u03b2-chain",
}


def _extract_chain_ids(pdb_text: str) -> list[str]:
    """Extract unique chain IDs from PDB ATOM records, in order of appearance."""
    chains: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
            cid = line[21]
            if cid not in chains:
                chains.append(cid)
    return chains


def _build_chain_legend(chains: list[str], peptide: str, allele: str) -> str:
    """Build an HTML table mapping chain IDs to biological components."""
    rows = []
    for cid in chains:
        name = _VIEWER_CHAIN_NAMES.get(cid, f"Chain {cid}")
        detail = ""
        if "MHC" in name:
            detail = allele
        elif name == "Peptide":
            detail = f"<code>{peptide}</code>"
        rows.append(
            f"<tr><td><strong>{cid}</strong></td><td>{name}</td><td>{detail}</td></tr>"
        )
    if not rows:
        return ""
    return (
        '<table class="chain-legend">'
        "<thead><tr><th>Chain</th><th>Component</th><th>Details</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _build_compnd_records(pdb_text: str, peptide: str, allele: str) -> str:
    """Build PDB COMPND records so Mol* shows meaningful chain names.

    PDB format: cols 1-6 = "COMPND", col 7 = blank, cols 8-10 = continuation
    (blank for the first line, right-justified integer from line 2 onwards),
    col 11+ = compound text.  Lines are padded to 80 characters.
    PDB is ASCII-only so Greek letters are transliterated.
    """
    compnd_texts: list[str] = []
    for idx, cid in enumerate(_extract_chain_ids(pdb_text)):
        name = _VIEWER_CHAIN_NAMES.get(cid, f"Chain {cid}")
        name = name.replace("\u03b1", "alpha").replace("\u03b2", "beta")
        if "MHC" in name and allele:
            name = f"{name} ({allele})"
        elif "Peptide" in name and peptide:
            name = f"Peptide ({peptide})"
        mol_num = idx + 1
        compnd_texts.append(f"MOL_ID: {mol_num};")
        compnd_texts.append(f"MOLECULE: {name};")
        compnd_texts.append(f"CHAIN: {cid};")
    if not compnd_texts:
        return pdb_text
    compnd_lines = []
    for i, text in enumerate(compnd_texts):
        if i == 0:
            compnd_lines.append(f"COMPND    {text}".ljust(80))
        else:
            compnd_lines.append(f"COMPND {i + 1:>3d} {text}".ljust(80))
    return "\n".join(compnd_lines) + "\n" + pdb_text

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Splice Neoepitope Report</title>
  <!-- Mermaid — pipeline diagram (requires network access to render).
       Version and integrity hash are a matched pair — update both together. -->
  <script src="https://cdn.jsdelivr.net/npm/mermaid@10.9.0/dist/mermaid.min.js"
          integrity="sha384-6F4Ibv/ylL12O35KFWTeGTHuBKDz5L6yjKsgv3QHQ8s4NTqlDXq7kMlYXGs7MHFc"
          crossorigin="anonymous"></script>
  <script>mermaid.initialize({{startOnLoad: true, theme: 'default'}});</script>
  {molstar_assets}
  <style>
    body  {{ font-family: Arial, sans-serif; margin: 2em; color: #333; }}
    h1    {{ color: #2c3e50; }}
    h2    {{ color: #34495e; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
    th    {{ background-color: #f2f2f2; }}
    tr:nth-child(even) {{ background-color: #f9f9f9; }}
    .mermaid {{ max-width: 560px; margin: 1.5em auto; }}

    /* Contig visualisation */
    .contig  {{ font-family: monospace; font-size: 0.85em; white-space: nowrap; }}
    .nt-up   {{ color: #555; }}
    .nt-down {{ color: #555; }}
    .nt-pep-up   {{ background: #d6eaf8; color: #1a5276; font-weight: bold; }}
    .nt-pep-down {{ background: #d5f5e3; color: #1e8449; font-weight: bold; }}
    .junction-mark {{ color: #e74c3c; font-weight: bold; }}

    /* Chain legend */
    .chain-legend {{
      margin: 1em 0; border-collapse: collapse; width: auto;
      font-size: 0.92em;
    }}
    .chain-legend th {{ background: #f2f2f2; padding: 5px 14px; }}
    .chain-legend td {{ padding: 5px 14px; }}
    .chain-legend code {{ background: #eee; padding: 1px 4px; border-radius: 3px; }}
  </style>
</head>
<body>
  <h1>Splice Neoepitope Pipeline Report</h1>
  <p>Generated by the splice neoepitope pipeline
     (reimplementation of Lee et al., 2015).</p>

  <h2>Pipeline overview</h2>
  {pipeline_diagram}

  {filtering_funnel}

  <h2>Junction origin summary</h2>
  <p>Splice junctions detected in the tumor that are absent from the GENCODE
  reference annotation, grouped by whether they also appear in the matched normal
  tissue. Only <strong class="keep">tumor_exclusive</strong> junctions — absent
  from both the reference and the normal sample — are used for neoepitope
  prediction. <strong class="discard">normal_shared</strong> junctions are
  present in normal tissue and are therefore not tumor-derived.</p>
  {origin_table}

  {hla_section}

  <h2>Neoepitope prediction summary</h2>
  {presenter_table}

  <h2>Top strong presentations (presentation_percentile &le; {presentation_percentile_strong}%)</h2>
  {strong_table}

  {structure_section}

  <hr/>
  <p><small>Pipeline source:
  <a href="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline">
  splice-neoepitope-pipeline</a>.
  Original work: Jin-Ho Lee, Seoul National University, 2015.</small></p>
</body>
</html>
"""


def _build_hla_section(hla_qc_tsv: str) -> str:
    """Build the HLA typing results section from hla_qc.tsv.

    Shows per-locus alleles, the source of each call (tumor / serology /
    normal / fallback), HLA read counts, normal/tumor discrepancy, and —
    when serology columns are present — a validation column comparing
    OptiType against the known clinical alleles.
    """
    import pandas as pd

    try:
        df = pd.read_csv(hla_qc_tsv, sep="\t")
    except Exception as exc:
        log.warning("Could not read HLA QC TSV %s: %s", hla_qc_tsv, exc)
        return ""

    if df.empty:
        return ""

    source_colours = {
        "tumor":    "#e67e22",
        "serology": "#8e44ad",
        "normal":   "#27ae60",
        "fallback": "#c0392b",
    }

    # Show serology column only when at least one locus has known alleles
    has_serology = (
        "serology_allele1" in df.columns
        and df["serology_allele1"].notna().any()
        and (df["serology_allele1"].astype(str).str.strip() != "").any()
    )

    rows = []
    has_discrepancy = False
    for _, row in df.iterrows():
        source = str(row.get("source", ""))
        colour = source_colours.get(source, "#888")
        source_badge = (
            f'<span style="color:{colour};font-weight:bold">'
            f'{html_mod.escape(source)}</span>'
        )
        _disc_raw = row.get("discrepancy", "")
        discrepancy = "" if pd.isna(_disc_raw) else str(_disc_raw).strip()
        if discrepancy:
            has_discrepancy = True
        disc_cell = (
            f'<td style="background:#fef3cd">{html_mod.escape(discrepancy)}</td>'
            if discrepancy else '<td style="color:#27ae60">&#10003; concordant</td>'
        )

        serology_cell = ""
        if has_serology:
            _ser1 = row.get("serology_allele1", "")
            _ser2 = row.get("serology_allele2", "")
            ser_allele1 = "" if pd.isna(_ser1) else str(_ser1).strip()
            ser_allele2 = "" if pd.isna(_ser2) else str(_ser2).strip()
            _val_raw = row.get("serology_validation", "")
            validation = "" if pd.isna(_val_raw) else str(_val_raw).strip()

            if ser_allele1:
                alleles_str = (
                    f"<code>{html_mod.escape(ser_allele1)}</code> / "
                    f"<code>{html_mod.escape(ser_allele2)}</code>"
                    if ser_allele2 and ser_allele2 != ser_allele1
                    else f"<code>{html_mod.escape(ser_allele1)}</code>"
                )
                if validation.startswith("match"):
                    val_badge = f' <span style="color:#27ae60">&#10003; {html_mod.escape(validation)}</span>'
                elif validation.startswith("mismatch"):
                    val_badge = f' <span style="color:#c0392b">&#10007; {html_mod.escape(validation)}</span>'
                else:
                    val_badge = ""
                serology_cell = f"<td>{alleles_str}{val_badge}</td>"
            else:
                serology_cell = "<td></td>"

        rows.append(
            f"<tr>"
            f"<td><strong>{html_mod.escape(str(row['locus']))}</strong></td>"
            f"<td><code>{html_mod.escape(str(row['allele1']))}</code></td>"
            f"<td><code>{html_mod.escape(str(row['allele2']))}</code></td>"
            f"<td>{source_badge}</td>"
            f"<td>{html_mod.escape(str(row.get('reads', '')))}</td>"
            f"{disc_cell}"
            f"{serology_cell}"
            f"</tr>"
        )

    serology_header = "<th>Known alleles / validation</th>" if has_serology else ""
    table = (
        "<table><thead><tr>"
        "<th>Locus</th><th>Allele 1</th><th>Allele 2</th>"
        "<th>Source</th><th>HLA reads</th><th>Normal/tumor discrepancy</th>"
        f"{serology_header}"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )

    note = ""
    if has_discrepancy:
        note = (
            "<p><em><strong>Note:</strong> Normal/tumor discrepancy detected in one or "
            "more loci. Tumor sample calls are used for prediction (tumor-first policy). "
            "Discrepancies may indicate loss of heterozygosity at the HLA locus in the "
            "tumor, or low read depth in the test subset.</em></p>"
        )

    return (
        "<h2>HLA typing (OptiType)</h2>"
        "<p>Patient HLA-A/B/C alleles typed by OptiType from RNA-seq reads. "
        "Priority order: <strong>tumor OptiType → serology → normal OptiType → fallback</strong>. "
        "Tumor calls are preferred as they reflect what the tumor cell actually presents, "
        "including any HLA loss-of-heterozygosity. "
        "These alleles were used for patient-specific neoepitope prediction.</p>"
        + note + table
    )


def _load_contigs(contigs_fasta: str | Path) -> dict[str, str]:
    """Parse contigs FASTA into {contig_key: sequence} dict.

    The contig key matches the contig_key column in the predictions TSV.
    """
    from Bio import SeqIO

    contigs: dict[str, str] = {}
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        contigs[record.description] = str(record.seq).upper()
    return contigs


def _render_contig_peek(peek: str) -> str:
    """Render a bracketed contig peek (from report_top_candidates.tsv) as styled HTML.

    Parses the plain-text format produced by ``_build_contig_peek``:
      - ``[`` and ``]`` enclose peptide nucleotides
      - ``|`` marks the splice junction
      - other chars are flanking nucleotides

    Output styling matches ``_render_contig`` so the artefact-driven render
    looks identical to the raw-driven one.
    """
    if not peek:
        return ""
    in_peptide = False
    seen_junction = False
    parts: list[str] = []
    for ch in peek:
        if ch == "[":
            in_peptide = True
        elif ch == "]":
            in_peptide = False
        elif ch == "|":
            parts.append('<span class="junction-mark">|</span>')
            seen_junction = True
        else:
            if in_peptide:
                cls = "nt-pep-down" if seen_junction else "nt-pep-up"
            else:
                cls = "nt-down" if seen_junction else "nt-up"
            parts.append(f'<span class="{cls}">{ch}</span>')
    return f'<span class="contig">{"".join(parts)}</span>'


def _render_contig(seq: str, start_nt: int, end_nt_incl: int,
                   upstream_nt: int = 26) -> str:
    """Render a 50 nt contig as HTML with junction marker and peptide highlighted."""
    html_parts = []
    for i, nt in enumerate(seq):
        if i == upstream_nt:
            html_parts.append('<span class="junction-mark">|</span>')
        in_pep = start_nt <= i <= end_nt_incl
        if in_pep and i < upstream_nt:
            html_parts.append(f'<span class="nt-pep-up">{nt}</span>')
        elif in_pep:
            html_parts.append(f'<span class="nt-pep-down">{nt}</span>')
        elif i < upstream_nt:
            html_parts.append(f'<span class="nt-up">{nt}</span>')
        else:
            html_parts.append(f'<span class="nt-down">{nt}</span>')
    return f'<span class="contig">{"".join(html_parts)}</span>'


def _build_strong_table_html(
    pred_df: pd.DataFrame,
    contigs: dict[str, str],
    upstream_nt: int = 26,
    max_rows: int = 50,
    presentation_percentile_weak: float = 2.0,
) -> str:
    """Build the top strong presentations table ranked by genotype_presentation_score."""
    strong_df = pred_df[pred_df["presentation_class"] == "strong"].copy()

    # Quality gate: at least one allele must reach weak-presenter threshold
    if "best_presentation_percentile" in strong_df.columns:
        strong_df = strong_df[strong_df["best_presentation_percentile"] <= presentation_percentile_weak]

    if strong_df.empty:
        return "<p><em>No strong presentations found.</em></p>"

    # Ranking: genotype_presentation_score desc → n_strong_alleles desc → best_presentation_percentile asc
    has_breadth_cols = "genotype_presentation_score" in strong_df.columns
    if has_breadth_cols:
        strong_df = strong_df.sort_values(
            ["genotype_presentation_score", "n_strong_alleles", "best_presentation_percentile"],
            ascending=[False, False, True],
        )
    else:
        strong_df = strong_df.sort_values("presentation_percentile")

    rows = []
    for _, row in strong_df.head(max_rows).iterrows():
        contig_key = row["contig_key"]
        start_nt = int(row["start_nt"])
        end_nt_incl = start_nt + len(row["peptide"]) * 3 - 1

        seq = contigs.get(contig_key, "")
        contig_html = _render_contig(seq, start_nt, end_nt_incl, upstream_nt) if seq else "<em>n/a</em>"

        gps_cell = (
            f"<td>{row['genotype_presentation_score']:.4f}</td>"
            if has_breadth_cols else ""
        )
        rows.append(
            f"<tr>"
            f"<td>{html_mod.escape(str(row['contig_key']))}</td>"
            f"<td>{html_mod.escape(str(row['peptide']))}</td>"
            f"<td>{html_mod.escape(str(row['best_allele']))}</td>"
            f"<td>{row['presentation_percentile']:.3f}</td>"
            f"{gps_cell}"
            f"<td>{row['ic50_nM']:.1f}</td>"
            f"<td>{contig_html}</td>"
            f"</tr>"
        )

    legend = (
        "upstream <span class='junction-mark'>|</span> downstream — "
        "<span class='nt-pep-up'>peptide (upstream)</span> "
        "<span class='nt-pep-down'>peptide (downstream)</span>"
    )
    gps_header = (
        "<th title='Probability ≥1 genotype allele presents peptide — primary rank'>GPS ▾</th>"
        if has_breadth_cols else ""
    )
    table = (
        f"<table><thead><tr>"
        f"<th>Source</th><th>Peptide</th><th>Best Allele</th>"
        f"<th title='Best-allele presentation percentile rank'>Best-allele %ile</th>"
        f"{gps_header}"
        f"<th title='Binding affinity in nM — informational'>IC50 (nM)</th>"
        f"<th>Contig ({legend})</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    )
    if len(strong_df) > max_rows:
        table += f"<p><em>Showing {max_rows} of {len(strong_df)} rows.</em></p>"
    return table


def _df_to_html(df: pd.DataFrame, max_rows: int = 100) -> str:
    if df.empty:
        return "<p><em>No data.</em></p>"
    out = df.head(max_rows).to_html(index=False, border=0, escape=True)
    if len(df) > max_rows:
        out += f"<p><em>Showing {max_rows} of {len(df)} rows.</em></p>"
    return out


def _build_filtering_funnel_html(filtering_stats_tsv: str | Path | None) -> str:
    """Render the cross-step filtering funnel section (Issue #215).

    Two tables:
      1. Per-sample junction-filter view (sample × funnel category counts +
         raw read distribution summary).
      2. Patient-level pipeline funnel for the four downstream steps
         (contig-assemble → mhc-affinity).

    Returns an empty string if the TSV is missing or unreadable —
    backward-compatible with pipelines that haven't run the aggregator.
    """
    if not filtering_stats_tsv or not Path(filtering_stats_tsv).exists():
        return ""
    import pandas as pd

    try:
        df = pd.read_csv(filtering_stats_tsv, sep="\t").fillna("")
    except Exception as exc:
        log.warning("Could not read filtering_stats TSV %s: %s", filtering_stats_tsv, exc)
        return ""
    if df.empty:
        return ""

    # ----- Per-sample junction-filter table -----
    junction_df = df[df["step"] == "junction-filter"].copy()
    funnel_cats = [
        "junctions_raw", "mean_reads_filtered", "annotated_discarded",
        "normal_shared", "tumor_exclusive",
    ]
    dist_cats = ["min_reads", "median_reads", "mean_reads", "max_reads"]

    if junction_df.empty:
        junction_html = "<p><em>No junction-filter stats available.</em></p>"
    else:
        funnel_pivot = (
            junction_df[junction_df["category"].isin(funnel_cats)]
            .pivot_table(
                index=["sample_id", "sample_type"],
                columns="category",
                values="count",
                aggfunc="first",
            )
            .reset_index()
            .reindex(columns=["sample_id", "sample_type"] + funnel_cats)
            .fillna(0)
        )
        dist_pivot = (
            junction_df[junction_df["category"].isin(dist_cats)]
            .pivot_table(
                index=["sample_id", "sample_type"],
                columns="category",
                values="count",
                aggfunc="first",
            )
            .reset_index()
            .reindex(columns=["sample_id", "sample_type"] + dist_cats)
            .fillna(0)
        )
        junction_html = (
            "<h3>Junction-level (per sample)</h3>"
            + _df_to_html(funnel_pivot)
            + "<h3>Raw read-count distribution (per sample, pre-filter)</h3>"
            + "<p><small>Mean is the threshold used by the per-file noise filter; "
            "min/median/max describe the input distribution.</small></p>"
            + _df_to_html(dist_pivot)
        )

    # ----- Patient-level pipeline funnel for downstream steps -----
    downstream = df[df["step"] != "junction-filter"][["step", "category", "count"]].copy()
    if downstream.empty:
        downstream_html = ""
    else:
        downstream_html = (
            "<h3>Pipeline funnel (patient-level)</h3>"
            + _df_to_html(downstream)
        )

    return (
        "<h2>Filtering funnel (audit trail)</h2>"
        "<p>Per-step counts of candidates retained vs filtered. "
        "Full long-format detail lives in <code>filtering_stats.tsv</code>.</p>"
        + junction_html
        + downstream_html
    )


def _build_strong_table_html_from_top_candidates(top_candidates_df: pd.DataFrame) -> str:
    """Build the top presenters table from the wide ``report_top_candidates.tsv``.

    The artefact is already quality-gated, sorted, capped, and carries a
    plain-text ``contig_peek`` that this function re-styles via
    ``_render_contig_peek``. No raw inputs needed.

    Parity gap with ``_build_strong_table_html``: this renderer cannot emit a
    "Showing N of M rows" notice because the writer caps the artefact at
    ``TOP_CANDIDATES_LIMIT`` and the pre-cap total is lost. Tracked in
    Issue #226 — fix is to surface the total in ``report.tsv``'s
    ``mhc_prediction`` stage and render the notice from there.
    """
    if top_candidates_df is None or top_candidates_df.empty:
        return "<p><em>No strong presentations found.</em></p>"

    import pandas as pd

    has_gps = "genotype_presentation_score" in top_candidates_df.columns

    rows = []
    for _, row in top_candidates_df.iterrows():
        contig_html = _render_contig_peek(str(row.get("contig_peek", ""))) or "<em>n/a</em>"

        gps_cell = ""
        if has_gps:
            gps_val = row.get("genotype_presentation_score", "")
            gps_cell = f"<td>{float(gps_val):.4f}</td>" if gps_val != "" and pd.notna(gps_val) else "<td></td>"

        pct_val = row.get("best_presentation_percentile", "")
        pct_cell = f"<td>{float(pct_val):.3f}</td>" if pct_val != "" and pd.notna(pct_val) else "<td></td>"

        ic50_val = row.get("ic50_nM", "")
        ic50_cell = f"<td>{float(ic50_val):.1f}</td>" if ic50_val != "" and pd.notna(ic50_val) else "<td></td>"

        rows.append(
            f"<tr>"
            f"<td>{html_mod.escape(str(row.get('contig_key', '')))}</td>"
            f"<td>{html_mod.escape(str(row.get('peptide', '')))}</td>"
            f"<td>{html_mod.escape(str(row.get('best_allele', '')))}</td>"
            f"{pct_cell}"
            f"{gps_cell}"
            f"{ic50_cell}"
            f"<td>{contig_html}</td>"
            f"</tr>"
        )

    legend = (
        "upstream <span class='junction-mark'>|</span> downstream — "
        "<span class='nt-pep-up'>peptide (upstream)</span> "
        "<span class='nt-pep-down'>peptide (downstream)</span>"
    )
    gps_header = (
        "<th title='Probability ≥1 genotype allele presents peptide — primary rank'>GPS ▾</th>"
        if has_gps else ""
    )
    return (
        f"<table><thead><tr>"
        f"<th>Source</th><th>Peptide</th><th>Best Allele</th>"
        f"<th title='Best-allele presentation percentile rank'>Best-allele %ile</th>"
        f"{gps_header}"
        f"<th title='Binding affinity in nM — informational'>IC50 (nM)</th>"
        f"<th>Contig ({legend})</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    )


def _presenter_counts_html(
    mp: dict[str, int],
    thresholds: dict[str, str] | None = None,
) -> str:
    """Render the presentation-class count table from the loaded report.tsv projection.

    ``mp`` is the ``mhc_prediction`` dict (from ``_load_report_tsv``); the
    ``total_predictions`` entry is excluded from the table and the remaining
    metrics become the per-class rows.
    """
    if not mp or all(k == "total_predictions" for k in mp):
        return "<p><em>No predictions available.</em></p>"

    rows = [(cls, cnt) for cls, cnt in mp.items() if cls != "total_predictions"]
    if not rows:
        return "<p><em>No predictions available.</em></p>"

    import pandas as pd

    df = pd.DataFrame(rows, columns=["presentation_class", "count"])
    return _df_to_html(df)


def _rank_presenters(df: pd.DataFrame) -> pd.DataFrame:
    """Sort presenters by genotype_presentation_score desc, then percentile asc."""
    if "genotype_presentation_score" in df.columns:
        return df.sort_values(
            ["genotype_presentation_score", "n_strong_alleles", "best_presentation_percentile"],
            ascending=[False, False, True],
        )
    return df.sort_values("presentation_percentile")


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def _resolve_top_candidate_for_structure(
    pred_df: pd.DataFrame,
    presentation_percentile_weak: float = 2.0,
) -> tuple[str, str]:
    """Resolve top candidate peptide + allele from raw predictions.

    Mirrors the candidate-selection ranking used by ``run_tcrdock.py`` and
    ``_build_report_top_candidates_tsv`` so all three surfaces stay aligned.
    Returns ``("NA", "NA")`` when no candidate passes the quality gate.
    """
    presenters = pred_df[pred_df["presentation_class"].isin(("strong", "weak"))].copy()
    if "best_presentation_percentile" in presenters.columns:
        presenters = presenters[presenters["best_presentation_percentile"] <= presentation_percentile_weak]
    presenters = _rank_presenters(presenters)
    top = presenters.head(1)
    if top.empty:
        return "NA", "NA"
    return str(top.iloc[0]["peptide"]), str(top.iloc[0]["best_allele"])


def _resolve_top_candidate_from_manifest(
    structure_manifest: pd.DataFrame,
) -> tuple[str, str]:
    """Read top candidate peptide + allele from the 3D structure manifest."""
    if structure_manifest is None or structure_manifest.empty:
        return "NA", "NA"
    top = structure_manifest.iloc[0]
    return str(top.get("peptide", "NA")), str(top.get("allele", "NA"))


def _build_structure_section(pdb_path: Path, peptide: str, allele: str) -> str:
    """Build the Mol* 3D viewer section for the top TCRdock candidate.

    The PDB is inlined as a JSON string so the HTML report is fully
    self-contained and does not depend on local file paths at view time.

    Args:
        pdb_path: Path to the TCRdock-predicted PDB file.
        peptide:  Top candidate peptide (already resolved by caller).
        allele:   Top candidate allele (already resolved by caller).

    Returns:
        HTML string containing the Mol* viewer and inlined PDB.
    """
    try:
        pdb_text = pdb_path.read_text()
    except OSError as exc:
        log.warning("Could not read PDB file %s: %s", pdb_path, exc)
        return "<p><em>Structure not available.</em></p>"

    peptide = html_mod.escape(peptide)
    allele = html_mod.escape(allele)

    pdb_text = _build_compnd_records(pdb_text, peptide, allele)

    # Inline PDB as a JSON string (safe for embedding in JS)
    pdb_json = json.dumps(pdb_text)

    chains = _extract_chain_ids(pdb_text)
    chain_legend = _build_chain_legend(chains, peptide, allele)

    viewer_html = _MOLSTAR_VIEWER.format(pdb_data=pdb_json)
    return (
        "<h2>TCR-peptide-MHC structure (TCRdock)</h2>"
        "<p>Predicted 3D structure of the TCR-peptide-MHC ternary complex for "
        f"the top neoepitope candidate (<strong>{peptide}</strong> / {allele}). "
        "Each chain is colored separately &mdash; hover over a region in the "
        "3D view to see chain details, or expand the sequence panel "
        "(top of the viewer) to browse per-chain sequences. "
        "Rendered with <a href='https://molstar.org'>Mol*</a>. "
        "TCR sequences: DMF5 fallback (see config).</p>"
        + chain_legend
        + viewer_html
    )


def _build_report_tsv(
    patient_id: str,
    origin_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    hla_qc_tsv: str | None,
    output_tsv: str | Path,
    presentation_percentile_strong: float,
    tcrdock_pdb: str | Path | None,
    presentation_percentile_weak: float = 2.0,
    junction_filter_stats_tsv: str | Path | None = None,
) -> None:
    """Write a structured summary TSV alongside the HTML report.

    Schema: patient_id | stage | metric | value | notes
    Stages: junction_filtering, mhc_prediction, top_candidate, hla_typing, tcrdock

    When ``junction_filter_stats_tsv`` is supplied, four patient-level funnel
    totals are added under stage=``junction_filtering`` (Issue #214):
    ``junctions_extracted_total``, ``junctions_mean_reads_filtered``,
    ``junctions_annotated_discarded``, ``junctions_unannotated_total``. These
    sum across tumor samples; normal samples are intentionally omitted from
    the funnel. The four rows reconcile arithmetically:
    ``extracted_total = mean_reads_filtered + annotated_discarded + unannotated_total``.
    """
    import pandas as pd

    rows: list[dict] = []

    # --- junction_filtering ---
    if not origin_df.empty:
        for _, r in origin_df.iterrows():
            sid = r.get("sample_id", "")
            stype = r.get("sample_type", "")
            note = f"{sid} ({stype})"
            for metric in ("unannotated", "tumor_exclusive", "normal_shared"):
                rows.append({
                    "patient_id": patient_id, "stage": "junction_filtering",
                    "metric": metric, "value": int(r.get(metric, 0)), "notes": note,
                })

    # --- junction_filtering: patient-level funnel totals (Issue #214) ---
    if junction_filter_stats_tsv is not None:
        stats_df = pd.read_csv(junction_filter_stats_tsv, sep="\t")
        if not stats_df.empty:
            by_cat = stats_df.groupby("category")["count"].sum()
            for metric, source_cats in (
                ("junctions_extracted_total", ("junctions_raw",)),
                ("junctions_mean_reads_filtered", ("mean_reads_filtered",)),
                ("junctions_annotated_discarded", ("annotated_discarded",)),
                ("junctions_unannotated_total", ("normal_shared", "tumor_exclusive")),
            ):
                value = sum(int(by_cat.get(c, 0)) for c in source_cats)
                rows.append({
                    "patient_id": patient_id, "stage": "junction_filtering",
                    "metric": metric, "value": value, "notes": "all tumor samples",
                })

    # --- mhc_prediction ---
    if not pred_df.empty:
        rows.append({
            "patient_id": patient_id, "stage": "mhc_prediction",
            "metric": "total_predictions", "value": len(pred_df), "notes": "",
        })
        presentation_notes = {
            "strong": f"presentation_percentile <= {presentation_percentile_strong}%",
            "weak": f"presentation_percentile <= {presentation_percentile_weak}%",
            "non": f"presentation_percentile > {presentation_percentile_weak}%",
        }
        for cls, cnt in pred_df["presentation_class"].value_counts().items():
            rows.append({
                "patient_id": patient_id, "stage": "mhc_prediction",
                "metric": str(cls), "value": int(cnt),
                "notes": presentation_notes.get(str(cls), ""),
            })

    # --- top_candidate ---
    if not pred_df.empty:
        top_candidates = pred_df[pred_df["presentation_class"].isin(["strong", "weak"])].copy()
        if "best_presentation_percentile" in top_candidates.columns:
            top_candidates = top_candidates[top_candidates["best_presentation_percentile"] <= presentation_percentile_weak]
        top_candidates = _rank_presenters(top_candidates)
        strong_df = top_candidates
    else:
        strong_df = pd.DataFrame()

    if not strong_df.empty:
        top = strong_df.iloc[0]
        numeric_metrics = ("presentation_percentile", "ic50_nM", "genotype_presentation_score")
        for metric, col in [
            ("peptide", "peptide"), ("allele", "best_allele"),
            ("presentation_percentile", "presentation_percentile"),
            ("genotype_presentation_score", "genotype_presentation_score"),
            ("ic50_nM", "ic50_nM"), ("presentation_class", "presentation_class"),
        ]:
            if col not in top.index:
                continue
            rows.append({
                "patient_id": patient_id, "stage": "top_candidate",
                "metric": metric,
                "value": round(float(top[col]), 4) if metric in numeric_metrics else str(top[col]),
                "notes": "top by genotype_presentation_score" if metric == "peptide" else "",
            })

    # --- hla_typing ---
    if hla_qc_tsv:
        try:
            hla_df = pd.read_csv(hla_qc_tsv, sep="\t")
            for _, r in hla_df.iterrows():
                locus = str(r.get("locus", ""))
                allele1 = str(r.get("allele1", ""))
                allele2 = str(r.get("allele2", ""))
                source = str(r.get("source", ""))
                reads = r.get("reads", 0)
                alleles = f"{allele1} / {allele2}" if allele1 != allele2 else allele1
                rows.append({
                    "patient_id": patient_id, "stage": "hla_typing",
                    "metric": f"HLA-{locus}", "value": alleles,
                    "notes": f"source: {source}, reads: {reads}",
                })
        except Exception as exc:
            log.warning("Could not read HLA QC TSV for report.tsv: %s", exc)

    # --- tcrdock ---
    pdb_available = (
        tcrdock_pdb is not None and Path(tcrdock_pdb).exists()
    )
    rows.append({
        "patient_id": patient_id, "stage": "tcrdock",
        "metric": "pdb_available", "value": str(pdb_available).lower(), "notes": "",
    })

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["patient_id", "stage", "metric", "value", "notes"]).to_csv(
        output_path, sep="\t", index=False
    )
    log.info("Report TSV written to %s", output_path)


# ---------------------------------------------------------------------------
# Machine-readable artefacts: top presenters and 3D structure manifest
# ---------------------------------------------------------------------------

# Cap on how many top-ranked presenters land in both report_top_candidates.tsv
# and the HTML rendering. Single source of truth: change here, both surfaces follow.
TOP_CANDIDATES_LIMIT = 10


def _allele_slot_map(predictions_columns: list[str]) -> dict[str, str]:
    """Map per-allele percentile column names to stable slot names.

    The mhc_presentation.tsv emits per-allele columns like
    ``HLA-A*02:01_presentation_percentile``. Slot mapping turns that into
    a stable wide schema (``hla_a1_pct``, ``hla_a2_pct``, ...) so
    ``report_top_candidates.tsv`` has a fixed column set across patients.

    Allocation: alleles are grouped by locus letter (A/B/C), sorted
    lexicographically within the locus, then assigned to slot 1, slot 2.
    Homozygous loci yield one slot mapped; the second slot stays empty.

    Returns:
        ``{allele_id: slot_letter_index}`` e.g. ``{"HLA-A*02:01": "a1"}``.
    """
    suffix = "_presentation_percentile"
    by_locus: dict[str, list[str]] = {"a": [], "b": [], "c": []}
    for col in predictions_columns:
        if not col.startswith("HLA-") or not col.endswith(suffix):
            continue
        allele = col[: -len(suffix)]
        # allele looks like "HLA-A*02:01"; locus letter is index 4
        locus_letter = allele[4].lower()
        if locus_letter in by_locus:
            by_locus[locus_letter].append(allele)

    slot_map: dict[str, str] = {}
    for locus, alleles in by_locus.items():
        for slot_idx, allele in enumerate(sorted(set(alleles))[:2], start=1):
            slot_map[allele] = f"{locus}{slot_idx}"
    return slot_map


def _build_contig_peek(seq: str, start_nt: int, end_nt_incl: int, upstream_nt: int = 26) -> str:
    """Render a 50 nt contig as a plain-text peek with junction and peptide markers.

    Format: brackets ``[ ... ]`` enclose the nucleotides that translate to the
    displayed peptide; ``|`` marks the splice junction inside the peptide.
    Peptides in this pipeline are junction-spanning by construction, so the
    bracketed region always contains the ``|``. Example (peptide spans junction):
    ``ATCGATCGATCG[CGAT|CGATCGAT]CGATCGATCGATCGATCGATCG`` — empty if seq is empty.
    """
    if not seq:
        return ""
    parts: list[str] = []
    for i, nt in enumerate(seq):
        if i == upstream_nt:
            parts.append("|")
        if i == start_nt:
            parts.append("[")
        parts.append(nt)
        if i == end_nt_incl:
            parts.append("]")
    return "".join(parts)


def _round_or_blank(value, ndigits: int) -> str | float:
    """Round numeric values; return empty string for missing/non-numeric."""
    import pandas as pd

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return ""


def _build_report_top_candidates_tsv(
    patient_id: str,
    pred_df: pd.DataFrame,
    contigs: dict[str, str],
    output_tsv: str | Path,
    presentation_percentile_weak: float = 2.0,
    upstream_nt: int = 26,
) -> None:
    """Write the top presenters as a wide-format TSV, ranked by genotype_presentation_score.

    Schema (stable, cross-patient comparable):
      patient_id | rank | peptide | best_allele
        | best_presentation_percentile | genotype_presentation_score | ic50_nM
        | n_strong_alleles | presentation_class
        | hla_a1_id | hla_a1_pct | hla_a2_id | hla_a2_pct
        | hla_b1_id | hla_b1_pct | hla_b2_id | hla_b2_pct
        | hla_c1_id | hla_c1_pct | hla_c2_id | hla_c2_pct
        | contig_key | contig_start_nt | contig_peek

    The HLA slot columns map the patient's alleles into stable positions
    (a1/a2/b1/b2/c1/c2) so downstream consumers can compare across patients
    without per-patient column lookups; allele identity is preserved in
    the ``hla_<slot>_id`` columns.

    Limited to ``TOP_CANDIDATES_LIMIT`` rows after the same quality gate the
    HTML render uses (``best_presentation_percentile <= presentation_percentile_weak``).
    """
    import pandas as pd

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    slot_columns = [
        "hla_a1_id", "hla_a1_pct", "hla_a2_id", "hla_a2_pct",
        "hla_b1_id", "hla_b1_pct", "hla_b2_id", "hla_b2_pct",
        "hla_c1_id", "hla_c1_pct", "hla_c2_id", "hla_c2_pct",
    ]
    schema = [
        "patient_id", "rank", "peptide", "best_allele",
        "best_presentation_percentile", "genotype_presentation_score", "ic50_nM",
        "n_strong_alleles", "presentation_class",
        *slot_columns,
        "contig_key", "contig_start_nt", "contig_peek",
    ]

    if pred_df.empty:
        pd.DataFrame(columns=schema).to_csv(output_path, sep="\t", index=False)
        log.info("Top presenters TSV written (empty) to %s", output_path)
        return

    df = pred_df[pred_df["presentation_class"].isin(["strong", "weak"])].copy()
    if "best_presentation_percentile" in df.columns:
        df = df[df["best_presentation_percentile"] <= presentation_percentile_weak]

    df = _rank_presenters(df)
    df = df.head(TOP_CANDIDATES_LIMIT).reset_index(drop=True)

    if df.empty:
        pd.DataFrame(columns=schema).to_csv(output_path, sep="\t", index=False)
        log.info("Top presenters TSV written (no qualifying presenters) to %s", output_path)
        return

    slot_map = _allele_slot_map(list(pred_df.columns))
    inverse_slot: dict[str, str] = {v: k for k, v in slot_map.items()}

    rows: list[dict] = []
    for rank, (_, row) in enumerate(df.iterrows(), start=1):
        contig_key = row.get("contig_key", "")
        start_nt = int(row.get("start_nt", 0)) if pd.notna(row.get("start_nt")) else 0
        peptide = str(row.get("peptide", ""))
        end_nt_incl = start_nt + len(peptide) * 3 - 1 if peptide else start_nt
        seq = contigs.get(contig_key, "")

        out: dict = {
            "patient_id": patient_id,
            "rank": rank,
            "peptide": peptide,
            "best_allele": row.get("best_allele", ""),
            "best_presentation_percentile": _round_or_blank(row.get("best_presentation_percentile"), 4),
            "genotype_presentation_score": _round_or_blank(row.get("genotype_presentation_score"), 4),
            "ic50_nM": _round_or_blank(row.get("ic50_nM"), 1),
            "n_strong_alleles": int(row["n_strong_alleles"]) if pd.notna(row.get("n_strong_alleles")) else "",
            "presentation_class": row.get("presentation_class", ""),
            "contig_key": contig_key,
            "contig_start_nt": start_nt,
            "contig_peek": _build_contig_peek(seq, start_nt, end_nt_incl, upstream_nt),
        }
        for slot in ("a1", "a2", "b1", "b2", "c1", "c2"):
            allele = inverse_slot.get(slot, "")
            out[f"hla_{slot}_id"] = allele
            if allele:
                pct_col = f"{allele}_presentation_percentile"
                out[f"hla_{slot}_pct"] = _round_or_blank(row.get(pct_col), 4)
            else:
                out[f"hla_{slot}_pct"] = ""
        rows.append(out)

    pd.DataFrame(rows, columns=schema).to_csv(output_path, sep="\t", index=False)
    log.info("Top presenters TSV (%d rows) written to %s", len(rows), output_path)


# Chain mapping for relabel_pdb_chains() in run_tcrdock.py — duplicated here to
# avoid a cross-script import (run_tcrdock imports torch/Bio at module load
# which we don't want pulled into report generation).
_TCRDOCK_CHAIN_NAMES = {"A": "MHC", "B": "peptide", "C": "TCR-alpha", "D": "TCR-beta"}


def _build_report_3d_structure_tsv(
    patient_id: str,
    docking_scores_tsv: str | Path,
    pdb_relative_path: str,
    output_tsv: str | Path,
) -> None:
    """Write the 3D structure manifest TSV.

    The PDB itself stays as the canonical machine-readable representation in
    its TCRdock output location; this manifest exposes the metadata + chain
    labelling so any downstream tool can locate and interpret the structure
    without re-running TCRdock or re-deriving the chain mapping.

    Schema:
      patient_id | rank | peptide | allele | pdb_path
        | chain_A | chain_B | chain_C | chain_D
        | <pass-through plddt/pae/ptm columns from docking_scores.tsv>

    ``pdb_path`` is recorded relative to the patient result directory
    (typically ``predictions/tcrdock/top_candidate.pdb``) so consumers can
    resolve it from the patient root regardless of how reports/ was copied.

    The manifest tracks the single PDB that's actually written today
    (n_candidates=1). If TCRdock starts emitting multiple PDBs per run, the
    rank column is the natural extension point.
    """
    import pandas as pd

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fixed_columns = [
        "patient_id", "rank", "peptide", "allele", "pdb_path",
        "chain_A", "chain_B", "chain_C", "chain_D",
    ]

    docking_path = Path(docking_scores_tsv)
    if not docking_path.exists():
        pd.DataFrame(columns=fixed_columns).to_csv(output_path, sep="\t", index=False)
        log.warning("docking_scores.tsv not found at %s — wrote empty manifest", docking_path)
        return

    scores_df = pd.read_csv(docking_path, sep="\t")
    if scores_df.empty:
        pd.DataFrame(columns=fixed_columns).to_csv(output_path, sep="\t", index=False)
        log.warning("docking_scores.tsv at %s is empty — wrote empty manifest", docking_path)
        return

    # Pass-through any plddt/pae/ptm columns that exist; preserves whatever
    # TCRdock currently emits without forcing a hard schema dependency.
    metric_cols = [c for c in scores_df.columns
                   if any(k in c.lower() for k in ("plddt", "pae", "ptm"))]

    top = scores_df.iloc[0]
    row: dict = {
        "patient_id": patient_id,
        "rank": 1,
        "peptide": top.get("peptide", ""),
        "allele": top.get("mhc") or top.get("allele", ""),
        "pdb_path": pdb_relative_path,
        "chain_A": _TCRDOCK_CHAIN_NAMES["A"],
        "chain_B": _TCRDOCK_CHAIN_NAMES["B"],
        "chain_C": _TCRDOCK_CHAIN_NAMES["C"],
        "chain_D": _TCRDOCK_CHAIN_NAMES["D"],
    }
    for col in metric_cols:
        row[col] = _round_or_blank(top.get(col), 4)

    pd.DataFrame([row], columns=fixed_columns + metric_cols).to_csv(
        output_path, sep="\t", index=False
    )
    log.info("3D structure manifest written to %s", output_path)


# Numeric metrics in the top_candidate stage — typed back to float on load.
# Mirrors the writer's set in _build_report_tsv. Keep in sync with the writer.
_TOP_CANDIDATE_NUMERIC_METRICS = frozenset({
    "presentation_percentile", "ic50_nM", "genotype_presentation_score",
})


def _load_report_tsv(path: str | Path) -> dict[str, Any]:
    """Load report.tsv into a stage-keyed structure for HTML rendering.

    Inverts the long-format ``patient_id | stage | metric | value | notes`` shape
    written by ``_build_report_tsv`` into per-stage projections. Non-lossless
    by design — only the projections HTML needs are exposed; consumers needing
    the raw long form should read the TSV directly.

    Returns:
        {
          'junction_filtering':        DataFrame[sample_id, sample_type,
                                                 unannotated, normal_shared, tumor_exclusive],
          'mhc_prediction':            {metric: int}      # total_predictions + class counts
          'mhc_prediction_thresholds': {metric: str}      # threshold notes per class
          'top_candidate':             {metric: str|float}# numerics typed; strings as-is
          'hla_typing':                {locus: {alleles, notes}}  # locus stripped of "HLA-" prefix
          'tcrdock':                   {metric: bool}     # pdb_available cast to bool
        }
    """
    import pandas as pd

    df = pd.read_csv(path, sep="\t")

    out: dict[str, Any] = {
        "junction_filtering": _empty_origin_df(),
        "mhc_prediction": {},
        "mhc_prediction_thresholds": {},
        "top_candidate": {},
        "hla_typing": {},
        "tcrdock": {},
    }

    # --- junction_filtering: long → wide, parsing "sid (stype)" from notes ---
    jf = df[df["stage"] == "junction_filtering"]
    if not jf.empty:
        rows: dict[tuple[str, str], dict] = {}
        for _, r in jf.iterrows():
            note = str(r.get("notes", "")).strip()
            # Format produced by writer: "<sample_id> (<sample_type>)"
            if " (" in note and note.endswith(")"):
                sid, rest = note.rsplit(" (", 1)
                stype = rest[:-1]
            else:
                sid, stype = note, ""
            key = (sid, stype)
            row = rows.setdefault(key, {"sample_id": sid, "sample_type": stype})
            row[str(r["metric"])] = int(r["value"])
        origin_df = pd.DataFrame(list(rows.values()))
        for col in ("unannotated", "normal_shared", "tumor_exclusive"):
            if col not in origin_df.columns:
                origin_df[col] = 0
        out["junction_filtering"] = origin_df[
            ["sample_id", "sample_type", "unannotated", "normal_shared", "tumor_exclusive"]
        ]

    # --- mhc_prediction: counts (int) + thresholds (str) ---
    mp = df[df["stage"] == "mhc_prediction"]
    for _, r in mp.iterrows():
        metric = str(r["metric"])
        out["mhc_prediction"][metric] = int(r["value"])
        if metric != "total_predictions":
            note = r.get("notes", "")
            if pd.notna(note) and str(note).strip():
                out["mhc_prediction_thresholds"][metric] = str(note)

    # --- top_candidate: numerics typed, strings as-is ---
    tc = df[df["stage"] == "top_candidate"]
    for _, r in tc.iterrows():
        metric = str(r["metric"])
        value = r["value"]
        if metric in _TOP_CANDIDATE_NUMERIC_METRICS:
            try:
                value = float(value)
            except (TypeError, ValueError):
                pass
        else:
            value = str(value)
        out["top_candidate"][metric] = value

    # --- hla_typing: keyed by locus (HLA- prefix stripped) ---
    ht = df[df["stage"] == "hla_typing"]
    for _, r in ht.iterrows():
        locus = str(r["metric"]).removeprefix("HLA-")
        out["hla_typing"][locus] = {
            "alleles": str(r["value"]),
            "notes": "" if pd.isna(r.get("notes")) else str(r["notes"]),
        }

    # --- tcrdock: pdb_available cast to bool ---
    td = df[df["stage"] == "tcrdock"]
    for _, r in td.iterrows():
        metric = str(r["metric"])
        if metric == "pdb_available":
            out["tcrdock"][metric] = (str(r["value"]).strip().lower() == "true")
        else:
            out["tcrdock"][metric] = r["value"]

    return out


def _empty_origin_df() -> pd.DataFrame:
    """Empty origin DataFrame with the canonical column set."""
    import pandas as pd

    return pd.DataFrame(columns=[
        "sample_id", "sample_type", "unannotated", "normal_shared", "tumor_exclusive",
    ])


def generate_report(
    novel_junctions_tsv: str | Path,
    predictions_tsv: str | Path,
    output_html: str | Path,
    contigs_fasta: str | Path | None = None,
    hla_qc_tsv: str | None = None,
    presentation_percentile_strong: float = 0.5,
    presentation_percentile_weak: float = 2.0,
    tcrdock_pdb: str | Path | None = None,
    docking_scores_tsv: str | Path | None = None,
    output_tsv: str | Path | None = None,
    output_top_candidates_tsv: str | Path | None = None,
    output_3d_structure_tsv: str | Path | None = None,
    junction_filter_stats_tsv: str | Path | None = None,
    filtering_stats_tsv: str | Path | None = None,
    patient_id: str = "",
) -> None:
    """Generate the summary HTML report and the machine-readable report artefacts.

    Args:
        novel_junctions_tsv:           Classified junctions TSV (with junction_origin column).
        predictions_tsv:               MHCflurry predictions TSV.
        output_html:                   Destination HTML file.
        contigs_fasta:                 Contig FASTA for junction visualisation (optional).
        hla_qc_tsv:                    HLA QC TSV from aggregate_hla_alleles (optional).
        presentation_percentile_strong: Strong-presentation percentile threshold.
        presentation_percentile_weak:  Weak-presentation percentile threshold (quality gate).
        tcrdock_pdb:                   TCRdock-predicted PDB file (optional).
        docking_scores_tsv:            TCRdock raw docking scores TSV; consumed as input
                                       for the 3D structure manifest (optional).
        output_tsv:                    Destination for the run-summary report.tsv (optional).
        output_top_candidates_tsv:     Destination for the wide top-presenters TSV (optional).
        output_3d_structure_tsv:       Destination for the 3D structure manifest TSV
                                       (optional; only written when TCRdock is enabled).
        patient_id:                    Patient identifier written into every report.tsv row.
    """
    import pandas as pd

    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    junc_df = pd.read_csv(novel_junctions_tsv, sep="\t")
    pred_df = pd.read_csv(predictions_tsv, sep="\t")
    contigs = _load_contigs(contigs_fasta) if contigs_fasta else {}

    # --- Build origin_df (intermediate; also written into report.tsv below) ---
    if junc_df.empty or "junction_origin" not in junc_df.columns:
        origin_df = _empty_origin_df()
    else:
        origin_rows = []
        for sample_id, grp in junc_df.groupby("sample_id"):
            sample_type = grp["sample_type"].iloc[0]
            counts = grp["junction_origin"].value_counts()
            origin_rows.append({
                "sample_id": sample_id,
                "sample_type": sample_type,
                "unannotated": len(grp),
                "normal_shared": counts.get("normal_shared", 0),
                "tumor_exclusive": counts.get("tumor_exclusive", 0),
            })
        origin_df = pd.DataFrame(origin_rows)

    if patient_id:
        effective_patient_id = patient_id
    else:
        effective_patient_id = output_html.parts[-3]
        log.warning(
            "patient_id not provided; falling back to output_html path component "
            "[-3]=%r. This assumes output is …/{patient_id}/reports/report.html — "
            "if the path layout changes, every artefact row will record the wrong "
            "patient_id silently. Pass patient_id explicitly to be safe.",
            effective_patient_id,
        )

    # === Phase 2: write artefacts FIRST, then render HTML from them ===

    if output_tsv:
        _build_report_tsv(
            patient_id=effective_patient_id,
            origin_df=origin_df,
            pred_df=pred_df,
            hla_qc_tsv=hla_qc_tsv,
            output_tsv=output_tsv,
            presentation_percentile_strong=presentation_percentile_strong,
            tcrdock_pdb=tcrdock_pdb,
            presentation_percentile_weak=presentation_percentile_weak,
            junction_filter_stats_tsv=junction_filter_stats_tsv,
        )

    if output_top_candidates_tsv:
        _build_report_top_candidates_tsv(
            patient_id=effective_patient_id,
            pred_df=pred_df,
            contigs=contigs,
            output_tsv=output_top_candidates_tsv,
            presentation_percentile_weak=presentation_percentile_weak,
        )

    if output_3d_structure_tsv and docking_scores_tsv:
        # PDB lives in predictions/tcrdock/ next to the docking scores.
        # reports/ and predictions/ are siblings under results/{patient_id}/
        # so this relative path resolves from any patient root.
        pdb_relative_path = "predictions/tcrdock/top_candidate.pdb"
        _build_report_3d_structure_tsv(
            patient_id=effective_patient_id,
            docking_scores_tsv=docking_scores_tsv,
            pdb_relative_path=pdb_relative_path,
            output_tsv=output_3d_structure_tsv,
        )

    # --- Reload artefacts to drive HTML rendering ---
    report_data = _load_report_tsv(output_tsv) if output_tsv else None
    top_candidates_df = (
        pd.read_csv(output_top_candidates_tsv, sep="\t")
        if output_top_candidates_tsv else None
    )
    structure_manifest = (
        pd.read_csv(output_3d_structure_tsv, sep="\t")
        if (output_3d_structure_tsv and docking_scores_tsv) else None
    )

    # --- Junction origin summary ---
    origin_html = _df_to_html(
        report_data["junction_filtering"] if report_data is not None else origin_df
    )

    # --- Prediction summary (presentation-class counts) ---
    if report_data is not None:
        presenter_html = _presenter_counts_html(report_data["mhc_prediction"])
    elif pred_df.empty:
        presenter_html = "<p><em>No predictions available.</em></p>"
    else:
        presenter_counts = pred_df["presentation_class"].value_counts().reset_index()
        presenter_counts.columns = ["presentation_class", "count"]
        presenter_html = _df_to_html(presenter_counts)

    # --- Top presenters table with contig visualisation ---
    if top_candidates_df is not None:
        strong_html = _build_strong_table_html_from_top_candidates(top_candidates_df)
    elif pred_df.empty:
        strong_html = "<p><em>No predictions available.</em></p>"
    else:
        strong_html = _build_strong_table_html(
            pred_df, contigs,
            presentation_percentile_weak=presentation_percentile_weak,
        )

    # --- HLA typing section (optional, reads upstream artefact directly) ---
    hla_section = _build_hla_section(hla_qc_tsv) if hla_qc_tsv else ""

    # --- Filtering funnel (Issue #215) ---
    filtering_funnel = _build_filtering_funnel_html(filtering_stats_tsv)

    # --- TCRdock 3D structure viewer (optional) ---
    if tcrdock_pdb is not None and Path(tcrdock_pdb).exists():
        if structure_manifest is not None:
            peptide, allele = _resolve_top_candidate_from_manifest(structure_manifest)
        else:
            peptide, allele = _resolve_top_candidate_for_structure(
                pred_df, presentation_percentile_weak=presentation_percentile_weak,
            )
        structure_section = _build_structure_section(Path(tcrdock_pdb), peptide, allele)
    else:
        structure_section = ""

    # Only load Mol* assets when a structure section is present
    if structure_section:
        molstar_assets = (
            '<!-- Mol* molecular viewer (MIT license) — loaded only when TCRdock\n'
            '       structural validation produced a PDB for the report. -->\n'
            '  <script src="https://unpkg.com/molstar@4.9.0/build/viewer/molstar.js"></script>\n'
            '  <link rel="stylesheet" href="https://unpkg.com/molstar@4.9.0/build/viewer/molstar.css"/>'
        )
    else:
        molstar_assets = ""

    report_html = _HTML_TEMPLATE.format(
        pipeline_diagram=_PIPELINE_DIAGRAM,
        filtering_funnel=filtering_funnel,
        origin_table=origin_html,
        hla_section=hla_section,
        presenter_table=presenter_html,
        strong_table=strong_html,
        structure_section=structure_section,
        presentation_percentile_strong=presentation_percentile_strong,
        molstar_assets=molstar_assets,
    )
    output_html.write_text(report_html, encoding="utf-8")
    log.info("Report written to %s", output_html)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    generate_report(
        novel_junctions_tsv=snakemake.input.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        predictions_tsv=snakemake.input.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_html=snakemake.output.report_html,  # type: ignore[name-defined]  # noqa: F821
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        hla_qc_tsv=getattr(snakemake.input, "hla_qc", None),  # type: ignore[name-defined]  # noqa: F821
        presentation_percentile_strong=float(snakemake.params.presentation_percentile_strong),  # type: ignore[name-defined]  # noqa: F821
        presentation_percentile_weak=float(snakemake.params.presentation_percentile_weak),  # type: ignore[name-defined]  # noqa: F821
        tcrdock_pdb=getattr(snakemake.input, "pdb", None),  # type: ignore[name-defined]  # noqa: F821
        docking_scores_tsv=getattr(snakemake.input, "scores_tsv", None),  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.report_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_top_candidates_tsv=snakemake.output.report_top_candidates_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_3d_structure_tsv=getattr(snakemake.output, "report_3d_structure_tsv", None),  # type: ignore[name-defined]  # noqa: F821
        junction_filter_stats_tsv=snakemake.input.junction_filter_stats,  # type: ignore[name-defined]  # noqa: F821
        filtering_stats_tsv=getattr(snakemake.input, "filtering_stats", None),  # type: ignore[name-defined]  # noqa: F821
        patient_id=snakemake.wildcards.patient_id,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate summary HTML report from pipeline results."
    )
    parser.add_argument("--novel-junctions", required=True, help="Classified junctions TSV")
    parser.add_argument("--predictions", required=True, help="MHCflurry predictions TSV")
    parser.add_argument("--output", required=True, help="Output HTML report")
    parser.add_argument("--contigs-fasta", default=None, help="Contigs FASTA for junction visualisation")
    parser.add_argument("--hla-qc-tsv", default=None, help="HLA QC TSV from aggregate_hla_alleles")
    parser.add_argument("--tcrdock-pdb", default=None, help="TCRdock PDB file for 3D viewer")
    parser.add_argument("--presentation-percentile-strong", type=float, default=0.5)
    parser.add_argument("--presentation-percentile-weak", type=float, default=2.0)
    parser.add_argument("--output-tsv", default=None, help="Output machine-readable summary TSV")
    parser.add_argument("--junction-filter-stats", default=None,
                        help="Per-tumor-sample junction funnel stats TSV (Issue #214)")
    parser.add_argument("--filtering-stats", default=None,
                        help="Aggregated cross-step filtering audit TSV (Issue #215)")
    parser.add_argument("--patient-id", default="", help="Patient identifier for report.tsv rows")
    args = parser.parse_args()

    generate_report(
        novel_junctions_tsv=args.novel_junctions,
        predictions_tsv=args.predictions,
        output_html=args.output,
        contigs_fasta=args.contigs_fasta,
        hla_qc_tsv=args.hla_qc_tsv,
        presentation_percentile_strong=args.presentation_percentile_strong,
        presentation_percentile_weak=args.presentation_percentile_weak,
        tcrdock_pdb=args.tcrdock_pdb,
        output_tsv=args.output_tsv,
        junction_filter_stats_tsv=args.junction_filter_stats,
        filtering_stats_tsv=args.filtering_stats,
        patient_id=args.patient_id,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
