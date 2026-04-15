#!/usr/bin/env python3
"""generate_report.py — Generate a summary HTML report for pipeline results.

Produces a self-contained HTML page showing:
  1. Junction origin summary — counts of tumor_specific vs patient_specific
     junctions per sample, so results can be quickly verified.
  2. Top strong binders — the highest-affinity predicted neoepitopes
     (IC50 < strong threshold).

Usage (standalone):
  python generate_report.py \\
      --novel-junctions results/junctions/local/novel_junctions.tsv \\
      --predictions results/predictions/local/predictions.tsv \\
      --output results/reports/local/report.html

Usage (Snakemake):
  Called automatically by the ``generate_report`` rule.
"""

import argparse
import html as html_mod
import json
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

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
<pre class="mermaid">
graph TD
    A["🔬 RNA-Seq reads<br/>(tumor + normal)"]
    B["HISAT2 alignment<br/>(regtools junction extraction)"]
    C["GENCODE reference filter<br/>(remove annotated junctions)"]
    D["Unannotated junctions"]
    E["Normal sample comparison"]
    F["patient_specific<br/>(also in normal) ✗"]
    G["tumor_specific<br/>(absent from normal) ✓"]
    H["Contig assembly<br/>(50 nt per junction)"]
    I["In-silico translation<br/>(3 reading frames → 9-mers)"]
    J["HLA typing<br/>(arcasHLA)"]
    K["Patient alleles<br/>(alleles.tsv)"]
    L["MHCflurry epitope prediction<br/>(IC50 affinity per 9-mer)"]
    M["📊 Neoepitope candidates<br/>(strong/weak binders)"]
    
    A --> B --> C --> D
    D --> E
    E --> F
    E --> G
    F -.->|excluded| E
    G --> H --> I
    I --> L
    B --> J --> K --> L
    L --> M
    
    style A fill:#eafaf1,stroke:#27ae60,stroke-width:2px
    style M fill:#eafaf1,stroke:#27ae60,stroke-width:2px
    style F fill:#fdecea,stroke:#c0392b,stroke-width:2px
    style G fill:#eafaf1,stroke:#27ae60,stroke-width:2px
    style K fill:#f5eeff,stroke:#8e44ad,stroke-width:2px
    style J fill:#f5eeff,stroke:#8e44ad,stroke-width:2px
</pre>
"""

_MOLSTAR_VIEWER = """\
<div id="molstar-container" style="width:100%;height:500px;position:relative;border:1px solid #ddd;border-radius:6px;overflow:hidden;margin:1em 0;"></div>
<script type="text/javascript">
  // Mol* viewer — loads the inlined PDB and renders the TCR-peptide-MHC complex.
  // Mol* is MIT-licensed and loaded from the official CDN (no registration required).
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
  {molstar_assets}
  <script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
  <style>
    body  {{ font-family: Arial, sans-serif; margin: 2em; color: #333; }}
    h1    {{ color: #2c3e50; }}
    h2    {{ color: #34495e; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
    th    {{ background-color: #f2f2f2; }}
    tr:nth-child(even) {{ background-color: #f9f9f9; }}

    /* Mermaid diagram */
    .mermaid {{ display: flex; justify-content: center; margin: 1.5em 0; }}

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

  <h2>Junction origin summary</h2>
  <p>Splice junctions detected in the tumor that are absent from the GENCODE
  reference annotation, grouped by whether they also appear in the matched normal
  tissue. Only <strong class="keep">tumor_specific</strong> junctions — absent
  from both the reference and the normal sample — are used for neoepitope
  prediction. <strong class="discard">patient_specific</strong> junctions are
  present in normal tissue and are therefore not tumor-derived.</p>
  {origin_table}

  {hla_section}

  <h2>Neoepitope prediction summary</h2>
  {binder_table}

  <h2>Top strong binders (IC50 &lt; {ic50_strong} nM)</h2>
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


def _load_contigs(contigs_fasta: str | Path) -> dict[str, str]:
    """Parse contigs FASTA into {contig_key: sequence} dict.

    The contig key matches the contig_key column in the predictions TSV.
    """
    contigs: dict[str, str] = {}
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        contigs[record.description] = str(record.seq).upper()
    return contigs


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
) -> str:
    """Build the top strong binders table with a contig visualisation column."""
    strong_df = pred_df[pred_df["binder_class"] == "strong"].sort_values("ic50_nM")
    if strong_df.empty:
        return "<p><em>No strong binders found.</em></p>"

    rows = []
    for _, row in strong_df.head(max_rows).iterrows():
        contig_key = row["contig_key"]
        start_nt = int(row["start_nt"])
        end_nt_incl = start_nt + 26  # 9 aa × 3 nt − 1

        seq = contigs.get(contig_key, "")
        contig_html = _render_contig(seq, start_nt, end_nt_incl, upstream_nt) if seq else "<em>n/a</em>"

        rows.append(
            f"<tr>"
            f"<td>{html_mod.escape(str(row['contig_key']))}</td>"
            f"<td>{html_mod.escape(str(row['peptide']))}</td>"
            f"<td>{html_mod.escape(str(row['allele']))}</td>"
            f"<td>{row['ic50_nM']:.1f}</td>"
            f"<td>{row['percentile_rank']:.3f}</td>"
            f"<td>{contig_html}</td>"
            f"</tr>"
        )

    legend = (
        "upstream <span class='junction-mark'>|</span> downstream — "
        "<span class='nt-pep-up'>peptide (upstream)</span> "
        "<span class='nt-pep-down'>peptide (downstream)</span>"
    )
    percentile_header = (
        "<th title='Percentage of random peptides that bind worse than this one. "
        "Lower is better — &lt;2% is the standard binder threshold.'>Percentile rank ▾</th>"
    )
    table = (
        f"<table><thead><tr>"
        f"<th>Source</th><th>9-mer</th><th>Allele</th><th>IC50 (nM)</th>"
        f"{percentile_header}"
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


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _build_hla_section(hla_qc_tsv: "str | Path | None") -> str:
    """Build the HLA allele typing section for the report.

    Args:
        hla_qc_tsv: Path to hla_qc.tsv from aggregate_hla_alleles, or None
                    when HLA typing was not run.

    Returns:
        HTML string for the HLA typing section.
    """
    heading = "<h2>HLA allele typing</h2>"
    if hla_qc_tsv is None:
        return (
            heading
            + "<p><em>HLA typing was not run (config.hla.enabled is false or "
            "data_source is not fastq). MHCflurry used the fallback alleles "
            "defined in config.</em></p>"
        )

    try:
        qc_df = pd.read_csv(hla_qc_tsv, sep="\t")
    except (FileNotFoundError, pd.errors.ParserError) as exc:
        log.warning("Could not read HLA QC file %s: %s", hla_qc_tsv, exc)
        return heading + "<p><em>HLA QC file not available.</em></p>"

    if qc_df.empty:
        return heading + "<p><em>No HLA allele data available.</em></p>"

    _SOURCE_STYLE = {
        "normal":   "color:#27ae60;font-weight:bold",
        "tumor":    "color:#e67e22;font-weight:bold",
        "fallback": "color:#c0392b;font-weight:bold",
    }
    rows = []
    has_discrepancy = False
    for _, row in qc_df.iterrows():
        source = str(row.get("source", ""))
        src_style = _SOURCE_STYLE.get(source, "")
        discrepancy = str(row.get("discrepancy", "") or "")
        if discrepancy:
            has_discrepancy = True
            disc_cell = (
                f'<td><span style="color:#c0392b">'
                f"{html_mod.escape(discrepancy)}</span></td>"
            )
        else:
            disc_cell = "<td>—</td>"
        reads = row.get("reads", 0)
        reads_str = f"{int(reads):,}" if reads else "—"
        rows.append(
            f"<tr>"
            f"<td><strong>HLA-{html_mod.escape(str(row['locus']))}</strong></td>"
            f"<td><code>{html_mod.escape(str(row.get('allele1', '')))}</code></td>"
            f"<td><code>{html_mod.escape(str(row.get('allele2', '')))}</code></td>"
            f"<td><span style='{src_style}'>{html_mod.escape(source)}</span></td>"
            f"<td>{reads_str}</td>"
            f"{disc_cell}"
            f"</tr>"
        )

    table = (
        "<table><thead><tr>"
        "<th>Locus</th><th>Allele 1</th><th>Allele 2</th>"
        "<th title='normal = from matched normal tissue (preferred); "
        "tumor = from tumor (normal unavailable); "
        "fallback = config default (no confident call)'>Source</th>"
        "<th title='Reads used for genotyping at this locus'>Reads</th>"
        "<th>Normal/tumor discrepancy</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )

    notes = ""
    if has_discrepancy:
        notes = (
            "<p style='color:#c0392b'><strong>Warning:</strong> normal/tumor "
            "discrepancy detected at one or more loci. The normal-tissue call "
            "was used (normal-first policy).</p>"
        )

    return (
        heading
        + "<p>Patient-specific HLA-A/B/C alleles typed from RNA-Seq reads "
        "using <a href='https://github.com/RabadanLab/arcasHLA'>arcasHLA</a>. "
        "Normal tissue is preferred over tumor (germline HLA alleles are less "
        "susceptible to LOH). Alleles below the minimum read threshold are "
        "replaced by config-defined fallback alleles.</p>"
        + table
        + notes
    )


def _build_structure_section(pdb_path: Path, pred_df: pd.DataFrame) -> str:
    """Build the Mol* 3D viewer section for the top TCRdock candidate.

    The PDB is inlined as a JSON string so the HTML report is fully
    self-contained and does not depend on local file paths at view time.

    Args:
        pdb_path: Path to the TCRdock-predicted PDB file.
        pred_df:  MHCflurry predictions DataFrame (used for annotation).

    Returns:
        HTML string containing the Mol* viewer and inlined PDB.
    """
    try:
        pdb_text = pdb_path.read_text()
    except OSError as exc:
        log.warning("Could not read PDB file %s: %s", pdb_path, exc)
        return "<p><em>Structure not available.</em></p>"

    # Annotate with the top candidate's peptide and allele.
    # Match run_tcrdock.py's candidate selection: best strong binder, falling
    # back to weak. Exclude non-binders to stay in sync.
    binders = pred_df[pred_df["binder_class"].isin(("strong", "weak"))]
    top = binders.sort_values("ic50_nM").head(1)
    if top.empty:
        peptide, allele = "NA", "NA"
    else:
        peptide = html_mod.escape(str(top.iloc[0]["peptide"]))
        allele  = html_mod.escape(str(top.iloc[0]["allele"]))

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


def generate_report(
    novel_junctions_tsv: str | Path,
    predictions_tsv: str | Path,
    output_html: str | Path,
    contigs_fasta: str | Path | None = None,
    ic50_strong: float = 50.0,
    tcrdock_pdb: str | Path | None = None,
    hla_qc_tsv: str | Path | None = None,
) -> None:
    """Generate the summary HTML report.

    Args:
        novel_junctions_tsv: Classified junctions TSV (with junction_origin column).
        predictions_tsv:     MHCflurry predictions TSV.
        output_html:         Destination HTML file.
        contigs_fasta:       Contig FASTA for junction visualisation (optional).
        ic50_strong:         Strong-binder IC50 threshold (nM).
        tcrdock_pdb:         TCRdock-predicted PDB file (optional). When provided,
                             an embedded Mol* 3D viewer is added to the report.
        hla_qc_tsv:          HLA QC TSV from aggregate_hla_alleles (optional).
                             When provided, an HLA allele section is added.
    """
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    junc_df = pd.read_csv(novel_junctions_tsv, sep="\t")
    pred_df = pd.read_csv(predictions_tsv, sep="\t")
    contigs = _load_contigs(contigs_fasta) if contigs_fasta else {}

    # --- Junction origin summary ---
    if junc_df.empty or "junction_origin" not in junc_df.columns:
        origin_df = pd.DataFrame(
            columns=["sample_id", "sample_type", "unannotated", "patient_specific", "tumor_specific"]
        )
    else:
        origin_rows = []
        for sample_id, grp in junc_df.groupby("sample_id"):
            sample_type = grp["sample_type"].iloc[0]
            counts = grp["junction_origin"].value_counts()
            origin_rows.append({
                "sample_id": sample_id,
                "sample_type": sample_type,
                "unannotated": len(grp),
                "patient_specific": counts.get("patient_specific", 0),
                "tumor_specific": counts.get("tumor_specific", 0),
            })
        origin_df = pd.DataFrame(origin_rows)

    origin_html = _df_to_html(origin_df)

    # --- Binder summary ---
    if pred_df.empty:
        binder_html = "<p><em>No predictions available.</em></p>"
    else:
        binder_counts = pred_df["binder_class"].value_counts().reset_index()
        binder_counts.columns = ["binder_class", "count"]
        binder_html = _df_to_html(binder_counts)

    # --- Top strong binders with contig visualisation ---
    if pred_df.empty:
        strong_html = "<p><em>No predictions available.</em></p>"
    else:
        strong_html = _build_strong_table_html(pred_df, contigs)

    # --- HLA allele typing section (optional) ---
    hla_section = _build_hla_section(hla_qc_tsv)

    # --- TCRdock 3D structure viewer (optional) ---
    if tcrdock_pdb is not None and Path(tcrdock_pdb).exists():
        structure_section = _build_structure_section(Path(tcrdock_pdb), pred_df)
    else:
        structure_section = ""

    # Only load Mol* assets when a structure section is present
    if structure_section:
        molstar_assets = (
            '<!-- Mol* molecular viewer (MIT license) — loaded only when TCRdock\n'
            '       structural validation produced a PDB for the report. -->\n'
            '  <script src="https://www.unpkg.com/molstar/build/viewer/molstar.js"></script>\n'
            '  <link rel="stylesheet" href="https://www.unpkg.com/molstar/build/viewer/molstar.css"/>'
        )
    else:
        molstar_assets = ""

    report_html = _HTML_TEMPLATE.format(
        pipeline_diagram=_PIPELINE_DIAGRAM,
        origin_table=origin_html,
        hla_section=hla_section,
        binder_table=binder_html,
        strong_table=strong_html,
        structure_section=structure_section,
        ic50_strong=int(ic50_strong),
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

    _hla_qc = getattr(snakemake.input, "hla_qc", None)  # type: ignore[name-defined]  # noqa: F821
    generate_report(
        novel_junctions_tsv=snakemake.input.novel_junctions,  # type: ignore[name-defined]  # noqa: F821
        predictions_tsv=snakemake.input.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_html=snakemake.output.report_html,  # type: ignore[name-defined]  # noqa: F821
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
        tcrdock_pdb=getattr(snakemake.input, "pdb", None),  # type: ignore[name-defined]  # noqa: F821
        hla_qc_tsv=_hla_qc if _hla_qc else None,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate summary HTML report from pipeline results."
    )
    parser.add_argument("--novel-junctions", required=True, help="Classified junctions TSV")
    parser.add_argument("--predictions", required=True, help="MHCflurry predictions TSV")
    parser.add_argument("--output", required=True, help="Output HTML report")
    parser.add_argument("--contigs-fasta", default=None, help="Contigs FASTA for junction visualisation")
    parser.add_argument("--tcrdock-pdb", default=None, help="TCRdock PDB file for 3D viewer")
    parser.add_argument("--hla-qc", default=None, help="HLA QC TSV from aggregate_hla_alleles")
    parser.add_argument("--ic50-strong", type=float, default=50.0)
    args = parser.parse_args()

    generate_report(
        novel_junctions_tsv=args.novel_junctions,
        predictions_tsv=args.predictions,
        output_html=args.output,
        contigs_fasta=args.contigs_fasta,
        ic50_strong=args.ic50_strong,
        tcrdock_pdb=args.tcrdock_pdb,
        hla_qc_tsv=args.hla_qc,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
