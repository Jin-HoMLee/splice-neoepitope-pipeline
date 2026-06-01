# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================

_HLA_QC_ENABLED = config.get("hla", {}).get("enabled", False)
_TCRDOCK_ENABLED = config.get("tcrdock", {}).get("enabled", False)
_PROTEOME_FILTER_ENABLED_REPORT = config.get("proteome_filter", {}).get("enabled", True)


def _aggregate_filtering_stats_input(wildcards):
    d = {
        "junction_filter": os.path.join(
            _RES, wildcards.patient_id, "junctions", "junction_filter_stats.tsv"
        ),
        "contig_assemble": os.path.join(
            _RES, wildcards.patient_id, "contigs", "contig_assemble_stats.tsv"
        ),
        "translate": os.path.join(
            _RES, wildcards.patient_id, "peptides", "translate_stats.tsv"
        ),
        "mhc": os.path.join(
            _RES, wildcards.patient_id, "predictions", "mhc_stats.tsv"
        ),
    }
    if _PROTEOME_FILTER_ENABLED_REPORT:
        d["proteome"] = os.path.join(
            _RES, wildcards.patient_id, "peptides", "proteome_stats.tsv"
        )
    return d


rule aggregate_filtering_stats:
    """Concatenate the per-step funnel stats TSVs into a single per-patient
    ``filtering_stats.tsv`` (Issue #215). Unified schema:
    ``patient_id, sample_id, sample_type, step, category, count``.
    The ``proteome`` input is omitted when ``proteome_filter.enabled: false``."""
    input:
        unpack(_aggregate_filtering_stats_input),
    output:
        filtering_stats=os.path.join(
            _RES, "{patient_id}", "reports", "filtering_stats.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "analysis", "aggregate_filtering_stats.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/aggregate_filtering_stats.py"


def _generate_report_input(wildcards):
    d = {
        "novel_junctions": os.path.join(
            _RES, wildcards.patient_id, "junctions", "novel_junctions.tsv"
        ),
        "junction_filter_stats": os.path.join(
            _RES, wildcards.patient_id, "junctions", "junction_filter_stats.tsv"
        ),
        "filtering_stats": os.path.join(
            _RES, wildcards.patient_id, "reports", "filtering_stats.tsv"
        ),
        "predictions_tsv": os.path.join(
            _RES, wildcards.patient_id, "predictions", "mhc_presentation.tsv"
        ),
        "contigs_fasta": os.path.join(
            _RES, wildcards.patient_id, "contigs", "contigs.fa"
        ),
    }
    if _HLA_QC_ENABLED:
        d["hla_qc"] = os.path.join(
            _RES, wildcards.patient_id, "hla_typing", "hla_qc.tsv",
        )
        # VDJdb reference panel (Issue #204/#206) — the fetch_vdjdb_panel rule
        # that produces these is gated on config[hla][enabled], same as hla_qc.
        d["vdjdb_panel"] = os.path.join(
            _RES, wildcards.patient_id, "tcr_panel", "vdjdb", "panel.tsv",
        )
        d["vdjdb_panel_qc"] = os.path.join(
            _RES, wildcards.patient_id, "tcr_panel", "vdjdb", "panel_qc.tsv",
        )
    if _TCRDOCK_ENABLED:
        d["pdb"] = rules.run_tcrdock.output.pdb.format(patient_id=wildcards.patient_id)
        d["scores_tsv"] = rules.run_tcrdock.output.scores_tsv.format(
            patient_id=wildcards.patient_id
        )
    return d


_generate_report_output = {
    "report_html": os.path.join(_RES, "{patient_id}", "reports", "report.html"),
    "report_tsv": os.path.join(_RES, "{patient_id}", "reports", "report.tsv"),
    "report_top_candidates_tsv": os.path.join(
        _RES, "{patient_id}", "reports", "report_top_candidates.tsv"
    ),
}
if _TCRDOCK_ENABLED:
    _generate_report_output["report_3d_structure_tsv"] = os.path.join(
        _RES, "{patient_id}", "reports", "report_3d_structure.tsv"
    )


rule generate_report:
    """Generate a summary HTML report showing:
    - Junction origin counts (tumor_exclusive vs normal_shared per sample)
    - HLA typing results with source and normal/tumor concordance (when enabled)
    - Neoepitope prediction summary (strong / weak / non presenter counts)
    - Top strong presenters table (presentation_percentile ≤ strong threshold)
    - VDJdb reference TCR panel + per-allele coverage stats (when HLA enabled; Issue #206)
    - Embedded Mol* 3D viewer for the top TCR-pMHC candidate, with the selected TCR's
      provenance (HLA-matched VDJdb or DMF5-fallback flag) + TCRdock confidence metrics
      (when TCRdock enabled; Issue #205/#206)"""
    input:
        unpack(_generate_report_input),
    output:
        **_generate_report_output,
    log:
        os.path.join(_LOGS, "{patient_id}", "analysis", "report.log"),
    params:
        presentation_percentile_strong=config["mhcflurry"]["presentation_percentile_strong"],
        presentation_percentile_weak=config["mhcflurry"]["presentation_percentile_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_report.py"
