# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================

_HLA_QC_ENABLED = config.get("hla", {}).get("enabled", False)
_TCRDOCK_ENABLED = config.get("tcrdock", {}).get("enabled", False)


def _generate_report_input(wildcards):
    d = {
        "novel_junctions": os.path.join(
            _RES, wildcards.patient_id, "junctions", "novel_junctions.tsv"
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
    if _TCRDOCK_ENABLED:
        d["pdb"] = os.path.join(
            _RES, wildcards.patient_id, "predictions", "tcrdock", "top_candidate.pdb"
        )
        d["scores_tsv"] = os.path.join(
            _RES, wildcards.patient_id, "predictions", "tcrdock", "docking_scores.tsv"
        )
    return d


rule generate_report:
    """Generate a summary HTML report showing:
    - Junction origin counts (tumor_exclusive vs normal_shared per sample)
    - HLA typing results with source and normal/tumor concordance (when enabled)
    - Neoepitope prediction summary (strong / weak / non binder counts)
    - Top strong binders table (presentation_percentile ≤ strong threshold)"""
    input:
        unpack(_generate_report_input),
    output:
        report_html=os.path.join(
            _RES, "{patient_id}", "reports", "report.html"
        ),
        report_tsv=os.path.join(
            _RES, "{patient_id}", "reports", "report.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "report", "report.log"),
    params:
        presentation_percentile_strong=config["mhcflurry"]["presentation_percentile_strong"],
        presentation_percentile_weak=config["mhcflurry"]["presentation_percentile_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_report.py"
