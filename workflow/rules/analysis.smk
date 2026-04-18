# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================

_HLA_QC_ENABLED = config.get("hla", {}).get("enabled", False)


def _generate_report_input(wildcards):
    d = {
        "novel_junctions": os.path.join(
            _RES, wildcards.patient_id, "junctions", "novel_junctions.tsv"
        ),
        "predictions_tsv": os.path.join(
            _RES, wildcards.patient_id, "predictions", "mhc_affinity.tsv"
        ),
        "contigs_fasta": os.path.join(
            _RES, wildcards.patient_id, "contigs", "contigs.fa"
        ),
    }
    if _HLA_QC_ENABLED:
        d["hla_qc"] = os.path.join(
            _RES, wildcards.patient_id, "hla_typing", "hla_qc.tsv",
        )
    return d


rule generate_report:
    """Generate a summary HTML report showing:
    - Junction origin counts (tumor_exclusive vs normal_shared per sample)
    - HLA typing results with source and normal/tumor concordance (when enabled)
    - Neoepitope prediction summary (strong / weak / non binder counts)
    - Top strong binders table (IC50 < ic50_strong threshold)"""
    input:
        unpack(_generate_report_input),
    output:
        report_html=os.path.join(
            _RES, "{patient_id}", "reports", "report.html"
        ),
    log:
        os.path.join(_LOGS, "report", "{patient_id}_report.log"),
    params:
        ic50_strong=config["mhcflurry"]["ic50_strong"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_report.py"
