# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================

_HLA_QC_ENABLED = config.get("hla", {}).get("enabled", False)


def _generate_report_input(wildcards):
    d = {
        "novel_junctions": os.path.join(
            OUT["junctions"], wildcards.patient_id, "novel_junctions.tsv"
        ),
        "predictions_tsv": os.path.join(
            OUT["predictions"], wildcards.patient_id, "mhc_affinity.tsv"
        ),
        "contigs_fasta": os.path.join(
            OUT["contigs"], wildcards.patient_id, "contigs.fa"
        ),
    }
    if _HLA_QC_ENABLED:
        d["hla_qc"] = os.path.join(
            os.path.dirname(OUT["raw_data"]), "hla_typing",
            wildcards.patient_id, "hla_qc.tsv",
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
            OUT["reports"], "{patient_id}", "report.html"
        ),
    log:
        os.path.join(OUT["logs"], "report", "{patient_id}_report.log"),
    params:
        ic50_strong=config["mhcflurry"]["ic50_strong"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_report.py"
