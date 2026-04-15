# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================


def _hla_qc_input(wildcards):
    """Return the hla_qc.tsv path when HLA typing is enabled, else []."""
    if (
        config.get("hla", {}).get("enabled", False)
        and config.get("data_source") == "fastq"
    ):
        hla_dir = os.path.join(os.path.dirname(OUT["raw_data"]), "hla_typing")
        return os.path.join(hla_dir, wildcards.patient_id, "hla_qc.tsv")
    return []


rule generate_report:
    """Generate a summary HTML report showing:
    - Junction origin counts (tumor_specific vs patient_specific per sample)
    - HLA allele typing summary (when HLA typing is enabled)
    - Neoepitope prediction summary (strong / weak / non binder counts)
    - Top strong binders table (IC50 < ic50_strong threshold)"""
    input:
        novel_junctions=rules.filter_junctions.output.novel_junctions,
        predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
        contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
        hla_qc=_hla_qc_input,
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
