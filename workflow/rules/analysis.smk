# =============================================================================
# Rule module: Step 6 — Report generation
# =============================================================================

rule generate_report:
    """Generate a summary HTML report showing:
    - Junction origin counts (tumor_specific vs patient_specific per sample)
    - Neoepitope prediction summary (strong / weak / non binder counts)
    - Top strong binders table (IC50 < ic50_strong threshold)"""
    input:
        novel_junctions=rules.filter_junctions.output.novel_junctions,
        predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
        contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
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
