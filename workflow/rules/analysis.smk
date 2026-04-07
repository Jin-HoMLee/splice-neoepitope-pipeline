# =============================================================================
# Rule module: Step 6 — Statistical analysis and report generation
# =============================================================================

rule statistical_analysis:
    """Count epitopes per sample; perform one-tailed Fisher's exact test for
    tumour-vs-normal enrichment (BRCA, LUAD); calculate normalised counts.
    Output: per-cancer-type statistics TSV."""
    input:
        predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
        manifest=os.path.join(OUT["raw_data"], "{cancer_type}", "manifest.tsv"),
    output:
        stats_tsv=os.path.join(
            OUT["analysis"], "{cancer_type}", "statistics.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "analysis", "{cancer_type}_stats.log"),
    params:
        cancer_types_with_normal=config["cancer_types_with_normal"],
        ic50_strong=config["mhcflurry"]["ic50_strong"],
        ic50_weak=config["mhcflurry"]["ic50_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/statistical_analysis.py"


rule generate_report:
    """Generate an HTML report with:
    - Epitope count distributions (box plots)
    - -log10(p) visualisations from Fisher's exact test
    - Summary tables of strong/weak binders
    - Annotation with junction coordinates and gene names"""
    input:
        stats_tsv=rules.statistical_analysis.output.stats_tsv,
        predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
        novel_junctions=rules.filter_junctions.output.novel_junctions,
    output:
        report_html=os.path.join(
            OUT["reports"], "{cancer_type}", "report.html"
        ),
    log:
        os.path.join(OUT["logs"], "report", "{cancer_type}_report.log"),
    params:
        cancer_type=lambda wc: wc.cancer_type,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_report.py"


rule summarise_all:
    """Aggregate per-cancer-type statistics into a single summary TSV."""
    input:
        stats=expand(
            os.path.join(OUT["analysis"], "{cancer_type}", "statistics.tsv"),
            cancer_type=get_data_types(),
        ),
    output:
        summary=os.path.join(OUT["reports"], "summary_table.tsv"),
    log:
        os.path.join(OUT["logs"], "report", "summary.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/statistical_analysis.py"
