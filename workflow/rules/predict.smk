# =============================================================================
# Rule module: Step 5 — Epitope prediction with MHCflurry 2.x
# =============================================================================

rule download_mhcflurry_models:
    """Download MHCflurry trained models (~1 GB) on first run.
    Models are cached in ~/.local/share/mhcflurry/ and reused across runs.
    A sentinel file is written to resources/ to prevent re-downloading."""
    output:
        sentinel=touch("resources/mhcflurry_models.done"),
    log:
        os.path.join(OUT["logs"], "predict", "mhcflurry_downloads.log"),
    conda:
        "../envs/python.yaml"
    shell:
        "mhcflurry-downloads fetch > {log} 2>&1"


rule run_mhcflurry:
    """Run MHCflurry 2.x on junction-spanning 9-mers.
    Output: TSV with columns: contig_key, start_nt, peptide, allele, IC50_nM,
            percentile_rank, binder_class (strong / weak / non)."""
    input:
        peptides_tsv=rules.translate_peptides.output.peptides_tsv,
        mhcflurry_models=rules.download_mhcflurry_models.output.sentinel,
    output:
        predictions_tsv=os.path.join(
            OUT["predictions"], "{cancer_type}", "predictions.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "predict", "{cancer_type}_predict.log"),
    params:
        hla_allele=config["mhcflurry"]["hla_allele"],
        ic50_strong=config["mhcflurry"]["ic50_strong"],
        ic50_weak=config["mhcflurry"]["ic50_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_mhcflurry.py"
