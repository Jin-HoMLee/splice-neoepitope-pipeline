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


def _mhcflurry_alleles_input(wildcards):
    """Return alleles.tsv path when HLA typing is enabled, else empty dict."""
    if config.get("hla", {}).get("enabled", False):
        hla_typing_dir = os.path.join(os.path.dirname(OUT["raw_data"]), "hla_typing")
        return {"alleles_tsv": os.path.join(hla_typing_dir, wildcards.patient_id, "alleles.tsv")}
    return {}


rule run_mhcflurry:
    """Run MHCflurry 2.x on junction-spanning 9-mers.

    When hla.enabled is true, reads patient-specific alleles from the
    alleles.tsv produced by aggregate_hla_alleles and runs predictions for
    all typed alleles. Otherwise uses the fallback allele from config.

    Output: TSV with columns: contig_key, start_nt, peptide, allele, IC50_nM,
            percentile_rank, binder_class (strong / weak / non)."""
    input:
        unpack(_mhcflurry_alleles_input),
        peptides_tsv=rules.translate_peptides.output.peptides_tsv,
        mhcflurry_models=rules.download_mhcflurry_models.output.sentinel,
    output:
        predictions_tsv=os.path.join(
            OUT["predictions"], "{patient_id}", "predictions.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "predict", "{patient_id}_predict.log"),
    params:
        fallback_alleles=config["hla"]["fallback_alleles"],
        ic50_strong=config["mhcflurry"]["ic50_strong"],
        ic50_weak=config["mhcflurry"]["ic50_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_mhcflurry.py"
