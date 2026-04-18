# =============================================================================
# Rule module: Step 5 — Epitope prediction with MHCflurry 2.x
# =============================================================================

_HLA_TYPING_ENABLED = config.get("hla", {}).get("enabled", False)

# When HLA typing is enabled, the patient-specific alleles.tsv produced by
# aggregate_hla_alleles is wired in as an input so Snakemake tracks the
# dependency and run_mhcflurry.py reads alleles from it.
# When disabled, the fallback alleles from mhcflurry.fallback_alleles are
# passed directly as a parameter.
def _run_mhcflurry_input(wildcards):
    d = {
        "peptides_tsv": os.path.join(_RES, wildcards.patient_id, "peptides", "peptides.tsv"),
        "mhcflurry_models": rules.download_mhcflurry_models.output.sentinel,
    }
    if _HLA_TYPING_ENABLED:
        d["alleles_tsv"] = os.path.join(
            _RES, wildcards.patient_id, "hla_typing", "alleles.tsv",
        )
    return d


rule download_mhcflurry_models:
    """Download MHCflurry trained models (~1 GB) on first run.
    Models are cached in ~/.local/share/mhcflurry/ and reused across runs.
    A sentinel file is written to resources/ to prevent re-downloading."""
    output:
        sentinel=touch("resources/mhcflurry_models.done"),
    log:
        os.path.join(_LOGS, "predict", "mhcflurry_downloads.log"),
    conda:
        "../envs/python.yaml"
    shell:
        "mhcflurry-downloads fetch > {log} 2>&1"


rule run_mhcflurry:
    """Run MHCflurry 2.x on junction-spanning 9-mers.
    Output: TSV with columns: contig_key, start_nt, peptide, allele, IC50_nM,
            percentile_rank, binder_class (strong / weak / non).
    When hla.enabled is true, predictions are run for all patient-specific
    alleles from alleles.tsv. Otherwise mhcflurry.fallback_alleles are used."""
    input:
        unpack(_run_mhcflurry_input),
    output:
        mhc_affinity_tsv=os.path.join(
            _RES, "{patient_id}", "predictions", "mhc_affinity.tsv"
        ),
    log:
        os.path.join(_LOGS, "predict", "{patient_id}_predict.log"),
    params:
        # Used only when hla.enabled is false (alleles_tsv input absent).
        fallback_alleles=list(config["mhcflurry"]["fallback_alleles"].values()),
        ic50_strong=config["mhcflurry"]["ic50_strong"],
        ic50_weak=config["mhcflurry"]["ic50_weak"],
        threads=config["mhcflurry"]["threads"],
    threads: config["mhcflurry"]["threads"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_mhcflurry.py"
