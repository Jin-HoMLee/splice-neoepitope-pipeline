# =============================================================================
# Rule module: Step 5 — Epitope prediction with MHCflurry 2.x
# =============================================================================

_HLA_TYPING_ENABLED  = config.get("hla", {}).get("enabled", False)
_PROTEOME_FILTER_ENABLED = config.get("proteome_filter", {}).get("enabled", True)

# nf-core layout (Issue #63): sentinel for "mhcflurry-downloads fetch ran" lives
# outside the repo (default ~/.mhcflurry/.download_done). The actual model cache
# remains in MHCflurry's own platformdirs location — wiring MHCFLURRY_DATA_DIR to
# fully co-locate the cache is a follow-up (touches run_mhcflurry.py too).
_MHCFLURRY_MODELS_DIR = os.path.expanduser(
    config.get("mhcflurry", {}).get("models_cache_dir", "~/.mhcflurry")
)
_MHCFLURRY_SENTINEL = os.path.join(_MHCFLURRY_MODELS_DIR, ".download_done")

# When proteome_filter is enabled, use the filtered novel peptides TSV.
# When disabled, fall back to the raw translated peptides TSV.
def _run_mhcflurry_input(wildcards):
    if _PROTEOME_FILTER_ENABLED:
        peptides_tsv = os.path.join(
            _RES, wildcards.patient_id, "peptides", "peptides_novel.tsv"
        ) # mirrors proteome_filter_peptides.output.novel_tsv
    else:
        peptides_tsv = os.path.join(
            _RES, wildcards.patient_id, "peptides", "peptides.tsv"
        ) # mirrors translate_junctions.output.peptides_tsv
    d = {
        "peptides_tsv": peptides_tsv,
        "mhcflurry_models": rules.download_mhcflurry_models.output.sentinel,
    }
    if _HLA_TYPING_ENABLED:
        d["alleles_tsv"] = os.path.join(
            _RES, wildcards.patient_id, "hla_typing", "alleles.tsv",
        ) # mirrors hla_typing.output.alleles_tsv
    return d


rule download_mhcflurry_models:
    """Download MHCflurry trained models (~1 GB) on first run.
    Models are cached in MHCflurry's default platformdirs location
    (~/.local/share/mhcflurry/) and reused across runs. The sentinel marker
    lives in mhcflurry.models_cache_dir (default ~/.mhcflurry/) — out of the
    repo per nf-core convention (Issue #63)."""
    output:
        sentinel=touch(_MHCFLURRY_SENTINEL),
    log:
        os.path.join(_LOGS, "mhc_affinity", "mhcflurry_downloads.log"),
    conda:
        "../envs/python.yaml"
    shell:
        "mhcflurry-downloads fetch > {log} 2>&1"


rule run_mhcflurry:
    """Run MHCflurry 2.x on junction-spanning peptides.
    Output: TSV with per-allele presentation scores plus genotype_presentation_score,
            n_strong_alleles, and best_presentation_percentile.
    When hla.enabled is true, predictions are run for all patient-specific
    alleles from alleles.tsv. Otherwise mhcflurry.fallback_alleles are used."""
    input:
        unpack(_run_mhcflurry_input),
    output:
        mhc_presentation_tsv=os.path.join(
            _RES, "{patient_id}", "mhcflurry", "presentation.tsv"
        ),
        stats=os.path.join(
            _RES, "{patient_id}", "mhcflurry", "stats.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "mhc_affinity", "predict.log"),
    params:
        # Used only when hla.enabled is false (alleles_tsv input absent).
        fallback_alleles=list(config["mhcflurry"]["fallback_alleles"].values()),
        presentation_percentile_strong=config["mhcflurry"]["presentation_percentile_strong"],
        presentation_percentile_weak=config["mhcflurry"]["presentation_percentile_weak"],
        hla_c_weight=config["mhcflurry"]["hla_c_weight"],
    threads: config["mhcflurry"]["threads"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_mhcflurry.py"
