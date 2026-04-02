# =============================================================================
# Rule module: Step 5 — Epitope prediction with MHCflurry 2.x
# =============================================================================

rule run_mhcflurry:
    """Run MHCflurry 2.x on the 16-mer peptides FASTA.
    MHCflurry slides a 9-mer window across each peptide internally.
    Output: parsed TSV with columns: peptide, allele, IC50_nM, percentile_rank,
            binder_class (strong / weak / non)."""
    input:
        peptides_fasta=rules.translate_peptides.output.peptides_fasta,
    output:
        predictions_tsv=os.path.join(
            OUT["predictions"], "{cancer_type}", "predictions.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "predict", "{cancer_type}_predict.log"),
    params:
        hla_allele=config["mhcflurry"]["hla_allele"],
        peptide_length=config["mhcflurry"]["peptide_length"],
        ic50_strong=config["mhcflurry"]["ic50_strong"],
        ic50_weak=config["mhcflurry"]["ic50_weak"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_mhcflurry.py"
