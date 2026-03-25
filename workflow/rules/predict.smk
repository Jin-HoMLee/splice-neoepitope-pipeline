# =============================================================================
# Rule module: Step 5 — Epitope prediction with NetMHCPan 4.1
# =============================================================================

rule run_netmhcpan:
    """Run NetMHCPan 4.1 on the 16-mer peptides FASTA.
    NetMHCPan internally slides a 9-mer window across each peptide.
    Output: parsed TSV with columns: peptide, allele, IC50_nM, rank,
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
        executable=config["netmhcpan"]["executable"],
        hla_allele=config["netmhcpan"]["hla_allele"],
        peptide_length=config["netmhcpan"]["peptide_length"],
        ic50_strong=config["netmhcpan"]["ic50_strong"],
        ic50_weak=config["netmhcpan"]["ic50_weak"],
        output_format=config["netmhcpan"]["output_format"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/run_netmhcpan.py"
