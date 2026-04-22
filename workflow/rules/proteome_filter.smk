# =============================================================================
# Rule module: Step 4b — proteome filter peptides against human proteome
# =============================================================================
#
# Removes any translated peptide that already exists in the canonical human
# proteome (UniProt Swiss-Prot). Such peptides are self-peptides the immune
# system is tolerized to and cannot be neoantigens.
#
# Runs between translate_peptides (Step 4) and run_mhcflurry (Step 5):
#
#   translate_peptides → proteome_filter_peptides → run_mhcflurry
#
#
# =============================================================================

_PROTEOME_FILTER_ENABLED = config.get("proteome_filter", {}).get("enabled", True)
_PROTEOME_FASTA = config.get("proteome_filter", {}).get(
    "proteome_fasta", "resources/human_proteome.fasta"
)

if _PROTEOME_FILTER_ENABLED:

    rule download_human_proteome:
        """Download UniProt Swiss-Prot human proteome FASTA.

        Source: UniProt REST API — reviewed human entries only (~20k canonical
        proteins). Downloaded once to resources/ and reused across all patients.
        """
        output:
            fasta=_PROTEOME_FASTA,
        log:
            os.path.join(_LOGS, "proteome_filter", "download_proteome.log"),
        params:
            url="https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true",
        conda:
            "../envs/python.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.fasta})

            curl --fail -L --no-progress-meter \
                "{params.url}" 2>> {log} \
                | gunzip -c > {output.fasta}

            echo "Downloaded $(grep -c '^>' {output.fasta}) Swiss-Prot entries" >> {log}
            """


    rule proteome_filter_peptides:
        """Filter translated peptides against the human proteome via k-mer index.

        Builds a {k_mer: [accessions]} index from Swiss-Prot in one FASTA pass,
        then checks each query peptide in O(1).  Excludes any peptide present
        in the canonical proteome; retains all source accessions per match for
        downstream inspection.
        """
        input:
            peptides_tsv=rules.translate_peptides.output.peptides_tsv,
            proteome_fasta=_PROTEOME_FASTA,
        output:
            novel_tsv=os.path.join(
                _RES, "{patient_id}", "peptides", "peptides_novel.tsv"
            ),
            excluded_tsv=os.path.join(
                _RES, "{patient_id}", "peptides", "peptides_excluded.tsv"
            ),
        log:
            os.path.join(_LOGS, "{patient_id}", "proteome_filter", "proteome_filter.log"),
        params:
            peptide_lengths=config["translation"]["peptide_lengths"],
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/proteome_filter.py"
