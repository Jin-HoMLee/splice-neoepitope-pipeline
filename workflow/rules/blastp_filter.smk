# =============================================================================
# Rule module: Step 4b — blastp filter peptides against human proteome
# =============================================================================
#
# Removes any translated 9-mer that already exists in the canonical human
# proteome (UniProt Swiss-Prot). Such peptides are self-peptides the immune
# system is tolerized to and cannot be neoantigens.
#
# Runs between translate_peptides (Step 4) and run_mhcflurry (Step 5):
#
#   translate_peptides → blastp_filter_peptides → run_mhcflurry
#
# Disable for speed on test runs: set blastp_filter.enabled: false in config.
#
# =============================================================================

_BLASTP_ENABLED   = config.get("blastp_filter", {}).get("enabled", True)
_PROTEOME_FASTA   = config.get("blastp_filter", {}).get(
    "proteome_fasta", "resources/human_proteome.fasta"
)
_BLASTP_DB_DIR    = config.get("blastp_filter", {}).get(
    "db_dir", "resources/blastp_db"
)
_BLASTP_DB_PREFIX = os.path.join(_BLASTP_DB_DIR, "human_proteome")

if _BLASTP_ENABLED:

    rule download_human_proteome:
        """Download UniProt Swiss-Prot human proteome and build a blastp database.

        Source: UniProt REST API — reviewed human entries only (~20k canonical
        proteins). Downloaded once to resources/ and reused across all patients.
        """
        output:
            fasta=_PROTEOME_FASTA,
            db_done=touch(os.path.join(_BLASTP_DB_DIR, "makeblastdb.done")),
        log:
            os.path.join(_LOGS, "blastp_filter", "download_proteome.log"),
        params:
            db_prefix=_BLASTP_DB_PREFIX,
            db_dir=_BLASTP_DB_DIR,
            url="https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true",
        conda:
            "../envs/blast.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.db_dir}

            curl --fail -L --no-progress-meter \
                "{params.url}" \
                | gunzip -c > {output.fasta} \
                2>> {log}

            echo "Downloaded $(grep -c '^>' {output.fasta}) Swiss-Prot entries" >> {log}

            makeblastdb \
                -in {output.fasta} \
                -dbtype prot \
                -out {params.db_prefix} \
                >> {log} 2>&1
            """


    rule blastp_filter_peptides:
        """Filter translated peptides against the human proteome via blastp.

        Excludes any peptide with a perfect full-length match (100% identity,
        full query coverage) to a UniProt Swiss-Prot human protein. Retains
        an audit TSV of excluded peptides with matching accessions.
        """
        input:
            peptides_tsv=rules.translate_peptides.output.peptides_tsv,
            db_done=os.path.join(_BLASTP_DB_DIR, "makeblastdb.done"),
        output:
            novel_tsv=os.path.join(
                _RES, "{patient_id}", "peptides", "peptides_novel.tsv"
            ),
            excluded_tsv=os.path.join(
                _RES, "{patient_id}", "peptides", "peptides_excluded.tsv"
            ),
            blastp_hits=os.path.join(
                _RES, "{patient_id}", "peptides", "blastp_hits.tsv"
            ),
        log:
            os.path.join(_LOGS, "{patient_id}", "blastp_filter", "blastp_filter.log"),
        params:
            db_prefix=_BLASTP_DB_PREFIX,
        threads: config.get("blastp_filter", {}).get("threads", 4)
        conda:
            "../envs/blast.yaml"
        script:
            "../scripts/blastp_filter.py"
