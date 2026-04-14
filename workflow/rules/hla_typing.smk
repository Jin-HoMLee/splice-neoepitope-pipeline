# =============================================================================
# Rule module: HLA typing with arcasHLA
# =============================================================================
#
# Runs arcasHLA on per-sample BAMs (HISAT2 or STAR) to predict patient-specific
# HLA-A/B/C alleles, then aggregates per patient using a normal-first policy.
# The resulting alleles.tsv is consumed by run_mhcflurry to make predictions
# patient-specific instead of relying on a hardcoded default allele.
#
# Enabled via config.hla.enabled: true. If false, these rules are not defined
# and run_mhcflurry uses the fallback alleles from config.
#
# =============================================================================

import os


if (
    config.get("hla", {}).get("enabled", False)
    and config.get("data_source") == "fastq"
):

    # Derive HLA typing output directory from the raw_data path so that test
    # mode (results/test/...) and production (results/...) stay separated
    # without hardcoding either path.
    _HLA_TYPING_DIR = os.path.join(os.path.dirname(OUT["raw_data"]), "hla_typing")

    _ARCASHLA_REF_DONE = "resources/arcashla_reference.done"


    def _sample_is_paired(sample_id):
        """Return True if the sample has a non-empty fastq2 in samples_tsv."""
        for s in _read_samples_tsv(config["samples_tsv"]):
            if s["sample_id"] == sample_id:
                return bool((s.get("fastq2") or "").strip())
        return False


    rule download_arcashla_reference:
        """Fetch the arcasHLA IMGT/HLA reference database (one-time, ~200 MB).

        Mirrors the download_mhcflurry_models pattern: a sentinel file prevents
        re-downloading on subsequent runs.
        """
        output:
            sentinel=touch(_ARCASHLA_REF_DONE),
        log:
            os.path.join(OUT["logs"], "hla_typing", "download_reference.log"),
        conda:
            "../envs/arcashla.yaml"
        shell:
            """
            mkdir -p $(dirname {log})
            arcasHLA reference --update > {log} 2>&1
            """


    rule arcashla_genotype:
        """Run arcasHLA extract + genotype on a single sample BAM.

        Produces {sample}.genotype.json (allele calls) and
        {sample}.genotype.log (per-locus read counts used by the aggregator's
        fallback logic). If either arcasHLA stage fails, an empty json is
        written so the aggregator can apply fallback alleles downstream.
        """
        input:
            bam=_BAM_OUTPATH,
            bai=_BAI_OUTPATH,
            ref_done=_ARCASHLA_REF_DONE,
        output:
            json=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "{sample}", "{sample}.genotype.json"),
            genolog=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "{sample}", "{sample}.genotype.log"),
            done=touch(os.path.join(_HLA_TYPING_DIR, "{patient_id}", "{sample}", "{sample}.done")),
        log:
            os.path.join(OUT["logs"], "hla_typing", "{patient_id}_{sample}.log"),
        params:
            outdir=lambda w: os.path.join(_HLA_TYPING_DIR, w.patient_id, w.sample),
            single_flag=lambda w: "" if _sample_is_paired(w.sample) else "--single",
        threads: 4
        conda:
            "../envs/arcashla.yaml"
        shell:
            r"""
            set -uo pipefail
            mkdir -p {params.outdir}
            mkdir -p $(dirname {log})

            # arcasHLA requires UCSC-style `chr6` contig naming. Fail fast with
            # a clear message if the BAM uses Ensembl-style naming instead.
            if ! samtools view {input.bam} chr6 2>/dev/null | head -n1 | grep -q .; then
                echo "[arcashla_genotype] ERROR: no reads aligned to 'chr6' in {input.bam}" | tee -a {log}
                echo "[arcashla_genotype] arcasHLA requires UCSC-style 'chr6' contig naming." | tee -a {log}
                exit 1
            fi

            # Stage 1: extract chr6/HLA reads into FASTQ(s)
            arcasHLA extract {params.single_flag} \
                -t {threads} \
                -o {params.outdir} \
                {input.bam} \
                >> {log} 2>&1 || echo "[arcashla_genotype] WARN: extract failed, fallback will be used" | tee -a {log}

            # Locate extracted FASTQ(s) — single-end vs paired-end
            FQ_SE="{params.outdir}/{wildcards.sample}.extracted.fq.gz"
            FQ_PE1="{params.outdir}/{wildcards.sample}.extracted.1.fq.gz"
            FQ_PE2="{params.outdir}/{wildcards.sample}.extracted.2.fq.gz"
            FQ_INPUT=""
            if [ -f "$FQ_SE" ]; then
                FQ_INPUT="$FQ_SE"
            elif [ -f "$FQ_PE1" ] && [ -f "$FQ_PE2" ]; then
                FQ_INPUT="$FQ_PE1 $FQ_PE2"
            fi

            # Stage 2: genotype from the extracted FASTQ(s)
            if [ -n "$FQ_INPUT" ]; then
                arcasHLA genotype {params.single_flag} \
                    -g A,B,C \
                    -t {threads} \
                    -o {params.outdir} \
                    $FQ_INPUT \
                    >> {log} 2>&1 || echo "[arcashla_genotype] WARN: genotype failed, fallback will be used" | tee -a {log}
            else
                echo "[arcashla_genotype] WARN: no extracted reads found, fallback will be used" | tee -a {log}
            fi

            # Guarantee the declared outputs exist so Snakemake is happy.
            # The aggregator handles empty/missing files via the fallback path.
            if [ ! -f "{output.json}" ]; then
                echo '{{}}' > {output.json}
            fi
            if [ ! -f "{output.genolog}" ]; then
                touch {output.genolog}
            fi
            """


    def _aggregate_hla_inputs(wildcards):
        sample_ids = [s["sample_id"] for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id)]
        return {
            "jsons": expand(
                os.path.join(_HLA_TYPING_DIR, wildcards.patient_id, "{sample}", "{sample}.genotype.json"),
                sample=sample_ids,
            ),
            "genologs": expand(
                os.path.join(_HLA_TYPING_DIR, wildcards.patient_id, "{sample}", "{sample}.genotype.log"),
                sample=sample_ids,
            ),
        }


    rule aggregate_hla_alleles:
        """Aggregate per-sample arcasHLA calls into a per-patient alleles TSV.

        Applies the normal-first policy (prefer Solid Tissue Normal / Blood
        Derived Normal samples over tumor samples, since HLA is germline) and
        substitutes fallback alleles from config for missing or low-confidence
        loci. Writes a side QC file flagging tumor/normal discrepancies.
        """
        input:
            unpack(_aggregate_hla_inputs),
        output:
            alleles=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "alleles.tsv"),
            qc=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "hla_qc.tsv"),
        log:
            os.path.join(OUT["logs"], "hla_typing", "aggregate_{patient_id}.log"),
        params:
            samples_tsv=config["samples_tsv"],
            min_reads=config.get("hla", {}).get("min_reads_per_locus", 30),
            fallback=config.get("hla", {}).get("fallback_alleles", {}),
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/aggregate_hla_alleles.py"
