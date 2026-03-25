# Resources Directory

This directory stores large reference files that are **not tracked in git** due
to their size.  Before running the pipeline you must download or generate the
following files and place them here (or update the paths in
`config/config.yaml`).

---

## Required Files

### 1. GRCh38 Reference Genome FASTA

**Filename**: `GRCh38.primary_assembly.genome.fa`
(and its index `GRCh38.primary_assembly.genome.fa.fai`)

**Download** from GENCODE:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
```

The index (`.fai`) is required by `bedtools getfasta` and will be created
automatically by the `samtools faidx` command above.

---

### 2. GENCODE GTF Annotation

**Filename**: `gencode.v47.annotation.gtf.gz`

**Download** from GENCODE:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```

This file is used by `build_reference_junctions.py` to extract all known
splice junctions, which are then used to filter out annotated junctions from
the TCGA data.

---

### 3. Reference Junction BED (auto-generated)

**Filename**: `reference_junctions.bed`

This file is **automatically generated** by the Snakemake pipeline from the
GENCODE GTF.  You do not need to create it manually.

---

## After Placing Files

Update the paths in `config/config.yaml` if you use different filenames or
locations:

```yaml
reference:
  genome_fasta: "resources/GRCh38.primary_assembly.genome.fa"
  gencode_gtf:  "resources/gencode.v47.annotation.gtf.gz"
  junction_bed: "resources/reference_junctions.bed"
```

---

## Disk Space Estimates

| File | Approximate Size |
|------|-----------------|
| GRCh38 genome FASTA (compressed) | ~840 MB |
| GRCh38 genome FASTA (uncompressed) | ~3.1 GB |
| GENCODE v47 GTF (compressed) | ~60 MB |
| TCGA-BRCA junction files | ~5 GB |
| TCGA-LUAD junction files | ~3 GB |
| TCGA-LAML junction files | ~1 GB |

Total: approximately **12–15 GB** for all input data.
