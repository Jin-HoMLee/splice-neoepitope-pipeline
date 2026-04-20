# Conclusions

*Living document — intended to be refined into a publication conclusions section.*

---

## Summary

We present a Snakemake-based pipeline for the discovery of splice-junction-derived
neoepitope candidates from RNA-seq data. The pipeline integrates splice junction
extraction (regtools), junction-level normal filtering, HLA typing (OptiType), peptide
translation with a junction-spanning filter, MHC class I binding prediction (MHCflurry
2.x), and structural validation (TCRdock), producing a ranked list of neoepitope
candidates with predicted TCR–peptide–MHC ternary complex structures.

---

## Key Findings

Applied to a matched gastric cancer tumor/normal RNA-seq pair (patient_001,
SRR9143066/SRR9143065), the pipeline:

1. **Identified 27,347 tumor-exclusive splice junctions** from 146,644 total extracted
   junctions (18.6%). Junction-level filtering against the matched normal eliminated
   79.5% annotated junctions and an additional 1.8% unannotated but normal-shared
   junctions.

2. **Translated 433,129 junction-spanning 9-mers** (424,133 unique sequences) from
   tumor-exclusive contigs across all three reading frames, applying a conservative
   complete-codon rule at the junction boundary.

3. **Predicted 14,990 strong binders** (IC50 ≤ 50 nM, 0.58% of 2,598,774 predictions)
   and 99,444 weak binders (IC50 ≤ 500 nM, 3.83%) across all six patient-specific HLA
   class I alleles. The HLA-B locus showed a normal/tumor typing discrepancy (normal:
   B\*15:01/B\*18:02; tumor: B\*15:63/B\*18:01), consistent with noise in tumor RNA
   typing rather than true LOH; normal calls were used per the normal-first policy.

4. **Identified EVAEYNASF / HLA-A\*26:01 (IC50 = 16.5 nM)** as the top neoepitope
   candidate. TCRdock structural prediction of the TCR–peptide–MHC ternary complex was
   successfully completed on GCP GPU infrastructure (NVIDIA P100, 16 GB VRAM).

---

## Significance

Aberrant splicing is a widespread but underexplored source of tumor neoantigens.
The junction-spanning filter and matched-normal subtraction strategy employed here
provide a computationally efficient means of prioritising peptides that are (i) derived
from tumor-specific splicing events and (ii) present at the MHC-bound peptide level.
The conservative design of the filter minimises false positives, which is critical
given the cost of downstream experimental validation.

The pipeline is fully automated, reproducible (Snakemake + conda), and cloud-ready
(GCP Compute Engine). Parallelisation of MHC binding predictions (ProcessPoolExecutor)
and GPU-accelerated structural prediction (TCRdock on P100) make it feasible to
process a full patient dataset within a single cloud session.

---

## Limitations

- **No proteome filter:** predicted peptides are not cross-referenced against the full
  human proteome. The junction-spanning filter removes most exonic false positives, but
  a cross-proteome BLAST check would provide a further layer of confidence.
- **HISAT2 sensitivity:** STAR has been shown to detect more novel splice junctions in
  benchmarks. Future runs will compare junction yield between the two aligners (issue #17).
- **HLA validation:** no ground-truth HLA alleles are available for patient_001.
  Patient_002 (osteosarcoma, known Red Cross serology) will provide the first direct
  validation of OptiType accuracy in this pipeline.
- **No matched normal for patient_002:** all unannotated junctions are treated as
  tumor-exclusive. Based on patient_001 statistics (~8.9% normal-shared rate), a
  corresponding fraction of candidates will be false positives from patient-specific
  but non-tumor splicing.

---

## Future Directions

- **Patient_002 (osteosarcoma IPISRC044):** apply the pipeline to a tumor-only case and
  validate OptiType HLA typing against known serology. Longitudinal samples (T0/T1/T2)
  enable future analysis of neoepitope evolution during treatment.
- **STAR aligner comparison:** benchmark STAR vs. HISAT2 junction sensitivity on the
  same dataset (issue #17).
- **Chimeric codon recovery:** optionally retain junction-boundary 9-mers where a
  chimeric codon introduces an amino acid change at an anchor position (P2 or P9),
  representing a small high-confidence supplementary candidate set.
- **Proteome filter:** add a BLAST or exact-match step against the reference proteome
  to catch any peptides matching a normal protein at a different genomic locus.
- **Pre-built HISAT2 index:** cache the GRCh38 index to reduce VM startup time from
  ~60–90 min to near-instant (issue #16).
