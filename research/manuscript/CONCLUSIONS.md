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

1. **Identified 27,348 tumor-exclusive splice junctions** from 146,647 total extracted
   junctions (18.6%). Junction-level filtering against the matched normal eliminated
   79.5% annotated junctions and an additional 1.8% unannotated but normal-shared
   junctions.

2. **Translated 1,286,492 junction-spanning peptides** (1,260,074 unique sequences,
   97.9%) from tumor-exclusive contigs across 8/9/10-mer lengths and all three reading
   frames, applying a conservative complete-codon rule at the junction boundary
   (`peptide_lengths: [8, 9, 10]`, PR #99).

3. **Predicted 44,916 strong presenters** (presentation percentile ≤ 0.5%, 3.5% of
   1,286,492 predictions) and 125,775 weak presenters (≤ 2.0%, 9.8%) across all six
   patient-specific HLA class I alleles, using MHCflurry 2.x `Class1PresentationPredictor`
   (genotype-level call returning one best-allele prediction per peptide). The HLA-B
   locus showed a normal/tumor typing discrepancy (normal: B\*15:01/B\*18:02; tumor:
   B\*18:01/B\*15:63), consistent with noise in tumor RNA typing rather than true LOH;
   the run used the tumor calls per `report.tsv`.

4. **Identified SQIPRTHSY / HLA-C\*07:01 (IC50 = 33.6 nM, presentation percentile
   = 0.0052%, GPS = 0.9999)** as the top neoepitope candidate, ranked by Genotype
   Presentation Score and presented as strong by 5 of 6 patient alleles. TCRdock
   structural prediction of the TCR–peptide–MHC ternary complex was successfully
   completed on GCP GPU infrastructure (NVIDIA P100, 16 GB VRAM).

Applied to an osteosarcoma tumor-only sample (patient_002, BG003082 T0) without a
matched RNA-seq normal, the pipeline:

5. **Identified 55,912 tumor-exclusive splice junctions** from 347,046 total extracted
   junctions (16.1%). Without an RNA-seq normal, junction-level subtraction was limited
   to a WES DNA proxy (3 junctions removed), so a fraction of candidates may reflect
   patient-specific non-tumor splicing.

6. **Translated 781,424 junction-spanning 9-mers** (775,440 unique sequences).

7. **Predicted 12,430 strong binders** (IC50 ≤ 50 nM, 0.32% of 3,907,120 predictions)
   across all six patient-specific HLA alleles.

8. **Validated OptiType HLA typing against Red Cross serology** (A\*01:01/A\*01:11N,
   B\*08:01/B\*27:05, C\*01:02/C\*07:01): all six alleles were an exact match, providing
   the first direct confirmation of OptiType accuracy in this pipeline. A\*01:11N is a
   null allele; OptiType correctly identified A\*01:01 as the expressed allele.

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
  Patient_002 (osteosarcoma) confirmed OptiType accuracy against Red Cross serology
  (A\*01:01/A\*01:11N, B\*08:01/B\*27:05, C\*01:02/C\*07:01): all six alleles were
  an exact match. A\*01:11N is a null allele; OptiType correctly identified A\*01:01
  as the expressed allele and called it homozygous.
- **No matched normal for patient_002:** all unannotated junctions are treated as
  tumor-exclusive. Based on patient_001 statistics (~8.9% normal-shared rate), a
  corresponding fraction of candidates will be false positives from patient-specific
  but non-tumor splicing.

---

## Future Directions

- **Patient_002 longitudinal (T1/T2):** T0 run is complete (55,912 tumor-exclusive
  junctions, 12,430 strong MHC binders, HLA typing validated against serology).
  Longitudinal samples T1 and T2 are available for neoepitope evolution analysis
  during treatment.
- **STAR aligner comparison:** benchmark STAR vs. HISAT2 junction sensitivity on the
  same dataset (issue #17).
- **Chimeric codon recovery:** optionally retain junction-boundary 9-mers where a
  chimeric codon introduces an amino acid change at an anchor position (P2 or P9),
  representing a small high-confidence supplementary candidate set.
- **Proteome filter:** add a BLAST or exact-match step against the reference proteome
  to catch any peptides matching a normal protein at a different genomic locus.
- **Pre-built HISAT2 index:** cache the GRCh38 index to reduce VM startup time from
  ~60–90 min to near-instant (issue #16).
