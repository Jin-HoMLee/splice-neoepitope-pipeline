# Conclusions

*Living document — intended to be refined into a publication conclusions section.*

---

## Summary

We present a Snakemake-based pipeline for the discovery of splice-junction-derived
neoepitope candidates from RNA-seq data. The pipeline integrates splice junction
extraction (STAR or regtools), junction-level normal filtering, HLA typing (OptiType), peptide
translation with a junction-spanning filter, MHC class I binding prediction (MHCflurry
2.x), and structural validation (TCRdock), producing a ranked list of neoepitope
candidates with predicted TCR–peptide–MHC ternary complex structures.

---

## Key Findings

Applied to a matched gastric cancer tumor/normal RNA-seq pair (patient_001,
SRR9143066/SRR9143065), the pipeline:

1. **Identified 8 tumor-exclusive splice junctions** from 131,321 extracted junctions.
   After a read-support filter and discarding GENCODE-annotated junctions, 141 unannotated
   junctions remained; matched-normal subtraction (94) and the GTEx pan-tissue population
   filter (39) left 8 tumor-exclusive candidates. (This corrects a pre-#370 draft that
   reported 27,348 - an artifact of the anchor-outer coordinate bug fixed in
   [#370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370).)

2. **Translated 395 novel junction-spanning peptides** (all unique) from tumor-exclusive
   contigs across 8/9/10-mer lengths and all three reading frames, applying a conservative
   complete-codon rule at the junction boundary (`peptide_lengths: [8, 9, 10]`, PR #99).

3. **Predicted 17 strong and 22 weak presenters** (of 395 predictions) across all six
   patient-specific HLA class I alleles, using MHCflurry 2.x `Class1PresentationPredictor`
   (genotype-level call returning one best-allele prediction per peptide); HLA-C alleles
   account for 12 of 17 strong presenters. The HLA-B locus showed a normal/tumor typing
   discrepancy (normal: B\*15:01/B\*18:02; tumor: B\*18:01/B\*15:63), consistent with noise
   in tumor RNA typing rather than true LOH; the run used the tumor calls.

4. **Identified SQVTRGLAM / HLA-B\*15:63 (IC50 = 64.7 nM, presentation percentile
   = 0.257%, GPS = 0.942)** as the top neoepitope candidate, ranked by Genotype
   Presentation Score and presented as strong by 3 of 6 patient alleles. TCRdock
   structural prediction of the TCR-peptide-MHC ternary complex was completed (model
   pLDDT 92.4), using a fallback TCR (no cognate TCR available) - so the structure is a
   geometry check rather than evidence of patient-specific recognition.

Applied to an osteosarcoma tumor sample (patient_002, BG003082 T0), which lacks a
tissue-matched normal, the pipeline:

5. **Demonstrated a matched-normal boundary condition.** The STAR run yielded 42 nominal
   tumor-exclusive junctions (of 442 unannotated: 154 normal-shared, 246 removed by the GTEx
   pan-tissue filter), but the only available normal is a CD3+ T-cell PBMC transcriptome, not
   tissue-matched to a bone tumor. Because a T-cell normal cannot subtract bone/mesenchymal
   splicing, this candidate set is not tissue-specificity-controlled and is not comparable to
   patient_001's; we report it as a limitation rather than a validated result. Downstream
   figures (peptides, presenters, GPS, and the top candidate VPQVRVTVL / HLA-B\*08:01) are
   regenerated in the patient_002 notebook for the record but not tabulated as results. The
   code/config inconsistency that let the blood normal silently drive subtraction is tracked in
   [#940](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/940).

6. **Validated OptiType HLA typing against Red Cross serology** (A\*01:01/A\*01:11N,
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
(GCP Compute Engine). MHC presentation prediction runs as a single genotype-level
`Class1PresentationPredictor.predict()` call with all peptides batched at once
(GPU parallelism applies within the call); GPU-accelerated structural prediction
(TCRdock on P100) handles the downstream ternary-complex modelling. Together these
make it feasible to process a full patient dataset within a single cloud session.

---

## Limitations

- **No proteome filter:** predicted peptides are not cross-referenced against the full
  human proteome. The junction-spanning filter removes most exonic false positives, but
  a cross-proteome BLAST check would provide a further layer of confidence.
- **HLA validation:** no ground-truth HLA alleles are available for patient_001.
  Patient_002 (osteosarcoma) confirmed OptiType accuracy against Red Cross serology
  (A\*01:01/A\*01:11N, B\*08:01/B\*27:05, C\*01:02/C\*07:01): all six alleles were
  an exact match. A\*01:11N is a null allele; OptiType correctly identified A\*01:01
  as the expressed allele and called it homozygous.
- **No tissue-matched normal for patient_002:** its only normal is a CD3+ T-cell PBMC,
  which cannot subtract bone/mesenchymal splicing, so its tumor-exclusive junctions are not
  tissue-specificity-controlled (see Key Findings 5). A tissue-appropriate normal, or an
  explicit population-only mode, is the prerequisite for a defensible candidate set
  ([#940](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/940)).

---

## Future Directions

- **Patient_002 longitudinal (T1/T2):** T0 STAR run is complete (HLA typing validated
  against serology; junction-level candidates pending a tissue-appropriate normal, see
  Key Findings 5). Longitudinal samples T1 and T2 are available for neoepitope evolution
  analysis during treatment.
- **Chimeric codon recovery:** optionally retain junction-boundary 9-mers where a
  chimeric codon introduces an amino acid change at an anchor position (P2 or P9),
  representing a small high-confidence supplementary candidate set.
- **Proteome filter:** add a BLAST or exact-match step against the reference proteome
  to catch any peptides matching a normal protein at a different genomic locus.
