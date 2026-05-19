# Glossary

Project-relevant abbreviations and acronyms. The pipeline mixes biology, ML, bioinformatics, and cloud infrastructure — entries here aim to make onboarding (or returning to the project after time away) easier than re-asking each time.

**Format:** `**ABBREV** — Expansion. One-line context. *Domain: bio | ml | cloud | pipeline | stats.*`

**Scope:** project-relevant terms only. Skip generic web/programming acronyms (HTTP, REST, etc.).

---

## A

**AC** — Acceptance Criteria. The checklist in an issue body defining what "done" means for that issue. Per the closure ritual, all AC boxes must be ticked (`- [x]`) or comment-deferred before closing. *Domain: project.*

**AF-Multimer** — AlphaFold-Multimer. DeepMind extension of AlphaFold (Evans et al. 2021) trained on protein complexes; predicts both monomeric structure and inter-chain interactions. Used by TCRdock as the structural backbone. *Domain: bio.*

**AF_confidence** — AlphaFold confidence score. Combines pTM (predicted Template Modeling) and ipTM (interface predicted TM); typically `0.8·ipTM + 0.2·pTM`. Standard global-quality metric for AF-Multimer model selection. *Domain: bio.*

**AFDB** — AlphaFold Database (alphafold.ebi.ac.uk). DeepMind/EBI's public repo of ~200M predicted protein structures. *Domain: bio.*

**AP** — Access Point. In HTCondor topology, the machine where jobs are submitted, files staged, and conda envs built. Distinct from the EP (execution point) when there is no shared FS; the AP→EP env-deployment gap is the same problem cloud executors (Google Batch) face. *Domain: cloud.*

## C

**CAPRI** — Critical Assessment of PRedicted Interactions. Community blind-prediction challenge for protein-protein / peptide-MHC docking (since 2001); CAPRI quality bands (high / medium / acceptable / incorrect) are the standard yardstick for docking models — DockQ is the continuous reformulation of these bands. *Domain: bio.*

**CAR-T** — Chimeric Antigen Receptor T-cell therapy. Adoptive cell therapy where patient T cells are transduced with a synthetic chimeric receptor (scFv + CD3ζ + costim domains); recognizes surface antigen directly, MHC-independent. FDA-approved CD19/BCMA-targeting products for B-cell malignancies and myeloma. Distinct from TCR-T (native MHC-restricted TCR). *Domain: bio.*

**CDR3** — Complementarity-Determining Region 3. Hypervariable loop on each TCR chain (α and β) that makes direct contact with the peptide in pMHC; the primary specificity-encoding region of the TCR. Average `pLDDT` across CDR3 residues used as a TCR-pMHC docking quality signal (Lu et al. 2025, [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316)). *Domain: bio.*

**CHTC** — Center for High Throughput Computing (UW-Madison). Develops and maintains HTCondor; took over maintenance of Snakemake's HTCondor executor plugin at the 2026 Snakemake Hackathon. *Domain: cloud.*

**CoREST** — Corepressor of REST (RE1-Silencing Transcription factor). Chromatin co-repressor complex (LSD1–HDAC1–CoREST core scaffold); originally characterised as a neuronal-gene silencer, recently shown to also regulate RNA-processing genes. Pharmacological inhibition by **corin** in melanoma cell lines disrupts spliceosome activity → induces immunogenic splice-derived neoantigens, reactivates anti-PD-1 in cold tumors ([Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). *Domain: bio.*

**CR** — Complete Remission (a.k.a. Complete Response). Oncologic response criterion: full disappearance of detectable disease following therapy; distinguished from PR (partial response), SD (stable disease), PD (progressive disease) in RECIST 1.1 solid-tumor criteria. The weakest TNBC vaccine responder in Sahin et al. 2026 achieved CR on subsequent anti-PD-1 rescue. *Domain: bio.*

**CVE** — Common Vulnerabilities and Exposures. Standardised identifier system for publicly disclosed security flaws (CVE-YYYY-NNNN format); maintained by MITRE, used industry-wide for tracking known vulns. *Domain: security.*

## D

**DAG** — Directed Acyclic Graph. Snakemake's compiled job-dependency graph; computed from rule wildcards + I/O files, traversed for scheduling. Rendered via `scripts/visualize_dag.sh` after a dry-run. *Domain: pipeline.*

**DDA** — Data-Dependent Acquisition. Mass-spec acquisition mode where the instrument picks the top-N most intense MS1 precursors per cycle and fragments only those. Stochastic and intensity-biased — low-abundance peptides (including most neoantigens) are routinely missed across replicate runs. Historical default; superseded by DIA-MS for immunopeptidomics. *Domain: bio.*

**DIA-MS** — Data-Independent Acquisition Mass Spectrometry. Acquisition mode where the instrument fragments every precursor in a stepped m/z window, regardless of intensity — no peptide is "missed", but resulting MS2 spectra are multiplexed and require a reference spectral library for deconvolution. Public libraries cover canonical proteomes; patient-specific neoantigens need custom libraries (e.g. Pepyrus, [Manakongtreecheep et al. 2026](https://www.nature.com/articles/s41587-026-03003-9)). *Domain: bio.*

**DLT** — Dose-Limiting Toxicity. Pre-defined adverse-event severity threshold (typically CTCAE grade ≥3) in Phase 1 dose-escalation trials; reaching DLT caps further escalation. "No DLTs" is the standard Phase 1 safety endpoint reported in neoantigen-vaccine trials (e.g. GNOS-PV01, [Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**DockQ** — Docking Quality. Continuous 0–1 score combining Fnat (fraction of native contacts), iRMSD (interface RMSD), and LRMSD (ligand RMSD); maps onto CAPRI bands; used to evaluate predicted protein-protein / TCR-pMHC complex structures against experimental ground truth. *Domain: bio.*

## E

**ENCODE** — Encyclopedia of DNA Elements. NIH consortium (since 2003) cataloging functional genome features via standardized high-throughput assays (RNA-seq, ChIP-seq, ATAC-seq, Hi-C, CAGE, bisulfite-seq, etc.). "ENCODE-style assays" = this standard set of functional genomics readouts, regardless of which project produced the data. *Domain: bio.*

**EP** — Execution Point. In HTCondor topology, the machine where a job actually runs. In no-shared-FS topologies the conda env built on the AP must be transported to (or rebuilt on) the EP — the same constraint cloud workers face. *Domain: cloud.*

**ESM-IF1** — Evolutionary Scale Modeling Inverse Folding (v1). Meta AI model (Hsu et al. 2022) trained on 12M protein structures; given a backbone, infers compatible sequences. Used as a structure-aware embedding source by structure-informed TCR-pMHC scorers (e.g. NetTCR-struc). *Domain: ml.*

## G

**GATK** — Genome Analysis Toolkit (Broad Institute). Canonical somatic/germline variant calling suite; reference for the Panel of Normals (PoN) concept reused by neoepitope filtering. *Domain: bio.*

**GBM** — Glioblastoma (Multiforme). WHO grade 4 astrocytoma; most aggressive primary adult brain tumor (~15-mo median OS post standard-of-care surgery + radiation + temozolomide). MGMT-unmethylated subset has the worst prognosis (resistant to temozolomide); subject of the GNOS-PV01 personalized DNA vaccine Phase 1 ([Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**GCP** — Google Cloud Platform. The cloud provider hosting this pipeline's VMs and storage bucket (zone `europe-west1-b`). *Domain: cloud.*

**GCS** — Google Cloud Storage. GCP's object store; pipeline results, logs, and reference data live in `gs://splice-neoepitope-project`. *Domain: cloud.*

**GKE** — Google Kubernetes Engine. Managed Kubernetes on GCP. *Domain: cloud.*

**GNOS-PV01** — Personalized DNA neoantigen vaccine (Geneos Therapeutics + WashU Medicine). Encodes up to 40 patient-specific neoantigens on a DNA backbone, delivered intramuscularly with electroporation; Phase 1 in MGMT-unmethylated GBM ([Garfinkle et al., Nat Cancer 2026](https://www.nature.com/articles/s43018-026-01163-w)). Adds **DNA** to the platform diversity census alongside peptide / mRNA / DC. *Domain: bio.*

**GPS** — Genotype Presentation Score. Project-specific scoring formula combining per-allele MHCflurry presentation percentiles into a single rank for the patient's genotype; primary ranking key for top neoepitope candidates. *Domain: pipeline.*

**GTEx** — Genotype-Tissue Expression project. NIH consortium with healthy-tissue RNA-seq across ~50 tissues from ~1,000 donors; canonical reference for tissue-matched normal expression panels. *Domain: bio.*

**GVP-GNN** — Geometric Vector Perceptron Graph Neural Network. Architecture (Jing et al. 2021) for learning on protein 3D structures; respects rotational/translational equivariance natively. Used in NetTCR-struc and similar structure-aware TCR-pMHC scorers. *Domain: ml.*

## H

**HLA** — Human Leukocyte Antigen. Human MHC class I/II proteins (HLA-A/B/C class I; HLA-DR/DP/DQ class II); patient typing via OptiType drives MHCflurry's per-allele predictions. *Domain: bio.*

**HLA-LOH** — HLA Loss of Heterozygosity. Tumor immune-escape mechanism: somatic deletion or copy-neutral loss of one HLA allele, narrowing the peptide-presentation repertoire and shielding the tumor from neoantigen-specific T cells. Detected from WES (e.g. LOHHLA); candidates predicted on a lost allele are no longer presented. *Domain: bio.*

**HPC** — High-Performance Computing. Cluster compute systems with job schedulers (SLURM, LSF, HTCondor); distinct from cloud "on-demand" compute models (GCP, AWS Batch). Pipeline currently runs on cloud, not HPC. *Domain: cloud.*

**HTCondor** — High-Throughput Condor. Distributed batch scheduler from UW-Madison's CHTC; dominant in physics/astronomy (e.g. CERN ATLAS). Supports no-shared-FS topologies natively — relevant precedent for cloud-executor design. *Domain: cloud.*

## I

**IC50** — half-maximal Inhibitory Concentration. Concentration (in nM) at which a peptide displaces 50% of a reference ligand from MHC; classic affinity metric, retained as informational column alongside the percentile-based presentation score. *Domain: bio.*

**ICI** — Immune Checkpoint Inhibitor. Class of monoclonal antibodies blocking inhibitory T-cell receptors (anti-CTLA-4, anti-PD-1, anti-PD-L1) to release brakes on anti-tumor immunity; standard-of-care across many solid tumors. "Immune cold" tumors (low T-cell infiltrate) often fail to respond — motivating combination strategies including neoantigen vaccines and pharmacologic mis-splicing inducers (corin; [Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). *Domain: bio.*

**IV** — Intravenous. Drug delivery route via direct venous injection; chosen for personalized mRNA cancer vaccines to deliver LNP-mRNA payloads systemically (Sahin et al. 2026; Rojas et al. 2023). *Domain: bio.*

## L

**LNP-mRNA** — Lipid Nanoparticle-formulated mRNA. Delivery format: ionizable-lipid encapsulation protects mRNA from RNases and mediates dendritic-cell uptake. Foundation of COVID-19 mRNA vaccines (BioNTech/Moderna) and current personalized cancer-vaccine trials (autogene cevumeran; Sahin et al. 2026). *Domain: bio.*

**LSD1–HDAC1–CoREST** — Core chromatin co-repressor complex: lysine-specific demethylase 1 (LSD1, KDM1A) + histone deacetylase 1 (HDAC1) + CoREST scaffold. Coordinates H3K4 demethylation with H3 deacetylation for transcriptional repression; the molecular target of the small-molecule inhibitor **corin** ([Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). Disruption reroutes the complex away from RNA-processing genes → mis-splicing → immunogenic neopeptides. *Domain: bio.*

**LSF** — Load Sharing Facility (IBM Spectrum). Commercial HPC batch scheduler; used at DKFZ Heidelberg among others. Snakemake has a dedicated LSF executor plugin. *Domain: cloud.*

## M

**MCP** — Model Context Protocol. Anthropic's open spec (late 2024) for AI agents to discover and invoke external tools / data sources via per-purpose MCP servers; the standardized way to expose resources to LLMs without bespoke per-tool integration. *"USB-C for AI agents."* Used by Claude Code, the GitHub MCP server, filesystem servers, etc. *Domain: ml.*

**MGMT** — O⁶-methylguanine-DNA methyltransferase. DNA repair enzyme; promoter methylation silences the gene and sensitizes glioma to alkylating chemotherapy (temozolomide). **MGMT-unmethylated** GBM has worse prognosis (temozolomide-resistant) — the poor-response subset targeted by the GNOS-PV01 vaccine trial ([Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**MHC** — Major Histocompatibility Complex. Cell-surface proteins that present peptides to T cells; class I (all nucleated cells, presents endogenous peptides) drives this pipeline. Human MHC = HLA. *Domain: bio.*

## N

**NJ** — Neojunction. Tumor-specific or recurrent splice junction absent from canonical normal-tissue annotation; the upstream substrate from which splice neoepitopes are translated. Notation introduced by Nejo et al. (*Nat. Med.* 2023). Kwok et al. additionally use **NEJ** (neoepitope-encoding junction) for the validated subset yielding a presented peptide (e.g. NEJ<sub>GNAS</sub>, NEJ<sub>RPL22</sub>). *Domain: bio.*

## O

**OKR** — Objectives and Key Results. Goal-setting framework popularised by Intel/Google: one qualitative *Objective* + 3–5 quantitative *Key Results*, usually quarterly. Common in software shops; we use S1–S7 lifecycle stages + iteration capacity instead. *Domain: project.*

**OOD** — Out-of-distribution. Inputs whose features fall in low-density regions of the model's training distribution; predictions on OOD inputs are typically over-confident yet unreliable since the model never learned to constrain them. Detected operationally via density estimators, embedding-space distance to training samples, or ensemble disagreement. Relevant for splice-junction-spanning peptides — systematically under-represented in canonical-proteome MHC training data. *Domain: ml.*

## P

**PAE** — Predicted Aligned Error. AlphaFold/TCRdock per-residue-pair error estimate (in Å); summarised globally to gauge inter-domain confidence (e.g. relative MHC↔TCR docking confidence). *Domain: bio.*

**PDB** — Protein Data Bank ([rcsb.org](https://www.rcsb.org)). Public archive of experimentally determined 3D structures (~220k entries); training source for AlphaFold / Boltz / HERMES; reference ground truth for DockQ scoring. *Domain: bio.*

**pLDDT** — predicted Local Distance Difference Test. AlphaFold/TCRdock per-residue confidence score (0–100); per-residue version of the global `AF_confidence`. >90 = high, >70 = confident, <50 = low confidence. CDR3-region average used as a TCR-pMHC docking quality flag (Lu et al. 2025, [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316)). *Domain: bio.*

**pMHC** — peptide-MHC. The complex of a peptide bound in the MHC groove; the molecular target a TCR recognises. "TCR-pMHC interaction" = the central event in adaptive immune recognition. *Domain: bio.*

**PoN** — Panel of Normals. Aggregated reference set of normal samples (originally GATK's somatic variant calling concept); used to filter recurrent germline variation / artefacts that single matched normals miss. *Domain: bio.*

**PSI** — Percent Spliced In. Fraction of transcripts that include a given alternative exon/junction (vs skip); the standard quantitative readout of splicing in RNA-seq analyses. *Domain: bio.*

**PSR** — Positive Sample Rate. Percentage of samples in a cohort expressing a splice junction at a read frequency ≥1% relative to the canonical junction (Kwok et al., *Nature* 2025); thresholded as `PSR_TCGA ≥ 10%` (tumor-recurrence) and `PSR_GTEx < 1%` (normal-exclusion). Distinct from **PSI** (within-sample inclusion fraction): PSR is a binary cross-sample detection rate at a fixed relative-frequency threshold; PSI is a continuous within-sample metric. *Domain: bio.*

## R

**RCE** — Remote Code Execution. Vulnerability class where an attacker runs arbitrary code on a remote system, typically via input-validation or memory-safety bugs (e.g. CVE-2026-3854's crafted git push trigger). *Domain: security.*

**RMSD** — Root Mean Square Deviation. Average distance (in Å) between corresponding atoms after optimal superposition; the standard structural-similarity metric. Variants used in DockQ: backbone-only, all-atom, interface-only (iRMSD), ligand-only (LRMSD). *Domain: bio.*

## S

**SLURM** — Simple Linux Utility for Resource Management. The dominant HPC job scheduler in academia. *Domain: cloud.*

## T

**TCR** — T-Cell Receptor. Heterodimeric (α/β) receptor on T cells that recognises peptide-MHC complexes; modelled here via TCRdock for the top neoepitope candidate. *Domain: bio.*

**TCR-T** — TCR-engineered T-cell therapy. Adoptive cell therapy where autologous patient T cells are transduced with a tumor-reactive TCR (typically discovered in HLA-matched healthy donors to bypass tumor-induced tolerance). MHC-restricted; distinct from CAR-T (synthetic chimeric receptor) and TIL therapy (non-engineered expanded patient T cells). *Domain: bio.*

**TF** — Transcription Factor. DNA-binding regulatory protein controlling gene expression; binding sites are one of the modalities ENCODE-style assays (ChIP-seq) and predictors like AlphaGenome map. *Domain: bio.*

**TIL** — Tumor-Infiltrating Lymphocyte. T cell isolated from a tumor biopsy; expanded ex vivo and re-infused as TIL therapy (FDA-approved 2024 for advanced melanoma as lifileucel/Amtagvi). Distinct from TCR-T (engineered TCR) and CAR-T (synthetic chimeric receptor). *Domain: bio.*

**TNBC** — Triple-Negative Breast Cancer. Breast cancer subtype lacking ER, PR, and HER2 expression; lacks targeted hormonal / HER2 therapy options. ~10–15% of breast cancers; relatively higher tumor mutation burden makes it a recurrent personalized-vaccine trial target (Sahin et al., *Nature* 2026). *Domain: bio.*

**TPM** — Transcripts Per Million. Normalised RNA expression unit; reads per kb of transcript scaled so per-sample values sum to 10⁶, making cross-sample comparisons direct (unlike RPKM/FPKM). Conventional thresholds: ~1 TPM = expressed, <0.5 TPM = below detection floor (functionally silent in that tissue), >10 TPM = moderately–highly expressed. Used as a tier-1 normal-tissue QC for public-neoantigen prioritisation (e.g. Zhang et al. 2026: <0.5 TPM in normal tissues). *Domain: bio.*

## W

**WES** — Whole Exome Sequencing. Sequencing restricted to protein-coding regions (~1–2% of the genome); cheaper than WGS, sufficient for variant calls in coding regions. *Domain: bio.*

**WGS** — Whole Genome Sequencing. Sequencing of the entire genome (~3 Gbp human); alternative input for upstream filtering when matched-normal RNA-seq is missing. *Domain: bio.*
