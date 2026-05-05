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

## C

**CAPRI** — Critical Assessment of PRedicted Interactions. Community blind-prediction challenge for protein-protein / peptide-MHC docking (since 2001); CAPRI quality bands (high / medium / acceptable / incorrect) are the standard yardstick for docking models — DockQ is the continuous reformulation of these bands. *Domain: bio.*

**CVE** — Common Vulnerabilities and Exposures. Standardised identifier system for publicly disclosed security flaws (CVE-YYYY-NNNN format); maintained by MITRE, used industry-wide for tracking known vulns. *Domain: security.*

## D

**DockQ** — Docking Quality. Continuous 0–1 score combining Fnat (fraction of native contacts), iRMSD (interface RMSD), and LRMSD (ligand RMSD); maps onto CAPRI bands; used to evaluate predicted protein-protein / TCR-pMHC complex structures against experimental ground truth. *Domain: bio.*

## E

**ENCODE** — Encyclopedia of DNA Elements. NIH consortium (since 2003) cataloging functional genome features via standardized high-throughput assays (RNA-seq, ChIP-seq, ATAC-seq, Hi-C, CAGE, bisulfite-seq, etc.). "ENCODE-style assays" = this standard set of functional genomics readouts, regardless of which project produced the data. *Domain: bio.*

**ESM-IF1** — Evolutionary Scale Modeling Inverse Folding (v1). Meta AI model (Hsu et al. 2022) trained on 12M protein structures; given a backbone, infers compatible sequences. Used as a structure-aware embedding source by structure-informed TCR-pMHC scorers (e.g. NetTCR-struc). *Domain: ml.*

## G

**GATK** — Genome Analysis Toolkit (Broad Institute). Canonical somatic/germline variant calling suite; reference for the Panel of Normals (PoN) concept reused by neoepitope filtering. *Domain: bio.*

**GKE** — Google Kubernetes Engine. Managed Kubernetes on GCP. *Domain: cloud.*

**GPS** — Genotype Presentation Score. Project-specific scoring formula combining per-allele MHCflurry presentation percentiles into a single rank for the patient's genotype; primary ranking key for top neoepitope candidates. *Domain: pipeline.*

**GTEx** — Genotype-Tissue Expression project. NIH consortium with healthy-tissue RNA-seq across ~50 tissues from ~1,000 donors; canonical reference for tissue-matched normal expression panels. *Domain: bio.*

**GVP-GNN** — Geometric Vector Perceptron Graph Neural Network. Architecture (Jing et al. 2021) for learning on protein 3D structures; respects rotational/translational equivariance natively. Used in NetTCR-struc and similar structure-aware TCR-pMHC scorers. *Domain: ml.*

## H

**HLA** — Human Leukocyte Antigen. Human MHC class I/II proteins (HLA-A/B/C class I; HLA-DR/DP/DQ class II); patient typing via OptiType drives MHCflurry's per-allele predictions. *Domain: bio.*

## I

**IC50** — half-maximal Inhibitory Concentration. Concentration (in nM) at which a peptide displaces 50% of a reference ligand from MHC; classic affinity metric, retained as informational column alongside the percentile-based presentation score. *Domain: bio.*

## M

**MCP** — Model Context Protocol. Anthropic's open spec (late 2024) for AI agents to discover and invoke external tools / data sources via per-purpose MCP servers; the standardized way to expose resources to LLMs without bespoke per-tool integration. *"USB-C for AI agents."* Used by Claude Code, the GitHub MCP server, filesystem servers, etc. *Domain: ml.*

**MHC** — Major Histocompatibility Complex. Cell-surface proteins that present peptides to T cells; class I (all nucleated cells, presents endogenous peptides) drives this pipeline. Human MHC = HLA. *Domain: bio.*

## O

**OKR** — Objectives and Key Results. Goal-setting framework popularised by Intel/Google: one qualitative *Objective* + 3–5 quantitative *Key Results*, usually quarterly. Common in software shops; we use S1–S7 lifecycle stages + iteration capacity instead. *Domain: project.*

## P

**PAE** — Predicted Aligned Error. AlphaFold/TCRdock per-residue-pair error estimate (in Å); summarised globally to gauge inter-domain confidence (e.g. relative MHC↔TCR docking confidence). *Domain: bio.*

**PDB** — Protein Data Bank ([rcsb.org](https://www.rcsb.org)). Public archive of experimentally determined 3D structures (~220k entries); training source for AlphaFold / Boltz / HERMES; reference ground truth for DockQ scoring. *Domain: bio.*

**pLDDT** — predicted Local Distance Difference Test. Per-residue confidence score (0–100) emitted by AlphaFold/TCRdock; high pLDDT = the model is confident in that region's local geometry. *Domain: bio.*

**pMHC** — peptide-MHC. The complex of a peptide bound in the MHC groove; the molecular target a TCR recognises. "TCR-pMHC interaction" = the central event in adaptive immune recognition. *Domain: bio.*

**PoN** — Panel of Normals. Aggregated reference set of normal samples (originally GATK's somatic variant calling concept); used to filter recurrent germline variation / artefacts that single matched normals miss. *Domain: bio.*

**PSI** — Percent Spliced In. Fraction of transcripts that include a given alternative exon/junction (vs skip); the standard quantitative readout of splicing in RNA-seq analyses. *Domain: bio.*

## R

**RCE** — Remote Code Execution. Vulnerability class where an attacker runs arbitrary code on a remote system, typically via input-validation or memory-safety bugs (e.g. CVE-2026-3854's crafted git push trigger). *Domain: security.*

**RMSD** — Root Mean Square Deviation. Average distance (in Å) between corresponding atoms after optimal superposition; the standard structural-similarity metric. Variants used in DockQ: backbone-only, all-atom, interface-only (iRMSD), ligand-only (LRMSD). *Domain: bio.*

## S

**SLURM** — Simple Linux Utility for Resource Management. The dominant HPC job scheduler in academia. *Domain: cloud.*

## T

**TCR** — T-Cell Receptor. Heterodimeric (α/β) receptor on T cells that recognises peptide-MHC complexes; modelled here via TCRdock for the top neoepitope candidate. *Domain: bio.*

**TF** — Transcription Factor. DNA-binding regulatory protein controlling gene expression; binding sites are one of the modalities ENCODE-style assays (ChIP-seq) and predictors like AlphaGenome map. *Domain: bio.*

## W

**WES** — Whole Exome Sequencing. Sequencing restricted to protein-coding regions (~1–2% of the genome); cheaper than WGS, sufficient for variant calls in coding regions. *Domain: bio.*

**WGS** — Whole Genome Sequencing. Sequencing of the entire genome (~3 Gbp human); alternative input for upstream filtering when matched-normal RNA-seq is missing. *Domain: bio.*
