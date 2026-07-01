# #736 baseline-predictor landscape - build brief

**Deep-research run #3 (predictor landscape for the #736 scoring harness), 2026-07-01.**
Source: 56 predictors catalogued -> 26 deep-dived -> 19 distinct after dedup, + 5 benchmark papers.
The workflow hit the usage limit mid-run; this brief was synthesized by hand from the 19 banked dossiers, run #2's scoring-design conclusions, and targeted re-verification of the baseline shortlist.

**Verification legend:** ✅ verified this pass (primary source) · ◻︎ dossier-asserted (run #3 deep-dive, not independently re-verified) · ⚠︎ verification deferred (fetch failed).

---

## 1. Recommended baseline set (ranked)

The harness ranks the A*02:01 registry positives by each score and compares enrichment against our `calibrated_immunogenicity_log_odds` (#907).
Because there is n=1 hard true-negative, everything below is a **ranking / enrichment** comparator, never a specificity comparator (run #2: specificity is structurally unmeasurable at n=1).

Two layers matter: a **presentation floor** (TCR-blind - what raw presentation already buys) and **immunogenicity comparators** (the post-presentation layer our score competes in).

### Tier 0 - mandatory controls (near-zero integration cost)
1. **MHCflurry 2.0 `presentation_score`** ◻︎ (it's our own upstream) - the decisive control. If `calibrated_immunogenicity_log_odds` does not out-rank raw `presentation_score` on the registry, the immunogenicity layer adds nothing. Already emitted in `mhc_presentation.tsv` - pull the column, no new dependency. Direction: higher = more presented (use `-presentation_percentile` to co-rank).
2. **NetMHCpan-4.1 EL (`%Rank_EL`)** ✅ - the field-standard presentation floor; the reference reviewers will expect. Drop-in, CPU, arbitrary standard-AA 8-14mers. Direction: `%Rank_EL` lower = better (negate to co-rank). Non-redistributable (DTU academic license) - pin version + command, do not vendor.

### Tier 1 - primary immunogenicity comparators (the interesting contrast)
3. **BigMHC_IM** ✅ - **the primary external comparator.** Open-source (KarchinLab/bigmhc), CPU-OK (no GPU), scores arbitrary peptide+HLA (even nonsense inputs), CSV in, both EL+IM shipped (`-m`). Same direction as ours (higher = more immunogenic); rank-based metric makes the probability-vs-log-odds scale mismatch irrelevant. ~5 GB repo. Confirm exact LICENSE on clone (README references a LICENSE file; type not quoted).
4. **PRIME 2.1** ✅ - second immunogenicity comparator, and the closest **training-domain analogue** to our registry (trained on ELISpot-validated point-mutation neoantigens). Open-source, academic-only, CPU. **Two caveats:** (a) requires **MixMHCpred v3.0+ installed in PATH** (so standing up PRIME also gives us MixMHCpred as a presentation baseline for free); (b) output is `%Rank`, lower = more immunogenic -> **invert** to co-rank. Repo ships **2.1** (dossier said 2.0 - superseded).
5. **IEDB Class-I Immunogenicity (Calis 2013)** ✅◻︎ - the canonical, most-cited reference immunogenicity predictor; sequence-only, higher = more immunogenic, allele used only for anchor-position masking. Confirmed available via both legacy `tools.iedb.org/immunogenicity/` and the NG IEDB tools (T-cell class-I). Include as the "classic baseline everyone reports." **Channel resolved (2026-07-01):** a standalone package (tar.gz) IS offered for offline use, academic-free (commercial via DTU). **BUT a licensing flag surfaced:** the download page states the bundle is governed by the **NetMHC 3.0 Academic License**, whose terms include single-site/internal-use-only, no redistribution, and - critically - **"benchmark results cannot be published without prior written consent."** For a *published* open benchmark this is a live constraint (see §5 open questions). The web / NG-IEDB-tools route may carry different terms - worth checking which channel #736 uses.

### Tier 2 - optional additional comparators (include if cheap / for breadth)
- **MixMHCpred 2.2/3.0** ✅ - presentation baseline; comes free with PRIME (v3.0, academic-free, needs MAFFT). Open-source.
- **DeepHLApan** ✅(GPU-caveat) - open-source (GPL-2.0) + docker, emits an explicit immunogenicity score. **But the shipped docker image is built on CUDA 9 / cuDNN 7 (GPU); CPU execution is unconfirmed.** Given our CPU-only posture, treat as GPU-gated until a CPU run is proven - do not assume drop-in.
- **ImmuneApp / ImmuneApp-Neo** ✅ - open-source (MIT), emits both EL presentation and a neo immunogenicity score (ImmuneApp-Neo module).
- **NeoTImmuML** ✅◻︎ - peptide-only weighted ensemble (RF+LightGBM+XGBoost, no allele needed), 8-13mers, reported AUC 0.8865 (Frontiers Immunol 2025). **Obtainability caveat:** access is via the TumorAgDB2.0 portal (tumoragdb.com.cn); no public GitHub repo / one-click weights download was confirmable - softer than a clean open-source dependency.
- **Repitope** ✅ - **orthogonal**: allele-agnostic (no HLA input; class-specific MHC-I vs -II models), TCR-repertoire contact-potential signal. R package (bridge from our Python harness). Different failure mode from the presentation-based scorers - useful as an independent axis.
- **DeepImmuno-CNN** ✅ - MIT, CPU/GPU. **Coverage caveat: accepts only 9- and 10-mers, REJECTS 11-mers** - it would silently drop our 11-mer positives (e.g. MAN2C1 `LLLGIAKLLKV`). Include only with an explicit "11mers excluded" footnote, or exclude.

---

## 2. Rejected as baselines (and why)

| Predictor | Why not a #736 baseline |
|---|---|
| **IRIS** ◻︎ | A discovery *pipeline* (rMATS -> translate -> IEDB -> tumor-specificity tiers), not a per-peptide scorer. Its only peptide number is IEDB median IC50, which we can compute directly. Heavyweight to stand up. |
| **NMD-escape / frameshift-indel evidence (Litchfield/Turajlic)** ◻︎ | Emits a *categorical* per-variant NMD-escape flag + a between-class enrichment - no continuous per-peptide value, cannot rank A*02:01 candidates. |
| **Roudko MS-frameshift poly-epitopes (Cell 2020)** ◻︎ | Not a released predictor; its per-peptide "score" is just NetMHC 4.0 binding, redundant with (and weaker than) our MHCflurry presentation ranker. |
| **HLApollo** ◻︎ | Mechanically near-drop-in, but presentation-only (redundant with NetMHCpan-EL) and under a restrictive academic-binary license. Reference-only. |
| **HLAthena** ◻︎ | Primarily web-server-only (research use); presentation-only. Usable but lower priority than the open-source presentation baselines. |

---

## 3. Splice-admissibility + availability table (all 19)

| Predictor | Layer | Admits our splice/frameshift 9-11mers | Obtainable | Verify |
|---|---|---|---|---|
| NetMHCpan-4.1 EL | presentation | yes (standard-AA arbitrary) | standalone, academic, non-redist | ✅ |
| NetMHCpan-4.2 EL | presentation | yes-caveat (canonical-ligand prior) | standalone, academic | ◻︎ |
| MHCflurry 2.0 presentation | presentation | yes (arbitrary) | open-source (ours) | ◻︎ |
| MixMHCpred 2.2/3.0 | presentation | yes (arbitrary, pan-allele) | open-source (v3.0), academic-free; needs MAFFT | ✅ |
| HLAthena | presentation | yes-caveat | web-server (docker/Terra) | ◻︎ |
| ImmuneApp / -Neo | presentation + immuno | yes (arbitrary) | open-source (MIT) | ✅ |
| HLApollo | presentation | yes-caveat | binary, restrictive license | ◻︎ |
| **BigMHC_EL / _IM** | presentation + **immuno** | **yes (arbitrary, even nonsense)** | **open-source, CPU-OK** | ✅ |
| **PRIME 2.1** | **immuno** | yes (8-14mer arbitrary) | open-source, academic; needs MixMHCpred | ✅ |
| **IEDB-Calis 2013** | **immuno** | yes (arbitrary, 9mer validated) | standalone tar.gz confirmed, academic-free; **NetMHC-3.0-license "no publishing benchmark results w/o consent"** | ✅ |
| **DeepImmuno-CNN** | **immuno** | **9-10mers only - REJECTS 11mers** | open-source (MIT) | ✅ |
| DeepHLApan | binding + immuno | yes (arbitrary) | open-source (GPL-2.0) + docker; **docker image is CUDA9/cuDNN7 GPU - CPU uncertain** | ✅ |
| DeepNeo-v2 | presentation + tcr | yes-caveat (9mer fixed MHC-I) | web-server-only | ◻︎ |
| Repitope | immuno (allele-agnostic, class-specific) | yes (8-11mer arbitrary) | open-source (R, MIT) | ✅ |
| NeoTImmuML | immuno (peptide-only) | yes (8-13mer) | **TumorAgDB2.0 portal-gated; no confirmed public GitHub repo** | ✅◻︎ |
| IRIS | pipeline | n/a (not a scorer) | open-source, heavyweight | ◻︎ |
| NMD-escape | variant-level | n/a (categorical) | open-source (annotation only) | ◻︎ |
| Roudko MSI-frameshift | study | n/a (no scorer) | reference-only | ◻︎ |

---

## 4. Methodology implications for #736 (from benchmark papers + run #2)

1. **Metric: enrichment/ranking, not calibrated accuracy.** Immunogenic base rate is ~5% (Frontiers Immunol 2023, `10.3389/fimmu.2023.1094236`) and there is no standardized negative set field-wide. With our n≈81 positives / n=1 hard negative, report positive-rank ECDF, top-k recall, and AUROC/AUPRC of positives-vs-decoys - **not** specificity or calibrated PPV. This matches run #2's power-ledger red-light conclusion exactly.
2. **Direct methodological precedent exists.** NAR Cancer 2024 (`zcae002`) benchmarked *exactly* NetMHCpanEL-4.1, MHCflurry 2.0.6, MixMHCpred 2.2, and PRIME 2.0 on CEDAR neo-epitopes - i.e. our Tier-0/Tier-1 shortlist is the field-standard comparator set. Mirror their metric definitions so #736 is directly comparable.
3. **Leakage is the sharpest threat, and it is worse for us than usual.** Many of our registry positives were *sourced from IEDB deposits* (Kwok GNAS/RPL22, POSTN, IRIS rows - see project memory). The immunogenicity baselines **IEDB-Calis, PRIME, and BigMHC_IM are trained on IEDB/CEDAR**, and the EL baselines (NetMHCpan, MHCflurry, BigMHC_EL) are trained on MS eluted-ligand sets that may contain our positives. So a naive comparison **structurally inflates the baselines** via train/test circularity. Required mitigation: annotate each registry positive with its membership in each baseline's training corpus, and prefer **leave-one-study-out** (run #2's LOSO) so a baseline is never scored on a peptide it trained on. Treat any un-annotated baseline lift as a **ceiling, not a fair number.**
4. **TESLA (Wells 2020) is the prior, but diverges.** TESLA defined the immunogenicity parameters and negative handling the field uses - but explicitly **excluded splice isoforms**, was multi-allele, and was far larger. #736 borrows its metric vocabulary and diverges on scope (splice-specific, tiny, mono-allele A*02:01). Cite it as the precedent we extend into the splice regime.

---

## 5. Open questions / human-decision points

1. **Leakage annotation is real work.** Building the per-baseline training-set-membership map (which of our 81 positives each baseline has seen) is a prerequisite for an honest comparison, not an optional extra. Scope it as its own task. Decision: do we build it, or ship #736 with an explicit "baselines are leakage-inflated ceilings" caveat and defer the clean LOSO comparison?
2. **How many baselines to actually stand up?** Recommend the **core 5** (MHCflurry [free], NetMHCpan-4.1 [pin], BigMHC_IM, PRIME 2.1 [+MixMHCpred], IEDB-Calis). Tier-2 adds breadth at real integration cost - worth it only if we want a predictor-panel figure.
3. **CPU-only is not a blocker.** Every recommended baseline runs on CPU (BigMHC, PRIME, MHCflurry, IEDB, DeepImmuno all confirmed CPU-OK) - good, given the GCP/GPU decommission. No revival needed for this benchmark.
4. **11-mer coverage.** DeepImmuno drops 11-mers; if we include it, footnote the reduced-N. Confirm how many registry positives are 11-mers before deciding.
5. **Redistribution.** NetMHCpan, PRIME, HLApollo are non-redistributable -> the open-benchmark repo pins version + exact command, never vendors the binary. The open-source ones (MHCflurry, BigMHC, MixMHCpred, DeepImmuno-MIT, ImmuneApp-MIT, DeepHLApan-GPL) can be scripted/vendored.
6. **Publishing-license blocker (NEW, verified 2026-07-01).** The **IEDB-Calis standalone** carries the **NetMHC 3.0 Academic License**, which states *benchmark results cannot be published without prior written consent* - and the same DTU clause plausibly covers **NetMHCpan-4.1** (our Tier-0 floor). #736 is a *published* open benchmark, so this is a real fork: (a) obtain DTU/IEDB written consent (the NAR Cancer 2024 benchmark published NetMHCpan/PRIME scores, so consent or non-enforcement is the norm - but confirm), or (b) check whether the **web / NG-IEDB-tools** route carries lighter terms, or (c) lean the *publishable* headline on the fully-open baselines (MHCflurry, BigMHC, MixMHCpred - all OSI/permissive, no publish restriction) and treat NetMHCpan/IEDB-Calis as internal reference points. Decide before the results go in a manuscript.
7. **Which score is "ours"?** Confirm the harness compares `calibrated_immunogenicity_log_odds` (#907) as the headline, with `presentation_score` as the Tier-0 control - not the other way around.

---

## 6. What NOT to do

- **Do not report specificity / calibrated PPV / ROC-specificity at n=1 hard negative** - structurally unmeasurable (run #2). Ranking/enrichment only.
- **Do not treat IEDB-Calis / PRIME / BigMHC_IM lift as a clean win without leakage annotation** - our positives partly come from their training corpus. This is the single easiest way to publish an inflated result.
- **Do not wire IRIS / NMD-escape / Roudko in as per-peptide scorers** - they are pipelines / categorical / no-scorer.
- **Do not add a predictor needing inputs we lack** (structure, patient expression beyond pipeline output) as a *core* baseline.
- **Do not vendor non-redistributable binaries** into the open repo - pin + document instead.
- **Do not silently drop 11-mers** by adding DeepImmuno without a footnote.

---

## Provenance / caveats
- 10 predictors verified against primary source (NetMHCpan-4.1 via workflow; BigMHC, PRIME 2.1, DeepImmuno, MixMHCpred, ImmuneApp, Repitope, DeepHLApan, NeoTImmuML, IEDB-Calis via targeted WebFetch/WebSearch). The remaining rows (IRIS, NMD-escape, Roudko, HLApollo, HLAthena, NetMHCpan-4.2, DeepNeo-v2; MHCflurry is our own) are run #3 deep-dive assertions, ◻︎ - low value since they are rejected-as-baseline or already-known.
- Corrections made during verification: PRIME repo ships **2.1** (not 2.0) and needs MixMHCpred in PATH; DeepImmuno **rejects 11-mers**; Repitope is **allele-agnostic but class-specific** (not flatly "MHC-agnostic"); **DeepHLApan's docker is GPU/CUDA9** (CPU unproven - not a clean CPU baseline); **NeoTImmuML is TumorAgDB2.0-portal-gated** (no confirmed public repo).
- IEDB channel resolved 2026-07-01 (earlier fetch failures were a client-side VPN block, not IEDB): standalone tar.gz confirmed + the NetMHC-3.0-license publish-restriction surfaced (see Tier-1 #5 / §5 item 6).
- Raw dossiers: `scratchpad/run3_partial_result.json`. Run status: `scratchpad/run3_STATUS.md`.
- Dedup note: BigMHC, DeepImmuno, DeepNeo, PRIME appeared under variant names in the raw catalog; collapsed here.
