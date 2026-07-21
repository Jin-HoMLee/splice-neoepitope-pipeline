# #1049 class-II presentation-predictor landscape - build brief

**Deep-research run #5 (open MHC class-II presentation-predictor landscape for the comprehensive-toolkit architecture), 2026-07-20.**
Fan-out: 6 search angles -> 23 sources fetched -> 98 falsifiable claims extracted -> top 25 adversarially verified (3-vote, need 2/3 refutes to kill).
Result: 21 claims confirmed (every one a **3-0** unanimous vote), 4 refuted, 0 left unverified.
Run: `wf_f1a9b274-0bc` (106 agents, ~6.85M tokens; provenance only - the per-agent outputs are session-ephemeral).
Feeds: [Issue #1045](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1045) Layer B (presentation), the class-II frontier of the open/RNA-native/validated thesis vs pVACtools.

**Verification legend:** ✅ verified this pass (primary source, 3-0) · ❌ confirmed gate-failure (verified and disqualifying) · ⚠︎ genuinely unresolved or refuted claim (do not rely on without re-check).

---

## 1. Headline

There are **two clean winners** that clear all three gates - permissive+commercial license, MS-eluted-ligand *presentation* training, and an open, retrainable architecture (CAPTAn additionally ships trained weights; AEGIS's weights-release is unconfirmed this pass - see §3):

1. **AEGIS** (Novartis) - BSD-3-Clause, compact ~3.5M-param transformer, `pytorch-lightning`, trained on MHC-II MS eluted-ligand (EL) + binding-affinity (BA) data from IEDB. Repo: [`github.com/Novartis/AEGIS`](https://github.com/Novartis/AEGIS). Paper: Bioinformatics 2023, [btad469](https://academic.oup.com/bioinformatics/article/39/8/btad469/7234610) (CC-BY).
2. **CAPTAn** (Broad / Xavier lab) - BSD-3-Clause, convolutional motif-detector core, trained on a large MS-eluted-ligand HLA-II immunopeptidome (1.4M peptides / 237,905 non-overlapping ligands / 87 alleles), ships model implementations **and** trained parameters. Repo: [`gitlab.com/xavier-lab-computation/public/captan`](https://gitlab.com/xavier-lab-computation/public/captan) (via `broad.io/captan`). Paper: Immunity 2023, [PMC10519123](https://www.sciencedirect.com/science/article/pii/S1074761323002261).

**Recommendation: retrain / extend an open architecture** (AEGIS or CAPTAn as the base), **not** build-our-own from scratch, and **not** wrap a blocked tool.
Build-our-own is unwarranted when two permissive, presentation-trained, retrainable bases already exist.
Wrapping a blocked tool (NetMHCIIpan / MixMHC2pred / BERTMHC / MHCnuggets-II) is barred by gate 1 or gate 2.

This is the direct class-II analogue of the class-I decision: for class I we adopt MHCflurry (Apache-2.0); for class II, AEGIS/CAPTAn (BSD-3) are the license-clean presentation bases, which is exactly the gap the 2026-07-06 note flagged ("no license-clean class-II presentation predictor exists") - and the survey now says the frame was too pessimistic: it is not that none exist, it is that the *field-standard, most-cited* ones (NetMHCIIpan, MixMHC2pred) are the blocked ones, while two lesser-known BSD-3 options clear the bar.

---

## 2. Three-gate scored table (all candidates)

| Candidate | Gate 1 License | Gate 2 Presentation vs Binding | Gate 3 Trainability | Verdict |
|---|---|---|---|---|
| **AEGIS** (Novartis) | ✅ BSD-3-Clause (permissive, commercial-OK) | ✅ MS-EL presentation (+ BA) from IEDB | ✅ open ~3.5M-param transformer, pytorch-lightning, public repo | **GO - build base** |
| **CAPTAn** (Broad/Xavier) | ✅ BSD-3-Clause | ✅ MS-EL immunopeptidome (87 alleles) | ✅ open CNN + **released trained weights** | **GO - build base** |
| **BERTMHC** | ❌ NEC noncommercial ("commercial usage not granted") | ✅ dedicated MS-EL presentation model | ✅ TAPE-BERT + MIL, open code | **NO-GO (gate 1)** |
| **MHCnuggets-II** | ❌ JHU Academic (non-commercial; Genentech holds separate commercial) | ❌ **class-II trained on affinity only** (no allele-specific class-II HLAp) | ✅ disclosed LSTM, retrainable | **NO-GO (gate 1 + gate 2)** |
| **NetMHCIIpan-4.0** | ❌ DTU non-commercial **+ no redistribution / no toolkit inclusion** | ✅ MS-EL + BA presentation | (n/a - blocked) | **NO-GO (gate 1, hard)** |
| **MixMHC2pred** | ❌ academic-only; separate paid for-profit license | ✅ MS-EL presentation | (n/a - blocked) | **NO-GO (gate 1)** |
| **arXiv:2512.14011** (multi-scale) | (n/a - not a shipped predictor) | ✅ models BA + peptide-EL + novel antigen-level EL | dataset+benchmark framework | **NOT A PRODUCT** - but its curated IEDB dataset is a reusable training resource ⚠︎ (data license unverified) |
| **PIA-M** (ikmb) | ⚠︎ permissive-license claim **refuted** - gate 1 unresolved | ✅ MS immunopeptidome (multimodal transformer) | code public | **HOLD** - verify license before considering |
| **NetMHCIIphosPan** | (DTU family - non-commercial) | presentation (class-II **phospho**-specific) | (blocked) | off-axis (phospho niche); already shelved `3FGTW6VK` |

---

## 3. Why the two winners, in detail

### AEGIS ✅ (all three gates, 3-0 each)
- **License:** paper code-availability points to `github.com/Novartis/AEGIS`, verified live as a public Novartis repo under **BSD-3-Clause** (permissive, commercial use and redistribution allowed).
- **Presentation:** trained on "human (H) and mouse (M) peptide data from MHCII binding affinity (BA) and mass spectrometry (MS)-eluted ligand (EL) assays from IEDB" - the MS-EL component is the presentation signal.
- **Trainability:** ~3.5M params, 4-layer transformer encoder, embedding dim 128, 2 attention heads, explicitly built in `pytorch-lightning` for portability. Small enough to retrain on our own immunopeptidome / registry data cheaply.
- **Weights caveat:** the verified evidence covers the open architecture + training code + public repo; it does **not** confirm AEGIS ships usable pretrained weights. So treat AEGIS as the *train-on-open-architecture* option and CAPTAn as the *adopt-pretrained-weights* option until the repo is checked at integration time.

### CAPTAn ✅ (all three gates, 3-0 each)
- **License:** repo displays a "BSD 3-Clause New/Revised License" badge.
- **Presentation:** training substrate is unambiguously MS immunopeptidome (naturally processed / presented peptides), not IC50 affinity. Caveat: CAPTAn-core frames its *output* as an HLA-II binding probability even though its *training data* is eluted-ligand - so read the label, trust the data.
- **Trainability:** the paper states "implementations of CAPTAn-core, CAPTAn-context and CAPTAn models as well as their corresponding trained parameters can be obtained in open source" - open architecture **and** open weights. `CAPTAn-context` adds a context-aware variant useful for antigen-level scoring.

**AEGIS vs CAPTAn:** AEGIS is tinier and trivially retrainable; CAPTAn ships pretrained weights over a much larger allele panel (87). These are complementary, so ensembling both is on the table (see open questions).

---

## 4. Why the rest are out

- **BERTMHC** ⚠︎ - a *genuine* MS-EL presentation model (passes gate 2), disqualified purely on license: NEC Laboratories Europe grants "noncommercial research purposes" only and states "Commercial usage/research is not granted." A clean gate-1 fail.
- **MHCnuggets-II** ⚠︎ - fails **two** gates. Its own Methods: "Due to the lack of allelic-specific HLAp training data for class II, all MHC class II networks were trained only on binding affinity measurements" (gate 2 fail). And the LICENSE is a "JHU Academic Software License Agreement" restricted to non-commercial use (gate 1 fail; Genentech holds a separately negotiated commercial license). Two web/metadata claims that MHCnuggets is permissively open-source and incorporates class-II presentation signal were both **refuted 0-3** - the PyPI "Apache-2.0" tag is wrong versus the actual LICENSE file.
- **NetMHCIIpan-4.0 / MixMHC2pred** ⚠︎ - the confirmed license-blocked baselines. DTU: "Any use of the software which results in any form of commercialization is not allowed" and "not give the program to third parties or grant licenses on software, which include the Software" (a hard no-redistribution / no-toolkit-inclusion clause on top of non-commercial). MixMHC2pred README: for-profit users "are required to obtain a separate license." Keep both as *reference EL baselines* only, pinned by version + command, never vendored.
- **arXiv:2512.14011** (the shelf item `E8QZE2VA`, "Accelerating MHC-II Epitope Discovery via Multi-Scale Prediction") - not a single shippable predictor to wrap; it is a curated IEDB dataset + a modular benchmark framework spanning BA, peptide-EL, and a novel antigen-level EL task. The dataset is potentially a reusable training resource for the retrain path, **but** the specific claim that it merely aggregates the NetMHCIIpan-4 / MixMHC2pred-2 MS sets was **refuted 0-3**, so its data provenance and reuse license are genuinely open questions - verify before training on it.
- **PIA-M** (ikmb/PIA-inference) - credible MS-immunopeptidome presentation model, but its permissive-license (OmLiT MIT) claim was **refuted 0-3**, and the bioRxiv PDF 403'd during verification, so gate 1 is unresolved. Worth a targeted license check as a possible third base if AEGIS/CAPTAn allele coverage falls short.

---

## 5. Caveats (from the run)

- **Re-verify LICENSE at integration time.** Repo licenses can change; confirm the AEGIS and CAPTAn `LICENSE` files at the moment we vendor/depend on them.
- **CAPTAn weights not byte-enumerated.** The paper asserts trained parameters are released; the repo tree was not walked file-by-file this pass.
- **PIA-M gate 1 is genuinely open** (refuted license claim + 403 on the primary), so PIA-M's presentation classification rests on title + search digests, not a verbatim fetch.
- **The arXiv:2512.14011 dataset license is unverified** - do not assume free commercial reuse.
- All 21 confirmed claims were unanimous 3-0, so the two-winner conclusion itself is high-confidence; the caveats above are integration-time checks, not doubts about the headline.

---

## 6. Open questions (carry to #1045 Layer B execution)

1. What is the actual **data license** on the arXiv:2512.14011 curated IEDB peptide+antigen MHC-II dataset - can it retrain a commercial model?
2. What license do **PIA-M** and its OmLiT preprocessing actually carry (given the refuted permissive claim)?
3. On held-out MS immunopeptidome benchmarks, how do **AEGIS vs CAPTAn** compare head-to-head (and against the blocked NetMHCIIpan-4.0 EL baseline) - pick a primary base, or ensemble both?
4. Do **CAPTAn's released weights** (and AEGIS's, if the repo ships usable pretrained weights) cover the patient HLA-II alleles the toolkit targets, or is retraining/extension needed regardless of license?

---

## 7. What feeds the board

- **-> [Issue #1045](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1045) Layer B:** the recommendation (retrain-an-open-architecture on AEGIS or CAPTAn) and the four open questions above. This resolves the "class II is the build frontier, no license-clean predictor exists" note into a concrete, actionable base choice.
- **Zotero folder 3** (MHC presentation & immunogenicity) - proposed adds (surface-then-add, pending confirmation): **AEGIS** (btad469), **CAPTAn** (Immunity 2023), the **Frontiers Immunol 2024 class-II presentation benchmark** ([10.3389/fimmu.2024.1293706](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1293706/full), the best single gate-2 entry point), and optionally **PIA-M** and **BERTMHC** as blocked-but-citable references. `E8QZE2VA` (arXiv:2512.14011) and `3FGTW6VK` (NetMHCIIphosPan) are already shelved.

## Strategic through-line

The class-II frontier is **not** a build-from-scratch problem and **not** a dead end.
The blocked incumbents (NetMHCIIpan, MixMHC2pred) are the reason the frontier looked closed, but two BSD-3 presentation-trained, retrainable bases (AEGIS, CAPTAn) clear all three gates, so the open/commercial-usable comprehensive toolkit can have a class-II presentation layer on the same terms as its class-I MHCflurry layer.
The remaining work is a base-selection benchmark (AEGIS vs CAPTAn vs ensemble) and an allele-coverage check, not a licensing dead-end.
