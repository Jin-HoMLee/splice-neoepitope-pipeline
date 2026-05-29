# ASNEO eval — design spec (Issue #546)

**Date:** 2026-05-29
**Author:** Scientist
**Parent Issue:** [#546](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/546) — *research: evaluate ASNEO as splice-neoepitope candidate generator (peer to our HISAT2/regtools + STAR/SJ.out.tab paths)*
**Source eval:** [#258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/258) (NeoGuider) — ASNEO surfaced as the actual splice tool NeoGuider delegates to.

> ⚠️ **DOI corrected during the eval:** the Issue #546 body and this spec's first draft cited `10.18632/aging.103581`; the verified DOI is **`10.18632/aging.103516`** (PubMed PMID 32697765 / PMC7425491 / Aging / GitHub). The references below are corrected.

## Goal

Decide whether **ASNEO** (Zhang et al., *Aging* 2020, DOI `10.18632/aging.103516`) — our closest published peer (personalized alternative-splicing neoepitopes from RNA-seq) — should change our pipeline. The decision is framed as the five pipeline-fit modes from the Issue body: **Triage / Replacement / Cross-check / Component reuse / Reject**.

This is a **desk eval** (decision-only), consistent with the sibling evals [#218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) (HERMES), [#201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) (ImmSET), [#188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) (Boltz-2), [#258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/258) (NeoGuider). **No code execution / no ASNEO run** — the multi-day Python2-era dependency archaeology is not justified for a P3 informational eval.

## Success criteria (maps to #546 ACs)

- Verified Zotero entry on DOI `10.18632/aging.103516` with a 3-section HTML note (**Findings / Methods / vs. our pipeline**), added after a dedup pre-check.
- ASNEO repo state verified: maintenance status, language (Py2 vs Py3), code + weights availability, license.
- Head-to-head comparison of ASNEO's junction-calling backbone vs. ours documented across five axes (below).
- A 5-mode pipeline-fit decision recorded with rationale on #546.
- Tool-primer Quarto deck shipped at `research/evals/issue_546_asneo/slides.qmd`, visually verified.
- #258 verdict comment updated *if* the ASNEO eval changes the NeoGuider component-reuse calculus.
- Lab-notebook entry capturing the decision and rationale.

## Execution — Workflow (approach C)

Multi-agent Workflow because (i) ultracode is on and (ii) two memory rules bite hardest on tool evals: *never quote numerical findings from unreachable sources* (caught on PR #518) and the head-to-head needs **our-side** facts a general research harness won't have.

**Phase 1 — Gather (parallel, schema'd):**
1. ASNEO *paper claims* — junction backbone, normal-junction filtering, peptide assembly + frame logic, MHC predictor, validation cohort.
2. ASNEO *repo state* — last commit, Py2 vs Py3, code + weights availability, license.
3. *Our backbone* facts read from source — `workflow/scripts/bed12_to_junctions.py`, `star_sj_to_junctions.py`, `filter_junctions.py`, `assemble_contigs.py`, `translate_peptides.py`, `run_mhcflurry.py`.

**Phase 2 — Adversarial verify:** each load-bearing claim (validation cohort + N, any reported AUROC/performance, MHC-predictor identity, license) refuted-by-default against ≥2 independent sources. Anything unverifiable is hedged or dropped — no quoting from unreachable sources.

**Phase 3 — Synthesize (main loop, not an agent):** Scientist writes the 5-mode decision, Zotero note, and deck from verified findings only. Repo writes stay under direct control.

## Head-to-head comparison axes (slide 5 + Zotero "vs. our pipeline")

| Axis | ASNEO | Our pipeline |
|---|---|---|
| Junction calling | STAR → `SJ.out.tab` (+ any ASNEO filtering/motif selection — to verify) | HISAT2/regtools → `bed12_to_junctions.py`; STAR → `star_sj_to_junctions.py` (strand rescue from intron motif) |
| Normal-junction filtering | matched-normal / GTEx panel / both — to verify | matched-normal at junction level (annotated → normal_shared → tumor_exclusive); GTEx pan-tissue planned (#212) |
| Peptide assembly + frame | flank handling + frameshift detection — to verify | `assemble_contigs.py` (`bedtools getfasta -s`) → `translate_peptides.py` |
| MHC predictor | netMHCpan / MHCflurry / in-house — to verify | MHCflurry `Class1PresentationPredictor` (presentation likelihood) |
| Validation cohort | TESLA-like / self-curated — to verify | n/a (we are the pipeline under comparison) |

## Deliverables

Deck-only family at `research/evals/issue_546_asneo/`:
- `slides.qmd` + `refs.bib` (theme `[simple, ../../slides/custom.scss]`; `csl: ../../slides/nature.csl`); rendered `slides.html` + `slides_files/` kept in the working tree.
- Zotero 3-section HTML note (post dedup-check).
- 5-mode decision comment on #546; conditional #258 verdict update.
- Lab-notebook entry (written *after* review, before merge).

## Deck slide map (title + 9, mirrors HERMES `issue_218_hermes/`)

1. Title
2. The question — *"Should ASNEO change how we generate splice neoepitopes?"*
3. ASNEO at a glance — Zhang 2020, license, repo state, key facts
4. How it works — STAR → `SJ.out.tab` → filtering → peptide assembly → MHC predictor → ranking
5. Head-to-head — ASNEO vs. our backbone across the five axes above
6. Integration map (mermaid) — where ASNEO or a component would plug into our DAG
7. Reasons by mode — Triage / Replacement / Cross-check / Component reuse / Reject table
8. Decision — verdict + rationale (+ sub-issue if integrate / decline-note if not)
9. Open questions & caveats
10. References

## Lifecycle

1. #546 → **In progress** (`47fc9ee4`) — done before branching.
2. Branch `research/scientist/issue-546-asneo-eval` via `gh issue develop` — done.
3. Build deliverables (Workflow → synthesis → deck + Zotero + decision).
4. `quarto render` + **visually verify every slide** (reveal.js silently overflows; headless screenshots per slide).
5. Commit (separate from push); surface SHA + diff and wait for OK before pushing.
6. PR → flip PR **Ready for review** (`8bf9192f`) — note: a PostToolUse hook (#558) may auto-board the PR + set Status; verify, don't double-set.
7. Offer `@claude review`.
8. Lab-notebook entry (after review, before merge).
9. Slide-surfaced actionables (open questions, follow-ups) filed on the board *before* the PR merges.
10. Tick ACs via `scripts/audit_and_merge.sh <PR>`.

## Out of scope

- Any code execution / ASNEO run (desk eval only).
- #547 (KDE + centered isotonic calibration) — already spun off from #258; not re-litigated here.
- MHC-II support.
- Re-running sibling evals; cite them for cohesion only.

## Risks / open questions

- ASNEO repo may be Python2-era and unmaintained — affects the "code + weights available" AC and the Replacement/Cross-check feasibility verdict.
- Some paper claims (specific PSR/normal-filter thresholds, per-cohort AUROC) may be unverifiable from accessible sources — hedge or defer rather than quote, per the unreachable-source rule.
- If ASNEO's normal-junction filtering meaningfully overlaps the planned GTEx pan-tissue filter (#212), the decision may surface a cross-reference to #212 rather than a standalone integration.
