# Junction Mapping + Documented Labeling Scheme — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the 79-row open splice-immunogenicity registry a documented, citable positive/negative labeling scheme plus per-peptide junction annotation, so it can ground the #680 benchmark reproducibly.

**Architecture:** Treat the registry TSV as the data artifact and a standalone Python validator as the test spine. Write the scheme doc first (source of truth), encode it in a validator that fails on the un-annotated registry, then add four columns (`evidence_strength`, `label_rationale`, `junction_id`, `junction_mapping_grade`) and a materialized #681 decoy seed until the validator passes. `evidence_strength` is derived mechanically from `readout`; junction annotation is a per-source procedure with a strict no-inference rule.

**Tech Stack:** Python 3.14 via `research/.venv`; pandas; tab-separated `registry.tsv`. No new dependencies.

## Global Constraints

- **No junction inference.** Never translate or reconstruct a junction coordinate a source did not publish. Annotate only what the source states; flag the rest. (Spec §Scope guard, §Junction mapping.)
- **Existing verified rows are ground truth.** Do not re-label outside the 14 boundary rows; any boundary change requires a cited scheme rule and a logged old→new entry. (Spec §Approach fork 1.)
- **Scope guard:** do NOT add `assay_context` (that is #823), do NOT build the FASTA+label file (that is #736), do NOT generate Tier-3 synthetic decoys (harness-time). (Spec §Scope guard.)
- **Prose style:** American spelling (`tumor`, `recognize`); plain dash `-`, never `-`; one sentence per physical line in long Markdown. (User global instructions.)
- **Terminology:** "presentation/presenter", not "binding/binder"; "splice-junction-derived neoepitope", never bare "splice peptide". (Project memory.)
- **All paths** are under `research/experiments/issue_680_splice_immunogenicity_registry/` unless stated; run Python as `research/.venv/bin/python`.

---

## File Structure

| Path | Responsibility | Task |
|---|---|---|
| `LABELING_SCHEME.md` | The documented, citable scheme + decoy-construction approach | 1 |
| `validate_registry.py` | Schema + consistency validator (the test spine) | 2 |
| `registry.tsv` | + `evidence_strength`, `label_rationale` columns; boundary re-audit | 3 |
| `registry.tsv` | + `junction_id`, `junction_mapping_grade` columns | 4 |
| `junction_evidence_by_source.md` | Working note: each source's published junction unit → grade | 4 |
| `decoy_negatives/presented_decoys_681.tsv` | Materialized Tier-2 presented-decoy seed (13 rows) | 5 |
| `README.md`, `PROVENANCE.md` | New-column docs, per-grade tally, re-audit changelog | 6 |

---

## Task 1: Write `LABELING_SCHEME.md`

The scheme is the source of truth the validator (Task 2) encodes and the data (Tasks 3-5) conform to. No code; the deliverable is the document.

**Files:**
- Create: `research/experiments/issue_680_splice_immunogenicity_registry/LABELING_SCHEME.md`

**Interfaces:**
- Produces: the controlled vocabularies `evidence_strength ∈ {strong, weak, hard, soft, na}` and `junction_mapping_grade ∈ {coords, event-id, gene-mechanism, none}`, and the label rules that Tasks 2-5 depend on verbatim.

- [ ] **Step 1: Write the scheme document.** Sections, in order:

  1. **Purpose** - one paragraph: a citable rule set making the registry's labels reproducible for the #680 benchmark.
  2. **Positive labels** (`label=positive`):
     - `evidence_strength=strong` - a measured effector readout on the peptide: IFN-γ, cytotoxicity/killing (LDH, Incucyte, caspase), ELISpot, granzyme-B/GZMB, CD107a/CD137, degranulation, cloned/engineered-TCR activation, TCR-T, or MHC-stabilization paired with function.
     - `evidence_strength=weak` - antigen-specific detection only: tetramer / dextramer / multimer binding with no effector readout. (Mass-spec `MS` is presentation, not T-cell function - it does NOT upgrade a detection-only row to strong.)
  3. **Negative labels** (`label=negative`):
     - `hard` - splice-derived, MHC-presented, functional assay performed, no T-cell response.
     - `soft` - tested-negative in healthy-donor in-vitro sensitization (failed-to-prime ≠ intrinsically non-immunogenic).
  4. **Excluded from the scorable set** (retained in TSV, never deleted):
     - `presentation-prevalence` → `label=untested`, `evidence_strength=na` (no functional assay).
     - `functional-nonscorable` → exact sequence unpublished; cannot be keyed.
     - `candidate` → positive whose 2nd adversarial-verify pass did not confirm; held until re-verified.
     - `negative-control-not-splice` → fails the splice gate; control only, `evidence_strength=na`.
  5. **Confidence de-conflation** - state that the existing `confidence=high/medium` is `f(evidence_strength, assay_context)`, where `assay_context` (patient vs healthy-donor-IVS vs engineered) is formalized later in #823. This task introduces `evidence_strength` only.
  6. **Junction mapping grade ladder** - the four grades with exact definitions (copy from Task 4 Step 1), and the no-inference rule.
  7. **Decoy-negative construction** - the 3-tier hierarchy (copy from Task 5 Step 1), with the hard rule: tiers are never pooled as equal; any score reports which negative tier(s) it used.

- [ ] **Step 2: Self-check the doc against the spec.** Confirm every vocabulary value and rule matches `docs/superpowers/specs/2026-06-26-issue-735-junction-labeling-scheme-design.md`. Fix drift inline.

- [ ] **Step 3: Commit.**

```bash
git add research/experiments/issue_680_splice_immunogenicity_registry/LABELING_SCHEME.md
git commit -m "docs(registry): documented pos/neg labeling scheme for #680 benchmark (#735)"
```

---

## Task 2: Write the registry validator (the failing test)

**Files:**
- Create: `research/experiments/issue_680_splice_immunogenicity_registry/validate_registry.py`

**Interfaces:**
- Consumes: `registry.tsv`; the vocabularies from Task 1.
- Produces: `validate_registry.py` exits 0 on a clean registry, 1 with a printed violation list otherwise. Tasks 3-5 re-run it as their pass condition.

- [ ] **Step 1: Write the validator.** It checks the *target* schema, so it fails today (columns absent).

```python
#!/usr/bin/env python3
"""Validate registry.tsv against the documented labeling scheme (Issue #735).

Run: research/.venv/bin/python validate_registry.py
Exit 0 = clean; exit 1 = violations printed to stderr.
"""
import sys
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.tsv"

GRADES = {"coords", "event-id", "gene-mechanism", "none"}
STRENGTHS = {"strong", "weak", "hard", "soft", "na"}
REQUIRED_NEW_COLS = ["evidence_strength", "label_rationale", "junction_id", "junction_mapping_grade"]

# Effector vs detection-only keyword sets used to cross-check positives.
EFFECTOR = ("ifn", "elispot", "cytotox", "granzyme", "gzmb", "cd107", "cd137",
            "degranul", "killing", "ldh", "caspase", "incucyte", "tcr", "tnf",
            "stabiliz", "activation", "in vivo")
DETECTION = ("tetramer", "dextramer", "multimer")


def violations(df: pd.DataFrame) -> list[str]:
    out = []
    for col in REQUIRED_NEW_COLS:
        if col not in df.columns:
            out.append(f"missing required column: {col}")
    if out:
        return out  # schema not yet present; stop here

    for i, r in df.iterrows():
        rid = f"row {i} ({r['peptide']})"
        if r["junction_mapping_grade"] not in GRADES:
            out.append(f"{rid}: bad junction_mapping_grade {r['junction_mapping_grade']!r}")
        if r["evidence_strength"] not in STRENGTHS:
            out.append(f"{rid}: bad evidence_strength {r['evidence_strength']!r}")
        if not str(r["label_rationale"]).strip():
            out.append(f"{rid}: empty label_rationale")
        # consistency: strong positives must name an effector readout
        ro = str(r["readout"]).lower()
        if r["label"] == "positive" and r["evidence_strength"] == "strong":
            if not any(k in ro for k in EFFECTOR):
                out.append(f"{rid}: strong positive but readout has no effector term")
        if r["label"] == "positive" and r["evidence_strength"] == "weak":
            if any(k in ro for k in EFFECTOR):
                out.append(f"{rid}: weak positive but readout names an effector term")
            if not any(k in ro for k in DETECTION):
                out.append(f"{rid}: weak positive but readout has no detection term")
        # grade 'none' must record a reason in notes
        if r["junction_mapping_grade"] == "none" and not str(r["notes"]).strip():
            out.append(f"{rid}: grade 'none' with no reason in notes")
    return out


def main() -> int:
    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    v = violations(df)
    if v:
        print(f"FAIL: {len(v)} violation(s):", file=sys.stderr)
        for line in v:
            print("  -", line, file=sys.stderr)
        return 1
    print(f"PASS: {len(df)} rows valid against the labeling scheme.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run it; confirm it fails on the missing columns.**

Run: `research/.venv/bin/python research/experiments/issue_680_splice_immunogenicity_registry/validate_registry.py`
Expected: FAIL, lists `missing required column: evidence_strength` (and the other three).

- [ ] **Step 3: Commit.**

```bash
git add research/experiments/issue_680_splice_immunogenicity_registry/validate_registry.py
git commit -m "test(registry): scheme validator (fails pre-annotation) (#735)"
```

---

## Task 3: Add `evidence_strength` + `label_rationale`; re-audit boundary rows

**Files:**
- Modify: `research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv`
- Create (transient, committed): `research/experiments/issue_680_splice_immunogenicity_registry/derive_evidence_strength.py`

**Interfaces:**
- Consumes: the validator from Task 2; the `EFFECTOR`/`DETECTION` rule from Task 1/2.
- Produces: `registry.tsv` with two new columns populated for all 79 rows.

- [ ] **Step 1: Write the derivation script.** It assigns `evidence_strength` mechanically and seeds `label_rationale` per tier; special rows are handled explicitly.

```python
#!/usr/bin/env python3
"""Derive evidence_strength + seed label_rationale (Issue #735, Task 3)."""
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"
EFFECTOR = ("ifn", "elispot", "cytotox", "granzyme", "gzmb", "cd107", "cd137",
            "degranul", "killing", "ldh", "caspase", "incucyte", "tcr", "tnf",
            "stabiliz", "activation", "in vivo")
DETECTION = ("tetramer", "dextramer", "multimer")


def strength(r):
    ro = str(r["readout"]).lower()
    tier = r["tier"]
    if r["label"] == "untested" or tier in ("presentation-prevalence", "negative-control-not-splice"):
        return "na"
    if r["label"] == "negative":
        return "soft" if "healthy-donor ivs" in ro or "ivs" in ro else "hard"
    # positive
    if any(k in ro for k in EFFECTOR):
        return "strong"
    if any(k in ro for k in DETECTION):
        return "weak"
    return "na"  # validator will catch; investigate by hand


RATIONALE = {
    "candidate-negative": "tested-negative in healthy-donor IVS; failed-to-prime != intrinsically non-immunogenic (soft negative).",
    "hard-negative-true-splice": "splice-derived, MS-presented, functional assay performed, no T-cell response (hard negative).",
    "presentation-prevalence": "presented + prevalence/survival evidence but no functional T-cell assay performed; untested.",
    "negative-control-not-splice": "constitutive same-locus control; fails the splice gate; control only.",
    "functional-nonscorable": "functionally validated but exact sequence unpublished; cannot be keyed.",
    "candidate": "positive whose 2nd adversarial-verify pass did not confirm; held below scorable pending re-verify.",
}


def rationale(r):
    if r["tier"] in RATIONALE:
        return RATIONALE[r["tier"]]
    es = r["evidence_strength"]
    if es == "strong":
        return f"effector readout ({r['readout']}) -> strong positive."
    if es == "weak":
        return f"antigen-specific detection only ({r['readout']}) -> weak positive."
    return "REVIEW: rule did not assign a rationale."


df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")
df["evidence_strength"] = df.apply(strength, axis=1)
df["label_rationale"] = df.apply(rationale, axis=1)
df.to_csv(REG, sep="\t", index=False)
print(df["evidence_strength"].value_counts().to_string())
print("rows needing manual review:",
      int((df["label_rationale"].str.startswith("REVIEW")).sum()))
```

- [ ] **Step 2: Run the derivation.**

Run: `research/.venv/bin/python research/experiments/issue_680_splice_immunogenicity_registry/derive_evidence_strength.py`
Expected: a value-count table (≈ `weak` 41, `strong` 21, `soft` 8, `hard` 1, `na` 8); `rows needing manual review: 0`. If any row prints `REVIEW`, inspect its `readout` and fix the rule or the row by hand before proceeding.

- [ ] **Step 3: Re-audit the 14 boundary rows.** For each, confirm the scheme-derived `label`/`evidence_strength` matches the existing `tier`; the expected dispositions are: 2 `candidate` (stay positive, held non-scorable), 8 `candidate-negative` (negative/soft), 1 `presentation-prevalence` (untested/na), 2 `functional-nonscorable` (positive/strong, non-scorable), 1 `negative-control-not-splice` (negative/na). Any genuine change gets an old→new line in a scratch list for Task 6's changelog. Expect zero label changes (the scheme codifies the existing curation); log the confirmation either way.

- [ ] **Step 4: Run the validator; the label/evidence checks now pass** (junction columns still missing → still FAIL on those).

Run: `research/.venv/bin/python .../validate_registry.py`
Expected: FAIL, but only `missing required column: junction_id` / `junction_mapping_grade` remain - no `evidence_strength`/`label_rationale`/consistency violations.

- [ ] **Step 5: Commit.**

```bash
git add research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv research/experiments/issue_680_splice_immunogenicity_registry/derive_evidence_strength.py
git commit -m "data(registry): apply evidence_strength + label_rationale per scheme (#735)"
```

---

## Task 4: Add `junction_id` + `junction_mapping_grade`

**Files:**
- Modify: `research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv`
- Create: `research/experiments/issue_680_splice_immunogenicity_registry/junction_evidence_by_source.md`

**Interfaces:**
- Consumes: `PROVENANCE.md` + each source's published junction evidence; the grade ladder from Task 1.
- Produces: `registry.tsv` with `junction_id` + `junction_mapping_grade` for all 79 rows.

- [ ] **Step 1: Record the grade ladder (authoritative definitions).**

| Grade | Assign when the source published... |
|---|---|
| `coords` | genomic donor/acceptor for that peptide's junction (or a coordinate that maps 1:1 to it) |
| `event-id` | a transcript/splice-event identifier as the unit (e.g. SF3B1 S8 junction, IRIS event ID) but not genomic coords |
| `gene-mechanism` | only gene + canonical mechanism (no junction-level identifier) |
| `none` | no recoverable junction; record the reason in `notes` |

No-inference rule: if a peptide's junction is not stated by the source, it is `gene-mechanism` or `none` - never reconstruct coordinates from the peptide or transcript.

- [ ] **Step 2: Build `junction_evidence_by_source.md`.** One row per source (10 sources: Bigot, Manoharan, Merlotti, SNAF, Kim, IRIS, Fisher, Kwok, Xiong, POSTN). For each, read its `PROVENANCE.md` entry + the registry README notes and classify the junction unit it actually published into the ladder. Known anchors from the README (verify against PROVENANCE, do not assume): Bigot's 5 S8-junction-validated rows = `event-id`, Bigot's 30 tetramer-only = `gene-mechanism` (panel-design origin, no per-peptide junction); IRIS = `event-id` (event IDs); SNAF = whichever unit its supplements give (verify - coords vs event). Each source line cites the PROVENANCE location backing the grade.

- [ ] **Step 3: Apply per-row.** Set `junction_mapping_grade` from the source classification; set `junction_id` to the published coordinate/event identifier where grade is `coords`/`event-id`, else blank. For any `none`, ensure `notes` states why (extend `notes`, do not overwrite existing text). Hand-apply or script from the source table - either is fine; the validator is the gate.

- [ ] **Step 4: Run the validator; expect PASS.**

Run: `research/.venv/bin/python .../validate_registry.py`
Expected: `PASS: 79 rows valid against the labeling scheme.`

- [ ] **Step 5: Produce the per-grade tally** (for Task 6).

Run: `research/.venv/bin/python -c "import pandas as pd; d=pd.read_csv('research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv',sep='\t'); print(d['junction_mapping_grade'].value_counts().to_string())"`
Expected: a distribution dominated by `gene-mechanism`/`event-id` (record the exact numbers).

- [ ] **Step 6: Commit.**

```bash
git add research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv research/experiments/issue_680_splice_immunogenicity_registry/junction_evidence_by_source.md
git commit -m "data(registry): junction annotation + per-source grade ladder, no inference (#735)"
```

---

## Task 5: Materialize the #681 presented-decoy seed (Tier 2)

**Files:**
- Create: `research/experiments/issue_680_splice_immunogenicity_registry/decoy_negatives/presented_decoys_681.tsv`

**Interfaces:**
- Consumes: the 13 MS-presented-untested peptides enumerated in `README.md` §"Deferred - MS-presented / immunogenicity-untested tier".
- Produces: a 13-row TSV usable by the harness (#736) as the Tier-2 negative seed.

- [ ] **Step 1: Create the seed file.** Columns: `peptide	gene	hla	source	decoy_tier	rationale`. Populate the 13 rows from the README enumeration (SNAF Supp Fig 4: `NQDEDPLEV`/C6orf52, `KGPWYPLSL`/C20orf204, `VAPGEAKNL`/RASA3, `YALANIKWI`/DYNLT5, `KEKLDQLVY`/FBXO7, `TELQRTLSL`/NGLY1; Supp Fig 7: `SQTPKSRAL`/PSMF1, `RKLEAPYLL`/MCF2L, `LSWPRSTPM`/CPN1, `VSTGCAVVL`/SERPINE2, `RRLPNPPAV`/RGS12, `IVKRPRSEL`/EXO1, `AVPLLQTNR`/ETV4). Set `decoy_tier=presented-decoy`, `source=SNAF (Li 2024)`, `hla` per the README where given else blank, `rationale="MHC-presented splice peptide, no functional T-cell assay (#681 feeder)."`

- [ ] **Step 2: Write a check that the seed matches the README enumeration.** Add to `validate_registry.py` a second entry point or a small standalone assert: the seed has exactly 13 rows, all `decoy_tier=presented-decoy`, peptides are unique and 8-11 aa.

```python
# decoy seed check (append to validate_registry.py main, gated on file existence)
DECOY = HERE / "decoy_negatives" / "presented_decoys_681.tsv"
if DECOY.exists():
    dd = pd.read_csv(DECOY, sep="\t", dtype=str).fillna("")
    assert len(dd) == 13, f"expected 13 presented decoys, got {len(dd)}"
    assert (dd["decoy_tier"] == "presented-decoy").all(), "bad decoy_tier value"
    assert dd["peptide"].is_unique, "duplicate decoy peptide"
    assert dd["peptide"].str.len().between(8, 11).all(), "decoy peptide length out of range"
    print("PASS: 13 presented decoys valid.")
```

- [ ] **Step 3: Run the check.**

Run: `research/.venv/bin/python .../validate_registry.py`
Expected: `PASS: 79 rows valid...` and `PASS: 13 presented decoys valid.`

- [ ] **Step 4: Commit.**

```bash
git add research/experiments/issue_680_splice_immunogenicity_registry/decoy_negatives/presented_decoys_681.tsv research/experiments/issue_680_splice_immunogenicity_registry/validate_registry.py
git commit -m "data(registry): materialize #681 presented-decoy Tier-2 seed (#735)"
```

---

## Task 6: Update README + PROVENANCE (docs + changelog)

**Files:**
- Modify: `research/experiments/issue_680_splice_immunogenicity_registry/README.md`
- Modify: `research/experiments/issue_680_splice_immunogenicity_registry/PROVENANCE.md`

**Interfaces:**
- Consumes: the per-grade tally (Task 4 Step 5), the re-audit log (Task 3 Step 3).

- [ ] **Step 1: Document the four new columns** in the README schema section, pointing to `LABELING_SCHEME.md` as the authoritative rule set.

- [ ] **Step 2: Add a "Junction-mapping coverage" subsection** with the per-grade tally and the one-line finding: coordinate-level grounding is a minority; most rows rest at event-id/gene-mechanism, mirroring the thin functional base.

- [ ] **Step 3: Add a "Labeling re-audit (#735)" note** recording the boundary-row review outcome (expected: 0 label changes; the scheme codifies existing curation) and any logged old→new changes.

- [ ] **Step 4: Document the decoy tiers** and link `decoy_negatives/presented_decoys_681.tsv` as the materialized Tier-2 seed; note Tier-3 generation is deferred to #736.

- [ ] **Step 5: Final validator run + commit.**

```bash
research/.venv/bin/python research/experiments/issue_680_splice_immunogenicity_registry/validate_registry.py
git add research/experiments/issue_680_splice_immunogenicity_registry/README.md research/experiments/issue_680_splice_immunogenicity_registry/PROVENANCE.md
git commit -m "docs(registry): document new columns, junction coverage, re-audit (#735)"
```

---

## Acceptance-criteria coverage

| AC | Task(s) |
|---|---|
| #1 peptide→junction where annotated; unmappable flagged | 4 (columns + `none` reasons) |
| #2 documented labeling scheme committed | 1 (`LABELING_SCHEME.md`) |
| #3 scheme applied → labeled set + per-peptide rationale | 3 (`evidence_strength` + `label_rationale`) |
| #4 decoy-negative construction approach specified | 1 (§Decoy) + 5 (materialized Tier-2 seed) |

## Self-review

- **Spec coverage:** every spec deliverable maps to a task (scheme→1, validator-as-test→2, evidence/rationale→3, junction→4, decoy seed→5, docs→6). Confidence de-conflation is documented (Task 1 Step 1.5) with `assay_context` correctly deferred to #823. No gap.
- **Placeholder scan:** the only non-pre-filled content is the per-source junction grades (Task 4) - this is source-research output, not a placeholder; the procedure, ladder, and no-inference rule are fully specified and the validator enforces the schema. Evidence_strength, rationale, decoy seed, and validator are concrete code.
- **Type consistency:** `evidence_strength` vocab `{strong,weak,hard,soft,na}` and grade vocab `{coords,event-id,gene-mechanism,none}` are identical in the validator (Task 2), derivation (Task 3), and scheme (Task 1). The `EFFECTOR`/`DETECTION` tuples are identical in Tasks 2 and 3.
