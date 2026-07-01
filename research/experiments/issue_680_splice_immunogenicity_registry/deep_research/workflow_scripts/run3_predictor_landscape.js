export const meta = {
  name: 'issue736-predictor-landscape',
  description: 'Map the immunogenicity/presentation-predictor landscape as baseline candidates + validation-methodology priors for the #736 scoring harness',
  phases: [
    { title: 'Catalog', detail: '7 predictor families swept in parallel' },
    { title: 'Deep-dive', detail: 'per-predictor extraction + adversarial capability verify' },
    { title: 'Methodology', detail: 'field validation practice: metrics, negatives, leakage' },
    { title: 'Synthesize', detail: 'build-oriented brief for #736' },
  ],
}

// ---- Shared context every agent gets ----
const CTX = `
CONTEXT — you are helping the Scientist of the splice-neoepitope-pipeline project.
We are designing Issue #736: an OPEN benchmark that scores splice-junction-derived neoepitopes
against MEASURED T-cell immunogenicity. The registry is tiny (81 functional-scorable positives,
89% HLA-A*02:01, n=1 hard true-negative) so #736 will ship as a NO-SCALAR RANKING / ENRICHMENT
benchmark on A*02:01, NOT a powered specificity benchmark. The pipeline ALREADY emits a
'calibrated_immunogenicity_log_odds' score (post-MHCflurry Class1PresentationPredictor, pre-TCRdock)
and uses MHCflurry 2.x Class1PresentationPredictor (presentation, not affinity) as its primary ranker;
TCRdock runs as optional structural validation.
Our epitopes are NON-CANONICAL: 9-11mers from intron-retention / poison-exon / novel-junction /
frameshift splice events. So "does this predictor admit non-canonical / frameshift / arbitrary peptides,
or only proteome-annotated ones?" is a FIRST-CLASS question.
Use WebSearch + WebFetch (load them via ToolSearch first). Cite a URL/DOI for every capability claim.
Do NOT invent scores, licenses, or availability — if unverified, say "unverified".
`.trim()

const CATALOG_SCHEMA = {
  type: 'object',
  properties: {
    predictors: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: { type: 'string' },
          family: { type: 'string' },
          one_liner: { type: 'string' },
          paper_or_url: { type: 'string' },
          score_type: { type: 'string', description: 'binding-affinity | eluted-ligand/presentation | immunogenicity | TCR-recognition | processing/cleavage | other' },
          admits_noncanonical: { type: 'string', enum: ['yes', 'no', 'unknown'] },
          why_relevant_to_736: { type: 'string' },
        },
        required: ['name', 'family', 'one_liner', 'score_type', 'admits_noncanonical'],
      },
    },
  },
  required: ['predictors'],
}

const DOSSIER_SCHEMA = {
  type: 'object',
  properties: {
    name: { type: 'string' },
    inputs: { type: 'string', description: 'exact required inputs: peptide, MHC allele, expression, flanks, etc.' },
    output_score_semantics: { type: 'string', description: 'what the number means; %rank vs raw vs probability; direction' },
    admits_splice_frameshift: { type: 'string', enum: ['yes-arbitrary-peptide', 'yes-with-caveat', 'no-proteome-locked', 'unknown'] },
    admits_evidence: { type: 'string', description: 'the concrete basis for the admits verdict (API/docs/paper quote) + URL' },
    availability: { type: 'string', description: 'open-source (lang/repo) | web-server-only | license-gated | commercial | unavailable' },
    availability_url: { type: 'string' },
    validation_summary: { type: 'string', description: 'how the authors validated it: dataset, metric, held-out design' },
    baseline_fit_for_736: { type: 'string', enum: ['drop-in', 'adaptable', 'reference-only', 'unfit'] },
    integration_notes: { type: 'string', description: 'concrete notes for wiring it into our A*02:01 ranking harness next to calibrated_immunogenicity_log_odds' },
  },
  required: ['name', 'admits_splice_frameshift', 'availability', 'baseline_fit_for_736', 'output_score_semantics'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    admits_verdict: { type: 'string', enum: ['CONFIRMED', 'REFUTED', 'UNCERTAIN'] },
    availability_verdict: { type: 'string', enum: ['CONFIRMED', 'REFUTED', 'UNCERTAIN'] },
    corrections: { type: 'string', description: 'any corrected value + the source that corrects it; empty if none' },
    residual_risk: { type: 'string', description: 'what a builder should still check before trusting this row' },
  },
  required: ['admits_verdict', 'availability_verdict'],
}

const METHOD_SCHEMA = {
  type: 'object',
  properties: {
    lens: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          point: { type: 'string' },
          implication_for_736: { type: 'string' },
          source: { type: 'string' },
        },
        required: ['point', 'implication_for_736'],
      },
    },
  },
  required: ['lens', 'findings'],
}

// ---- Phase 1: Catalog (parallel families) ----
phase('Catalog')
const FAMILIES = [
  { key: 'presentation-EL', prompt: `List the major PRESENTATION / eluted-ligand MHC-I predictors (e.g. NetMHCpan-4.1 EL, MHCflurry 2.0 presentation, MixMHCpred 2.x/3.x, HLAthena, BigMHC-EL). For each give the catalog fields. Emphasize whether each accepts an ARBITRARY peptide+allele pair (non-canonical friendly) vs. only annotated proteome.` },
  { key: 'immunogenicity', prompt: `List predictors that specifically target T-cell IMMUNOGENICITY (not just binding/presentation): PRIME / PRIME2.0, DeepImmuno, IEDB "Immunogenicity" (Calis 2013), Repitope, iPred, Gao/DeepHLApan, DeepNeo, any 2023-2026 immunogenicity models. Catalog fields for each.` },
  { key: 'deep-foundation', prompt: `List deep-learning / foundation-model MHC-I and pMHC models usable for scoring peptide immunogenicity or presentation: BigMHC (IM + EL), TransPHLA, MHCnuggets, ImmunoBERT, ESM/protein-LM-based scorers, capsule/transformer pMHC models. Catalog fields; flag which take arbitrary peptides.` },
  { key: 'tcr-recognition', prompt: `List pMHC-TCR RECOGNITION predictors relevant to whether a presented splice peptide is actually T-cell-recognized: NetTCR-2.x, ERGO-II, pMTnet, PISTE, epiTCR, TCRdock, AlphaFold-based pMHC-TCR. We already run TCRdock as optional structural validation — note relevance/overlap. Catalog fields.` },
  { key: 'noncanonical-splice', prompt: `Find scoring approaches built for NON-CANONICAL neoantigens specifically — splice / frameshift / fusion / intron-retention derived: SNAF's own immunogenicity scoring, IRIS scoring, frameshift-neoantigen immunogenicity work, ScanNeo/NeoFuse scoring layers, any 2023-2026 splice-neoantigen ranking method. This family is the MOST relevant to us — be thorough. Catalog fields.` },
  { key: 'processing-upstream', prompt: `List antigen-PROCESSING / proteasomal-cleavage / TAP predictors that could add an orthogonal feature to a splice-neoepitope score: NetChop, NetCTL/NetCTLpan, MHCflurry processing predictor, proteasomal cleavage models. Catalog fields; note whether each is combinable with a presentation score.` },
  { key: 'recent-sota-benchmarks', prompt: `Sweep 2024-2026 for (a) NEW immunogenicity/presentation predictors not covered by the classic list, and (b) BENCHMARK / review papers that compare neoantigen immunogenicity predictors head-to-head (TESLA and successors, immunogenicity-prediction reviews). Catalog the predictors; for benchmark papers put the paper in paper_or_url and use family="benchmark-paper".` },
]
const catalogResults = await parallel(
  FAMILIES.map(f => () => agent(`${CTX}\n\nTASK (${f.key}): ${f.prompt}`, {
    label: `catalog:${f.key}`, phase: 'Catalog', schema: CATALOG_SCHEMA, effort: 'low',
  }))
)

// dedupe by normalized name (barrier: needs all families before deep-dive)
const norm = s => (s || '').toLowerCase().replace(/[^a-z0-9]/g, '')
const seen = new Map()
const benchmarkPapers = []
for (const r of catalogResults.filter(Boolean)) {
  for (const p of (r.predictors || [])) {
    if ((p.family || '').includes('benchmark')) { benchmarkPapers.push(p); continue }
    const k = norm(p.name)
    if (!k) continue
    if (!seen.has(k)) seen.set(k, p)
  }
}
let predictors = [...seen.values()]
// prioritise splice/immunogenicity-relevant + arbitrary-peptide-friendly; cap with a logged cut
const rank = p => (p.admits_noncanonical === 'yes' ? 2 : p.admits_noncanonical === 'unknown' ? 1 : 0)
  + (/immunogen|splice|frameshift|presentation|eluted/i.test(p.score_type || '') ? 1 : 0)
predictors.sort((a, b) => rank(b) - rank(a))
const CAP = 26
if (predictors.length > CAP) {
  log(`Catalog surfaced ${predictors.length} distinct predictors; deep-diving top ${CAP} by 736-relevance, dropping ${predictors.length - CAP}: ${predictors.slice(CAP).map(p => p.name).join(', ')}`)
  predictors = predictors.slice(0, CAP)
} else {
  log(`Catalog surfaced ${predictors.length} distinct predictors + ${benchmarkPapers.length} benchmark papers; deep-diving all.`)
}

// ---- Phase 2+3: deep-dive pipeline (extract -> adversarial verify) AND methodology, concurrently ----
const METHOD_LENSES = [
  { key: 'metrics', prompt: `How is neoantigen/T-cell IMMUNOGENICITY prediction benchmarked METRICALLY? Compare AUROC vs AUPRC vs positive-rank / top-k enrichment / ECDF under heavy class imbalance and tiny positive sets. What metric is honest when you have ~80 positives, ~1 hard negative, and one dominant allele? Give implications for #736.` },
  { key: 'negative-sets', prompt: `How do immunogenicity benchmarks construct NEGATIVE sets? Random-peptide decoys, DUD-E-style property-matched decoys, MS-presented-but-non-immunogenic ligands, allele-matched vs mismatched. What are the documented biases of each? This is our core #681/#736 problem (we have n=1 hard negative). Give concrete negative-construction recommendations for a splice-neoepitope ranking benchmark.` },
  { key: 'leakage', prompt: `What EVALUATION PITFALLS / DATA-LEAKAGE modes bite immunogenicity-predictor benchmarks? Allele leakage, peptide-similarity/homology leakage between train and test, source-study leakage, publication bias toward A*02:01. What held-out designs mitigate them (LOSO/leave-one-study-out, cluster-restricted splits, allele-holdout)? Map each to our #736 LOSO plan.` },
  { key: 'tesla-prior', prompt: `Summarize what the TESLA consortium benchmark (and any successor neoantigen-immunogenicity benchmarks) actually did, measured, and concluded — especially their handling of presentation vs immunogenicity, their negative sets, and their explicit EXCLUSION of splice isoforms. What can #736 borrow, and where must it diverge because it is splice-specific and tiny?` },
]

const [dossiers, methodology] = await Promise.all([
  pipeline(
    predictors,
    (p) => agent(`${CTX}\n\nProduce a structured dossier for predictor "${p.name}" (${p.one_liner || p.family}). Verify inputs, exact score semantics + direction, whether it admits arbitrary/non-canonical peptides (splice/frameshift), real availability + license, and how the authors validated it. Judge its fit as a BASELINE next to our calibrated_immunogenicity_log_odds on an A*02:01 ranking harness. Reference: ${p.paper_or_url || '(find it)'}.`, {
      label: `dive:${p.name}`, phase: 'Deep-dive', schema: DOSSIER_SCHEMA, effort: 'medium',
    }),
    (dossier, p) => agent(`${CTX}\n\nADVERSARIALLY VERIFY the two error-prone capability claims in this dossier for "${p.name}". These are "what can my tooling do" facts — check them against the primary source (repo/docs/paper), do not trust the dossier.\n(1) admits_splice_frameshift = "${dossier?.admits_splice_frameshift}" — can it REALLY score an arbitrary non-proteome peptide? Try to REFUTE it.\n(2) availability = "${dossier?.availability}" (${dossier?.availability_url || 'no url'}) — is it actually obtainable/runnable as claimed? Try to REFUTE it.\nDefault to UNCERTAIN if you cannot confirm from a primary source.\nDossier: ${JSON.stringify(dossier)}`, {
      label: `verify:${p.name}`, phase: 'Deep-dive', schema: VERDICT_SCHEMA, effort: 'high',
    }).then(v => ({ ...dossier, _verdict: v })).catch(() => ({ ...dossier, _verdict: null }))
  ),
  parallel(METHOD_LENSES.map(m => () => agent(`${CTX}\n\nMETHODOLOGY LENS (${m.key}): ${m.prompt}`, {
    label: `method:${m.key}`, phase: 'Methodology', schema: METHOD_SCHEMA, effort: 'medium',
  }))),
])

const cleanDossiers = dossiers.filter(Boolean)
const cleanMethod = methodology.filter(Boolean)

// ---- Phase 4: Synthesize ----
phase('Synthesize')
const brief = await agent(`${CTX}

You are writing the DECISION BRIEF that unblocks the #736 scoring-harness build. You have:
- ${cleanDossiers.length} verified predictor dossiers: ${JSON.stringify(cleanDossiers)}
- ${cleanMethod.length} methodology lenses: ${JSON.stringify(cleanMethod)}
- benchmark papers surfaced: ${JSON.stringify(benchmarkPapers)}

Write a tight, build-oriented markdown brief with these sections:
1. **Recommended baseline set** — a RANKED shortlist of predictors to score against calibrated_immunogenicity_log_odds in the A*02:01 ranking harness. Only include predictors whose admits_splice_frameshift survived verification as yes/yes-with-caveat AND that are actually obtainable. For each: score semantics, how to wire it, and the one caveat that matters. Explicitly list the predictors REJECTED as baselines and why (proteome-locked, unavailable, refuted claim).
2. **Splice-admissibility table** — every predictor, its verified admits verdict, availability verdict, and residual risk. This is the reusable artifact.
3. **Methodology implications for #736** — what the field's benchmarking practice says our harness MUST do (metric choice under n≈80/1, negative-set construction, leakage-safe held-out design) and where TESLA's precedent applies vs. must diverge (splice-specific, tiny, mono-allele).
4. **Open questions / human-decision points** — concrete forks a human must resolve before build (e.g. which predictors are worth the integration cost, licensing blockers, whether any needs GPU we no longer have).
5. **What NOT to do** — honest anti-recommendations grounded in the findings.
Forbid vague claims — every capability statement must trace to a dossier verdict. Be honest about what is unverified.`, {
  label: 'synthesize:brief', phase: 'Synthesize', effort: 'high',
})

return {
  brief,
  n_predictors_deepdived: cleanDossiers.length,
  n_benchmark_papers: benchmarkPapers.length,
  dossiers: cleanDossiers,
  methodology: cleanMethod,
  benchmarkPapers,
}
