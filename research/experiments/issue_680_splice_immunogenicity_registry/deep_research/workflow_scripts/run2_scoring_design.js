export const meta = {
  name: 'splice-scoring-harness-design',
  description: 'Generate + adversarially judge candidate scoring methodologies for the sparse #680 splice-immunogenicity registry (#736 keystone); synthesize a recommended design + claims-it-supports/forbids',
  phases: [
    { title: 'Frame', detail: 'compute live registry constraints (counts by tier/label/allele/negative-tier)' },
    { title: 'Design', detail: '5 independent candidate scoring methodologies, distinct methodological stances' },
    { title: 'Judge', detail: 'adversarial multi-lens critique of each design (small-n validity, reviewer-proofness, implementability, honesty)' },
    { title: 'Synthesize', detail: 'winning design + grafted best-of-runners-up + claims supported/forbidden + open questions for morning review' },
  ],
}

const CONTEXT = `
PROJECT CONTEXT - the #680 splice-immunogenicity registry and the #736 scoring harness it feeds.
The registry (research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv) is the first field-wide
collection of splice-junction-derived neoepitopes with MEASURED T-cell immunogenicity. It exists to ground-truth a
splice-neoantigen predictor. #736 is the scoring harness that will score a predictor against it.

THE DOCUMENTED CONSTRAINTS (from our own #737 sparsity analysis - these are the hard reality any scoring design must confront):
- The POSITIVE base is small (~81 functional-scorable) and an HLA-A*02:01 MONOCULTURE (~89%). Any allele-stratified
  result collapses to an A2 result; effective allele diversity ~1.26.
- Source clustering: ~10 studies but effective ~3.8 independent (top study ~43%). Not i.i.d.
- The NEGATIVE set is the binding constraint: effectively ONE hard true-negative field-wide (VELEDHVML, and it is A*11:01,
  not even the dominant allele) + 8 soft (failed-to-prime IVS) negatives. A powered AUC needs ~19-31 negatives. Specificity
  is, today, essentially unmeasurable.
- The negatives come in THREE NON-POOLABLE TIERS (LABELING_SCHEME.md section 7): Tier 1 experimental true-negatives
  (hard n=1 + soft n=8), Tier 2 presented-but-untested decoys (#681, ~13 seed), Tier 3 matched-synthetic shuffled/decoy
  (abundant, weakest claim). HARD RULE: tiers must NEVER be pooled into one undifferentiated negative set; any metric must
  report which tier(s) it used. A result on Tier 1 only vs Tier 1+2 vs Tier 1+2+3 are different, non-interchangeable claims.
- The scorable positive set is 'label==positive AND tier==functional-scorable' (NOT label alone - candidate and
  functional-nonscorable rows are excluded). evidence_strength splits positives into 'strong' (effector readout) vs
  'weak' (tetramer/dextramer detection only).
- TESLA (the SNV/indel gold standard) explicitly EXCLUDED splice isoforms, so there is no prior art to copy - this design
  is establishing the first open splice-immunogenicity scoring methodology.
`

phase('Frame')

const FRAME_SCHEMA = {
  type: 'object',
  required: ['counts_summary', 'binding_constraints', 'available_signal'],
  properties: {
    counts_summary: { type: 'string', description: 'live registry counts: total rows, scorable positives, by allele, by evidence_strength, hard vs soft negatives, by source' },
    binding_constraints: { type: 'array', items: { type: 'string' }, description: 'the statistical facts that constrain any scoring design' },
    available_signal: { type: 'string', description: 'what CAN legitimately be measured given the data, and what cannot' },
  },
}

const frame = await agent(
  `You are framing the data reality for a scoring-methodology design exercise. Compute the LIVE constraints from the actual registry.\n` +
  `${CONTEXT}\n` +
  `Tasks (use Bash + read the files; do NOT rely on remembered numbers):\n` +
  `1. Load research/experiments/issue_680_splice_immunogenicity_registry/registry.tsv. Compute: total rows; scorable positives (label==positive AND tier==functional-scorable); positives split by evidence_strength (strong vs weak); breakdown by HLA allele (count A*02:01 vs each other); count of hard negatives (tier==hard-negative-true-splice), soft negatives (tier==candidate-negative), and any presentation-prevalence / non-splice-control rows; breakdown by source.\n` +
  `2. Read LABELING_SCHEME.md section 7 (the three negative tiers) and decoy_negatives/presented_decoys_681.tsv (count Tier-2 seed rows).\n` +
  `3. NOTE: an enrichment deep-research (run #1) may have just produced confirmed new rows not yet folded into registry.tsv. If a file research/experiments/issue_680_splice_immunogenicity_registry/run1_confirmed_adds* exists, fold its counts in; otherwise score against the current 96-row registry and say so.\n` +
  `Return the live counts, the binding statistical constraints, and an honest statement of what can vs cannot be measured.`,
  { label: 'frame:live-constraints', phase: 'Frame', schema: FRAME_SCHEMA, effort: 'high' }
)

log(`Frame: ${frame?.counts_summary?.slice(0, 200) || 'no frame'}`)

phase('Design')

const DESIGN_SCHEMA = {
  type: 'object',
  required: ['stance', 'one_line', 'metrics', 'negative_handling', 'allele_handling', 'uncertainty_handling', 'claims_supported', 'claims_forbidden', 'implementation_sketch', 'weaknesses'],
  properties: {
    stance: { type: 'string' },
    one_line: { type: 'string' },
    metrics: { type: 'array', items: { type: 'string' }, description: 'the concrete metrics this design reports' },
    negative_handling: { type: 'string', description: 'how it uses the 3 negative tiers without pooling them' },
    allele_handling: { type: 'string', description: 'how it deals with the A*02:01 monoculture' },
    uncertainty_handling: { type: 'string', description: 'CIs, bootstrap, permutation, effect-size - how small-n uncertainty is quantified' },
    claims_supported: { type: 'array', items: { type: 'string' } },
    claims_forbidden: { type: 'array', items: { type: 'string' }, description: 'what this design explicitly REFUSES to claim' },
    implementation_sketch: { type: 'string', description: 'how it runs against the registry.tsv schema concretely' },
    weaknesses: { type: 'array', items: { type: 'string' } },
  },
}

const STANCES = [
  { key: 'ranking-only', prompt: `STANCE A - RANKING-ONLY. Accept that absolute calibrated probabilities are impossible at this n; design a purely rank-based evaluation. A predictor is scored by how well it ranks the scorable positives ABOVE the negative tiers. Use AUC-ROC and AUPRC computed SEPARATELY per negative tier (never pooled), plus top-k recall / enrichment at the top of the ranked list. Treat the abundant Tier-3 synthetic negatives as the ranking-denominator floor. Be explicit that this measures discrimination/ranking, not specificity or probability.` },
  { key: 'tiered-stratified', prompt: `STANCE B - TIERED + STRATIFIED REPORTING. The headline IS the stratification. Report every metric as a small table: (negative-tier) x (allele stratum: A2-only vs all) x (evidence_strength: strong-only vs all positives). No single number is ever the result; the result is the table + an explicit statement that the A2-stratified cell is the only well-powered one. Bakes the LABELING_SCHEME non-pooling rule into the output format itself.` },
  { key: 'uncertainty-first', prompt: `STANCE C - UNCERTAINTY-FIRST. Given tiny n, point estimates are dishonest. Lead with uncertainty: bootstrap confidence intervals on every metric, permutation/randomization tests for whether ranking beats chance, effect sizes with intervals rather than p-values. The deliverable for any metric is '[estimate, 95% CI]' and a powered/underpowered verdict. Explicitly compute and report the achievable power given the current negative count (cf. the ~19-31 needed).` },
  { key: 'negative-sensitivity', prompt: `STANCE D - NEGATIVE-SAMPLING SENSITIVITY AS THE HEADLINE (cf. field_gap_map section 3.6). Reported AUCs swing wildly with how negatives are constructed. So make the SWING the result: re-score the same positives under multiple standardized negative-sampling schemes (Tier-1 only; Tier-1+2; Tier-1+2+3; within-allele shuffle; background-repertoire decoys; hard-negative-only) and report the AUC range/variance. The claim becomes 'discrimination is X under scheme 1...N, swing = delta', which is more honest and more reviewer-proof than any single number.` },
  { key: 'positive-only-enrichment', prompt: `STANCE E - POSITIVE-ONLY / PRESENTATION-ANCHORED. With negatives nearly absent, lean on the positives + the presentation background. Score whether a predictor enriches the validated positives within the full presented-peptide background (the Tier-2/Tier-3 presented set as a realistic decoy universe), using precision@k and enrichment-factor metrics that need only positives + a presented background, not labeled true-negatives. Frame specificity as a known gap with a concrete data-acquisition ask (#911), not a computed number.` },
]

const designs = await parallel(
  STANCES.map((s) => () =>
    agent(
      `You are designing ONE candidate scoring methodology for the #736 harness. Commit fully to your assigned stance - produce the strongest version of it.\n` +
      `${CONTEXT}\n` +
      `LIVE DATA CONSTRAINTS (from the Frame phase):\n${JSON.stringify(frame)}\n\n` +
      `${s.prompt}\n\n` +
      `Be concrete and implementable against the registry.tsv columns (peptide, gene, hla, length, splice_mechanism, source, readout, label, tier, evidence_strength, assay_context, ...). State exactly which metrics, how negatives are used WITHOUT pooling tiers, how the A2 monoculture is handled, how uncertainty is quantified, and - critically - the explicit list of claims your design SUPPORTS and the claims it FORBIDS. Name your own weaknesses honestly.`,
      { label: `design:${s.key}`, phase: 'Design', schema: DESIGN_SCHEMA, effort: 'high' }
    ).then((d) => ({ ...d, key: s.key }))
  )
)

const validDesigns = designs.filter(Boolean)
log(`Design: ${validDesigns.length} candidate methodologies generated.`)

phase('Judge')

const JUDGE_LENSES = ['statistical-validity', 'reviewer-proofness', 'implementability', 'honesty']
const VERDICT_SCHEMA = {
  type: 'object',
  required: ['lens', 'score', 'rationale', 'fatal_flaws', 'salvageable_ideas'],
  properties: {
    lens: { type: 'string' },
    score: { type: 'number', description: '0-10' },
    rationale: { type: 'string' },
    fatal_flaws: { type: 'array', items: { type: 'string' } },
    salvageable_ideas: { type: 'array', items: { type: 'string' }, description: 'ideas worth grafting into the winning design even if this design loses' },
  },
}

// pipeline: each design is judged by all 4 lenses as soon as it is designed
const judged = await pipeline(
  validDesigns,
  (d) =>
    parallel(
      JUDGE_LENSES.map((lens) => () =>
        agent(
          `You are an ADVERSARIAL judge evaluating a candidate #736 scoring methodology through the ${lens.toUpperCase()} lens ONLY. Be skeptical - your job is to find what breaks.\n` +
          `${CONTEXT}\n` +
          `CANDIDATE DESIGN (JSON):\n${JSON.stringify(d)}\n\n` +
          `For the ${lens} lens specifically:\n` +
          (lens === 'statistical-validity' ? `Is this design statistically sound at n≈1 hard negative, ~81 clustered positives, 89% A2? Does it avoid claiming power it doesn't have? Are the metrics appropriate for the sample sizes? Would a statistician reject it?` :
           lens === 'reviewer-proofness' ? `Would a hostile peer reviewer of a methods paper accept the resulting claims? Where would they attack? Does it overclaim?` :
           lens === 'implementability' ? `Can this actually be coded against registry.tsv + the Tier-2 seed + a synthetic Tier-3 generator, without data we don't have? Is anything hand-wavy?` :
           `Is the design HONEST about what it can and cannot claim - does it forbid the right claims (no calibrated probability, no specificity number from n=1), or does it quietly smuggle in an unsupported claim?`) +
          `\nScore 0-10, give rationale, name fatal flaws, and list any salvageable ideas worth keeping even if this design loses.`,
          { label: `judge:${d.key}:${lens}`, phase: 'Judge', schema: VERDICT_SCHEMA }
        )
      )
    ).then((verdicts) => ({ design: d, verdicts: verdicts.filter(Boolean) }))
)

const scored = judged.filter(Boolean).map((j) => {
  const vs = j.verdicts
  const avg = vs.length ? vs.reduce((a, v) => a + (v.score || 0), 0) / vs.length : 0
  return { key: j.design.key, design: j.design, verdicts: vs, avgScore: avg }
}).sort((a, b) => b.avgScore - a.avgScore)
log(`Judge: ranked ${scored.length} designs. Top: ${scored[0]?.key} (${scored[0]?.avgScore?.toFixed(1)}).`)

phase('Synthesize')

const FINAL_SCHEMA = {
  type: 'object',
  required: ['executive_summary', 'recommended_design', 'metric_list', 'claims_supported', 'claims_forbidden', 'grafted_ideas', 'open_questions', 'implementation_next_steps'],
  properties: {
    executive_summary: { type: 'string', description: 'the recommended scoring methodology in a paragraph + why it won' },
    recommended_design: { type: 'string', description: 'the full recommended #736 methodology, synthesized from the winner + grafts' },
    metric_list: { type: 'array', items: { type: 'string' }, description: 'the concrete metrics to compute, each with its negative tier + allele stratum' },
    claims_supported: { type: 'array', items: { type: 'string' } },
    claims_forbidden: { type: 'array', items: { type: 'string' } },
    grafted_ideas: { type: 'array', items: { type: 'string' }, description: 'best ideas taken from non-winning designs' },
    open_questions: { type: 'array', items: { type: 'string' }, description: 'decisions left for Jin-Ho in the morning' },
    implementation_next_steps: { type: 'array', items: { type: 'string' }, description: 'concrete #736 build steps' },
  },
}

const final = await agent(
  `You are writing the morning-review brief: a recommended scoring methodology for the #736 harness, synthesized from a judged panel of 5 candidate designs.\n` +
  `${CONTEXT}\n` +
  `LIVE CONSTRAINTS:\n${JSON.stringify(frame)}\n\n` +
  `RANKED + JUDGED DESIGNS (JSON, best first):\n${JSON.stringify(scored)}\n\n` +
  `Synthesize the RECOMMENDED design: start from the top-ranked design, GRAFT the salvageable_ideas the judges flagged from the runners-up, and drop every fatal flaw the judges found. The result must obey the LABELING_SCHEME non-pooling rule and be honest about the n≈1-negative / A2-monoculture reality.\n` +
  `Produce: an executive summary; the full recommended methodology; a concrete metric list (each metric tagged with its negative tier + allele stratum); the explicit claims-it-supports and claims-it-forbids; the grafted ideas; the open questions for Jin-Ho to decide in the morning; and concrete next implementation steps for #736. Write it so Jin-Ho can read it over coffee and make the convergent call.`,
  { label: 'synthesize:morning-brief', phase: 'Synthesize', schema: FINAL_SCHEMA, effort: 'high' }
)

return { frame, ranked: scored.map((s) => ({ key: s.key, score: s.avgScore })), brief: final }
