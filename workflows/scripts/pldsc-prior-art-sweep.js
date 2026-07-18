export const meta = {
  name: 'pldsc-prior-art-sweep',
  description: 'Full prior-art landscape sweep: is the PLDSC note\'s niche (decomposition + unification + aliasing corollary) unoccupied?',
  phases: [
    { title: 'Search', detail: '5 parallel cluster searchers, each trying to REFUTE novelty' },
    { title: 'Synthesize', detail: 'per-contribution occupancy verdict' },
  ],
}

const CLAIMS = `
PLDSC theory-note contributions to check for prior art (the note is a
clarifying/unifying methods note, NOT a new estimator):

CLAIM-DECOMP: the genome-wide LD score is EXACTLY the GRM spectral decomposition
  sum_k r_jk^2 = (M/N) sum_p lambda_p corr(x_j,u_p)^2, used as an ORGANIZING
  local(linkage)/non-local(structure) split; the classic "PC score" (eigenvalue-
  weighted squared PC loadings) is its leading-eigenvector truncation.
CLAIM-UNIF: cov-LDSC (Luo 2021), the Cholesky-whitened "WhiLD" score, the PC
  score, and a LOCO mixed-model residual-precision operator
  P_c = V_c^{-1} - V_c^{-1}C(C^T V_c^{-1}C)^+ C^T V_c^{-1} (V_c = sigma_g^2 K + sigma_e^2 I)
  are all instances of removing/measuring the non-local component; the published
  structure-corrected-LD operators are HARD rank-K PC projections (Mangin 2012,
  Bercovich 2025 on pairwise r^2; cov-LDSC on LD scores) and P_c is their
  eigenvalue-weighted GLS GENERALIZATION.
CLAIM-ALIAS: in the marginal E[z_j^2], directional ancestry-aligned GENETIC
  variance and environmental stratification enter as ONE regressor
  corr(x_j,u_p)^2 (exact aliasing); consequently every structure-only LD
  correction (cov-LDSC/P_c/PC-score/Bercovich) conservatively CO-ATTENUATES
  directional ancestry-aligned h^2. (We frame this as an LD-score formalization of
  Vilhjalmsson & Nordborg 2013.)
CLAIM-GLS-OP: P_c as the LD-score whitener MATCHED to a LOCO linear-mixed-model
  score statistic (random-effect-adjusted LD, which cov-LDSC left as future work).
`

phase('Search')

const CLUSTERS = [
  { key: 'ldscore-decomposition', prompt: `Search the literature for work that DECOMPOSES the LD score or views it through the GRM/eigen-spectrum. Target CLAIM-DECOMP and CLAIM-GLS-OP. Look for: "LD score decomposition", spectral/eigen decomposition of LD or the GRM tied to LD scores, local vs long-range/structure LD-score partition, GRM-spectrum views of heritability/LD, the "PC score" or eigenvalue-weighted PC-loading scores. Is the exact spectral identity as an ORGANIZING device already published? Is a mixed-model/GLS-whitened LD score matched to an LMM score statistic published?` },
  { key: 'structure-corrected-LD', prompt: `Search for structure/relatedness-corrected LD and heritability methods. Target CLAIM-UNIF. Cover: cov-LDSC (Luo 2021), Bercovich 2025, Mangin 2012, Price 2008 long-range LD, admixture LD, any 2022-2026 structure-corrected LD-SCORE method, and whether anyone has generalized the hard PC-projection correction to an eigenvalue-weighted / GLS / mixed-model operator on LD scores. Find the closest competitor to P_c as an LD-score operator.` },
  { key: 'enrichment-competitors', prompt: `Map the active heritability-partitioning / enrichment landscape (context, not direct claim). Cover S-LDSC + baseline-LD (Finucane 2015, Gazal 2017), graphREML (2024), FGLS functional enrichment (2024), LDSC++ (2025), and any PC-loading-as-annotation work. For each, note whether it touches the local/non-local decomposition or the structure-aliasing conservatism. Verify DOIs/identifiers.` },
  { key: 'ldsc-variants-gls', prompt: `Search for LDSC variants and GLS/likelihood/eigenvalue reformulations. Target CLAIM-UNIF and CLAIM-GLS-OP. Cover HEELS, LDER, RHE-mc, GRE, any "generalized least squares LD score regression", whitened/precision-weighted LD-score regression, and unifying-perspective/review papers on LDSC vs REML. Does any of them frame LDSC as a GLS operation with a V_c-precision metric, or unify cov-LDSC/WhiLD/PC-score under one decomposition?` },
  { key: 'aliasing-conservatism', prompt: `Target CLAIM-ALIAS. Search the polygenic-signal-vs-stratification literature (Vilhjalmsson & Nordborg 2013; Berg & Coop 2014; Berg 2019; Sohail 2019; Zaidi & Mathieson 2020; Young 2018 RDR; Kong 2018) AND, crucially, whether anyone has EXPLICITLY stated the LD-score-specific corollary: that cov-LDSC / PC-correction / PC-annotation conservatively ATTENUATE genuine ancestry-aligned heritability along with confounding (i.e. de-confounding has a real-signal cost, stated for the LD-score family). This corollary is our most defensible novel slice -- try hard to find it already stated.` },
]

const found = await parallel(CLUSTERS.map(c => () =>
  agent(`You are a prior-art searcher. Use live web search and the bio-research literature tools (pubmed, biorxiv, consensus) -- do NOT rely on memory; actually search and cite DOIs/URLs. Your job is ADVERSARIAL: try to find published work that ALREADY OCCUPIES the PLDSC claim(s) below, not to confirm novelty.\n\nYour cluster: ${c.prompt}\n\n${CLAIMS}\n\nReturn the papers you found (with identifiers) and, for each PLDSC claim your cluster touches, whether it is OCCUPIED (someone published essentially this), PARTIAL (adjacent/weaker), or UNOCCUPIED, with the specific evidence.`,
    { label: `search:${c.key}`, phase: 'Search',
      schema: {
        type: 'object', additionalProperties: false,
        required: ['cluster', 'papers', 'claim_assessments', 'closest_competitor'],
        properties: {
          cluster: { type: 'string' },
          papers: { type: 'array', items: { type: 'object', additionalProperties: false,
            required: ['title', 'id', 'relevance'],
            properties: { title: { type: 'string' }, id: { type: 'string', description: 'DOI/PMID/arXiv/URL' }, relevance: { type: 'string' } } } },
          claim_assessments: { type: 'array', items: { type: 'object', additionalProperties: false,
            required: ['claim', 'status', 'evidence'],
            properties: { claim: { type: 'string', enum: ['CLAIM-DECOMP', 'CLAIM-UNIF', 'CLAIM-ALIAS', 'CLAIM-GLS-OP'] },
              status: { type: 'string', enum: ['occupied', 'partial', 'unoccupied'] },
              evidence: { type: 'string' } } } },
          closest_competitor: { type: 'string' },
        } } })
))

phase('Synthesize')

const synthesis = await agent(
  `You are synthesizing a prior-art landscape sweep for the PLDSC theory note. Below are ${CLUSTERS.length} cluster-searcher reports (JSON). For EACH of the four claims (CLAIM-DECOMP, CLAIM-UNIF, CLAIM-ALIAS, CLAIM-GLS-OP) render a final occupancy verdict (occupied / partial / unoccupied) with the single closest competitor and its identifier. Then give an OVERALL verdict: is there a defensible unoccupied niche for a clarifying/unifying note, and if so what is the sharpest honest positioning? Be skeptical and honest -- if the niche is thin, say so.\n\n${CLAIMS}\n\nCLUSTER REPORTS:\n${JSON.stringify(found.filter(Boolean), null, 1)}`,
  { label: 'synthesize:verdict', phase: 'Synthesize',
    schema: {
      type: 'object', additionalProperties: false,
      required: ['per_claim', 'overall_verdict', 'sharpest_positioning', 'must_cite'],
      properties: {
        per_claim: { type: 'array', items: { type: 'object', additionalProperties: false,
          required: ['claim', 'verdict', 'closest_competitor', 'note'],
          properties: { claim: { type: 'string' }, verdict: { type: 'string', enum: ['occupied', 'partial', 'unoccupied'] },
            closest_competitor: { type: 'string' }, note: { type: 'string' } } } },
        overall_verdict: { type: 'string', enum: ['clear-niche', 'thin-but-defensible', 'largely-occupied'] },
        sharpest_positioning: { type: 'string' },
        must_cite: { type: 'array', items: { type: 'string' } },
      } } })

return { found: found.filter(Boolean), synthesis }
