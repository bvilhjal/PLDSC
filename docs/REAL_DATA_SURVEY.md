# Real-data survey: scoping

Date: 2026-07-18

The decomposition and diagnostics are validated in simulation
(`experiments/`). Their *value* — the "resulting knowledge" — can only be settled
on real data: how much of real LD / real GWAS inflation is structure, and does
correcting it change any published conclusion. This document scopes that study.

**Status (2026-07-19): question 1 now has a chr22 pilot and a 3-chromosome
benchmark suite** — `experiments/realdata_1000g_structure.py` +
`experiments/realdata_1000g_benchmarks.py` (1000G phase 3, EUR / AMR / AFR /
full panel). Headline: the framework reproduces (structure fraction tracks λ1;
the non-local remainder tracks the PC score at corr 0.9–1.0), and the LOCO-clean
structure share of the genome-wide LD score is ~90% in the multi-ancestry panel
but ~0–5% in single-ancestry cohorts (~17% vs ~0% in ±250 kb windows). A key
methodological finding: the operator matters — an in-chromosome GRM over-whitens
(removes local LD with structure); cov-LDSC needs only ~2–10 PCs to match the
GLS operator on real data (B3). The remaining questions are still unrun.

**Status (2026-07-23): the gating engine is built and validated.** The efficient
decomposition engine below is now `pldsc/genome.py` (numpy-only, 21 unit tests in
`tests/test_genome.py`), and a genome-wide LOCO driver with cM windows is
`experiments/realdata_genome_loco.py`. The enabling primitive is the exact
spectral genome-wide score (METHOD §14 Eq. D2), `O(M N²)` with no `M × M` matrix;
it reproduces the brute-force benchmark to `3e-13` on real chr22 and the driver
reproduces the B2 LOCO shares (FULL 17%/90%, EUR −0%/−4%). Genetic-distance (cM)
windows are supported via the Beagle GRCh37 map, surfacing a first new finding:
the *windowed* structure share depends on the window unit (FULL chr22: 17% at
±250 kb vs ~11% at ±0.25 cM), because a fixed physical window spans variable
genetic distance. **It also depends on the width, and at LDSC's own default
(±1 cM) the FULL-panel share is 29.3%** — so 29%, not 17%, is the figure any
practitioner-facing claim should quote. What remains is the full 22-autosome run
(data acquisition,
`realdata_genome_loco.py --plan`) and questions 2–4.

## Questions the survey answers (that simulations cannot)

1. **How large is the structure fraction of LD in real cohorts?** Compute
   `ell_struct / ell_total` (Section 14 of METHOD.md; `ld_decomposition.py`
   generalized) across cohorts spanning homogeneous → admixed → multi-ancestry.
   The sim gave ~50% genome-wide / ~5% within-window at F_st=0.1; real value is
   the empirical headline.
2. **Which published GWAS are inflated along which PCs?** Run the per-PC inflation
   spectrum (`detect_pc_inflation.py`) on a panel of public summary statistics
   (LD Hub / GWAS Catalog), using reference-panel LOCO loadings. Output: a per-
   trait × per-PC map of residual stratification — more actionable than the
   aggregate LDSC intercept.
3. **Do published functional enrichments survive structure?** For a set of
   S-LDSC enrichment results, add the PC-loading annotation
   (`sldsc_pc_annotation.py`) as a sensitivity test: which enrichments are
   structure-robust, which vanish (confounded OR genuinely ancestry-correlated —
   the identifiability limit applies).
4. **Does structure-corrected LD move any headline h²/r_g?** Compare LDSC on the
   total vs local (`ell_local`) LD score for a panel of traits/cohorts.

## Gating first step: an efficient decomposition engine — DONE (2026-07-23)

Built as `pldsc/genome.py` and validated (see the 2026-07-23 status above). The
route taken was the exact **spectral genome-wide score** (METHOD §14 Eq. D2):
the total per-SNP squared-correlation sum is the GRM quadratic form
`||Xtilde^T xtilde_j||²`, evaluable in `O(M N²)` from the `N × N` sample Gram
with no `M × M` matrix — which is what makes millions of SNPs tractable. The
windowed (local) score stays banded and now takes cM windows. The Woodbury /
matrix-free `P_c` options below remain the path for a biobank-scale `N`; at
`N ≈ 2504` the dense eigendecomposition is trivial and is what the driver uses.

Original scoping (still the plan for large `N`):

- `V_c = sigma_g^2 K + sigma_e^2 I` with **low-rank K** (from M_grm SNPs) →
  Woodbury: `V_c^{-1} = sigma_e^{-2}(I - Z(sigma_e^2/sigma_g^2 I + Z^T Z)^{-1}Z^T)`,
  cost `O(N M_grm + M_grm^3)`; or **matrix-free conjugate gradient** (BOLT-style)
  for `P_c x`; or a **sparse GRM** (close relatives only, fastGWA-style).
- Verify the operator reproduces the dense `pldsc/_lmm` oracle to
  `rtol=1e-8`-ish, then benchmark scaling (this is the deferred M6b, now
  motivated by the decomposition rather than a new estimator).
- The windowed cross-products `X^T P_c X` stay banded (local window), so the
  LD-score computation is `O(M · window · N)`.

Deliverable: `pldsc/_operator.py` (matrix-free `P_c`) + a scaling benchmark.
Build and verify on simulation first; it is the prerequisite for touching real
data.

## Data (in rough order of accessibility)

- **1000 Genomes** (public) — multi-ancestry; the structure-fraction-of-LD map
  and the reference-panel LOCO loadings for the per-PC spectrum.
- **Public summary statistics** — LD Hub / GWAS Catalog / the Neale-lab UKB
  round; the per-PC inflation spectrum and enrichment-robustness need only
  sumstats + reference loadings (no individual data).
- **The 2017 manuscripts' cohorts** — 1000GP five European pops + the Danish
  "Where Are You From?" sample (if still accessible) — the natural reproduction
  of the original WhiLD/PC-score analyses under the new framework.
- **UK Biobank** (application/access required) — the biobank-scale demonstration
  and the relatedness/founder-structure cases; needs the efficient engine.

## Stages

1. Build + verify the efficient `P_c` operator on simulation (gating).
2. 1000G: the structure-fraction-of-LD map + reference LOCO loadings.
3. Sumstats-only: per-PC inflation spectrum across a trait panel; enrichment
   robustness for a set of S-LDSC results.
4. If (2)-(3) show material, real effects: the biobank / manuscript-cohort
   reproduction and the "does it change a published conclusion" analysis.

## Honest caveats

- **Identifiability** (METHOD.md §14): the diagnostics flag, they do not prove;
  ancestry-aligned true signal and confounding are not separable from the PCs
  alone.
- **Bercovich 2025** already does structure-corrected LD via PC residuals
  (PRIOR_ART.md) — the survey's novelty is the LDSC-tied decomposition, the
  per-PC/enrichment diagnostics, and the empirical map, not structure-corrected
  LD per se.
- **Reference-panel matching** — the per-PC spectrum on public sumstats needs
  loadings from an ancestry-matched panel; mismatch is itself a confounder.
- **Compute / access** — UK Biobank needs an application; the biobank stage is
  contingent on the efficient engine.
