# PLDSC scientific idea review

Date: 2026-07-17

## 1. Verdict

The research problem is real: pooled genotypes can contain long-range
correlation induced by population structure, while GWAS statistics are usually
computed after covariate or mixed-model adjustment. LD scores should reflect
the same statistical geometry as the association statistics.

The 2017 `newldsc` implementation is not a sound starting implementation. It
contains useful ideas, but it combines an unproved GRM-whitening transform, an
ad hoc PC score, unweighted regressions, and hard-coded data handling. It also
contains correctness bugs and has no tests.

The recommended restart is therefore:

- **baseline track:** reproduce covariate-adjusted LDSC with explicit,
  independently tested array operations;
- **research track:** derive and falsify an LD score matched to
  leave-one-chromosome-out (LOCO) mixed-model association statistics;
- **legacy track:** retain the original whitening and PC-score constructions
  only as named comparators.

This is a **go** for a research project and a **no-go** for a mechanical port.

## 2. What the legacy prototype was trying to do

For each target chromosome, the prototype builds a relationship matrix from
the other 21 chromosomes, transforms target-chromosome genotypes with a
Cholesky factor of its pseudoinverse, and computes local LD on the transformed
genotypes. It also computes a per-SNP score from squared correlations with the
top PCs, weighted by their eigenvalues. Finally, it compares ordinary LD,
transformed LD, and PC scores against GWAS chi-square statistics.

That contains two separate hypotheses:

1. LD should be calculated after removing covariance attributable to sample
   structure.
2. SNP loadings on ancestry PCs may explain nonuniform confounding in GWAS
   statistics.

Both hypotheses are plausible. Neither is validated by the legacy code.

## 3. What is already established

Let \(X\in\mathbb{R}^{N\times M}\) be standardized genotypes and let
\(C\in\mathbb{R}^{N\times K}\) contain the centered covariates used by the
GWAS. Covariate-adjusted LDSC uses the residual-forming matrix

\[
P_C = I-C(C^T C)^+C^T, \qquad \widetilde X=P_CX. \tag{1}
\]

It then defines a covariate-adjusted LD score from correlations between columns
of \(\widetilde X\):

\[
\ell_j^{(C)}=\sum_{k\in W_j}
\operatorname{Cor}(\widetilde X_j,\widetilde X_k)^2. \tag{2}
\]

When the GWAS and LD calculation use the same covariates, the working LDSC
moment is

\[
E[\chi_j^2]\approx 1+Na+\frac{Nh_g^2}{M}\ell_j^{(C)}. \tag{3}
\]

Luo et al. derived and evaluated this approach as cov-LDSC in 2021. A separate
2025 study by Bercovich et al. gives a population-genetic justification for LD
based on genotype residuals and shows that it removes structure-induced LD in
pooled and admixed samples. Thus, residualized LD is the correct baseline, not
a novel PLDSC claim.

The legacy PC-score idea also overlaps PC-loading regression. The important
difference is that PC-loading regression uses **signed** SNP loadings and
signed GWAS effects. The legacy score squares and aggregates loadings using PC
eigenvalues. Those weights describe genotype variation, not the phenotype's
confounding direction, and cross-PC terms are lost.

## 4. The research opportunity

Covariate-adjusted LDSC does not directly cover summary statistics from a
linear mixed model. The cov-LDSC authors explicitly identify random-effect-
adjusted LD as future work. The legacy LOCO idea points at this gap, but its
specific transform is not derived from the mixed-model score statistic.

For a target chromosome \(c\), let

\[
V_c=\sigma_g^2K_{-c}+\sigma_e^2I \tag{4}
\]

be the null covariance using a LOCO relationship matrix. With fixed covariates
\(C\), the generalized residual precision is

\[
P_c=V_c^{-1}-V_c^{-1}C(C^TV_c^{-1}C)^+C^TV_c^{-1}. \tag{5}
\]

A natural correlation matched to the null mixed-model score test is then

\[
r_{jk}^{(P_c)}=
\frac{x_j^TP_cx_k}
{\sqrt{(x_j^TP_cx_j)(x_k^TP_cx_k)}}. \tag{6}
\]

Equation (6) is a **research hypothesis**, not yet a PLDSC estimator. The next
scientific task is to derive the corresponding expectation of \(\chi_j^2\),
including how the fitted variance components, LOCO construction, causal-effect
model, and finite-sample correction enter. Simulations must then try to break
that derivation.

## 5. Technical assessment

**Table 1. Legacy components and restart decisions**

| Component | Legacy approach | Restart decision |
|---|---|---|
| Structure-adjusted LD | Cholesky factor of a pseudoinverse LOCO GRM | Use covariate residualization as the baseline; derive Equation (6) for the mixed-model track |
| PC contribution | Eigenvalue-weighted sum of squared SNP-PC correlations | Export signed per-PC loadings; evaluate signed PC-loading regression separately |
| Regression | Binned and unbinned ordinary least squares | Use LDSC weights and block jackknife inference |
| Reference data | External 1000 Genomes panel with a fixed population mixture | Require study-matched in-sample data or a documented random subsample for adjusted scores |
| Windows | Fixed 1 cM and an adjusted-r-squared threshold | Make window convergence an empirical diagnostic; do not threshold by default |
| Software | Python 2 scripts with hard-coded cluster paths | Reimplement array-first in Python 3; add file formats only after numerical validation |

Important failure modes are:

1. **Mismatched adjustment.** LD scores are biased if the GWAS and LD
   calculations use different covariates or different mixed-model geometry.
2. **Reference mismatch.** A convenient external population mixture is not a
   substitute for the study's ancestry and relatedness distribution.
3. **Unstable whitening.** A pseudoinverse plus a nearly singular Cholesky
   transform is sensitive to regularization and factor orientation.
4. **Wrong heritability interpretation.** In admixed populations, directional
   LD can make LDSC estimate genic variance rather than total genetic variance;
   the gap can exist without environmental confounding.
5. **PC-score nonidentifiability.** A generic squared loading score cannot
   identify the phenotype-specific direction of stratification.
6. **Overstated novelty.** Covariate-adjusted LD scores and PC-loading
   regression already exist. Any new contribution must be the mixed-model
   derivation, a demonstrably better estimator, or substantially better
   software and validation.

## 6. Falsification-first validation plan

**Table 2. Minimum simulation gates before real-data claims**

| Scenario | Required result |
|---|---|
| Homogeneous population, ordinary least squares GWAS | Adjusted and ordinary LD scores agree within Monte Carlo error; heritability bias is below 5% |
| Structured population with no within-population LD | Adjusted cross-chromosome LD centers on its finite-sample null while pooled ordinary LD does not |
| Admixed polygenic trait with matched PC adjustment | Covariate-adjusted LDSC improves absolute bias over ordinary LDSC and remains below 10% bias |
| Environmental effect aligned with ancestry | The intercept absorbs null inflation without a material slope shift |
| LOCO mixed-model GWAS | Equation (6)-based scores outperform both ordinary and fixed-covariate scores across variance-component settings |
| Directional LD or ancestry-correlated effects | PLDSC reports the estimand gap and does not label genic variance as total heritability |
| Mismatched covariates or ancestry composition | Bias is quantified and the software emits a clear incompatibility diagnostic |

The first executable milestone should use tiny deterministic genotype matrices
and simulated phenotypes. Real 1000 Genomes I/O, summary-statistic munging, and
plotting come later. Otherwise the project risks producing polished answers to
an equation it never actually specified. A classic research-software hobby.

## 7. Recommended project scope

The first PLDSC release should contain:

1. covariate projection and residualized correlations;
2. standard and covariate-adjusted local LD scores;
3. signed per-PC SNP loadings;
4. a LOCO mixed-model score-correlation prototype based on Equation (6);
5. deterministic simulations and explicit estimand checks;
6. weighted regression and block-jackknife inference only after the numerical
   kernels pass those checks.

It should not initially contain genotype-format plumbing, phenotype-specific
dictionaries, plotting dashboards, or claims that the method removes all
population structure.

## 8. Primary references

1. Bulik-Sullivan et al. (2015), [LD Score regression distinguishes confounding
   from polygenicity](https://doi.org/10.1038/ng.3211).
2. Bhatia et al. (2016), [Correcting subtle stratification in summary
   association statistics](https://doi.org/10.1101/076133).
3. Luo et al. (2021), [Estimating heritability and its enrichment in
   tissue-specific gene sets in admixed populations](https://doi.org/10.1093/hmg/ddab130).
4. Bercovich et al. (2025), [Measuring linkage disequilibrium and improvement
   of pruning and clumping in structured populations](https://doi.org/10.1093/genetics/iyaf009).
5. Huang et al. (2025), [Interpreting SNP heritability in admixed
   populations](https://doi.org/10.1093/genetics/iyaf100).
