# PLDSC development plan

Date: 2026-07-17

## 1. Objective and decision rule

PLDSC will test whether LD scores matched to a specified leave-one-chromosome-
out (LOCO) linear mixed-model association statistic yield calibrated LD Score
regression for quantitative traits. It is not a general LDSC rewrite, and
covariate-residualized LD is a baseline rather than the novelty claim.

The project advances only when a milestone's scientific and software gates
pass. A negative result is acceptable: if mixed-model-matched scores provide no
reproducible advantage over covariate-adjusted LDSC, PLDSC becomes a validation
and interoperability project rather than a new estimator.

The initial scope is deliberately narrow:

- quantitative traits and autosomes;
- complete individual-level genotypes from the GWAS sample for the exact
  matched-score derivation;
- fixed covariates identical to those used by the association model;
- an exact Gaussian LOCO score test before approximate association engines;
- genic or random-effect variance unless a derivation supports a stronger
  heritability interpretation.

Binary traits, meta-analysis, annotations, external reference panels, and
in-sample subsampling are out of exact scope until the quantitative-trait
method passes its oracle gates. Subsampling changes the GRM, precision
operator, and score covariance; it is a separately tested approximation, not
an implementation detail.

## 2. Scientific contract

For target chromosome \(c\), define the null covariance
\(V_c=\sigma_g^2K_{-c}+\sigma_e^2I\). Here \(C\) is the complete fixed-effect
design, including the intercept. M0 will freeze genotype centering and scaling,
the normalization of \(K_{-c}\), and the modeled target-chromosome variant set
\(X_c=\{x_k\}\). The generalized residual precision is

\[
P_c=V_c^{-1}-V_c^{-1}C(C^TV_c^{-1}C)^+C^TV_c^{-1}. \tag{1}
\]

For variant \(j\), define

\[
q_j=x_j^TP_cx_j, \qquad
z_j=\frac{x_j^TP_cy}{\sqrt{q_j}}. \tag{2}
\]

If
\(y=C\alpha+X_c\beta_c+\varepsilon\),
\(\varepsilon\sim N(0,V_c)\), and independent effects satisfy
\(E[\beta_k]=0\), \(\operatorname{Var}(\beta_k)=\tau_k^2\), and
\(\beta_c\perp\varepsilon\mid X\), then the oracle conditional moment is

\[
E[z_j^2\mid X,V_c]
=1+\sum_k\tau_k^2\frac{(x_j^TP_cx_k)^2}{q_j}
=1+\sum_k\tau_k^2q_k\{r_{jk}^{(P_c)}\}^2. \tag{3}
\]

For a homoscedastic effect model, the candidate mixed-model LD score is

\[
L_j^{(P_c)}=\sum_{k\in W_j}q_k\{r_{jk}^{(P_c)}\}^2. \tag{4}
\]

Equation (4), not merely \(\sum_k r_{jk}^2\), is the first hypothesis to test.
The sum in Equation (3) covers all modeled target-chromosome effects; replacing
it by the local window \(W_j\) in Equation (4) is a truncation approximation
that must pass a window-convergence diagnostic. The \(q_k\) terms disappear
only under an explicitly justified normalization. Variants with zero or
numerically negligible \(q_j\) have no defined score and must be rejected or
excluded under a documented, scale-aware threshold.
Before implementation is promoted, the derivation must also state:

1. the effect-size scale and exact slope estimand;
2. whether modeled-effect, tested, and regression-variant sets coincide and,
   if not, how their filters enter the summation;
3. whether \(V_c\) is known or uses fitted variance components;
4. whether the target statistic is a score, Wald, likelihood-ratio, or
   engine-specific calibrated statistic;
5. whether a separate finite-sample correction is appropriate for a
   population/reference-LD estimand; no adjusted \(r^2\) is inserted into the
   exact conditional identity in Equation (3);
6. the regression variance, weighting, and block-resampling model rather than
   assuming conventional LDSC choices remain valid;
7. what the intercept represents under environmental covariance
   misspecification and directional LD.

This work is motivated by a specific gap. Covariate-adjusted LDSC computes LD
after applying the same fixed-covariate adjustment used in the GWAS, but its
authors state that mixed-model summary statistics require random-effect-
adjusted LD and leave that computation to future work
([Luo et al., 2021](https://doi.org/10.1093/hmg/ddab130)). Residual-genotype LD
has independent theoretical support in structured populations
([Bercovich et al., 2025](https://doi.org/10.1093/genetics/iyaf009)). LOCO is
established in association testing to avoid proximal contamination
([Yang et al., 2014](https://doi.org/10.1038/ng.2876)), but that precedent does
not prove Equation (4).

## 3. Milestone roadmap

**Table 1. Scientific and engineering milestones**

| Milestone | Deliverable | Gate before proceeding |
|---|---|---|
| M0. Freeze the theorem | `docs/METHOD.md` deriving Equations (1)-(4), the slope estimand, genotype/effect/GRM scaling, modeled/tested/regression variant sets, regression weights, uncertainty model, and limiting cases | Independent algebra review; with declared normalization, the moment reduces to fixed-covariate LDSC when \(V_c=\sigma^2I\) and to ordinary LDSC up to a known scale with only an intercept |
| M1. Harden the fixed-covariate baseline (`v0.1`) | Complete the array API, deterministic simulator, an internal full-projector oracle fixture, and a separately versioned external cov-LDSC comparison; numerical kernels reject missing values and defer imputation to I/O | Small well-conditioned internal fixtures agree at `atol=rtol=1e-12`; cross-software tolerance is justified from the frozen formula/version; invariance to full-rank covariate reparameterization and declared genotype shift/scale; wheels pass on Python 3.10/NumPy 1.24 and the latest supported pair |
| M2. Build a dense LMM oracle | Private small-matrix implementation of \(P_c\), \(q_j\), \(z_j\), and Equations (3)-(4), accepting supplied \(K_{-c}\) and variance components | Scale-normalized residuals for symmetry, \(P_cC=0\), and \(P_cV_cP_c=P_c\) meet condition-aware tolerances; score covariance agrees with direct calculation; indefinite GRMs, singular \(V_c\), ill-conditioned designs, and negligible \(q_j\) produce explicit diagnostics |
| M3. Falsify the oracle moment (`v0.2`) | Repeated draws of effects and noise conditional on fixed \(X,V_c\); exact full-chromosome moment and window-truncated approximation | Per-variant/bin Monte Carlo residuals support Equation (3) with covariance-aware uncertainty and no trend by MAF, \(q_j\), ancestry loading, or GRM leverage; the omitted-window remainder is reported |
| M4. Quantify fitted-component error | Run a predeclared GEMMA score-test protocol with a supplied LOCO kinship, import its fitted null variance components, construct \(\widehat P_c\), and compare with the known-\(V_c\) oracle; PLDSC does not implement REML here | Conditional-moment error meets its preregistered gate or a correction/scope limit is documented; GEMMA fitting mode, allele scaling, and score output field are reproduced exactly and are not treated as proof for a different statistic |
| M5. Add regression and inference (`v0.3`) | Implement the M0-derived weights and resampling scheme, explicit estimand labels, provenance/compatibility schema, and machine-readable results; compare ordinary LDSC, cov-LDSC, and optional experiment-only legacy scores | Synthetic slope/intercept bias and interval coverage pass; resampling implementation matches a direct oracle and achieves calibrated coverage; incompatible adjustment metadata fails before fitting |
| M6a. Block local products | Blocked evaluation of full and windowed candidate scores | Dense and blocked results agree at `rtol=1e-10`, `atol=1e-12`; memory grows with active blocks/windows rather than chromosome-wide pair count |
| M6b. Add scalable solves | Matrix-free \(P_cx\), iterative-solver diagnostics, and preconditioning | Operator results agree with the dense oracle under condition-aware tolerances; convergence failures are explicit and reproducible |
| M6c. Interoperate (`v0.4`) | PLINK 2 adapter, minimal CLI, provenance sidecars, and only scientifically justified LDSC-format output | One end-to-end PGEN fixture passes; optional dependencies stay isolated; output round-trips with all required count metadata |
| M7. Validate and release (`v1.0`) | Synthetic phenotypes on real genotypes, independent reproduction, stable schemas/API, archival release, and limited real-trait analyses | Every preregistered gate passes; conclusions survive ancestry, window, GRM, and approximation sensitivity analyses; all headline results reproduce from committed configurations |

Current status: the fixed-covariate kernel and five deterministic tests exist,
so M1 is partly complete. M0 is the next scientific checkpoint and blocks M2;
the remaining M1 packaging and baseline work may proceed independently.

GEMMA is the first external association benchmark because it implements an
exact LMM suitable for checking a dense implementation
([Zhou and Stephens, 2012](https://doi.org/10.1038/ng.2310)). BOLT-LMM and
REGENIE are later compatibility targets because their whole-genome prediction
and calibration are not identical to the oracle Gaussian score test
([Loh et al., 2015](https://doi.org/10.1038/ng.3190);
[Mbatchou et al., 2021](https://doi.org/10.1038/s41588-021-00870-7)).

## 4. Simulation protocol

**Table 2. Required simulation dimensions**

This is a targeted scenario matrix, not a Cartesian product. M3 starts with the
smallest scenarios that distinguish the candidate moment from its limits; M5
adds estimator-level stress cases.

| Dimension | Required settings | Purpose |
|---|---|---|
| Genotypes | Homogeneous haplotypes; two discrete populations; recent admixture; pedigrees or cryptic relatives | Separate local LD, structure-induced LD, and relatedness |
| Architecture | Infinitesimal Gaussian; sparse; MAF/LD-dependent; ancestry-correlated and directional effects | Test the random-effect assumptions and estimand interpretation |
| Environment | Independent; ancestry-aligned; family-shared; covariance misspecified relative to \(V_c\) | Test intercept behavior and residual confounding |
| Adjustment | Correct and mismatched covariates; correct and misspecified GRMs; oracle and fitted variance components | Measure geometry mismatch and plug-in error |
| Reference approximation | Full in-sample oracle first; only later, random in-sample subset, ancestry-matched external, and deliberately mismatched external | Quantify rather than assume when score reuse is defensible |
| Numerical stress | GRM condition number; sample size; marker density; LD-window truncation; low MAF | Expose instability and finite-sample failure |
| Outputs | Moment calibration by score bin; bias/RMSE; intercept bias; interval coverage; type-I calibration; runtime and peak memory | Prevent selection of one flattering metric |

Each oracle-moment simulation fixes \(X,V_c\), repeatedly samples independent
effects and noise, and computes a per-variant Monte Carlo mean. Correlated SNPs
are not treated as independent replicates. Calibration uses per-variant or
predeclared-bin standardized residuals with covariance-aware uncertainty.
M0b preregisters each trend diagnostic: the test or regression, minimum effect
size of interest, covariance estimator, and multiplicity correction for MAF,
\(q_j\), ancestry loading, and GRM leverage.

Each simulation record must contain its full configuration, seed, Git revision,
truth, estimator output, and Monte Carlo standard error. Ordinary tests remain
small and deterministic; larger stochastic gates run separately from unit
tests. Gate scenarios use a 200-replicate pilot followed by at least 1,000
replicates and continue until the Monte Carlo standard error of coverage is at
most one percentage point and that of nonzero-truth relative bias is at most
one percentage point, subject to a predeclared compute cap.

### Promotion and falsification rules

1. **M3 moment gate.** Empirical mean \(z_j^2\) must agree with Equation (3)
   within covariance-aware Monte Carlo uncertainty. A regression summary may
   report intercept 0 and slope 1, but it cannot replace standardized residual
   checks when variants are correlated or the predictor range is narrow.
2. **M3 truncation diagnostic.** Full-chromosome and windowed scores report the
   omitted remainder. A less than 1% change on doubling the window marks
   convergence; genuine long-range LD may instead require a whole-chromosome
   calculation or an explicit nonconverged result.
3. **M4 plug-in gate.** Persistent conditional-moment error from fitted
   components narrows the scope. It is not rebranded as an intercept problem.
4. **M5 bias gate.** For the named M0 slope estimand in predeclared nonzero-
   truth Gaussian scenarios, absolute relative bias must be below 5%. Zero-
   truth cases use absolute bias and type-I calibration instead.
5. **M5 coverage gate.** Nominal 95% intervals must have 92%-98% empirical
   coverage, and their binomial Monte Carlo interval must contain 95%.
6. **M5 limit gate.** In predeclared nonzero-truth homogeneous scenarios,
   PLDSC must not worsen absolute relative bias by more than 0.02 compared with
   ordinary LDSC.
7. **M5 comparison gate.** Baseline selection and the aggregate RMSE metric are
   fixed before held-out seeds are run. In predeclared structured-LOCO
   scenarios, the fitted method must improve paired aggregate RMSE; one
   favorable scenario is not evidence.
8. **M6 equivalence gate.** Scalable implementations must agree with the dense
   oracle before any performance claim is reported.

The heritability label is especially constrained in admixed populations:
random-effect estimators can omit directional-LD covariance and therefore
differ from total realized genetic variance, even without environmental
confounding
([Huang et al., 2025](https://doi.org/10.1093/genetics/iyaf100)). PLDSC will
therefore report the derived estimand by name instead of calling every slope
"heritability."

## 5. Software shape

Keep the numerical layer array-first and small. Public APIs are added only
after their corresponding scientific gate passes.

**Table 3. Planned modules and constraints**

| Component | Earliest milestone | Constraint |
|---|---|---|
| `pldsc/ld.py` | Existing / M1 | Sample-by-variant NumPy arrays; independent dense oracle tests |
| `pldsc/simulate.py` | M1 | Explicit random generator; configuration and truth returned with data |
| `pldsc/_lmm.py` | M2 | Private dense oracle first; supplied GRM and variance components; no public estimator claim |
| `pldsc/regression.py` | M5 | Derived weighting/resampling separated from score construction; conventional LDSC choices are not assumed |
| `pldsc/schema.py` | M5 | Versioned adjustment and provenance metadata required before regression |
| `pldsc/_blocks.py` | M6a | Blocked products must reproduce the dense reference path |
| `pldsc/_operator.py` | M6b | Matrix-free solves with convergence and conditioning diagnostics |
| `pldsc/io.py` and CLI | M6c | Optional adapters around kernels; no bespoke genotype format |
| `tests/` | Every milestone | Algebraic properties, deterministic fixtures, then end-to-end tests |
| `experiments/` | M3 onward | Committed configurations and summaries, not generated bulk data |

Engineering rules:

- use `float64` for the reference path and make tolerances explicit;
- never construct an \(M\times M\) LD matrix;
- permit dense \(N\times N\) matrices only in the reference implementation;
- require an operator interface before claiming biobank scale;
- postpone pandas/xarray models, GPU/JAX/Numba paths, distributed execution,
  dashboards, and multiple genotype readers;
- benchmark only after numerical agreement, recording machine, Python, NumPy,
  and BLAS details.

M1 records a same-machine baseline for \(N=500\), \(M=5{,}000\), and a median
active window of 101 variants. M6a may be at most 25% slower on that fixture and
must process \(N=2{,}000\), \(M=50{,}000\), median active window 201 in at most
15 minutes with no more than 512 MiB peak working memory beyond the input and
output arrays on the documented benchmark machine. M6b uses a predeclared
\(N=5{,}000\) positive-definite operator with condition number at most \(10^4\),
requires a relative linear-solve residual at most \(10^{-8}\) within 500
iterations, and reports failures rather than silently returning approximate
scores. These are engineering promotion gates, not portable speed claims.

M1 supports the minimum pair Python 3.10/NumPy 1.24 and the latest stable pair
in continuous integration. `pgenlib` becomes an optional dependency only when
the M6c adapter lands.

PLINK 2 PGEN/PVAR/PSAM is the sole planned genotype adapter
([PLINK 2 format documentation](https://www.cog-genomics.org/plink/2.0/formats)).
For interoperability, do not label Equation (4) as conventional `L2` unless M0
proves the required normalization. If that equivalence is established, emit
the conventional `CHR`, `SNP`, `BP`, and `L2` columns together with `.l2.M` and
`.l2.M_5_50` counts; otherwise use an explicitly versioned PLDSC schema
([official LDSC LD-score tutorial](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial)).
The JSON sidecar introduced in M5 records formula revision, package and Git
revision, sample and variant counts, adjustment source, genome build and map,
window, finite-sample treatment, input identifiers/checksums, and simulation
seed when applicable.

## 6. Commit and review discipline

Every commit checkpoint is small and green. Scientific derivations are reviewed
before their optimized implementations.

**Table 4. Planned commit checkpoints**

| Milestone | Commit intent | Required check before commit |
|---|---|---|
| M0a | `docs: derive the PLDSC mixed-model moment` | Algebra review and rendered equations |
| M0b | `docs: preregister PLDSC validation gates` | Scenario manifest and uncertainty audit |
| M1a | `feat: harden fixed-covariate LD scores` | Unit and cov-LDSC fixture tests |
| M1b | `build: add the supported Python matrix` | Clean builds and imports on both dependency pairs |
| M2a | `test: specify mixed-model oracle identities` | Scale-aware reference assertions |
| M2b | `feat: add the dense mixed-model oracle` | Identity, failure-mode, and direct-covariance tests |
| M3 | `test: validate the oracle conditional moment` | Seed replay, Monte Carlo uncertainty, and truncation summaries |
| M4a | `feat: import fitted LOCO components` | Known-\(V_c\)/fitted-component comparison |
| M4b | `test: reproduce the GEMMA score protocol` | Statistic, fitting, scaling, and output parity |
| M5a | `feat: add derived PLDSC regression` | Synthetic slope/intercept recovery |
| M5b | `feat: add inference and provenance checks` | Coverage, resampling-oracle, and mismatch tests |
| M6a | `feat: block local mixed-model scores` | Dense/blocked equivalence and memory benchmark |
| M6b | `feat: add matrix-free mixed-model solves` | Dense/operator equivalence and convergence tests |
| M6c | `feat: add PLINK interoperability` | PGEN integration and schema round-trip tests |
| M7 | `release: validate PLDSC 1.0` | Full validation manifest, clean build, independent reproduction |

No release milestone is awarded because the version number looks lonely.

## 7. Primary references

1. Bulik-Sullivan et al. (2015), [LD Score regression distinguishes
   confounding from polygenicity](https://doi.org/10.1038/ng.3211).
2. Zhou and Stephens (2012), [Genome-wide efficient mixed-model analysis for
   association studies](https://doi.org/10.1038/ng.2310).
3. Yang et al. (2014), [Advantages and pitfalls in the application of
   mixed-model association methods](https://doi.org/10.1038/ng.2876).
4. Loh et al. (2015), [Efficient Bayesian mixed-model analysis increases
   association power in large cohorts](https://doi.org/10.1038/ng.3190).
5. Luo et al. (2021), [Estimating heritability and its enrichment in
   tissue-specific gene sets in admixed populations](https://doi.org/10.1093/hmg/ddab130).
6. Mbatchou et al. (2021), [Computationally efficient whole-genome regression
   for quantitative and binary traits](https://doi.org/10.1038/s41588-021-00870-7).
7. Bercovich et al. (2025), [Measuring linkage disequilibrium and improvement
   of pruning and clumping in structured populations](https://doi.org/10.1093/genetics/iyaf009).
8. Huang et al. (2025), [Interpreting SNP heritability in admixed
   populations](https://doi.org/10.1093/genetics/iyaf100).
