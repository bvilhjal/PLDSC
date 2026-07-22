# PLDSC experiments

Committed experiment configurations and their summaries (per
`docs/DEVELOPMENT_PLAN.md` §5). Not generated bulk data.

## `m_spike_kill_or_continue.py` — the M-spike (kill-or-continue gate)

Throwaway dense oracle that decides whether the mixed-model-matched LD score is
worth building, before any packaging. Builds the dense generalized residual
precision `P_c` (Eq 1) on one tiny structured cohort and runs three readouts.
Run in the project env:

```
conda run -n pldsc python experiments/m_spike_kill_or_continue.py
```

(env `pldsc`: Python 3.14, NumPy 2.4, numba 0.66; RNG bit generator PCG64,
seed 20260717.)

### Verdict: CONTINUE

Two scenarios (structured two-population + sib pairs; homogeneous control),
N=400, 240 SNPs on 3 chromosomes, 6000 Monte-Carlo replicates.

**Readout B — Eq (3) is numerically correct.** The Monte-Carlo `E[z_j^2]` matches
the Eq (3) moment: per-SNP standardized residuals have RMS ~1.55 sigma (max ~4.3),
which is expected because the 80 SNP means share the same MC draws and are
correlated — the aggregate agreement rules out a systematic error. Window
truncation (1.5 cM) drops ~12% of the full sum (median), the expected Eq-4
approximation to quantify in M3. Separately, the slope-estimand check (this is the
*Readout A2* "correct regressor" row, not Readout B): regressing the calibrated
oracle `z^2` on the full-chromosome `L_j^(P_c)` recovers
`(intercept, slope) = (0.997, 6.246e-4)` against the true `tau^2 = 6.229e-4`
(homogeneous 0.983 / 6.111e-4 vs 5.86e-4) — the cleanest confirmation that the
slope is `tau^2`.

**Readout A1 — the calibration/confounding use-case, under a null target
chromosome (`tau^2 = 0`).** Mean `z^2` is the LDSC-intercept-style confounding:

| association model | structured | homogeneous |
|---|---|---|
| ordinary OLS | **12.29** | **9.07** |
| covariate-adjusted OLS (intercept + 2 LOCO PCs) | 1.04 | 0.94 |
| mixed-model `P_c` (PLDSC) | **1.00** | **1.00** |

The mixed-model score is exactly calibrated in both; ordinary OLS is massively
inflated because it ignores `V_c` (the polygenic-background + structure
covariance) — the confounding PLDSC removes.

**Readout A2 — LD-score choice matters.** Regressing the *same* calibrated
oracle `z^2` on each LD score, only the correct `L^(P_c)` gives intercept ~1
with a good fit (R^2 = 0.90); `L_ordinary` (a=1.65, R^2=0.46) and `L_cov-LDSC`
(a=0.90, R^2=0.38) are the wrong regressor and bias the fit.

**Long-range LD.** `P_c` shrinks mean chr0-chr2 `r^2` from 0.053 (raw) to 0.0015
(structured) — the structure-induced cross-chromosome LD the manuscripts were
chasing, removed by the residual-precision metric.

### Honest caveats / what this does NOT yet show

- With two *discrete* populations, covariate-adjusted LDSC (2 fixed PCs) is
  already ~calibrated (A1: 0.94–1.04), so this scenario does not separate PLDSC
  from cov-LDSC on the intercept. The PLDSC-over-cov-LDSC advantage should
  appear where fixed PCs cannot capture the covariance — cryptic relatedness and
  continuous admixture. That is the first M3 scenario to add.
- Everything here uses an **oracle** `V_c`. Fitted-`V_c` plug-in error (M4) is
  untested.
- This is a throwaway spike, not the `pldsc/_lmm.py` oracle (M2) or its unit
  tests.

## `m3_relatedness_admixture.py` — does PLDSC beat cov-LDSC?

The M-spike showed PLDSC beats *ordinary* LDSC and that discrete-population
cov-LDSC is already ~calibrated. This experiment tests the actual novelty claim:
does PLDSC beat **cov-LDSC**? cov-LDSC corrects structure with a few fixed PCs,
which handles a low-rank admixture axis but not high-rank cryptic relatedness
(each family is its own block). PLDSC uses the full GRM `V_c`.

Headline metric: **mean null chi^2** = `mean_j (x_j^T P V_c P x_j)/(x_j^T P x_j)`,
the LDSC-intercept-style confounding under a null target chromosome. It is
**exactly 1** for the oracle `P_c` (property P2, `docs/METHOD.md`) and inflated
for OLS operators. Computed analytically (no Monte Carlo); one MC pass
cross-checks Eq (3). N=400, `sigma_g^2=0.7`, `sigma_e^2=0.3`, PCG64 seed
20260718.

Mean null chi^2 (want 1.000); the script sweeps 2/5/10/20/50 PCs, table shows a
representative subset:

| operator | A. relatedness | B. admixture | C. both (biobank-like) |
|---|---|---|---|
| ordinary (0 PC) | 1.55 | 1.64 | 2.35 |
| cov-LDSC (2 PC) | 1.48 | **0.97** | 1.47 |
| cov-LDSC (10 PC) | 1.26 | 0.88 | 1.25 |
| cov-LDSC (20 PC) | 1.06 | 0.78 | 1.04 |
| cov-LDSC (50 PC) | 0.66 | 0.60 | 0.65 |
| **PLDSC (full GRM)** | **1.000** | **1.000** | **1.000** |

### Verdict: PLDSC's advantage over cov-LDSC is real and specific to relatedness

- **PLDSC is exactly calibrated (1.000) in all three, with no tuning** — the P2
  identity guarantees it.
- **Admixture (B):** cov-LDSC with the *right* small number of PCs (2) is also
  ~calibrated (0.97) — an honest tie, consistent with Luo et al. 2021. PCs
  handle a low-rank admixture axis.
- **Relatedness (A) and the biobank case (C):** cov-LDSC stays inflated for any
  feasible PC count (1.48 at 2 PC, still 1.04–1.06 at 20 PC) because family
  structure is high-rank — you cannot absorb it into a few fixed covariates.
- **cov-LDSC needs tuning and has no safe setting:** too few PCs → inflation,
  too many → deflation (0.60–0.66 at 50 PC, from over-projecting variance out).
  The sweet spot is scenario-dependent (≈2 PC for admixture, ≈20 PC for
  relatedness) and fragile. PLDSC removes the tuning problem entirely.
- Eq (3) MC cross-check (scenario A, 15000 reps): regressing `E[z^2]` on the
  full `L^(P_c)` gives (intercept, slope) = (1.021, within ~4% of `tau^2`),
  converging to the exact `(1, tau^2)` as replicates grow.

This is the result that justifies the mixed-model-matched score over cov-LDSC:
the win is **relatedness**, exactly the regime mixed models exist for. Still
oracle `V_c` only; fitted-`V_c` (M4) is tested below.

## `m4_fitted_vc.py` — does calibration survive a FITTED `V_c`?

M-spike and M3 used an oracle `V_c`. Real use estimates the variance
components, so this experiment quantifies the plug-in error on an
admixture+relatedness cohort (N=400, true `sigma_g^2=0.7`, `sigma_e^2=0.3`,
`h^2=0.70`). It uses the tested `pldsc._lmm` oracle and `pldsc.simulate`
generators; the REML fit is local to the experiment (the package ships no
REML — it stands in for the fitted null components M4 would import from GEMMA).
`P_c` depends on `V_c` through its scale and its `h^2` ratio, and REML fits both.

**M4a — deliberate misspecification sweep** (true `V_c` fixed, vary assumed `h^2`):

| assumed h^2 | mean null chi^2 |
|---|---|
| 0.00 (= cov-LDSC) | 2.341 |
| 0.50 | 0.727 |
| 0.60 | 0.813 |
| **0.70 (truth)** | **1.000** |
| 0.80 | 1.427 |
| 0.90 | 2.824 |

The intercept passes exactly through 1 at the true `h^2` and is a steep,
monotone function of it: a ±0.1 error in `h^2` moves the intercept ~20–40%.
This is the calibration-diagnostic curve (`docs/METHOD.md` M9) — and it says the
intercept genuinely depends on an accurate `V_c`.

**M4b — REML-fitted plug-in error** (300 null replicates, fit `V_c` from each):

- REML recovery is unbiased: `sigma_g^2 = 0.700 ± 0.105`, `h^2 = 0.697 ± 0.040`.
- **Fitted-`P_c` intercept: mean 1.010, sd 0.078, [5%,95%] = [0.89, 1.15]** — versus
  the oracle 1.000 and cov-LDSC 2.341.

### Verdict: calibration survives fitting, with a caveat to carry into M5

With REML-fitted components the intercept is nearly unbiased (1.01) and far
better than cov-LDSC (2.34), so the M3 advantage is not an oracle artifact.
**But** the intercept inherits the sampling uncertainty of `V_hat_c` (sd 0.078
here at N=400, driven by the ±0.04 `h^2` error through the steep M4a curve).
Consequences to act on:

- The intercept's standard error must **propagate variance-component
  uncertainty**, not just the LD-score regression residual — a concrete
  requirement for the M5 inference/resampling design.
- Plug-in error shrinks as REML precision improves with N, so the biobank regime
  is more favorable than this N=400 probe.
- Fitted-`V_c` parity against GEMMA's actual null components (rather than this
  local REML) remains the full M4 deliverable.

## `m5_inference.py` — LD-score regression inference and the fitted-`V_c` caveat

Uses the `pldsc.regression` OLS path with `n_blocks=M`, hence a
delete-one-SNP jackknife, and the `pldsc._lmm` oracle. Two arms per replicate on
an admixture+relatedness
cohort (N=400, 80 target SNPs, R=300): an oracle-`V_c` arm and a REML-fitted-`V_c`
arm. Coverage of the intercept (truth 1) and slope (truth `tau^2`):

| arm / SE | intercept coverage | slope coverage |
|---|---|---|
| oracle `V_c`, jackknife | 94.3% | 95.3% |
| fitted `V_c`, jackknife only | 95.0% | 94.3% |
| fitted `V_c`, jackknife + `V_hat_c` term | 95.0% | — |

### Verdict: the toy fixed-operator experiment covers; genomic inference remains open

- The delete-one-SNP jackknife gives nominal coverage in this weak-dependence
  toy. That does **not** validate chromosome/position-aware genomic blocks,
  correlated SNPs, recomputed two-step weights, or cross-chromosome pooling.
- **The nuance M4 had hidden:** M4's stable "mean inflation" and M5's
  extrapolated "regression intercept" are different estimators. At only M=80
  target SNPs the regression intercept is noise-dominated (SE ≈ 1.63), so the
  `V̂_c` term (0.145) is negligible and coverage is fine with or without it.
- **But `V̂_c` uncertainty is a *shared systematic*** — one fitted `V_c` shifts
  every SNP's statistic together — so it does **not** shrink with the number of
  SNPs, whereas the regression SE does (~`1/sqrt(M)`). The decomposition
  (regression SE 1.63 at M=80 vs `V̂_c` SE 0.145) puts the crossover at
  **~10,000 SNPs**: beyond that the extrapolation predicts the jackknife-only
  intercept SE may under-cover, so
  a genome-scale M5 must propagate variance-component uncertainty into the
  intercept SE. The slope, far less sensitive to `V_c`, stays calibrated
  throughout.

The fitted-`V_c` concern is real, but the displayed crossover is illustrative:
the derivative and variance are estimated across Monte Carlo replicates, other
variance-component directions and same-phenotype cross-covariance are omitted,
and adding one term in quadrature is not a deployable correction. An outer
refit/bootstrap or joint sandwich is the M0 requirement.

## `map_ldsc_miscalibration.py` — when is LDSC miscalibrated on LMM sumstats?

Investigates the paper angle *"LDSC is miscalibrated for mixed-model GWAS
summary statistics — here's the corrected LD score, and a map of when it
matters."* An LMM GWAS produces the exact moment `E[z_j^2] = 1 + tau^2 L_j^(P_c)`;
a practitioner regresses `chi^2` on a conventional LD score (ordinary /
cov-LDSC, in-sample / external). We regress the *analytic* moment (no MC) on each
candidate LD score, holding N=400 and true `tau^2` fixed, averaging over 6
genotype realizations (single realizations are too noisy — an early run showed a
spurious relatedness "effect" that averaged away). Two readouts: intercept
(target 1) and h² drift = slope(regime)/slope(clean) (target 1).

**Structure sweep (unrelated, vary fst): a clean, strong effect.**

| fst | ordinary in-sample | cov-LDSC 2PC | PLDSC oracle |
|---|---|---|---|
| 0.05 | drift 1.08, int 1.05 | 1.08 / 1.05 | 1.01 / 1.06 |
| 0.10 | drift 0.88, int 1.09 | 1.01 / 1.06 | 0.97 / 1.07 |
| 0.15 | **drift 0.69**, int **1.13** | **1.00** / 1.06 | **1.02** / 1.05 |

**Relatedness sweep (homogeneous, vary sib-pairs at fixed N): no effect.** Even
with 180/200 sib-pairs, h² drift and intercepts stay ~1.0 for *every* method,
ordinary LDSC included.

### Verdict: the angle does NOT differentiate PLDSC

1. **The miscalibration is a structure effect, not a relatedness effect.**
   Ordinary LDSC underestimates h² by ~30% and inflates the intercept as `fst`
   grows; relatedness — already handled by the LMM in the GWAS — does not
   propagate into the *local* (windowed) LD-score regression (Wall 1, confirmed).
   So PLDSC's one edge over cov-LDSC (relatedness) is invisible here; it lives in
   the association step, not the LD-reference step.
2. **cov-LDSC already fixes the structure effect** (Luo 2021, by design). PLDSC
   matches it but adds nothing in the LDSC-application setting.
3. The structure result largely **reproduces Luo 2021** (standard LDSC
   underestimates h² in structured samples; cov-LDSC corrects it), plus a *null*
   result on relatedness — thin for a standalone paper and not a new-method one.

Conclusion: this paper angle reduces to "use cov-LDSC," which exists. PLDSC's
genuine niche is the **association / within-family step** (see the M3 relatedness
result and the applications analysis), not LDSC application. A useful negative
result that redirects effort away from a paper whose "corrected LD score" is
already the standard tool.

## `sldsc_pc_annotation.py` — PC loadings as S-LDSC annotations

Tests using the squared PC loading (the 2017 "PC score", LOCO) as an annotation
in stratified LDSC. Rationale: structure inflates `chi^2` by `~ b^2 gamma_{j,1}^2`
(Bhatia et al.), which is *linear in the squared loading* — so a squared-loading
annotation is the natural chi^2-scale sponge for structure confounding. Two-
population cohort (N=500, fst=0.12), a functional annotation that is enriched for
heritability, ordinary (OLS) GWAS, WLS S-LDSC, 40 paired replicates. Compares
standard S-LDSC (base+func) vs +PC-loading, with the functional annotation either
correlated with or independent of ancestry.

**Functional annotation correlated with ancestry:** standard S-LDSC inflates the
functional coefficient up to **×3.5** at b=1 (structure misread as enrichment);
adding the PC-loading annotation removes the inflation entirely (stable ×1.0),
absorbing it into the PC coefficient. **But** at b=0 the PC annotation already
cuts the functional coefficient to **0.51×** — it is collinear with the
ancestry-correlated functional annotation and cannibalizes genuine enrichment.

**Functional annotation independent of ancestry (control):** standard S-LDSC is
*not* inflated (structure hurts power, deflating it) and the PC annotation costs
nothing (B/A at b=0 = 1.03).

### Verdict: a conservative structure de-confounder for enrichment

The PC-loading annotation cleanly removes structure-driven *false* enrichment for
ancestry-correlated categories, but **cannot separate genuine ancestry-correlated
enrichment from confounding** (they are collinear), so it attenuates real
enrichment too. Best used as a **sensitivity analysis**: enrichment that survives
adding a PC-loading annotation is robust to structure; enrichment that vanishes is
either confounded or genuinely ancestry-correlated (indistinguishable here). As a
primary estimator it is biased toward the null for ancestry-correlated categories.

Why it matters: this is an *enrichment / robustness* application, not h²/intercept
— so it sidesteps the "reduces to cov-LDSC/REML" problem that sank the other
angles, and it gives the 2017 "PC score" a concrete home. Caveats: overlaps
Bhatia et al. (signed PC-loading regression on beta-hat; this is the squared,
chi^2-scale, S-LDSC-native analog) so **prior art must be checked** before
claiming novelty; and the S-LDSC here is crude (small M, negative intercepts show
the fit needs proper constrained/weighted S-LDSC over many SNPs), so numbers are
directional.

## `detect_pc_inflation.py` — mapping GWAS structure inflation to specific PCs

Regresses the GWAS per-SNP statistic on per-SNP LOCO PC loadings, one coefficient
per PC, to get an inflation *spectrum* (which PCs, not just an aggregate lambda_GC
/ LDSC intercept). Signed (beta-hat on loadings) and chi^2-scale (squared
loadings) versions; per-PC z from a block jackknife; 4-population cohort.

- **S1 (confound along PC1 & PC3):** cleanly mapped — signed |z| 96/96 at PC1/PC3,
  others ~0-2; chi^2 8.5/8.2. Detection works with PC resolution and specificity.
- **S2 (genuine signal aligned with PC2, no confound):** falsely flags PC2
  (|z|=56) — the **identifiability limit**: ancestry-aligned true signal is
  indistinguishable from confounding.
- **S3 (null):** chi^2 version flat (|z|<=1.6); signed version shows residual
  |z|~4 at structure PCs (structured genotypes project polygenic signal onto
  structure axes -> naive jackknife SE is anti-conservative). Signal/null ~24x,
  but a proper permutation/theoretical null is needed for calibration.

### Verdict: works, but largely Bhatia + incremental

A PC-resolved structure-QC diagnostic, more informative than an aggregate
inflation number. But the signed version is essentially Bhatia et al. 2016; the
value-adds are the chi^2-scale magnitude spectrum (cleaner null), the PC
resolution as a QC output, and the summary-statistics-only framing. Same
identifiability wall (flag, not proof) and calibration caveat as the S-LDSC
annotation. Natural companion to `sldsc_pc_annotation.py`: together a coherent
(if modest) "PC-loading structure diagnostics for GWAS/LDSC" methods note —
detect which PCs inflate a GWAS, and test whether an enrichment is
structure-robust — building on Bhatia and the 2017 PC score.

## `ld_decomposition.py` — the local/non-local decomposition (framework linchpin)

Tests the central identity of the "local vs non-local LD" framing:
`ell_j (ordinary LD score) = ell_j^local + ell_j^struct`, with
`ell_j^local = sum_k (r_jk^(P_c))^2` (structure-whitened linkage) and the claim
that the non-local part equals the 2017 PC score,
`ell_j^struct ~ s_j = sum_p (lambda_p/lambda_t) corr(g_j, PC_p)^2`. Multi-
population cohort, windowed and genome-wide, 10 replicates.

| scope | corr(ell_struct, PC score) | struct/total | ell_ord = ell_local + ell_struct |
|---|---|---|---|
| windowed | 0.70 | 0.05 | 1.41 = 1.34 + 0.07 |
| genome-wide | **0.988** | 0.50 | 4.20 = 2.12 + 2.08 |

### Verdict: the decomposition and unification hold

Genome-wide, the non-local component of the LD score correlates **0.99** with the
PC score and the split is cleanly additive. So **the non-local part of the LD
score IS the PC score** — one identity that unifies the P_c operator, cov-LDSC,
the legacy WhiLD score, and the 2017 PC score, with standard LDSC's structure
bias = conflating `ell^local` (h²) with `ell^struct` (confounding). The structure
fraction is 50% genome-wide but only 5% within a local window, which explains why
LDSC's windowed slope is mostly protected while structure surfaces in the
intercept and long-range LD. This is the validated backbone for reframing the
project around the decomposition + its structure diagnostics, with the
identifiability limit (genuine ancestry-aligned signal loads on the same PCs) as
a stated limit rather than a hidden flaw.

## `misspecified_background.py` — the one non-circular robustness test

Every other experiment generates data under exactly the assumed model. This one
deliberately violates `epsilon ~ (0, V_c)` (METHOD.md §6's flagged axis): the
off-chromosome background is *sparse* or *directional*, so the true covariance
differs from the `V_c` used to build `P_c`. Exact null-target-chromosome intercept
= `mean_j null_inflation(P_c, X, V_true)`; 6 realizations; control must return
1.000 (the P2 identity).

| off-chromosome background | mean null intercept |
|---|---|
| infinitesimal (control) | 1.000 |
| sparse — 10% of SNPs causal | 0.996 |
| sparse — 2% | 0.996 |
| sparse — 0.5% | 0.993 |
| directional — top-axis extra variance α=0.5 | 1.001 |
| directional — α=1 | 1.002 |
| directional — α=2 | 1.004 |

### Verdict: robust, with honest scope

The oracle intercept is **remarkably robust** to off-chromosome architecture
misspecification — even a 0.5%-sparse or doubled-directional background moves it
<1%. Mechanism: both misspecifications live mostly in the low-rank structure
subspace `P_c` already whitens (structure is robustly estimated from any SNP
subset; extra variance along the top structure axis is suppressed by `V_c^{-1}`).
This is the **first non-circular experiment** — it violates the model assumption
rather than assuming it — and the answer is reassuring within its scope. Scope
limits (stated, not hidden): it tests the null *intercept* only; calibration IS
sensitive to variance-component/scale misspecification (M4a); and it does **not**
test a background structure axis *absent* from the GRM used to build `V_c` (a
misspecification orthogonal to the whitened subspace), which is the harder,
still-untested case.


## `within_family_killtest.py` — the within-family kill-test (fair version)

The one candidate niche after the summary-statistics angle collapsed was
within-family / sib-GWAS, where relatedness is exactly known by pedigree so `P_c`
is exact with no fitted-`V_c` fragility. Claim under test: the within-family
genetic covariance that survives family-mean-centering (`M_B = I − B(BᵀB)⁻¹Bᵀ`)
biases the standard sib-GWAS LD-score regression, and `P_c` removes it. The design
was hardened by a 3-lens adversarial review *before* coding (all three lenses
returned *fundamentally-flawed* on the first draft, which used the circular
P2-identity metric and a strawman baseline — see the note in the script).

- **Part 1 (the crux).** For full sibs `R = 0.5(I + B)` and `M_B B = 0`, so
  `M_B R M_B = 0.5·M_B` **exactly** — measured `c_hat ≈ 0.47–0.55`. The surviving
  covariance is a homoscedastic *scale*, not LD-correlated confounding; a real
  sib-GWAS's cluster-robust / effective-N SE (Howe et al. 2022) divides it out and
  the null χ² is 1 with no modelling of `V`.
- **Part 2 (MC).** Given each estimator a consistent (cluster-robust-equivalent)
  SE, within-sib FE, pedigree-`P_c`, and oracle-`P_c` **all** recover intercept ~1
  and slope ~τ² — equivalent. Pedigree-`P_c` is the most off under admixture.
- **Part 3 (Wall-1).** `corr(FE surviving-leverage, windowed L) ≈ 0` — the
  structured remainder is off-target/cross-chromosome, orthogonal to the windowed
  LD score, so it cannot bias the slope either.

### Verdict: KILLED for standard full-sib designs

The standard sib-GWAS pipeline is already calibrated; PLDSC's exact `P_c`
reproduces it, it does not beat it, and pedigree-only `P_c` is *more* fragile
under admixture (it forfeits the "no fitted-`V_c`" premise if it uses the realized
GRM instead). The only unruled-out regime is **heterogeneous pedigrees** (mixed
full/half sibs, cousins) where `M_B K M_B` is not proportional to `M_B`, evaluated
on a **genome-wide** (not windowed) intercept — but that estimand is already owned
by Sib-Regression (Visscher 2006) and RDR (Young 2018), and testing it honestly
needs a per-meiosis-recombination simulator this repo lacks. This closes candidate
option (a) of DEVELOPMENT_PLAN §1.1.


## `identifiability_aliasing.py` — the structure-aliasing limit, illustrated

Backs the METHOD §14 "structure aliasing" remark (adversarially verified). Marginal
uncorrected `z_j = x_jᵀy/√N` under three phenotype models that all load variance on
the same top structure axis `u₁`; regress `z_j²` on `[1, corr(x_j,u₁)²,
corr(x_j,u₁)²·(ℓ_j^local − mean)]` (centered local factor to orthogonalize).

| model | mean χ² | b[corr²] | t[corr²×local] |
|---|---|---|---|
| (A) environmental axis shift | 2.49 | 28.73 | 5.29 |
| (B) directional genetic (rank-1) | 2.49 | 28.73 | 5.29 |
| (C) heteroscedastic-along-axis | 2.50 | 20.26 | 8.48 |

### Verdict: aliasing is exact for the directional component, not the LD-borne one

(A) and (B) are **identical to 2 dp** — a genuine *genetic* effect producing a
phenotype along the ancestry axis is indistinguishable from an *environmental* shift
along it (the shared nonzero interaction is a finite-sample orthogonalization
residual, identical in both, so not a distinguishing feature). (C) heteroscedastic
genetic signal loads less on pure `corr²` and leaves a larger `corr²×ℓ^local`
interaction — it rides on local LD, so that part *is* identifiable via the local
channel. This is the exact scope of the limit: the rank-1 directional component is
aliased (every structure-only correction co-attenuates it); the LD-borne component
is not.



## `local_adaptation_channel.py` — kill-test for the one surviving direction

The aliasing limit (METHOD §14) leaves one identifiable channel for ancestry-aligned
signal: a *heteroscedastic-along-axis* architecture (per-SNP h² elevated for
high-loading SNPs — what polygenic local adaptation along an ancestry cline should
look like) leaves a `corr²×ℓ^local` interaction an environmental axis shift lacks.
This kill-test asks whether that channel has **power** at wild-population-study
scale (N = 600/1200, discrete populations, in-sample PCs and LD scores — no
reference panel). Three amplitude-matched arms per realization: ENV (rank-1
environmental shift; the null), DIR (rank-1 directional genetic; the aliased
control), HET (the local-adaptation architecture; the alternative). Thresholds are
empirically calibrated on the ENV arm (95th percentile, one-sided). Pre-registered
verdict rule (fixed before running): PASS = power ≥ 80% at N=600, χ²=2.5, κ=1;
KILL = power < 50% at χ²=4.0, κ=2; else MARGINAL. PCG64 seed 20260719, 16–24
genotype realizations per setting, 2000 phenotype replicates per arm.

| setting | med t ENV | med t DIR | med t HET | thr95 | power |
|---|---|---|---|---|---|
| N600 χ²=2.5 κ=0.5 | 3.77 | 3.77 | 4.06 | 6.71 | 0.17 |
| N600 χ²=2.5 κ=1 | 3.77 | 3.77 | 6.18 | 6.71 | 0.38 |
| N600 χ²=2.5 κ=2 | 3.77 | 3.77 | 15.11 | 6.71 | **1.00** |
| N600 χ²=1.5 κ=1 | 4.03 | 4.03 | 5.93 | 6.53 | 0.33 |
| N600 χ²=4.0 κ=1 | 3.46 | 3.46 | 6.18 | 6.65 | 0.38 |
| N1200 χ²=2.5 κ=1 | 3.94 | 3.94 | 6.97 | 6.59 | 0.75 |
| N600 χ²=4.0 κ=2 | 3.46 | 3.46 | 15.57 | 6.65 | 1.00 |

### Verdict: MARGINAL (by the pre-registered rule) — and a null-calibration warning

- **The channel has power only for strong architectures or larger samples**: power
  1.00 at κ=2, 0.75 at N=1200, but only 0.38 at the moderate setting (N=600, κ=1).
  Doubling N clearly helps (0.38 → 0.75); the trend says the channel is real but
  weak at realistic wild-study sizes and gradient strengths.
- **The naive interaction test is badly anti-conservative**: the ENV null itself
  has median t ≈ 3.8 and exceeds t>2 in 88–100% of realizations — centering
  `ℓ^local` does not remove the collinearity-driven leakage of the large `corr²`
  term into the interaction coefficient. Any real test on this channel needs a
  calibrated null (permutation / parametric bootstrap), which will cost further
  power. (This mirrors the signed-version anti-conservatism in
  `detect_pc_inflation.py`.)
- **DIR sits exactly on ENV in every setting** (identical medians) — the aliasing
  limit holds at the decision level, not just in the moment.
- Consequence: a within-sample local-adaptation test on this channel is **not
  validated** for moderate architectures at N≲1000. It is not dead (strong
  gradients, or N≥1200 with a calibrated null), but it cannot be the next
  milestone without a better-powered statistic or a calibrated-null design.


## `misspecified_orthogonal.py` — axes absent from the GRM (the untested case)

`misspecified_background.py` found the oracle intercept robust (<1%) to sparse and
directional off-chromosome backgrounds — but both live in the structure subspace
`P_c` already whitens. This tests the harder, previously untested case: true-
covariance structure **absent from the GRM** used to build `V_c`. Exact (no MC),
6 realizations, PCG64 seed 20260719.

| arm | setting | mean null intercept |
|---|---|---|
| (a) batch axis orthogonal to ancestry (shared microenvironment; projection onto top-10 GRM eigenvectors = 0.016) | v/σ_e² = 0.25 → 2.0 | 1.000 → 1.003 |
| (b) GRM from every 10th off-chr SNP (LD-thinning practice); truth uses full GRM | — | **3.185** |
| (c) family-shared environment + FULL GRM (100 sibships ×4) | v/σ_e² = 0.25 → 2.0 | 1.091 → 1.731 |
| (c) family-shared environment + PRUNED GRM | v/σ_e² = 0.25 → 2.0 | **4.383 → 7.687** |
| (c) oracle (`P_c` from `V_true`) | all | 1.000 ✓ |

### Verdict: fragility lives in GRM coverage, not in orthogonal axes

- **(a) An axis orthogonal to the genotype structure axes does nothing**, even at
  v = 2σ_e²: `P_c` does not down-weight it, but genotypes have no loading on it, so
  `x_j^T P_c u ≈ 0`. Confirms the precise mechanism of the robustness found in
  `misspecified_background.py`: inflation requires the missed axis to be correlated
  with *genotype* structure.
- **(b) A thin GRM is the realistic failure mode.** Building `V_c` from every 10th
  SNP — standard LD-thinning practice — triples the null intercept (3.19), because
  the unmodeled SNPs carry real covariance the assumed `V_c` cannot see. The
  intercept's oracle calibration depends on the GRM *spanning* the true covariance,
  not just on getting two variance-component scalars right.
- **(c) Family-shared environment** is partly whitened when the full GRM captures
  the sibship axes (1.09–1.73 with a full GRM) but severe under a pruned GRM
  (4.4–7.7). Together with (b): the intercept is only as good as the GRM's coverage
  of the sample's actual covariance structure.
- Practical upshot for METHOD §6: the "structure axis absent from the GRM" case is
  now tested. It is benign when genuinely orthogonal to genotypes, and damaging
  exactly when the GRM under-covers the genetic/family covariance — i.e. it folds
  into the known fitted-`V_c` sensitivity (M4a) rather than adding a new axis of
  concern, but it prices in GRM density as a first-order requirement.

## `realdata_1000g_structure.py` — the structure fraction of LD on real data

REAL_DATA_SURVEY question 1, smallest version: chr22 of 1000 Genomes phase 3
(20130502; biallelic SNPs, full-panel MAF ≥ 0.05, 1-in-25 thinning, ~3.6–3.9k
cohort-common SNPs at ~13 kb spacing; `.scratch_1000g/`, not committed). Ordinary
vs `P_c`-whitened LD scores (in-chromosome GRM — NOT LOCO, chr22 is the whole
dataset; window = ±50 SNP index, no genetic map). `ell_struct = ell_ord −
ell_local` by construction; the empirical content is its correlation with the
top-10 PC score.

| cohort | N | λ1 | struct share, genome-wide | struct share, windowed | corr(ℓ^struct, PC score) full / windowed |
|---|---|---|---|---|---|
| EUR | 503 | 5.3 | 28% | 19% | 0.75 / 0.80 |
| AFR | 661 | 13.0 | 35% | 20% | 0.93 / 0.49 |
| FULL panel | 2504 | 283.1 | **96%** | **62%** | **1.00** / 0.79 |

> **Correction (2026-07-19, see `realdata_1000g_benchmarks.py` B2 below):** these
> shares used an **in-chromosome GRM**, which removes local LD along with
> structure and thereby overstates the structure share. The LOCO-clean numbers
> are: genome-wide ~0–5% (EUR) / ~90% (FULL); windowed ~0% (EUR) / ~17% (FULL).
> The corr(ℓ^struct, PC score) findings are unaffected.

### Verdict: the framework is real; the sim's quantitative "~5% windowed" is not

- **The qualitative story reproduces on real data**: the LD score decomposes, the
  non-local remainder tracks the leading-PC score (corr 0.93–1.00 genome-wide),
  and the structure fraction grows steeply with actual structure strength — 28% →
  35% → **96%** of the genome-wide score from EUR to the multi-ancestry panel,
  tracking λ1 (5.3 → 13 → 283). In a real multi-ancestry sample the genome-wide
  LD score is almost entirely structure, and whitening removes 96% of it.
- **But the windowed share is NOT ~5% on real data at this density**: 19–20% in
  single-ancestry cohorts and 62% in the full panel. The sim's 5% came from dense
  simulated SNPs with strong local LD; here thinning leaves ~13 kb spacing where
  pairwise local LD is weak, so structure dominates even within a ±50-SNP window.
  The honest statement: the windowed structure share is **density- and
  structure-dependent**, and in real multi-ancestry panels it is large. (This
  cuts both ways for the project: it strengthens "structure matters in windows
  too" for real panels, while the sim's "windows are safe" comfort is gone —
  consistent with cov-LDSC mattering most in exactly such panels.)
- Caveats: one chromosome; in-chromosome GRM (not LOCO); SNP-index window; the
  `P_c` metric assumes h²=0.5; GT parsing treats any non-0/1 genotype as
  homozygous reference (phase-3 calls are essentially complete, so immaterial).
- Next steps if pursued: all chromosomes with a LOCO GRM and a recombination map
  (cM window), denser common-SNP sets (less thinning), and the EUR→admixed→
  multi-ancestry gradient at matched density. Compute is linear in chromosomes.

## `realdata_1000g_benchmarks.py` — thorough decomposition benchmarks (B1–B6)

Follow-up to the chr22 pilot, closing its caveats: chr22 (1-in-10 thin) + chr21/10
(1-in-25, GRM fodder), cached parses; **normalized correlations throughout**
(the cov-LDSC convention — an early version of this script used unnormalized
cross-products, which inflated the whitened scores; fixed and re-run). Windows
are physical (kb). All numbers below are on real 1000G phase-3 data.

**B1 — density × window (in-chromosome operator).** Windowed structure share is
robust to SNP density (EUR: 26% dense vs 27% sparse at ±250 kb; FULL: 63% vs 52%)
— the pilot's windowed numbers were not a thinning artifact.

**B2 — LOCO vs in-chromosome GRM: the operator matters more than the data.**
The pilot whitened with an in-chromosome GRM; this benchmark adds a proper LOCO
GRM (chr21+10). Genome-wide / ±250 kb structure shares:

| cohort | in-chrom full / win | LOCO full / win |
|---|---|---|
| EUR | 33% / 26% | **−4% / −0%** |
| FULL | 96% / 63% | **90% / 17%** |

The in-chromosome operator removes **structure plus local LD** (its middle
eigenspace spans the target chromosome's own haplotype background — for EUR it
whitens the genome-wide score down to the M/N noise floor), so the pilot's
shares overstated structure, especially in windows. The LOCO operator cannot
touch local LD, so it isolates structure — but slightly *inflates* residual
scores (up-weighting of sub-1 eigenvalues + the q-renormalization of METHOD §5),
visible as the small negative EUR shares. The two operators bracket the truth:
**the honest structure share of the genome-wide LD score is ~90% in a
multi-ancestry panel and only ~0–5% in single-ancestry cohorts; in ±250 kb
windows it is ~17% (FULL) vs ~0% (EUR).** The sim's "5% in windows" was the
right order for single-ancestry cohorts, wrong for multi-ancestry panels — and
the pilot's "19–62% in windows" overstated both.

**B3 — PC-projection vs P_c (FULL): cov-LDSC needs very few PCs on real data.**
Exact spectral capture of the genome-wide score: **top-1 eigenvector 79%,
top-2 92%, top-10 95%** (50% in top-1, 90% in top-2). Hard-projecting top-2 PCs
leaves only 8% of what P_c removes (top-10: 5%; top-50: 5% genome-wide but 41%
in windows — GLS down-weights the mid-spectrum that hard projection leaves
behind). Real-data support for the project's central negative result: on real
multi-ancestry genotypes the GLS operator ≈ cov-LDSC with a handful of PCs.

**B4 — metric sensitivity.** Structure share vs the assumed h² in P_c: EUR
15%/33%/25% (full) at h²=0.1/0.5/0.9 — NOT insensitive (the pilot's claim was
wrong for mild-structure cohorts); FULL is stable (93–96%).

**B5 — MAF gradient (in-chrom convention; relative gradient robust).** Structure
share rises with MAF: EUR 22% → 27% → 30% across [0.01–0.05, 0.05–0.2, 0.2–0.5];
FULL 88% → 90%. Common SNPs carry proportionally more structure signal (they are
older and more ancestry-informative).

**B6 — ancestry gradient (matched SNP set; in-chrom convention).**

| cohort | N | λ1 | struct% full | struct% ±250kb | corr(ℓ^struct, PC10) |
|---|---|---|---|---|---|
| EUR | 503 | 5.3 | 33% | 26% | 0.885 |
| AMR (CLM/MXL/PUR) | 262 | 8.4 | 34% | 19% | 0.976 |
| AFR | 661 | 13.2 | 38% | 27% | 0.959 |
| FULL | 2504 | 283.0 | 96% | 63% | 0.994 |

The gradient tracks λ1 as expected; the corr(ℓ^struct, PC score) ≈ 0.9–1.0
finding of the pilot is robust on real data.

**Caveats.** Three chromosomes, thinned SNP sets; in-chromosome vs LOCO bracket
the truth (B2) — future whole-genome runs should use LOCO throughout; MAF
gradient (B5) is in-chrom convention; windows are physical (no recombination map
yet); the P_c metric is intercept-only.

## `relatedness_effective_df.py` — the residual-df correction under relatedness

Post-archive (2026-07-23). Section 13 says relatedness "does not propagate into
the local, windowed LD-score regression." True of the *slope* — but there is one
channel it does reach, and the implementation had it wrong. The Eq. (7)
adjustment `r² − (1−r²)/(ν−1)` subtracts the noise floor two *unlinked* variants
generate by chance, and that floor is `1/(N−1)` only for independent rows. With
relationship matrix `A` and LD covariates `C`, it is `tr(Ã²)/tr(Ã)² = 1/N_eff`
for `Ã = P_C A P_C` (METHOD §10.1, Eq. M22).

Real Mendelian transmission (founder haplotypes, one transmitted per parent per
SNP), unlinked SNPs by construction, N=400 throughout, 8 replicates, PCG64 seed
20260723. Predictions come from the shipped `pldsc.effective_pair_df`, so this
validates the released code path.

| design | N_eff | observed floor | naive 1/(ν−1) | tr-ratio | obs/naive | obs/tr |
|---|---|---|---|---|---|---|
| unrelated (control) | 399.0 | 0.002504 | 0.002513 | 0.002506 | 1.00 | **1.00** |
| 100 fam × 4 sibs | 227.8 | 0.004384 | 0.002513 | 0.004391 | 1.74 | **1.00** |
| 50 fam × 8 sibs | 144.9 | 0.006812 | 0.002513 | 0.006903 | 2.71 | **0.99** |
| 25 fam × 16 sibs | 83.9 | 0.011752 | 0.002513 | 0.011925 | 4.68 | **0.99** |

### Verdict: correct, and an intercept effect — not a heritability result

- **(M22) is accurate to ≤1%** across designs spanning a **4.7× inflation** of
  the noise floor, while `ν−1` is wrong by that same 4.7×. The unrelated control
  reproduces `1/(N−1)`.
- Accumulated over a window of `S` variants the omitted floor reaches **≈19
  LD-score units** at S=2000 in the 16-sib design — large against a real local
  score.
- **But it is near-uniform across variants at fixed S**, so to first order it
  shifts every `ℓ_j` by a constant: it moves the LDSC **intercept**, not the
  slope. h² is protected; an intercept read as a confounding barometer on a
  related in-sample panel is not. Residual non-uniformity enters only through
  variation in variant density per window.
- **Does not revive the relatedness case for `P_c`** (§13). This channel is an
  intercept effect whose standard remedy is to compute LD from unrelated
  individuals (which sets `N_eff = N−1` and makes the correction moot); the
  biases `P_c` was built for are slope effects. (M22) matters only when dropping
  relatives is not an option — family or founder cohorts with in-sample LD.

Exposed as `pldsc.effective_pair_df(A, covariates)` and
`local_ld_scores(..., pair_df=...)`; 9 unit tests in `tests/test_ld.py`.

## `realdata_genome_loco.py` — genome-wide LOCO run (the #2 engine, end to end)

Post-archive (2026-07-23). Scales the B2/B6 benchmarks from three thinned
chromosomes to a genuine leave-one-chromosome-out sweep, on the tested
`pldsc.genome` engine, with genetic-distance (cM) windows and the exact
`O(M N²)` spectral genome-wide score (METHOD §14 Eq. D2 — no `M × M` matrix, so
it scales to millions of SNPs). This is the shared infrastructure for
`docs/WHEN_STRUCTURE_MATTERS.md` (#6) and the per-axis study (#4).

Engine validation (against the brute-force benchmark, real chr22): the spectral
full score matches to **3e-13**, the cM windowed score matches the physical
benchmark exactly, and the LOCO shares reproduce **B2** (EUR −0.2%/−4.4%, FULL
17.3%/89.6% vs the published −0%/−4% and 17%/90%). 21 unit tests in
`tests/test_genome.py`.

Rotating-LOCO demonstration on the cached three chromosomes (target 22 ← 21+10
fodder, 21 ← 22+10, 10 ← 22+21), 1.0 cM windows, `P_c` at h²=0.5, SNP-weighted
genome-wide aggregate:

| cohort | M | struct% (1 cM window) | struct% (genome-wide) |
|---|---|---|---|
| EUR | 25 308 | −0.7% | −9.2% |
| AMR | 27 395 | 1.2% | −0.2% |
| AFR | 26 362 | 1.5% | −2.0% |
| FULL | 27 600 | **29.3%** | **90.0%** |

### Verdict: the machinery scales; the genome-wide headline needs all 22 autosomes

- The multi-ancestry genome-wide share is **90.0%**, matching B2/B6; the
  windowed cM share (**29%** at 1 cM) is higher than the 250 kb physical share
  (17%) because a 1 cM window is wider than 250 kb on average — the window unit
  matters, and cM is the recombination-aware choice.
- Single-ancestry shares are ~0% windowed and slightly **negative** genome-wide
  — the documented B2 LOCO over-whitening / q-renormalization inflation for
  low-structure cohorts, faithfully reproduced.
- **Three chromosomes is a demonstration, not the headline.** The genome-wide
  run needs all 22 autosomes: `python experiments/realdata_genome_loco.py --plan`
  prints the manifest (~15 GB phase-3 VCFs; maps and panel already local). Gated
  on that download — a large resource commitment, so it is deliberately not run
  automatically.
