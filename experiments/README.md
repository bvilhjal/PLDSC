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

**Readout B — Eq (3) is numerically correct.** Regressing the Monte-Carlo
`E[z_j^2]` on the full-chromosome `L_j^(P_c)` recovers the oracle line
`(intercept, slope) = (0.997, 6.246e-4)` against the true `tau^2 = 6.229e-4`
(structured; homogeneous 0.983 / 6.111e-4 vs 5.86e-4). Per-SNP standardized
residuals have RMS ~1.55 sigma — expected, because the 80 SNP means share the
same MC draws and are correlated; the aggregate intercept pinned to 1 rules out
a systematic error. Window truncation (1.5 cM) drops ~12% of the full sum
(median), the expected Eq-4 approximation to quantify in M3.

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

Mean null chi^2 (want 1.000):

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

Uses the tested `pldsc.regression` (weighted LS + delete-a-block jackknife) and
the `pldsc._lmm` oracle. Two arms per replicate on an admixture+relatedness
cohort (N=400, 80 target SNPs, R=300): an oracle-`V_c` arm and a REML-fitted-`V_c`
arm. Coverage of the intercept (truth 1) and slope (truth `tau^2`):

| arm / SE | intercept coverage | slope coverage |
|---|---|---|
| oracle `V_c`, jackknife | 94.3% | 95.3% |
| fitted `V_c`, jackknife only | 95.0% | 94.3% |
| fitted `V_c`, jackknife + `V_hat_c` term | 95.0% | — |

### Verdict: jackknife is calibrated at a fixed operator; the `V̂_c` caveat is a scale effect

- The block jackknife gives nominal coverage in both arms, so the regression
  inference is sound at a fixed operator.
- **The nuance M4 had hidden:** M4's stable "mean inflation" and M5's
  extrapolated "regression intercept" are different estimators. At only M=80
  target SNPs the regression intercept is noise-dominated (SE ≈ 1.63), so the
  `V̂_c` term (0.145) is negligible and coverage is fine with or without it.
- **But `V̂_c` uncertainty is a *shared systematic*** — one fitted `V_c` shifts
  every SNP's statistic together — so it does **not** shrink with the number of
  SNPs, whereas the regression SE does (~`1/sqrt(M)`). The decomposition
  (regression SE 1.63 at M=80 vs `V̂_c` SE 0.145) puts the crossover at
  **~10,000 SNPs**: beyond that the jackknife-only intercept SE under-covers, so
  a genome-scale M5 must propagate variance-component uncertainty into the
  intercept SE. The slope, far less sensitive to `V_c`, stays calibrated
  throughout.

This corrects the earlier framing: the fitted-`V_c` intercept concern is real
but is a large-`M`/at-scale effect, and the fix (propagate `V̂_c` uncertainty,
e.g. via `dIntercept/dh^2 x SE(h_hat^2)` added in quadrature) is now a concrete,
quantified requirement for the M5 estimator.


