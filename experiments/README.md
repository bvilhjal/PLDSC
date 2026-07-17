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
