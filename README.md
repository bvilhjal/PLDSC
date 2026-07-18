# PLDSC

PLDSC is a clean restart of the 2017 `newldsc` research prototype. The project
studies LD Score regression when GWAS samples contain population structure or
when association statistics come from covariate-adjusted or mixed models.

The core theory is derived, validated, and tested; the open question is whether
it has a differentiated application. It does not provide a production-ready
estimator.

## Status (2026-07-18)

**Theory: done and validated.** The LD score matched to a leave-one-chromosome-
out mixed-model score statistic is derived in [docs/METHOD.md](docs/METHOD.md)
(the `P_c` operator and the exact moment `E[z_j^2] = 1 + Σ_k τ_k² q_k
r_jk^(P_c)²`), unit-tested (`pldsc/_lmm.py`), and numerically confirmed across a
full experiment arc (`experiments/`): it reduces to cov-LDSC and ordinary LDSC,
calibrates exactly, beats cov-LDSC on relatedness *in the association step*, and
survives a REML-fitted `V_c`.

**Application: the LDSC-application niche reduces to cov-LDSC.** A direct test
([map_ldsc_miscalibration.py](experiments/map_ldsc_miscalibration.py)) showed
the miscalibration of LDSC on mixed-model summary statistics is a *structure*
effect that cov-LDSC (Luo et al. 2021) already corrects; relatedness — PLDSC's
one edge — does not propagate into the local LD-score regression. So PLDSC offers
no differentiated advantage in LDSC application. Its remaining candidate niche is
**within-family / sib-GWAS designs** (the association step), which is proposed
but not yet validated. See [docs/METHOD.md](docs/METHOD.md) §13 and
[docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md) §1.1 for the honest scope.

The legacy GRM-whitening (WhiLD) and PC-score methods are treated as named
comparators, not established estimators; see
[the scientific idea review](docs/IDEA_REVIEW.md).

## First numerical slice

The initial API calculates local LD scores after projecting the same fixed
covariates used by a GWAS out of a sample-by-variant genotype matrix:

```python
from pldsc import local_ld_scores

scores = local_ld_scores(
    genotypes,
    positions_cm,
    covariates=principal_components,
    window_cm=20.0,  # illustrative; select by an LD-score convergence check
)
```

It operates on one chromosome at a time, includes the LD Score self-term, and
retains negative finite-sample-adjusted terms. It is intentionally not an LMM
or LOCO implementation. Covariates are treated as fixed: when PCs are estimated
from the same genotypes, the residual-rank finite-sample correction is only
approximate. The appropriate window is population- and data-dependent; 20 cM
is not a universal default.

Run the tests with:

```text
python -m unittest discover -s tests -v
```
