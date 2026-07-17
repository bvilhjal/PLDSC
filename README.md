# PLDSC

PLDSC is a clean restart of the 2017 `newldsc` research prototype. The project
studies LD Score regression when GWAS samples contain population structure or
when association statistics come from covariate-adjusted or mixed models.

The project is currently in the specification and validation phase. It does
not yet provide a production-ready estimator.

## Current direction

1. Reproduce covariate-adjusted LD scores with a small, testable Python 3 core.
2. Derive and test LD scores matched to leave-one-chromosome-out mixed-model
   association statistics.
3. Treat the legacy GRM-whitening and aggregated PC-score methods as
   experimental comparators, not established estimators.

See [the scientific idea review](docs/IDEA_REVIEW.md) for the rationale,
failure modes, and validation gates.

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
