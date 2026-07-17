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
