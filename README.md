# PLDSC

PLDSC is a clean restart of the 2017 `newldsc` research prototype. The project
studies LD Score regression when GWAS samples contain population structure or
when association statistics come from covariate-adjusted or mixed models.

The central idea is a **local/non-local decomposition of the LD score**, with the
structure component and its diagnostics as the contribution. It is a
framework/synthesis project, not a production estimator.

## Status (2026-07-18)

**Central idea: the LD score decomposes into linkage + structure.** The LD score
splits as `ℓ_j = ℓ_j^local + ℓ_j^struct` — recombination-driven linkage plus a
population-structure component — and the non-local component **is the "PC score"**
(the eigenvalue-weighted sum of squared SNP–PC loadings), validated at
**corr = 0.99 genome-wide** ([experiments/ld_decomposition.py](experiments/ld_decomposition.py);
structure is ~50% of the genome-wide LD score but only ~5% within a local window).
One identity **unifies the existing structure-corrections**: covariate-adjusted
LDSC (cov-LDSC), the legacy whitened-LD (WhiLD) score, the PC score, and the
mixed-model residual-precision operator `P_c` are all ways of removing or
measuring `ℓ^struct`. Standard LDSC's structure bias is exactly the conflation of
`ℓ^local` (heritability) with `ℓ^struct` (confounding). See
[docs/METHOD.md](docs/METHOD.md).

**Applications — structure diagnostics** (all prototyped in `experiments/`):
a per-PC inflation spectrum that maps a GWAS's confounding to specific axes
([detect_pc_inflation.py](experiments/detect_pc_inflation.py)); a
structure-robust S-LDSC enrichment test via PC-loading annotations
([sldsc_pc_annotation.py](experiments/sldsc_pc_annotation.py)); and a map of when
LDSC on mixed-model sumstats is affected
([map_ldsc_miscalibration.py](experiments/map_ldsc_miscalibration.py)). All share
one **stated limit**: structure confounding and genuine ancestry-aligned signal
load on the same PCs, so no structure-only method can separate them.

**Honest positioning.** Framework/synthesis, not a new estimator. The qualitative
concept is established (cov-LDSC, Bhatia 2016, Price 2008, Bercovich 2025); the
novel slices are the explicit decomposition identity, the unification, and the
diagnostics + identifiability limit. See [docs/PRIOR_ART.md](docs/PRIOR_ART.md).

**Validated machinery.** The `P_c` operator, the exact moment, the reductions to
cov-LDSC/ordinary LDSC, and the diagnostics are derived in
[docs/METHOD.md](docs/METHOD.md), unit-tested (`pldsc/_lmm.py`), and confirmed
across the `experiments/` arc. The earlier "LMM-matched LD score as a new
*estimator*" angle was tested and reduces to cov-LDSC in LDSC application
([docs/METHOD.md](docs/METHOD.md) §13); the decomposition reframing is what
survives and now organizes the project. The legacy WhiLD and PC-score methods are
special cases of the decomposition, not standalone estimators; see
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
