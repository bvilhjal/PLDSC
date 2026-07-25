# PLDSC

PLDSC is a clean restart of the 2017 `newldsc` research prototype. The project
studies LD Score regression when GWAS samples contain population structure or
when association statistics come from covariate-adjusted or mixed models.

## Status: CONCLUDED / archived (2026-07-19)

**Outcome: the machinery is correct and well-tested, but an exhaustive search
found no differentiated contribution, so the project is archived as an honest
record rather than pursued to publication.** What the search covered and found:

- The **LMM-matched LD score as a new estimator** reduces to cov-LDSC in the
  summary-statistics LDSC setting (structure is what biases LDSC, and cov-LDSC
  already fixes it; relatedness does not reach the windowed regression). See
  [docs/METHOD.md](docs/METHOD.md) §13, [map_ldsc_miscalibration.py](experiments/map_ldsc_miscalibration.py).
- The **within-family / sib-GWAS** angle was kill-tested and failed: the surviving
  within-family genetic covariance is a homoscedastic scale (`M_B R M_B = 0.5 M_B`)
  that a standard sib-GWAS's cluster-robust / effective-N SE already removes
  ([within_family_killtest.py](experiments/within_family_killtest.py)).
- The **theory** (local/non-local decomposition, unification, identifiability
  limit) is real but textbook-adjacent and largely occupied: the identity is a GRM
  quadratic form; the GLS operator is Mangin 2012; the aliasing corollary is
  Huang 2025 / Sotoudeh 2026; unification-of-estimators is MQS/GEMMA. Full
  adversarial prior-art sweep in [docs/PRIOR_ART.md](docs/PRIOR_ART.md).
- A structured **contribution-ideation** pass (15 candidates, adversarially
  vetted) produced one high-death-prior long-shot (a within-sample
  local-adaptation test for wild-population studies) and otherwise confirmed the
  wall. The decision was to stop.

**Durable value of this repo:** (1) a small, tested, CI-green implementation of the
mixed-model residual-precision operator `P_c` and the LD-score decomposition
(`pldsc/`, 59 tests); (2) a verified theory write-up ([docs/METHOD.md](docs/METHOD.md),
including the adversarially-checked structure-aliasing remark in §14); (3) an
unusually candid **negative-results and prior-art record** (`experiments/`,
[docs/PRIOR_ART.md](docs/PRIOR_ART.md), [docs/EXTERNAL_REVIEW.md](docs/EXTERNAL_REVIEW.md),
[docs/NOTE_OUTLINE.md](docs/NOTE_OUTLINE.md)) documenting precisely where a natural
idea does and does not buy anything. If revived, the honest ceiling is a short
clarifying/negative-results note ([docs/NOTE_OUTLINE.md](docs/NOTE_OUTLINE.md)),
crediting Mangin 2012 / Huang 2025 / cov-LDSC / de Vlaming 2017 throughout.

---

## The science (record)

**Central idea: compare operator-defined local/non-local LD-score components.**
The total genome-wide correlation score has an exact spectral GRM identity,
`Σ_k r_jk² = (M/N) Σ_p λ_p corr(x_j,u_p)²` (textbook algebra, not a new
finding). Calling a component "structure" additionally requires a declared
window, operator, or justified PC subspace; the split is not unique. The legacy
"PC score" is a leading-eigenvector truncation — which correlates 0.99 with the
declared non-local remainder but captures only
~37% of its magnitude ([experiments/ld_decomposition.py](experiments/ld_decomposition.py)).
The honest headline (simulation): within the **windowed** LD score LDSC actually
uses, structure is only **~5%** (it is ~50% genome-wide, but LDSC never regresses
on that) — which is precisely *why* cov-LDSC suffices and the new-estimator angle
died. **On real 1000G data** the LOCO-clean windowed share is ~0% (EUR) to ~17%
(multi-ancestry panel) and the genome-wide share ~90% in the panel
([realdata_1000g_benchmarks.py](experiments/realdata_1000g_benchmarks.py)). The value of the identity is as an **organizing device**: cov-LDSC, the
legacy whitened-LD (WhiLD) score, the PC score, and the mixed-model
residual-precision operator `P_c` can be compared on one metric/eigenspace map.
They do not all remove or measure one unique `ℓ^struct`, and ancestry-aligned
genetic signal can be aliased with environmental confounding. See
[docs/METHOD.md](docs/METHOD.md).

**Post-archive theory addendum (2026-07-22).** [Section 15 of the method
note](docs/METHOD.md#15-theory-addendum-assortative-mating-and-genetic-nurture)
now records two complementary limits. Assortative mating creates signed,
trait-aligned cross-autosome GPD that a squared LD score cannot identify as mate
choice; under random mating and the usual isotropic-effect assumptions, pure
parental genetic nurture can produce a textbook LDSC slope even when direct SNP
heritability is zero. Offspring GWAS plus offspring LD alone cannot causally
identify mate choice or parental nurture. This synthesis extends the
interpretation/identification record; it does not revive the archived \(P_c\)
estimator claim.

**Post-archive correctness fix (2026-07-23).** [METHOD §10.1](docs/METHOD.md#101-relatedness-in-the-ld-panel-the-residual-df-correction)
records one channel through which relatedness *does* reach the LD score, and
which the implementation had wrong: the Eq. (7) noise-floor subtraction
`(ν−1)` assumes independent rows and under-subtracts by up to 4.7× on a
related panel. The correct denominator is `tr(Ã)²/tr(Ã²)`
(`pldsc.effective_pair_df`, verified in
[relatedness_effective_df.py](experiments/relatedness_effective_df.py)). The
resulting LD-score bias is large but near-uniform, so it moves the LDSC
**intercept**, not the slope: h² is protected, an intercept-based confounding
readout on a related in-sample panel is not. This does **not** revive the
relatedness case for `P_c` (§13) — it is an intercept effect whose standard
remedy is to compute LD from unrelateds.

**Project-direction review (2026-07-25).** [docs/PROJECT_DIRECTIONS.md](docs/PROJECT_DIRECTIONS.md)
reviews the whole record and asks what a *differentiated* revival would look like.
Diagnosis: every framing that died was a variant of "a better estimator of `h²`
under structure," which is bounded by the windowed structure share this repo itself
measured; the claims that survived four prior-art sweeps were all *measurements and
corrections*, never constructions. It proposes four directions with kill-tests, two
of which cost a day each and have not been run. The highest-ranked one is backed by
a new experiment ([raw_loading_control.py](experiments/raw_loading_control.py)):
S-LDSC's annotation LD scores include the `r_jj² = 1` self term, which makes the raw
ancestry loading and its LD-smeared score 0.82-collinear — so standard S-LDSC
returns a `τ` **~9× larger under pure environmental confounding than under a genuine
4× enrichment** on an ancestry-correlated annotation, and adding the *raw*
(un-smeared) loading removes it at ~no cost to genuine enrichment elsewhere. The
§14 identifiability limit is unchanged: ancestry-correlated real enrichment is still
co-attenuated.

**Proposed successor note (2026-07-23).** [docs/WHEN_STRUCTURE_MATTERS.md](docs/WHEN_STRUCTURE_MATTERS.md)
outlines a practice-facing alternative to the archived theory note: lead with
the empirical magnitudes (windowed structure share ~0% single-ancestry vs ~17%
multi-ancestry) and a two-moment decision rule (relatedness → intercept,
structure → slope) rather than with theory claims the prior-art sweep graded
partial-or-occupied. Gating task is a prior-art check on that framing.

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
surviving value is the explicit synthesis, diagnostics, and identifiability
discussion—not novelty of the textbook spectral identity. See
[docs/PRIOR_ART.md](docs/PRIOR_ART.md).

**Validated machinery.** The `P_c` operator and the fixed-operator conditional
mean are exact. The cov-LDSC correlation geometry is an exact limit, while
equality of full cov-LDSC/ordinary regression moments requires the scale/prior
conditions now stated in
[docs/METHOD.md](docs/METHOD.md), unit-tested (`pldsc/_lmm.py`), and confirmed
across the `experiments/` arc. The earlier "LMM-matched LD score as a new
*estimator*" angle was tested and reduces to cov-LDSC in LDSC application
([docs/METHOD.md](docs/METHOD.md) §13); the decomposition reframing is what
survives and now organizes the project. The legacy WhiLD and PC-score methods are
related metric/eigenspace constructions, not literal special cases under the
same effect prior and not standalone estimators; see
[the scientific idea review](docs/IDEA_REVIEW.md).

**Real data (1000 Genomes, post-archive).** The decomposition was tested on real
genotypes (phase-3 chr22/21/10; EUR/AMR/AFR/full panel): the structure share of
the genome-wide LD score tracks λ1 (~0–5% single-ancestry, ~90% multi-ancestry,
LOCO-clean; an in-chromosome GRM over-whitens), ~2–10 hard-projected PCs match
the GLS operator (real-data support for the reduces-to-cov-LDSC verdict), and the
non-local component tracks the PC score at corr 0.9–1.0. A kill-test of the one
surviving application direction (a within-sample local-adaptation test) came back
MARGINAL, and the intercept's robustness boundary is GRM coverage. See
[experiments/README.md](experiments/README.md) and
[docs/REAL_DATA_SURVEY.md](docs/REAL_DATA_SURVEY.md).

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
