# PLDSC development plan

Date: 2026-07-17 (revised 2026-07-18 after multi-lens review and the
application investigation)

## 1. Objective, use-case, and decision rule

PLDSC tests whether an LD score matched to a specified leave-one-chromosome-out
(LOCO) linear mixed-model (LMM) association statistic yields *calibrated* LD
Score regression for quantitative traits. It is not a general LDSC rewrite, and
covariate-residualized LD is a baseline rather than the novelty claim. The
central identity (Section 2) is already proved and independently checked, so the
project's risk was always delivery and use-case, not the theorem. As of
2026-07-18 both the theorem and the full MVP arc (M-spike–M5) are validated, and
the use-case has been tested directly: the LDSC-application niche
**reduces to cov-LDSC** (Section 1.1), so PLDSC's differentiated value, if any,
is upstream. Read Section 1.1 before the roadmap.

**Reframing (2026-07-18): the central idea is the local/non-local LD-score
decomposition.** \(\ell = \ell^{\text{local}} + \ell^{\text{struct}}\). The
genome-wide LD score is *exactly* the GRM spectral decomposition
\(\sum_k r_{jk}^2=(M/N)\sum_p\lambda_p\operatorname{corr}(x_j,u_p)^2\) (textbook
algebra); the legacy PC score is its leading-eigenvector truncation (corr 0.99
with the non-local part but only ~37% of its magnitude; see [METHOD.md](METHOD.md)
§14). The honest number is that structure is only ~5% of the *windowed* LD score
LDSC uses (~50% genome-wide, which LDSC never regresses on) — which is why
cov-LDSC suffices. The identity's value is as an **organizing device** unifying
cov-LDSC, the legacy WhiLD score, the PC score, and the \(P_c\) operator as ways
of removing or measuring \(\ell^{\text{struct}}\), pinning standard-LDSC structure
bias to conflating \(\ell^{\text{local}}\) (\(h^2\)) with \(\ell^{\text{struct}}\)
(confounding). The project is now organized as a
**framework/synthesis + diagnostics** effort, not a new estimator: theory = the
decomposition + unification + identifiability limit; applications = the structure
diagnostics (per-PC inflation spectrum, structure-robust S-LDSC enrichment, the
miscalibration map). Honest prior-art positioning is in [PRIOR_ART.md](PRIOR_ART.md)
— the qualitative concept is established (cov-LDSC, Bhatia 2016, Price 2008,
Bercovich 2025); the novel slices are the explicit identity, the unification, and
the diagnostics.

### 1.1 Use-case: the intended deliverable and its empirical verdict

The intended deliverable was a calibration/confounding **diagnostic**. PLDSC's
exact estimator consumes individual-level genotypes **and** a fitted null
covariance \(V_c=\sigma_g^2K_{-c}+\sigma_e^2I\); a null REML fit already returns
\(\sigma_g^2/(\sigma_g^2+\sigma_e^2)\) as SNP heritability, so **PLDSC does not
claim a better \(h^2\) point estimate than REML on individual data** and the
non-redundant object was to be the intercept — its displacement off 1 under a
fitted/misspecified \(V_c\) as a mixed-model analog of the LDSC intercept. The
slope is reported as \(\tau^2\) on the per-\(P_c\) metric (an \(h_g^2/M\)-type
quantity), never as total heritability.

**Empirical verdict (2026-07-18): this use-case reduces to cov-LDSC.** A direct
test (`experiments/map_ldsc_miscalibration.py`; Section 5 of METHOD, Section 13)
regressed the exact LMM moment on conventional LD scores across regimes:

- On mixed-model summary statistics, ordinary LDSC underestimates \(h^2\) by
  \(\sim30\%\) at \(F_{st}=0.15\) and inflates the intercept — but **cov-LDSC
  (Luo et al. 2021) already corrects this**, reproducing its known result.
- **Cryptic relatedness — PLDSC's one edge over cov-LDSC — does not propagate
  into the local, windowed LD-score regression** (the LMM absorbs it upstream).

So in the LDSC-application setting PLDSC gives **no differentiated advantage over
cov-LDSC**. The theory (LMM-matched LD score) is correct and validated, but the
practical gap the cov-LDSC authors flagged as future work is, in practice,
closed by cov-LDSC.

**Where the differentiated value was thought to be — and the kill-test verdict.**
Upstream, at the association step — specifically **within-family / sib-GWAS
designs**, where relatedness is definitional and exactly known, so \(P_c\) is
exact (no fitted-\(V_c\) fragility). This was the one surviving candidate niche
(from an adversarially-vetted 31-application panel). It has now been **kill-tested
and comes back NEGATIVE** (2026-07-18; `experiments/within_family_killtest.py`,
design hardened by a 3-lens adversarial review before coding). The reason is
exact: for full sibs \(R=0.5(I+B)\) and the family-mean-centering projector
annihilates the family span (\(M_BB=0\)), so \(M_BRM_B=0.5\,M_B\) — the
within-family genetic covariance that survives mean-centering is a homoscedastic
*scale*, not LD-correlated confounding. A standard sib-GWAS's cluster-robust /
effective-N standard error (Howe et al. 2022) removes it, so the standard pipeline
is already calibrated; \(P_c\) reproduces it but does not beat it, and pedigree-only
\(P_c\) is *more* fragile under admixture. The structured remainder is off-target
and orthogonal to the windowed LD score (the same Wall-1 that killed the
summary-stat angle). The only unruled-out regime — heterogeneous pedigrees on a
genome-wide intercept — is already owned by Sib-Regression (Visscher 2006) and RDR
(Young et al. 2018).

**Decision (2026-07-19): (c) STOP / archive.** After the within-family kill-test
(option (a)) failed, a full adversarial prior-art sweep (PRIOR_ART.md) showed the
theory-note niche is thin-but-defensible-at-best (the GLS operator is Mangin 2012,
the aliasing corollary is Huang 2025), and a structured contribution-ideation pass
(15 candidates, adversarially vetted) surfaced only one high-death-prior long-shot.
The project is archived as an honest record (see README "Status: CONCLUDED"). The
exposition note (former option (b)) remains available if revived
([NOTE_OUTLINE.md](NOTE_OUTLINE.md)) but was not pursued. Do **not** re-open the
LDSC-application angles (\(h^2\), \(r_g\), the intercept), the within-family angle,
or the killed `q_j` / effective-N / portability ideas as PLDSC differentiators.

### 1.2 Minimum viable deliverable (status: the original MVP was tested)

The original MVP was a simulation-grounded methods letter — *"an LMM-matched LD
score that calibrates LDSC where ordinary and covariate-adjusted LD scores do
not."* That MVP was **built (M-spike–M5) and its differentiating claim tested**,
and the claim is **false as stated**: cov-LDSC *does* calibrate for the regime
(structure) that actually biases LDSC on LMM sumstats (Section 1.1). The letter
would reduce to "use cov-LDSC," which exists — so this MVP is retired.

The candidate MVP was the **within-family matched-covariance** result
(Section 1.1). Its kill-test has now been **run and failed**
(`experiments/within_family_killtest.py`): the intercept-1 identity does transfer,
but the exact within-family covariance does **not** beat the effective-\(N\) /
cluster-robust standard error (Howe et al. 2022) — the surviving within-family
genetic covariance is a homoscedastic scale (\(M_BRM_B=0.5\,M_B\)) that the
standard SE already removes, and its structured remainder is orthogonal to the
windowed LD score. So the fallback applies: a **theory/methods note** (the
decomposition, unification, and identifiability limit), or stop.

Software novelty is nil and is not claimed: modernized, tested LDSC code already
exists. What is established is the **theorem and its validation**; what is *not*
established is a differentiated summary-statistics application.

### 1.3 Scope

The exact-method scope is deliberately narrow:

- quantitative traits and autosomes;
- complete individual-level genotypes from the GWAS sample for the exact
  matched-score derivation;
- fixed covariates identical to those used by the association model;
- an exact Gaussian LOCO score test before approximate association engines;
- oracle and deliberately-misspecified \(V_c\) in simulation; **fitted** variance
  components (GEMMA) are Phase 2, not MVP;
- genic or random-effect variance unless a derivation supports a stronger
  heritability interpretation.

Binary traits, meta-analysis, annotations, external reference panels, and
in-sample subsampling are out of exact scope until the quantitative-trait
method passes its oracle gates. Subsampling changes the GRM, precision
operator, and score covariance; it is a separately tested approximation, not
an implementation detail.

### 1.4 Decision rule

The project advances only when a milestone's scientific and software gates
pass. The **M-spike is the first kill-or-continue decision** and runs before any
packaging or write-up. A negative result is acceptable: if mixed-model-matched
scores provide no reproducible calibration advantage over covariate-adjusted
LDSC, PLDSC becomes a validation and interoperability project rather than a new
diagnostic.

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
\(\varepsilon\sim(0,V_c)\), and effects satisfy
\(E[\beta_k]=0\), \(E[\beta_k\beta_l]=\tau_k^2\delta_{kl}\), and
\(\beta_c\perp\varepsilon\mid X\), then the oracle conditional moment is

\[
E[z_j^2\mid X,V_c]
=1+\sum_k\tau_k^2\frac{(x_j^TP_cx_k)^2}{q_j}
=1+\sum_k\tau_k^2q_k\{r_{jk}^{(P_c)}\}^2. \tag{3}
\]

For a homoscedastic effect model \(\tau_k^2\equiv\tau^2\), the candidate
mixed-model LD score is

\[
L_j^{(P_c)}=\sum_{k\in W_j}q_k\{r_{jk}^{(P_c)}\}^2,
\qquad E[z_j^2]\approx1+\tau^2L_j^{(P_c)}. \tag{4}
\]

Equation (3) is **proved exact** under the stated assumptions and has been
independently re-derived three times; \(P_cC=0\) annihilates the fixed-effect
term and \(P_cV_cP_c=P_c\) pins the intercept to exactly 1. Equation (4), not
merely \(\sum_k r_{jk}^2\), is the object to test empirically. The sum in
Equation (3) covers all modeled target-chromosome effects; replacing it by the
local window \(W_j\) in Equation (4) is a truncation approximation that must
pass a window-convergence diagnostic. The \(q_k\) terms disappear only under an
explicitly justified normalization (they cancel iff \(q_k\) is constant, e.g.
ordinary LDSC with standardized genotypes and \(P=I\)); under a non-scalar
\(V_c\) metric they genuinely vary and carry the novelty over cov-LDSC.
Variants with zero or numerically negligible \(q_j\) have no defined score and
must be rejected or excluded under a documented, scale-aware threshold.

### 2.1 Exactness assumptions (record verbatim in METHOD.md)

Equation (3) is exact **only** if all of the following hold; each is a
precondition, not bookkeeping:

1. **Oracle \(V_c\):** the covariance generating \(\varepsilon\) equals the
   \(V_c\) inside \(P_c\). This is what forces the intercept to 1; a fitted or
   misspecified \(V_c\) moves it, and that displacement is the diagnostic signal
   of Section 1.1.
2. **Second-moment effect model only:** just \(E[\beta_k]=0\) and
   \(E[\beta_k\beta_l]=\tau_k^2\delta_{kl}\) are used — full independence is not
   required, and **Gaussianity is not required**, so the M3 moment gate may draw
   non-Gaussian effects.
3. \(\beta_c\perp\varepsilon\).
4. **The same \(C\)** appears in the GWAS and in \(P_c\).
5. The Equation (3) sum ranges over the **causal** target-chromosome set; an
   omitted causal chromosome-\(c\) SNP biases the moment.
6. \(q_j>0\) (\(P_c\) is PSD; \(q_j=0\) iff \(x_j\in\operatorname{span}(C)\)).

Locality lives in the \(P_c\) metric, not in centimorgans: define the window via
a threshold on \(|r_{jk}^{(P_c)}|\), and treat window non-convergence as itself a
diagnostic of \(V_c\) misspecification (it is most fragile in exactly the
structured regime the method targets).

Before implementation is promoted, the derivation must also state:

1. the effect-size scale and exact slope estimand (\(\tau^2\) on the per-\(P_c\)
   metric, an \(h_g^2/M\)-type quantity, not total \(h^2\));
2. whether modeled-effect, tested, and regression-variant sets coincide and,
   if not, how their filters enter the summation;
3. whether \(V_c\) is known or uses fitted variance components (MVP: known /
   deliberately-misspecified; Phase 2: fitted);
4. whether the target statistic is a score, Wald, likelihood-ratio, or
   engine-specific calibrated statistic;
5. whether a separate finite-sample correction is appropriate for a
   population/reference-LD estimand; no adjusted \(r^2\) is inserted into the
   exact conditional identity in Equation (3), which conditions on \(X\);
6. the regression variance, weighting, and block-resampling model rather than
   assuming conventional LDSC choices remain valid;
7. what the intercept represents under environmental covariance
   misspecification and directional LD;
8. **the deployment regime.** The exact \(P_c\) requires individual genotypes
   and \(V_c\). A summary-statistics-only path (a reference-panel approximation
   of \(P_c\)) is the eventual product; its approximation error is a named
   Phase-2 research question, and until then the exact \(P_c\) is the oracle to
   validate against, never the shipped artifact.

### 2.2 Relation to the legacy WhiLD score

The 2017 WhiLD score computed ordinary LDSC on genotypes whitened by the
Cholesky factor of the inverse of a \(w\)-weighted hybrid GRM/covariance matrix.
On paper this is the Equation (4) object with metric \(M^{-1}\) instead of
\(V_c^{-1}\), \(q_k\equiv1\), and no covariate projection — so \(P_c\)
essentially supplies the derivation WhiLD never had. M0 nests WhiLD as this
named special case, yielding a WhiLD \(\rightarrow\) correct-metric \(\rightarrow\)
add-projection \(\rightarrow\) add-\(q_k\) ablation ladder that isolates which
correction buys calibration. The unstable pseudoinverse-plus-Cholesky
*estimator* is **not** implemented (see Section 4).

### 2.3 Motivation

This work is motivated by a specific gap. Covariate-adjusted LDSC computes LD
after applying the same fixed-covariate adjustment used in the GWAS, but its
authors state that mixed-model summary statistics require random-effect-
adjusted LD and leave that computation to future work
([Luo et al., 2021](https://doi.org/10.1093/hmg/ddab130)). Residual-genotype LD
has independent theoretical support in structured populations
([Bercovich et al., 2025](https://doi.org/10.1093/genetics/iyaf009)). LOCO is
established in association testing to avoid proximal contamination
([Yang et al., 2014](https://doi.org/10.1038/ng.2876)), but that precedent does
not prove Equation (4). A 2026 landscape scan found no published LD score matched
to an LMM score statistic, so the estimator is unclaimed; confirm this scan
before M0 (Section 7).

## 3. Milestone roadmap

The roadmap is split into an **MVP** track (the first result) and a **Phase 2**
track (deferred; not a blocker for the MVP). Nothing in Phase 2 is scheduled
until the M-spike and the M3 moment gate pass.

**Table 1a. MVP milestones (first result, `v0.x`)**

| Milestone | Deliverable | Gate before proceeding |
|---|---|---|
| M-spike. Kill-or-continue | Throwaway dense \(V_c,P_c,q_j,z_j\) on one tiny structured simulation. Readout A: size of ordinary- and cov-LDSC miscalibration on LOCO-LMM z-scores. Readout B: empirical \(E[z_j^2]\) vs Equation (3), and windowed \(L_j\) vs full-chromosome remainder | **Explicit go/no-go.** If Readout A shows negligible miscalibration, the estimator is unnecessary and the project pivots (Section 1.4). If Readout B falsifies Equation (3) empirically at a scale not explained by the assumptions in 2.1, stop. 200-replicate pilot suffices to see gross miscalibration |
| M0. Freeze the theorem (light) | `docs/METHOD.md` deriving Equations (1)-(4), the assumption list (2.1), slope estimand, genotype/effect/GRM scaling, variant sets, WhiLD-as-special-case (2.2), and limiting cases | **Timeboxed to days with a named owner** (the 2017 death was an unwritten, unowned derivation, not a hard proof). With declared normalization, the moment reduces to fixed-covariate LDSC when \(V_c=\sigma^2I\) and to ordinary LDSC up to a known scale with only an intercept. No months-long "independent algebra review" blocking gate — the algebra is already verified |
| M1. Harden the fixed-covariate baseline (`v0.1`) | Complete the array API, deterministic simulator, an internal full-projector oracle fixture, a cov-LDSC comparison, and minimal CI; numerical kernels reject missing values and defer imputation to I/O | Small well-conditioned fixtures agree at `rtol=0, atol=1e-12` (tests must set `rtol=0` explicitly); invariance to full-rank covariate reparameterization; a GitHub Actions run executes the suite on the pinned interpreter; missing-value and zero-variance guards are covered |
| M2. Build a dense LMM oracle | Private small-matrix implementation of \(P_c\), \(q_j\), \(z_j\), and Equations (3)-(4), accepting supplied \(K_{-c}\) and variance components (oracle and deliberately-misspecified) | Scale-normalized residuals for symmetry, \(P_cC=0\), and \(P_cV_cP_c=P_c\) meet condition-aware tolerances; score covariance agrees with direct calculation; the oracle path consumes **raw** \(x_j^TP_cx_k\) and does not route through the finite-sample-corrected `local_ld_scores` kernel; indefinite GRMs, singular \(V_c\), ill-conditioned designs, and negligible \(q_j\) produce explicit diagnostics |
| M3. Falsify the oracle moment and demonstrate calibration (`v0.2`) | Repeated draws of effects and noise conditional on fixed \(X,V_c\); exact full-chromosome moment and window-truncated approximation; the intercept diagnostic under oracle vs deliberately-misspecified \(V_c\) | **The pivotal gate.** Per-variant/bin Monte Carlo residuals support Equation (3) with covariance-aware uncertainty and no trend by MAF, \(q_j\), ancestry loading, or GRM leverage; the omitted-window remainder is reported; under a known-misspecified \(V_c\) the intercept moves as predicted while ordinary/cov-LDSC do not calibrate |

**Table 1b. Phase 2 milestones (deferred; a separate software paper, `v0.3`+)**

| Milestone | Deliverable | Gate before proceeding |
|---|---|---|
| M4. Quantify fitted-component error | Import fitted null variance components, construct \(\widehat P_c\), compare with the known-\(V_c\) oracle; PLDSC does not implement REML. Fitted-component *science* validated via a local REML (`experiments/m4_fitted_vc.py`, `m5_inference.py`); external GEMMA byte-parity specified in [GEMMA_PARITY.md](GEMMA_PARITY.md) and deferred to a Linux environment (no Windows GEMMA build) | Conditional-moment error meets its preregistered gate or a correction/scope limit is documented; GEMMA fitting mode, allele scaling, and score output field are reproduced exactly |
| M5. Add regression and inference (`v0.3`) | M0-derived weights and resampling scheme, explicit estimand labels, provenance/compatibility schema, machine-readable results; compare ordinary LDSC, cov-LDSC, and optional legacy scores | Synthetic slope/intercept bias and interval coverage pass; resampling matches a direct oracle and achieves calibrated coverage; incompatible adjustment metadata fails before fitting |
| M6a. Block local products | Blocked evaluation of full and windowed candidate scores | Dense and blocked results agree at `rtol=1e-10`, `atol=1e-12`; memory grows with active blocks/windows rather than chromosome-wide pair count |
| M6b. Add scalable solves | Matrix-free \(P_cx\), iterative-solver diagnostics, preconditioning | Operator results agree with the dense oracle under condition-aware tolerances; convergence failures are explicit and reproducible |
| M6c. Interoperate (`v0.4`) | PLINK 2 adapter, minimal CLI, provenance sidecars, only scientifically justified LDSC-format output | One end-to-end PGEN fixture passes; optional dependencies stay isolated; output round-trips with all required count metadata |
| M7. Validate and release (`v1.0`) | Synthetic phenotypes on real genotypes, independent reproduction, stable schemas/API, archival release, limited real-trait analyses | Every preregistered gate passes; conclusions survive ancestry, window, GRM, and approximation sensitivity analyses; all headline results reproduce from committed configurations |

Current status (2026-07-18): the MVP arc is built and validated end-to-end.
The **M-spike** passed (`experiments/m_spike_kill_or_continue.py`): Eq (3)
confirmed numerically and PLDSC calibrates where ordinary LDSC does not. **M0**
is written (`docs/METHOD.md`). **M1** is hardened: CI workflow, `rtol=0`
tolerance gate, dependency floor raised, invariance/guard tests added. **M2**
landed as the tested `pldsc/_lmm.py` oracle (with `pldsc/simulate.py`); **M3**
(`experiments/m3_relatedness_admixture.py`) shows PLDSC beats cov-LDSC
specifically under relatedness; **M4** fitted-component science is validated by a
local REML (intercept `1.010 ± 0.078`), with external GEMMA byte-parity
specified in [GEMMA_PARITY.md](GEMMA_PARITY.md) and deferred to Linux; **M5**
(`pldsc/regression.py` + `experiments/m5_inference.py`) provides block-jackknife
inference and quantifies the fitted-`V_c` intercept-SE requirement. 29 tests
pass. **Then the use-case was investigated directly** and angle #1 ("LDSC
miscalibrated for LMM sumstats") came back negative for PLDSC
(`experiments/map_ldsc_miscalibration.py`): the bias is a structure effect
cov-LDSC already fixes, and relatedness does not reach the local LD-score
regression (Section 1.1). Remaining, contingent on the Section 1.1 decision: the
within-family kill-test; else a theory note. Deferred regardless: the GEMMA
cross-check on Linux, at-scale `V_c`-uncertainty propagation, and Phase-2
scale/interop (M6). **After two external reviews (2026-07-18):** the decomposition
was restated *exactly* as the GRM spectral identity (not an empirical discovery)
and the honest windowed ~5% figure led; one **non-circular** robustness result was
added (`experiments/misspecified_background.py`: the null intercept is robust to a
sparse/directional off-chromosome background); and **Bercovich 2025 was read in
full** — it is a pairwise-\(r^2\) QC method (pruning/clumping/FST/PCA) with a hard
rank-K projection and no heritability content, so it does not collide with the
heritability-facing claims (see [PRIOR_ART.md](PRIOR_ART.md)).

### 3.1 Forward roadmap (decomposition-framed paper)

Tables 1a/1b document the **completed estimator-MVP foundation** — the proved
moment (Eq 3), the dense \(P_c\) oracle, the calibration/relatedness experiments,
and the inference layer. That foundation stands and is reused, but the framing has
moved (Section 1) from "a new LMM-matched estimator" to a **framework/synthesis +
diagnostics** paper built on the local/non-local decomposition. Forward
deliverables and status:

| Deliverable | Track | Status |
|---|---|---|
| D1. Theory: exact spectral decomposition \(\ell=\ell^{local}+\ell^{struct}\); PC score as its leading-eigenvector truncation; unification of cov-LDSC/WhiLD/PC-score/\(P_c\); identifiability-limit theorem | theory | METHOD §14 written; identity stated exactly; limit stated, not yet a formal proof |
| D2. Diagnostics: per-PC inflation spectrum; structure-robust S-LDSC enrichment (PC-loading annotation); the miscalibration map | applications | prototyped (`detect_pc_inflation`, `sldsc_pc_annotation`, `map_ldsc_miscalibration`) with honest caveats; not yet real-data |
| D3. Prior-art positioning | positioning | Bercovich 2025 gate **closed** (full read); broader landscape sweep still to re-run |
| D4. One non-circular robustness result | validation | **done** (`misspecified_background.py`): null intercept robust to sparse/directional off-chr background |
| D5. Real-data survey / efficient decomposition engine | efficiency + empirics | scoped ([REAL_DATA_SURVEY.md](REAL_DATA_SURVEY.md)); gated on an efficient \(P_c\) engine; not started |
| D6. Alternative differentiated application: within-family matched covariance | application | **kill-tested, NEGATIVE** (`within_family_killtest.py`); option (a) closed |

The honest gate for the paper: D6 (the last candidate differentiated application)
has now failed its kill-test, so the choice narrows to whether **D1 + D3 + D2/D4**
is a clarifying contribution on its own (a theory/methods note) or whether the
project stops. That is the central open decision (Section 1.1's (b) vs (c)).

GEMMA is the first external association benchmark (Phase 2) because it implements
an exact LMM suitable for checking a dense implementation
([Zhou and Stephens, 2012](https://doi.org/10.1038/ng.2310)). BOLT-LMM,
fastGWA, and REGENIE are later compatibility targets because their whole-genome
prediction and calibration are not identical to the oracle Gaussian score test
([Loh et al., 2015](https://doi.org/10.1038/ng.3190);
[Jiang et al., 2019](https://doi.org/10.1038/s41588-019-0530-8);
[Mbatchou et al., 2021](https://doi.org/10.1038/s41588-021-00870-7)).

## 4. Simulation protocol

**Table 2. Required simulation dimensions**

This is a targeted scenario matrix, not a Cartesian product. The M-spike and M3
start with the smallest scenarios that distinguish the candidate moment from its
limits; M5 (Phase 2) adds estimator-level stress cases.

| Dimension | Required settings | Purpose |
|---|---|---|
| Genotypes | Homogeneous haplotypes; two discrete populations; recent admixture; pedigrees or cryptic relatives | Separate local LD, structure-induced LD, and relatedness |
| Target-chromosome architecture | Infinitesimal Gaussian; sparse; MAF/LD-dependent; ancestry-correlated and directional effects | Test the random-effect assumptions and estimand interpretation |
| Other-chromosome architecture | Infinitesimal Gaussian background vs sparse/directional \(g_{-c}\) in structure-LD with target-chromosome SNPs | Test assumption 2.1(1): non-infinitesimal off-target effects violate \(\varepsilon\sim(0,V_c)\) and leak into \(z_j\) |
| Environment | Independent; ancestry-aligned; family-shared; covariance misspecified relative to \(V_c\) | Test intercept behavior and the calibration diagnostic |
| Adjustment | Correct and mismatched covariates; correct and misspecified GRMs; oracle and (Phase 2) fitted variance components | Measure geometry mismatch and plug-in error |
| Reference approximation | Full in-sample oracle first; only later (Phase 2), random in-sample subset, ancestry-matched external, and deliberately mismatched external | Quantify rather than assume when score reuse is defensible |
| Numerical stress | GRM condition number; sample size; marker density; LD-window truncation; low MAF | Expose instability and finite-sample failure |
| Outputs | Moment calibration by score bin; intercept displacement vs \(V_c\) misspecification; bias/RMSE; interval coverage; type-I calibration; runtime and peak memory | Prevent selection of one flattering metric |

Each oracle-moment simulation fixes \(X,V_c\), repeatedly samples independent
effects and noise, and computes a per-variant Monte Carlo mean. Correlated SNPs
are not treated as independent replicates. Calibration uses per-variant or
predeclared-bin standardized residuals with covariance-aware uncertainty.
M0b preregisters each trend diagnostic. Internal exploratory trend checks (MAF,
\(q_j\), ancestry loading, GRM leverage) are eyeballed, not multiplicity-
corrected; the moment-calibration and coverage gates are the two that carry
promotion decisions.

Each simulation record must contain its full configuration, seed, **RNG bit
generator (pin PCG64 alongside the seed — "seed" alone is not reproducible
across NumPy majors)**, Git revision, truth, estimator output, and Monte Carlo
standard error. Ordinary tests remain small and deterministic; larger stochastic
gates run separately from unit tests.

Replicate budget is right-sized to the decision: the M-spike and exploratory
work use a **200-replicate pilot**; **1,000+ replicates and the 1-percentage-
point Monte Carlo standard-error caps are reserved for the final confirmatory
figures**, subject to a predeclared compute cap.

The legacy WhiLD/PC-score constructions are **paper-only comparators** (Section
2.2): derived as special cases and used in the ablation ladder, but the unstable
pseudoinverse-plus-Cholesky estimator and the eigenvalue-weighted PC score are
not implemented as release code. Only ordinary LD and cov-LDSC are implemented
baselines.

### Promotion and falsification rules

1. **M-spike gate.** Readout A must show that ordinary and/or cov-LDSC
   materially miscalibrate LOCO-LMM z-scores; if not, the estimator is
   unnecessary. Readout B must show empirical mean \(z_j^2\) tracking Equation
   (3) within covariance-aware Monte Carlo uncertainty on one tiny scenario.
2. **M3 moment gate.** Empirical mean \(z_j^2\) must agree with Equation (3)
   within covariance-aware Monte Carlo uncertainty. A regression summary may
   report intercept 0 and slope 1, but it cannot replace standardized residual
   checks when variants are correlated or the predictor range is narrow.
3. **M3 truncation diagnostic.** Full-chromosome and windowed scores report the
   omitted remainder. A less than 1% change on doubling the window marks
   convergence; genuine long-range LD may instead require a whole-chromosome
   calculation or an explicit nonconverged result.
4. **M3 calibration diagnostic.** Under a known-misspecified \(V_c\), the
   PLDSC intercept must move in the predicted direction and magnitude while
   ordinary/cov-LDSC do not, demonstrating the Section 1.1 use-case in
   simulation.
5. **M4 plug-in gate (Phase 2).** Persistent conditional-moment error from
   fitted components narrows the scope. It is not rebranded as an intercept
   problem.
6. **M5 bias gate (Phase 2).** For the named M0 slope estimand in predeclared
   nonzero-truth Gaussian scenarios, absolute relative bias must be below 5%.
   Zero-truth cases use absolute bias and type-I calibration instead.
7. **M5 coverage gate (Phase 2).** Nominal 95% intervals must have 92%-98%
   empirical coverage, and their binomial Monte Carlo interval must contain 95%.
8. **M5 limit gate (Phase 2).** In predeclared nonzero-truth homogeneous
   scenarios, PLDSC must not worsen absolute relative bias by more than 0.02
   compared with ordinary LDSC.
9. **M5 comparison gate (Phase 2).** Baseline selection and the aggregate RMSE
   metric are fixed before held-out seeds are run. In predeclared structured-
   LOCO scenarios, the fitted method must improve paired aggregate RMSE; one
   favorable scenario is not evidence.
10. **M6 equivalence gate (Phase 2).** Scalable implementations must agree with
    the dense oracle before any performance claim is reported.

The heritability label is especially constrained in admixed populations:
random-effect estimators can omit directional-LD covariance and therefore
differ from total realized genetic variance, even without environmental
confounding
([Huang et al., 2025](https://doi.org/10.1093/genetics/iyaf100)). PLDSC will
therefore report the derived estimand by name instead of calling every slope
"heritability."

## 5. Software shape

Keep the numerical layer array-first and small. Public APIs are added only
after their corresponding scientific gate passes. Only the MVP components are on
the critical path; Phase 2 components are built after M3.

**Table 3. Planned modules and constraints**

| Component | Earliest milestone | Track | Constraint |
|---|---|---|---|
| `pldsc/ld.py` | Existing / M1 | MVP | Sample-by-variant NumPy arrays; independent dense oracle tests |
| `pldsc/simulate.py` | M1 / M-spike | MVP | Explicit random generator (pin bit generator + seed); configuration and truth returned with data |
| `pldsc/_lmm.py` | M2 | MVP | Private dense oracle first; supplied GRM and variance components; raw cross-products; no public estimator claim |
| `pldsc/regression.py` | M5 | Phase 2 | Derived weighting/resampling separated from score construction |
| `pldsc/schema.py` | M5 | Phase 2 | Versioned adjustment and provenance metadata |
| `pldsc/_blocks.py` | M6a | Phase 2 | Blocked products must reproduce the dense reference path |
| `pldsc/_operator.py` | M6b | Phase 2 | Matrix-free solves with convergence and conditioning diagnostics |
| `pldsc/io.py` and CLI | M6c | Phase 2 | Optional adapters around kernels; no bespoke genotype format |
| `tests/` | Every milestone | Both | Algebraic properties, deterministic fixtures, then end-to-end tests |
| `experiments/` | M-spike onward | Both | Committed configurations and summaries, not generated bulk data |

Engineering rules:

- use `float64` for the reference path and make tolerances explicit;
- never construct an \(M\times M\) LD matrix;
- permit dense \(N\times N\) matrices only in the reference implementation;
- **the M2/M3 oracle path consumes raw \(x_j^TP_cx_k\), not `local_ld_scores`,
  which defaults to `finite_sample_correction=True`**;
- postpone pandas/xarray models, GPU/JAX/Numba paths, distributed execution,
  dashboards, multiple genotype readers, and the operator/PGEN/scale work until
  Phase 2;
- benchmark only after numerical agreement, recording machine, Python, NumPy,
  and BLAS details.

The biobank-scale benchmark, the matrix-free operator solve gate, and the PGEN
adapter are Phase-2 engineering promotion gates and are **not** MVP work; their
detailed thresholds are recorded in the Phase-2 addendum when M4 begins, not
here.

### Dependencies and CI

The dependency floor is currently `requires-python = ">=3.10"` and
`numpy>=1.24` (`pyproject.toml`); both are near end-of-life and the development
machine has no 3.10 interpreter. Raise the floor to a currently-supported pair
(Python 3.11+/NumPy 2.0+) unless a concrete downstream consumer needs the older
pins, and document the reason if the old floor is kept.

There is currently **no `.github/workflows`**, so the M1 "wheels/tests pass"
gate has no infrastructure. Add a minimal GitHub Actions workflow that runs
`python -m unittest discover -s tests` on the pinned interpreter before the M1b
checkpoint is claimed. Cross-software `atol=rtol=1e-12` agreement across
different BLAS builds is not achievable and is not gated; only internal-fixture
agreement is.

PLINK 2 PGEN/PVAR/PSAM is the sole planned genotype adapter (Phase 2)
([PLINK 2 format documentation](https://www.cog-genomics.org/plink/2.0/formats)).
For interoperability, do not label Equation (4) as conventional `L2` unless M0
proves the required normalization. If that equivalence is established, emit
the conventional `CHR`, `SNP`, `BP`, and `L2` columns together with `.l2.M` and
`.l2.M_5_50` counts; otherwise use an explicitly versioned PLDSC schema
([official LDSC LD-score tutorial](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial)).
The JSON sidecar introduced in M5 records formula revision, package and Git
revision, sample and variant counts, adjustment source, genome build and map,
window, finite-sample treatment, input identifiers/checksums, and simulation
seed and bit generator when applicable.

## 6. Commit and review discipline

Every commit checkpoint is small and green. Scientific derivations are reviewed
before their optimized implementations. Only the MVP checkpoints are scheduled;
the Phase-2 tail (M4-M7) is regenerated once M3 passes.

**Table 4. Planned commit checkpoints (MVP)**

| Milestone | Commit intent | Required check before commit |
|---|---|---|
| C1 | `docs: fix the PLDSC use-case, estimand, and M0 ownership` | Section 1.1 use-case, estimand labels, intercept-displacement = confounding, named M0 owner |
| C2 | `feat+test: dense P_c oracle and kill-or-continue spike` | Readout A miscalibration + Readout B Monte Carlo moment; raw cross-products; seed + bit generator replay |
| C3 | `docs: legacy triage, WhiLD-as-special-case, authorship note` | Every 2017 number marked non-citable-until-reproduced; WhiLD derived as Equation (4) special case; co-author/scope paragraph |
| M0 | `docs: derive the PLDSC mixed-model moment` | Rendered equations, assumption list, reductions; named owner sign-off |
| M1a | `feat: harden fixed-covariate LD scores` | Unit and cov-LDSC fixture tests; `rtol=0, atol=1e-12` |
| M1b | `build+ci: add supported Python pin and Actions` | Clean build and green Actions run on the pinned pair |
| M2a | `test: specify mixed-model oracle identities` | Scale-aware reference assertions |
| M2b | `feat: add the dense mixed-model oracle` | Identity, failure-mode, and direct-covariance tests; raw-cross-product path |
| M3 | `test: validate the oracle moment and calibration diagnostic` | Seed replay, Monte Carlo uncertainty, truncation summaries, intercept-under-misspecification |

Phase-2 checkpoints (M4-M7) follow the same discipline and are listed in the
Phase-2 addendum after the M3 gate passes. No release milestone is awarded
because the version number looks lonely.

## 7. Primary references

1. Bulik-Sullivan et al. (2015), [LD Score regression distinguishes
   confounding from polygenicity](https://doi.org/10.1038/ng.3211).
2. Finucane et al. (2015), [Partitioning heritability by functional annotation
   using genome-wide association summary statistics](https://doi.org/10.1038/ng.3404).
3. Bhatia et al. (2016), [Subtle stratification confounds estimates of
   heritability from polygenic models](https://doi.org/10.1101/048181).
4. Zhou and Stephens (2012), [Genome-wide efficient mixed-model analysis for
   association studies](https://doi.org/10.1038/ng.2310).
5. Yang et al. (2014), [Advantages and pitfalls in the application of
   mixed-model association methods](https://doi.org/10.1038/ng.2876).
6. Loh et al. (2015), [Efficient Bayesian mixed-model analysis increases
   association power in large cohorts](https://doi.org/10.1038/ng.3190).
7. Jiang et al. (2019), [A resource-efficient tool for mixed-model association
   analysis of large-scale data](https://doi.org/10.1038/s41588-019-0530-8).
8. Luo et al. (2021), [Estimating heritability and its enrichment in
   tissue-specific gene sets in admixed populations](https://doi.org/10.1093/hmg/ddab130).
9. Mbatchou et al. (2021), [Computationally efficient whole-genome regression
   for quantitative and binary traits](https://doi.org/10.1038/s41588-021-00870-7).
10. Bercovich et al. (2025), [Measuring linkage disequilibrium and improvement
    of pruning and clumping in structured populations](https://doi.org/10.1093/genetics/iyaf009).
11. Huang et al. (2025), [Interpreting SNP heritability in admixed
    populations](https://doi.org/10.1093/genetics/iyaf100).

Additional comparators surfaced by the 2026 landscape scan, to position against
once confirmed (verify identifiers before citing): S-LDSC baseline-LD model
(Gazal et al. 2017); RHE-mc randomized heritability estimation; LDER and HEELS
likelihood/eigenvalue LDSC variants; and recent modernized Python LDSC
reimplementations. The DOIs above for Finucane 2015, Bhatia 2016, and Jiang 2019
were added during revision and should be confirmed against the canonical
citations.
