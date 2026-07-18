# PLDSC — independent review: idea, theory, implementation, and plan

Date: 2026-07-18
Reviewer: independent code/theory review (requested by the project owner).
Scope: the whole repository at commit `f82df69` — `README.md`, all of `docs/`,
all of `pldsc/`, all of `tests/`, all of `experiments/` (scripts, not just their
README), plus re-running the test suite and an independent numerical check of
the central identity.

This review was asked to be honest and careful. It is therefore direct: it
credits what is genuinely good (a lot, especially process and documentation
honesty) and states plainly what is weak (the circularity of the empirical
evidence and the thinness of the remaining novelty). Nothing here is based on
the project's own summaries alone: every load-bearing claim was checked against
the code or re-derived.

---

## 1. Executive summary

**What PLDSC is.** A research project (not a production tool) studying LD Score
regression when GWAS samples contain population structure or when summary
statistics come from covariate-adjusted or mixed-model association tests. It is
a disciplined restart of a 2017 prototype (`newldsc`), built around one exact
conditional-moment identity for the mixed-model score statistic and — after an
honest internal kill-test — reframed around a local/non-local decomposition of
the LD score.

**Verdict in one paragraph.** The theory is **correct but modest**: the central
theorem is the standard second moment of the GLS score statistic, proved
cleanly and unit-tested, and the project itself then demonstrated that its
practical payoff over published cov-LDSC is nil in the LDSC-application
setting. The software is **genuinely good**: small, readable, well-tested
(29/29 tests pass on this machine), honest docstrings. The empirical evidence
is **the weak point**: every "validation" experiment simulates data under
exactly the model the theorem assumes, so the experiments confirm algebra,
not robustness; and the headline decomposition identity is (a) additive by
construction in the experiment and (b) — as verified independently below — a
restatement of the spectral decomposition of the GRM, i.e. exact textbook
linear algebra rather than a new empirical finding. The plan and documentation
are **unusually honest** — negative results are recorded, scope is narrowed
when the data say so, prior art is self-assessed almost harshly — and that is
the project's greatest strength. The remaining novel slice (explicit identity +
unification + diagnostics) is real but thin, and the project knows it.

**Bottom line.** Sound, honest, well-executed research plumbing in search of a
contribution. As a *theory/methods note* ("one identity unifies the structure
corrections; here are two diagnostics and their identifiability limit") it is
publishable-quality thinking, provided the decomposition claim is restated as
the exact spectral identity it actually is and at least one non-circular
robustness experiment is added. As a *new estimator* it is finished — and the
project's own experiments are what finished it, which is to its credit.

---

## 2. The idea

### 2.1 What the idea is now

The current organizing idea (README, METHOD §14, DEVELOPMENT_PLAN §1) is:

> The LD score decomposes as ℓ = ℓ^local + ℓ^struct — recombination-driven
> linkage plus a population-structure component — and ℓ^struct "is" the
> eigenvalue-weighted sum of squared SNP–PC loadings (the "PC score"). One
> identity thereby unifies cov-LDSC, the legacy WhiLD score, the PC score, and
> the mixed-model residual-precision operator P_c. Standard LDSC's structure
> bias is the conflation of the two components. Applications are diagnostics:
> a per-PC inflation spectrum, a structure-robust S-LDSC annotation, and a
> miscalibration map — all bounded by a stated identifiability limit
> (ancestry-aligned signal and confounding load on the same PCs).

### 2.2 Assessment of the idea

- **The problem is real and well-chosen.** LDSC's assumptions are visibly
  violated by structure and by mixed-model association statistics, and the
  cov-LDSC authors did leave random-effect-adjusted LD as future work. A
  project that pins down exactly what that gap contains is worthwhile.
- **The idea's trajectory is a strength, not a weakness.** The project began
  with "an LMM-matched LD score as a new estimator," built a kill-or-continue
  spike first, let the spike and the follow-up map experiment decide, recorded
  the negative verdict ("reduces to cov-LDSC in practice," METHOD §13), and
  reframed around what survived. This is how methods research is supposed to
  work and rarely does. The reframed idea (framework + diagnostics +
  identifiability limit) is more honest and more durable than the original.
- **The idea's ceiling is limited by its own correct finding.** The genuinely
  discriminating result of the whole arc is: *relatedness does not propagate
  into the local, windowed LD-score regression; structure does, and cov-LDSC
  already handles structure.* That is a useful, clarifying, publishable-in-
  principle negative/clarifying result — but it caps the idea at
  "synthesis + diagnostics." The docs say this themselves (PLAN §1.1–1.2);
  this review independently agrees (§5 below).

---

## 3. Theory

### 3.1 The central theorem (METHOD §2–4, Eq. M6)

Claim: for the matched score z_j = x_jᵀP_c y / √(x_jᵀP_c x_j) with
P_c = V⁻¹ − V⁻¹C(CᵀV⁻¹C)⁺CᵀV⁻¹, under ε ~ (0, V_c) and second-moment
effect assumptions,

E[z_j² | X, V_c] = 1 + Σ_k τ_k² (x_jᵀP_c x_k)² / q_j.

**The proof is correct.** I re-derived it: P1 (P_cC = 0, including the
rank-deficient pseudoinverse case via B(BᵀB)⁺BᵀB = B), P2 (P_cV_cP_c = P_c,
which is what pins the intercept to exactly 1), symmetry and PSD-ness all
check out, and the moment calculation is standard. The unit tests verify the
algebraic identities numerically (tests/test_lmm.py), and I confirmed all 29
tests pass (§4.1). The careful statements — intercept exactly 1 only under an
oracle V_c, the self-term j = k present, q_k weights not cancelling in
general, no Gaussianity needed for the mean — are all right.

**But the theorem is elementary, and the docs slightly over-dress it.** Eq. M6
is the second moment of the GLS/score statistic under a random-effects model —
the kind of identity that underlies every LMM association paper's calibration
argument. "Derived independently three times" (METHOD.md header) is true but
the derivation is five lines of conditional-moment algebra; the risk the
project faced was never the theorem, as the plan itself admits (PLAN §1). The
theorem's actual scientific content is not the identity but the *decision* to
use Σ_k q_k {r_jk^(P_c)}² as an LD-score regressor — and that is precisely
the piece the project's own map experiment later showed adds nothing over
cov-LDSC in practice. To METHOD §13's credit, this is stated plainly.

### 3.2 The reductions and the WhiLD nesting (METHOD §7–8)

Correct. V = σ²I gives P_c = σ⁻²P_C exactly (cov-LDSC), and the intercept-only
case with standardized genotypes gives constant q_k and recovery of the
ordinary LD score up to scale. Both are unit-tested
(tests/test_lmm.py::ReductionTests). The ladder
"ordinary → correct metric → add projection → add q_k" is a genuinely
clarifying way to present the method family, and nesting the 2017 WhiLD score
as the q_k ≡ 1, no-projection special case finally gives that prototype the
derivation it lacked.

### 3.3 The decomposition identity (METHOD §14) — checked independently

This is the current centerpiece, so I checked it independently rather than
trusting `experiments/ld_decomposition.py`. With standardized genotypes
(N×M), r_jk = x_jᵀx_k/N and K = XXᵀ/M with eigenpairs (λ_p, u_p):

Σ_k r_jk² = (1/N²)‖Xᵀx_j‖² = (M/N²) x_jᵀK x_j = (M/N) Σ_p (λ_p/λ̄) corr(x_j, u_p)²,

with λ̄ = tr(K)/N ≈ 1. **The genome-wide LD score (sum over all k) is exactly
an eigenvalue-weighted sum of squared PC loadings over all eigenvectors.**
I verified this numerically on the project's own simulator
(`pldsc.simulate.simulate_admixture`, N=600, M=900, fst=0.10, PCG64 seed 123):

- max |ℓ_full − s_all| / max ℓ_full = **0.014** — the identity is exact up to
  the standardization slack (`standardize()` scales by √(2p(1−p)), expected
  rather than sample variance, so columns are not exactly unit-norm);
- truncating to the top 10 PCs: corr(ℓ_full, s₁₀) = **0.993** — reproducing the
  project's headline "0.99" — but s₁₀ captures only **~37% of the magnitude**
  (mean 1.42 vs 3.82);
- within a ±2 cM window the mean score is 1.40, i.e. the non-local share of the
  full score is ~63% here (their run: ~50%), while the windowed score that
  LDSC actually uses is mostly local.

Three consequences the docs should absorb:

1. **The "identity" ℓ^struct = s_j is not an empirical discovery; it is the
   spectral decomposition of the GRM** (equivalently, of the genome-wide r²
   matrix) restricted to its leading eigenvectors. Framed exactly, it is a
   nice, citable way to *state* the unification — but the current framing
   ("validated at corr = 0.99") presents a scale/shift-invariant correlation
   from a 10-PC truncation as validation of an equality that, in magnitude,
   misses ~63% of the full score at 10 PCs. README line 16 ("the non-local
   component **is the 'PC score'**") overstates it; what is exact is the
   all-eigenvector identity above.
2. **The additive split is by construction in the experiment.**
   `experiments/ld_decomposition.py:97` computes
   `ell_struct = ell_ord − ell_local`, so "4.20 = 2.12 + 2.08" cannot be
   evidence of additivity. The only real empirical content is
   corr(ℓ^ord − ℓ^local, s₁₀) = 0.988 across the 400 target-chromosome SNPs —
   cross-chromosome pairs never enter any score in that script.
3. **The "structure is ~50% of the genome-wide LD score" headline** uses the
   genome-wide (cross-chromosome-inclusive) score that LDSC never regresses
   on; within the windowed score LDSC actually uses, structure is ~5%. The
   docs do note this (METHOD §14), but the README leads with the 50% figure.
   The honest framing is the opposite ordering: structure barely touches the
   windowed slope (which is *why* cov-LDSC suffices and why the original
   estimator angle died).

### 3.4 The identifiability limit

Stated as "a stated result" (METHOD §14): genuine ancestry-aligned polygenic
signal and structure confounding both load on the same PCs, so no PC-only
method can separate them. This is correct and important — but it is a
well-known observation (it is the core caveat of Bhatia 2016, of PC-corrected
GWAS generally, and of the differential-selection literature), not a new
theorem. Calling it "a stated result" is fair; the experiments/README phrase
"as a stated theorem" (experiments/README.md:339) is not.

### 3.5 What is theoretically missing

- **No non-circular robustness result.** METHOD §6 itself lists the
  assumption that matters most — non-infinitesimal, sparse, or directional
  effects on the *other* chromosomes, which enter ε and violate
  ε ~ (0, V_c) — and calls it "a required simulation axis." It is implemented
  nowhere (§5.1). The theory's fragility is therefore asserted, never
  measured.
- **The window-locality claim** ("locality lives in the P_c metric,"
  METHOD §10) is asserted from one empirical finding, not derived; it is
  plausible and consistent with §14 but deserves the same care as the rest.

---

## 4. Implementation

### 4.1 Verification performed

- **Tests: 29/29 pass** (Python 3.14.6, NumPy 2.4.6, conda env `pldsc`,
  `python -m unittest discover -s tests -v`, 0.29 s). Note: on this Windows
  machine the env's `python.exe` hard-crashes (exit 127, no traceback) on any
  LAPACK call *unless the env is conda-activated* (DLL loading) — an
  environment quirk, not a code bug, but it will bite anyone running the
  interpreter directly.
- The package installs cleanly in CI terms (pure NumPy, `numpy>=2.0`,
  `python>=3.11`); `.github/workflows/ci.yml` runs the suite on 3.11 and 3.13.

### 4.2 Code quality — genuinely good

- `pldsc/ld.py` — covariate-adjusted local LD scores. Correct: least-squares
  residualization via `lstsq` (rank-aware), scale-aware zero-residual-variance
  guard, the adjusted-r² finite-sample correction with residual rank
  ν = N − rank(C), negative adjusted terms retained (right call — truncation
  at 0 would bias), self-term fixed to exactly 1, window via `searchsorted`.
  The per-variant Python loop is O(M·window) and will not scale to real
  chromosomes — fine for a reference path, and the plan correctly defers
  blocking to M6a. Input validation is thorough and tested.
- `pldsc/_lmm.py` — the dense P_c oracle. Clean, minimal, private by design,
  and its docstring correctly scopes it as a reference implementation.
  Symmetrization `0.5(P + Pᵀ)` is a reasonable round-off guard. Dense O(N³)
  inverse — appropriate for an oracle.
- `pldsc/regression.py` — WLS + delete-a-block jackknife. Correct and
  honestly documented. One overstatement: the `two_step` weighting
  1/(2·E[χ²]²) is called "the standard LDSC heteroscedasticity weighting";
  real LDSC additionally down-weights by the sampling variance of the LD
  score itself (Bulik-Sullivan 2015 eq. for w_j). For the simulation regime
  here (in-sample exact scores) the omission is harmless; label it as the
  simplified weighting it is.
- `pldsc/simulate.py` — deterministic, seeded, metadata-returning; Acklam's
  ndtri to avoid a scipy dependency is a reasonable choice. **But the
  genetics are crude in ways that matter for the conclusions** (§5.3):
  `_sib_genome` transmits whole parental haplotypes with **zero recombination
  within a chromosome** (pldsc/simulate.py:89-94), so sib IBD is quantized to
  whole-chromosome blocks; `simulate_admixture` gives every individual one
  genome-wide-constant ancestry fraction (no admixture tracts, hence no
  admixture LD by construction) and both ancestral populations share the same
  AR(1) LD process differing only in Balding–Nichols frequencies. Adequate
  for algebra checks; not adequate to support population-genetic claims like
  "relatedness does not propagate" with real-meiosis confidence.

### 4.3 Repo hygiene

Good: tests mirror the theory (algebraic identities, reductions, guards,
`rtol=0` discipline), CI exists, seeds and bit generators pinned in
experiments. Minor: `build/` and `pldsc.egg-info/` are committed artifacts
that should be gitignored; the `build/lib` copy of the package is stale
relative to `pldsc/` (it predates `regression.py`/`simulate.py`).

---

## 5. Experiments and evidence

I read all eight experiment scripts (not just experiments/README.md) and
checked them against the documented claims. Headline: **the scripts do
compute what the README says, with two exceptions noted below — but the
evidence is almost entirely circular or toy-scale.**

### 5.1 Circularity — the central evidentiary weakness

Every experiment that "confirms" the theory generates data under exactly the
model the theorem assumes:

- m_spike (experiments/m_spike_kill_or_continue.py:168-169, 209-211):
  ε = cholesky(V_c)·z with the *same* V_c used to build P_c, and
  β ~ N(0, τ²I) i.i.d. — assumption M3 verbatim. The "confirmation" of Eq. M6
  ((intercept, slope) ≈ (1, τ²)) is the law of total variance plus P2.
- m3 (experiments/m3_relatedness_admixture.py:177-191): the headline "PLDSC =
  1.000 in all three scenarios" is computed **analytically** as
  xᵀP_cV_cP_cx / xᵀP_cx with the data-generating V_c — it is the P2 identity
  evaluated, not an estimate. experiments/README.md:97 does admit this
  ("the P2 identity guarantees it"), but the number still sits in a results
  table dressed as a result.
- m4's REML is a genuine profile restricted-likelihood maximization
  (eigenbasis golden-section, logdet term included — correctly done), but the
  fitted family contains the truth exactly (same GRM, two unknown scalars),
  so "calibration survives fitting" is shown only for a correctly specified
  model at N=400.
- **No script anywhere misspecifies the covariance *structure*** — sparse,
  directional, or large off-chromosome background effects, the failure mode
  METHOD §6 flags as the required axis. The claim "PLDSC calibrates where
  ordinary LDSC does not" is thus established only where the assumed model
  holds by construction.

This does not make the experiments worthless — algebra checks catch
implementation bugs, and the m_spike did its kill-or-continue job. But the
docs occasionally drift from "numerically confirms the identity" toward
"validates," and the distinction matters for everything downstream.

### 5.2 The map experiment's key negative result is partly designed-in

`experiments/map_ldsc_miscalibration.py` regresses the *analytic oracle-LMM
moment* (relatedness already removed by P_c, line 141) on **local windowed**
conventional LD scores (non-local LD discarded by the ±2 cM window,
lines 49, 130). Relatedness can therefore enter only through two doors, both
locked. "Relatedness produces no measurable drift" is close to a tautology of
the design; an OLS-GWAS dependent variable (where relatedness *would* show
up — m_spike found intercepts of 9–12) is never included. The *mechanism* is
correctly disclosed in METHOD §10, and the conclusion is very likely right —
but it should be presented as a structural argument with an illustrative
computation, not as an empirical discovery. (Notably this is the load-bearing
result for the whole "reduces to cov-LDSC" verdict, so its status matters.)

### 5.3 Simulator fidelity

As detailed in §4.2: no recombination in meiosis, no admixture tracts,
identical AR(1) LD in both ancestral populations, no mutation/drift beyond
the Balding–Nichols draw, N = 400–600, 80–900 SNPs on 2–3 chromosomes. The
extrapolations built on top — e.g. m5's "~10,000-SNP crossover" for when
fitted-V_c uncertainty dominates (m5_inference.py:132), an extrapolation from
M=80 assuming exact 1/√M scaling — should be read as order-of-magnitude
arguments, not measurements.

### 5.4 Inference hygiene

Mostly good, with specific lapses:

- **Disclosed post-hoc re-design:** the map's averaging basis was raised to
  n_geno = 6 *after* an early run showed a spurious relatedness "effect"
  (experiments/README.md:207-209). Disclosed honestly — but it is
  outcome-driven, and 6 realizations is a small basis either way.
- **In-sample fix validation:** m5's fitted-V_c SE correction is fit
  (`polyfit` on ĥ² vs intercept) and its "restored coverage" reported on the
  same 300 replicates (m5_inference.py:106-108).
- **Provenance mismatch:** experiments/README.md:25-28 quotes Readout-B
  numbers "(0.997, 6.246e-4) vs τ² = 6.229e-4" that the script's Readout-B
  block does not print; they come from the Readout-A2 "correct regressor" row
  (m_spike lines 283-286). Mislabeled, not fabricated, but should be fixed.
- m3's README table shows PC counts 0/2/10/20/50; the script runs
  [2,5,10,20,50].
- The plan's own gate 9 (metrics fixed "before held-out seeds are run") has
  no held-out-seed protocol in any script.
- m_spike's per-SNP MCSE ignores cross-SNP correlation of shared Monte-Carlo
  draws (the plan's gate 2 demands "covariance-aware Monte Carlo
  uncertainty"); the RMS-1.55σ outcome is then rationalized post hoc in the
  README. The aggregate intercept ≈ 1 does rescue the conclusion.

### 5.5 What the experiments *do* establish soundly

- The algebra of the implementation (the tests and m_spike genuinely derisk
  coding errors in P_c, q_j, the moment).
- The M4 fitted-component plug-in error at N=400 (intercept 1.010 ± 0.078) —
  a real, non-circular-in-spirit measurement, albeit in a correctly specified
  family.
- The qualitative ordering of methods across structure regimes (ordinary <
  cov-LDSC < oracle) — unsurprising but cleanly demonstrated.
- The honest negative result that the LDSC-application angle reduces to
  cov-LDSC (with the §5.2 caveat that it is structural, not measured).

---

## 6. The plan

### 6.1 Strengths — this is the best part of the project

- **Falsification-first, kill-or-continue gated, and the gates were actually
  obeyed.** The M-spike ran before packaging; when the use-case investigation
  came back negative, the plan absorbed it (§1.1) instead of rationalizing.
  "No release milestone is awarded because the version number looks lonely"
  is exactly the right culture.
- **Preregistered promotion gates with explicit estimand labels** (the
  τ²-not-h² discipline, the Huang-2025 admixture caveat) show real
  statistical maturity.
- **Scope honesty:** "framework/synthesis, not a production estimator,"
  "proposed, not yet validated," "software novelty is nil and is not
  claimed." PRIOR_ART.md downgrades its own novelty claims after a scan
  ("thinner than the first pass suggested") — almost no one does this.
- The decision rule (within-family kill-test → else theory note → else stop)
  is concrete, cheap, and correctly ordered.

### 6.2 Weaknesses and risks

- **Churn.** The plan was restructured twice in two days (lean-MVP →
  LMM-connection → decomposition reframing). Each move was data-driven and
  honestly recorded, but the moving target means several plan sections no
  longer match reality: the milestone tables still describe the estimator
  arc as the spine, while the project now says the spine is the
  decomposition + diagnostics. A consolidating rewrite of PLAN §3 is due.
- **The remaining niche is unvalidated and may not survive contact.**
  Within-family/sib-GWAS is the one differentiated direction left, and the
  plan correctly gates it behind a cheap kill-test. Be aware the same
  structural argument that killed the LDSC angle may apply: within-family
  designs already remove relatedness by construction (sib-difference), so the
  exact-covariance advantage may again live upstream of anything LDSC-like.
  Run the kill-test before any further investment.
- **The prior-art scan is provisional by its own admission** (session-limited,
  Bercovich 2025 not read in full). Given that Bercovich is the closest
  competitor to the current framing, the novelty section of any write-up
  should not be finalized before that read.
- **The real-data survey is only scoped** — no real data has touched the
  framework. All "value" statements are simulation-based, and (§5.3) the
  simulations are toy-genetic. The survey's gating step (efficient P_c
  operator) is correctly identified as the prerequisite; nothing should be
  claimed about real GWAS until then.
- **GEMMA parity is deferred**, reasonably (no Windows build), and the
  protocol document is good — but until it runs, "reproduces an established
  exact-LMM engine" rests on the local REML, not on GEMMA.
- Plan-vs-repo drift: PLAN §5 says "no `.github/workflows`" (added since);
  the dependency-floor item is done (pyproject now 3.11+/NumPy 2.0+). Minor,
  but the plan advertises itself as current.

---

## 7. Novelty assessment (independent view)

Agreeing with, and slightly sharpening, PRIOR_ART.md:

- **Published:** structure inflates long-range LD and biases LDSC
  (Price 2008, Bhatia 2016); correcting LD with PC-residualized genotypes
  (cov-LDSC 2021; Bercovich 2025 — closest); per-PC stratification mapping
  (Bhatia 2016, signed, on β̂); the admixed-h² estimand caveat (Huang 2025).
- **Not new although presented as empirical:** the genome-wide LD score as an
  eigenvalue-weighted sum of squared PC loadings — that is the spectral
  decomposition of the GRM (§3.3), exact algebra. The *observation* that the
  leading-PC truncation of it correlates 0.99 with the non-local remainder is
  a nice illustration, not a finding.
- **Genuinely contributed (thin but real):** (i) the explicit statement and
  use of the identity as the *organizing* device tying cov-LDSC / WhiLD /
  PC-score / P_c onto one ladder; (ii) the q_k-weighted regressor as the
  exact LMM-matched object, together with the honest negative theorem of
  practice — relatedness does not propagate into windowed LDSC, so the LMM
  gap cov-LDSC flagged is closed by cov-LDSC; (iii) the chi²-scale per-PC
  spectrum and the squared-loading S-LDSC annotation as *sensitivity*
  diagnostics, with the identifiability limit stated as a first-class caveat.
- **Realistic ceiling:** a unifying-perspective methods note. The honest bar
  the project set itself — "is identity + unification + diagnostics enough of
  a clarifying contribution?" — is the right bar, and it is borderline-
  passable only if the decomposition is restated exactly (§3.3) and at least
  one diagnostic is shown on real data (the REAL_DATA_SURVEY plan).

---

## 8. Recommendations

Ordered by leverage:

1. **Restate the decomposition exactly.** Replace "ℓ^struct is the PC score
   (validated corr = 0.99)" with the exact identity
   Σ_k r_jk² = (M/N) Σ_p (λ_p/λ̄) corr(x_j, u_p)² and report magnitude
   agreement (not only correlation) for the leading-PC truncation. Lead with
   the windowed ~5% figure, not the genome-wide ~50% — the 5% is both the
   honest number and the actual explanation of the project's own central
   negative result.
2. **Run one non-circular experiment before any write-up.** METHOD §6's own
   required axis: sparse/directional background effects on the other
   chromosomes (ε ≠ (0, V_c)). This is the only experiment that can turn
   "the theory is correct where its assumptions hold" into "here is where
   the theory bends." The machinery already exists.
3. **Fix meiosis in the simulator** (per-meiosis crossovers at the chromosome's
   genetic length) before relying on any relatedness conclusion; add
   ancestry tracts to `simulate_admixture` if admixture LD is ever discussed.
4. **Run the within-family kill-test now** (it is cheap) and be prepared for
   the structural answer; otherwise converge on the theory note. Do not open
   the efficient-operator / real-data work (REAL_DATA_SURVEY stage 1) until
   that decision is made — it is the largest remaining time sink and is only
   justified by the framework track surviving.
5. **Inference hygiene:** implement the plan's held-out-seed protocol; move
   m5's SE-correction validation out-of-sample; fix the Readout-B provenance
   and the m3 PC-count mismatch in experiments/README.md; either implement
   covariance-aware MCSE or drop the plan's gate-2 wording.
6. **Presentation:** present the map experiment's relatedness result as a
   structural argument (with the OLS-statistic arm added so the contrast is
   visible); stop calling the identifiability limit a theorem; label the
   two-step weights as simplified-LDSC.
7. **Repo hygiene:** gitignore `build/` and `*.egg-info`; consolidate
   DEVELOPMENT_PLAN §3 to match the decomposition framing; complete the
   Bercovich 2025 read before finalizing novelty claims.

---

## Appendix A — verification performed for this review

- Read in full: README.md, all six docs, all four package modules, all three
  test files, experiments/README.md, pyproject.toml, CI workflow.
- Read in full (via a delegated static review, confirmed against the files):
  all eight experiment scripts, with file:line cross-checks of the claims
  quoted in §5.
- Ran the full test suite in the project's conda env (`pldsc`, Python 3.14.6,
  NumPy 2.4.6): **29 tests, all pass** (0.29 s).
- Re-derived Eq. M6 and properties P1–P4 on paper; agree with METHOD §2–4.
- Independently re-derived and numerically verified the spectral identity
  underlying the decomposition claim (§3.3) on the project's own simulator
  with an unrelated seed: exact to 1.4% (standardization slack); top-10-PC
  truncation reproduces corr = 0.993 at ~37% magnitude.
- Did **not** do: run the long Monte-Carlo experiments end-to-end (minutes
  each; static reading sufficed for the claims checked), read Bercovich 2025
  or the other cited papers in full, or audit the 2017 legacy code.
