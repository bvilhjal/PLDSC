# PLDSC method: the mixed-model-matched LD score

Date: 2026-07-17

This is the M0 "freeze the theorem" document (see
[the development plan](DEVELOPMENT_PLAN.md)). It states the model, derives the
oracle conditional moment on which PLDSC rests, records the exact assumptions,
and shows the reductions to covariate-adjusted and ordinary LD Score regression.
The identity in Section 4 has been derived independently three times and
confirmed numerically by the M-spike (`experiments/m_spike_kill_or_continue.py`;
Section 11).

**Organizing framework (reframed 2026-07-18).** These results are best read
through the local/non-local decomposition of the LD score in **Section 14**: the
LD score splits into a linkage component and a structure component, the latter
equal to the PC score, and every structure-correction below (cov-LDSC, WhiLD,
the PC score, \(P_c\)) is an operation on that decomposition.

## 1. Model and notation

Fix a target chromosome \(c\). Let

- \(y\in\mathbb{R}^{N}\) be the quantitative phenotype on \(N\) individuals;
- \(C\in\mathbb{R}^{N\times K}\) the complete fixed-effect design used by the
  GWAS, including the intercept (rank \(K'\le K\));
- \(X_c=[x_1,\dots,x_{M_c}]\in\mathbb{R}^{N\times M_c}\) the modeled
  target-chromosome genotypes, columns centered and scaled by a frozen
  convention;
- \(K_{-c}\in\mathbb{R}^{N\times N}\) a genetic relationship matrix built from
  the **other** chromosomes (leave-one-chromosome-out), symmetric PSD with a
  frozen normalization.

The null covariance for the target-chromosome score test is

\[
V_c=\sigma_g^2K_{-c}+\sigma_e^2I,\qquad \sigma_e^2>0, \tag{M1}
\]

so \(V_c\) is symmetric positive definite. The data-generating model is

\[
y=C\alpha+X_c\beta_c+\varepsilon,\qquad \varepsilon\sim(0,V_c), \tag{M2}
\]

with random target-chromosome effects satisfying

\[
E[\beta_k]=0,\quad E[\beta_k\beta_l]=\tau_k^2\,\delta_{kl},\quad
\beta_c\perp\varepsilon. \tag{M3}
\]

Only the first two moments in (M3) are used; **no Gaussianity of \(\beta_c\) or
\(\varepsilon\) is required** for the moment identity (it is required only if one
later wants exact sampling distributions rather than the mean of \(z^2\)).

## 2. The generalized residual precision \(P_c\)

Write \(W=V_c^{-1}\). Define

\[
P_c=W-WC(C^{\top}WC)^{+}C^{\top}W, \tag{M4}
\]

where \((\cdot)^{+}\) is the Moore–Penrose pseudoinverse (so a rank-deficient
\(C\) is handled). Equivalently \(P_c=W Q\) with the oblique projector
\(Q=I-C(C^{\top}WC)^{+}C^{\top}W\).

**Property P1 (annihilates the fixed design): \(P_cC=0\).**
For full-rank \(C\), \((C^{\top}WC)^{+}=(C^{\top}WC)^{-1}\) and
\(P_cC=WC-WC(C^{\top}WC)^{-1}(C^{\top}WC)=0\). For rank-deficient \(C\), put
\(B=W^{1/2}C\); then \(C^{\top}WC=B^{\top}B\) and
\(B(B^{\top}B)^{+}B^{\top}B=B\) (a standard pseudoinverse identity), so
\(WC(C^{\top}WC)^{+}(C^{\top}WC)=WC\) and again \(P_cC=0\). By symmetry (P3),
\(C^{\top}P_c=0\) as well.

**Property P2 (\(V_c\)-idempotence): \(P_cV_cP_c=P_c\).**
\(Q\) is idempotent: with \(M=(C^{\top}WC)^{+}\),
\(Q^2=I-2CMC^{\top}W+CM(C^{\top}WC)MC^{\top}W=I-CMC^{\top}W=Q\), using
\(M(C^{\top}WC)M=M\). Then
\(P_cV_cP_c=(WQ)V_c(WQ)=WQ(V_cW)Q=WQ\,Q=WQ^2=WQ=P_c\).

**Property P3 (symmetry): \(P_c^{\top}=P_c\).**
\(P_c^{\top}=W-WC(C^{\top}WC)^{+}C^{\top}W=P_c\) since \(W\) and \(M\) are
symmetric.

**Property P4 (positive semidefinite):** \(P_c\succeq0\), with
\(x^{\top}P_cx=0\iff x\in\operatorname{span}(C)\). Hence the self-weight
\(q_j\) below is nonnegative and vanishes only for variants explained by \(C\).

## 3. The matched score statistic

For variant \(j\) define the self-weight and standardized score

\[
q_j=x_j^{\top}P_cx_j,\qquad
z_j=\frac{x_j^{\top}P_cy}{\sqrt{q_j}}. \tag{M5}
\]

\(x_j^{\top}P_cy\) is the efficient (generalized least squares) score for the
target-chromosome effect of variant \(j\) under the null covariance \(V_c\) with
the fixed effects \(C\) projected out; \(q_j\) is its null variance, so \(z_j\)
is the matched score z-statistic. Variants with numerically negligible \(q_j\)
have no defined score and are excluded under a documented, scale-aware
threshold.

## 4. The oracle conditional moment (the theorem PLDSC rests on)

Define the \(P_c\)-metric correlation
\(r_{jk}^{(P_c)}=\dfrac{x_j^{\top}P_cx_k}{\sqrt{q_jq_k}}\).

**Theorem.** Under (M1)–(M3), conditional on \(X\) and \(V_c\),

\[
\boxed{\;
E\!\left[z_j^2\mid X,V_c\right]
=1+\sum_{k}\tau_k^2\,\frac{(x_j^{\top}P_cx_k)^2}{q_j}
=1+\sum_{k}\tau_k^2\,q_k\,\bigl\{r_{jk}^{(P_c)}\bigr\}^2 .
\;} \tag{M6}
\]

*Proof.* Insert (M2) into the numerator of (M5) and use P1
(\(x_j^{\top}P_cC=0\)):

\[
x_j^{\top}P_cy=\underbrace{x_j^{\top}P_cC\alpha}_{0}
+\sum_k\beta_k\,(x_j^{\top}P_cx_k)+x_j^{\top}P_c\varepsilon
=\sum_k a_{jk}\beta_k+u_j,
\]

with \(a_{jk}=x_j^{\top}P_cx_k\) and \(u_j=x_j^{\top}P_c\varepsilon\). Taking the
conditional second moment and using \(\beta\perp\varepsilon\),
\(E[\beta]=0\), \(E[u_j]=0\):

- signal term \(E\bigl[(\sum_k a_{jk}\beta_k)^2\bigr]
  =\sum_{k,l}a_{jk}a_{jl}E[\beta_k\beta_l]=\sum_k a_{jk}^2\tau_k^2\)
  by (M3);
- cross term \(2\,E[(\sum_k a_{jk}\beta_k)\,u_j]=0\);
- noise term \(E[u_j^2]=x_j^{\top}P_c\,E[\varepsilon\varepsilon^{\top}]\,P_cx_j
  =x_j^{\top}P_cV_cP_cx_j=x_j^{\top}P_cx_j=q_j\) by P2.

Hence \(E[(x_j^{\top}P_cy)^2\mid X,V_c]=q_j+\sum_k a_{jk}^2\tau_k^2\). Dividing
by \(q_j\) and substituting \(a_{jk}^2=q_jq_k\{r_{jk}^{(P_c)}\}^2\) gives (M6).
\(\square\)

Two structural facts follow immediately: the intercept is **exactly 1** (from
P2, i.e. the noise contributes exactly its null variance), and the \(j=k\) term
is present and contributes \(\tau_j^2q_j\) (since \(r_{jj}^{(P_c)}=1\)).

## 5. The candidate LD score and the role of \(q_k\)

Under a homoscedastic effect model \(\tau_k^2\equiv\tau^2\), (M6) becomes a
simple line in a single regressor,

\[
E[z_j^2\mid X,V_c]=1+\tau^2\,L_j^{(P_c)},\qquad
L_j^{(P_c)}=\sum_{k\in W_j}q_k\bigl\{r_{jk}^{(P_c)}\bigr\}^2, \tag{M7}
\]

where \(W_j\) is a locality window (Section 10); \(W_j=\){all target-chromosome
\(k\)} recovers the exact moment. **The \(q_k\) weights do not cancel in
general.** They cancel iff \(q_k\) is constant in \(k\), which happens under a
per-\(P_c\)-metric standardization \(x_k\mapsto x_k/\sqrt{q_k}\); but that same
rescaling turns homoscedastic \(\tau_k^2\) into heteroscedastic \(q_k\tau^2\),
so one cannot simultaneously have homoscedastic effects and constant \(q_k\)
unless \(q_k\) is already constant (e.g. ordinary LDSC, Section 7). The weights
therefore carry genuine information about the effect-size scale in the \(P_c\)
metric, and are the primary point of departure from covariate-adjusted LDSC.

## 6. Exactness assumptions

(M6) is exact only if all of the following hold; each is a precondition, not
bookkeeping:

1. **Oracle \(V_c\).** The covariance generating \(\varepsilon\) equals the
   \(V_c\) inside \(P_c\). This is what forces the intercept to 1 (P2). A fitted
   or misspecified \(V_c\) moves the intercept; that displacement is the
   diagnostic signal (Section 9).
2. **Second-moment effect model.** Only \(E[\beta_k]=0\) and
   \(E[\beta_k\beta_l]=\tau_k^2\delta_{kl}\) are used. Full independence and
   Gaussianity are not needed for (M6).
3. \(\beta_c\perp\varepsilon\).
4. **Same \(C\)** in the GWAS and in \(P_c\).
5. The sum in (M6) ranges over the **causal** target-chromosome set; an omitted
   causal target SNP biases the moment.
6. \(q_j>0\) (P4).

A violation of (1) or the infinitesimal-\(V_c\) structure by large,
sparse, or directional effects on the **other** chromosomes (which enter
\(\varepsilon\) and can be in structure-LD with \(x_j\)) breaks the noise term
\(E[u_j^2]=q_j\). This is now tested (`experiments/misspecified_background.py`):
the null intercept is **robust** to sparse (down to 0.5% causal → 0.993) and
directional (top-axis variance doubled → 1.004) off-chromosome backgrounds,
because both live mostly in the low-rank structure subspace \(P_c\) already
whitens. It remains sensitive to variance-component/scale misspecification (M4a),
and a background structure axis *absent* from the GRM used to build \(V_c\) — a
misspecification orthogonal to the whitened subspace — is the harder, still-untested
case.

## 7. Reductions

**Covariate-adjusted LDSC (cov-LDSC).** If \(V_c=\sigma^2I\) then
\(P_c=\sigma^{-2}\bigl(I-C(C^{\top}C)^{+}C^{\top}\bigr)=\sigma^{-2}P_C\), so
\(r_{jk}^{(P_c)}=\operatorname{Cor}(\widetilde x_j,\widetilde x_k)\) with
\(\widetilde X=P_CX\) — exactly the covariate-adjusted LD score of Luo et al.
(2021).

**Ordinary LDSC.** With only an intercept and column-standardized genotypes,
\(P_c\propto I-\tfrac1N\mathbf{1}\mathbf{1}^{\top}\) and
\(q_k=\sigma^{-2}(N-1)\) is constant, so it cancels into the slope and
\(L_j^{(P_c)}\to\sum_k r_{jk}^2\), the Bulik-Sullivan et al. (2015) LD score, up
to a known scale.

Both reductions are the M0 acceptance checks and both hold.

## 8. The legacy WhiLD score as a special case

The 2017 WhiLD score computed ordinary LDSC on genotypes whitened by the
Cholesky factor of the inverse of a \(w\)-weighted hybrid GRM/covariance matrix
\(M\). On paper this is (M7) with metric \(M^{-1}\) in place of \(V_c^{-1}\),
weights \(q_k\equiv1\), and no covariate projection. So \(P_c\) supplies the
derivation WhiLD lacked, and the two constructions sit on one ladder:

\[
\text{ordinary}\ \to\ \text{correct metric }V_c^{-1}\ \to\
\text{add projection }C\ \to\ \text{add }q_k\ =\ \text{PLDSC}.
\]

The ablation isolates which correction buys calibration. The unstable
pseudoinverse-plus-Cholesky **estimator** is not implemented; WhiLD is a
paper-only comparator.

## 9. Estimand and the calibration diagnostic

Under the homoscedastic model, the regression slope of \(z_j^2\) on
\(L_j^{(P_c)}\) estimates \(\tau^2\): the per-\(P_c\)-metric effect-size variance,
an \(h_g^2/M\)-type quantity on the declared genotype scale. It is **not** total
realized heritability, and PLDSC labels it as \(\tau^2\), not \(h^2\)
(directional-LD covariance can separate the two in admixed samples; Huang et
al. 2025).

The intercept is the deliverable. Under an oracle \(V_c\) it is exactly 1
(Section 4). Under a fitted or deliberately misspecified \(V_c\), the noise term
becomes \(x_j^{\top}P_c\Sigma_{\text{true}}P_cx_j/q_j\neq1\), and its
displacement from 1 is the residual confounding that survives the modeled
mixed-model geometry — a mixed-model analog of the LDSC intercept that a REML
variance-component fit does not provide. This was the *intended* primary
use-case; Section 13 records the empirical verdict that, in the
summary-statistics LDSC-application setting, it reduces to cov-LDSC and gives
PLDSC no differentiated edge.

## 10. Window locality

Locality lives in the \(P_c\) metric, not in centimorgans: \(W_j\) should be a
threshold on \(|r_{jk}^{(P_c)}|\), and the omitted-window remainder relative to
the whole-chromosome sum is itself a diagnostic of \(V_c\) misspecification. In
a structured sample the \(P_c\)-metric neighborhood need not coincide with the
physical-map neighborhood.

**Empirical locality finding.** Within a physical (cM) window, however,
\(r_{jk}^{(P_c)}\approx\) the raw local correlation: tight recombination-driven
LD is nearly orthogonal to the few genome-wide structure axes that \(V_c^{-1}\)
down-weights, so \(P_c\)'s action is concentrated in the *long-range /
between-locus* block — exactly the block a local LD-score window discards. The
map experiment (Section 11) confirms the consequence: cryptic relatedness, even
with 90% of the sample in sib-pairs, does **not** measurably change the windowed
LD-score regression on mixed-model summary statistics, because the LMM already
absorbs relatedness in the association step. The one regime that does bias the
windowed regression is population structure — which cov-LDSC (Section 7) already
corrects.

## 11. Numerical confirmation

The M-spike (`experiments/m_spike_kill_or_continue.py`, env `pldsc`, PCG64 seed
20260717, \(N=400\), 6000 Monte-Carlo replicates) confirms (M6): regressing the
empirical \(E[z_j^2]\) on the full-chromosome \(L_j^{(P_c)}\) recovers
\((\text{intercept},\text{slope})=(0.997,\ \tau^2)\); the mixed-model intercept
under a null target chromosome is \(1.00\) in both a structured and a
homogeneous cohort, while ordinary OLS inflates to \(9\text{–}12\); and
\(P_c\) shrinks structure-induced cross-chromosome \(r^2\) from \(0.053\) to
\(0.0015\). See `experiments/README.md`.

The miscalibration map (`experiments/map_ldsc_miscalibration.py`) regresses the
exact analytic moment on conventional LD scores across regimes (fixed \(N\),
fixed \(\tau^2\), averaged over genotype realizations). On mixed-model summary
statistics, ordinary LDSC underestimates the slope by \(\sim30\%\) at
\(F_{st}=0.15\) and inflates the intercept, while **both** cov-LDSC and PLDSC
recover it; cryptic relatedness produces no measurable drift for any method.
This is the empirical basis for the scope verdict in Section 13.

## 12. Deferred to later milestones

- Fitted-\(V_c\) plug-in error: validated by a local REML (M4b; intercept
  \(1.010\pm0.078\)); external GEMMA byte-parity awaits a Linux environment
  (protocol in [GEMMA_PARITY.md](GEMMA_PARITY.md)).
- At-scale propagation of variance-component uncertainty into the intercept SE
  (M5 quantified the requirement; not yet implemented in the estimator).
- The within-family / sib-GWAS matched-covariance direction (Section 13) — the
  one differentiated niche, proposed but not yet validated.

## 13. Scope: what the LMM connection does and does not buy

The *theoretical* connection is exact and validated (Sections 2–4, 11): the
\(P_c\)-matched score is the LOCO-LMM score statistic, and its LD score is the
random-effect-adjusted LD that cov-LDSC (Luo et al. 2021) explicitly left as
future work. The *empirical* verdict on that gap, for a new summary-statistics
estimator, is negative:

- The genuinely random-effect part of the LMM null — **relatedness** — is
  handled by the LMM in the association step and does not propagate into the
  local, windowed LD-score regression (Sections 10–11).
- The part that does bias LDSC on mixed-model summary statistics is **population
  structure**, which cov-LDSC already corrects.

Hence, in the LDSC-application setting, the random-effect-adjusted LD score
**reduces in practice to cov-LDSC**, and PLDSC provides no differentiated
advantage there. PLDSC's genuine differentiated value is *upstream* — at the
association step, and specifically in **within-family / sib-GWAS designs** where
relatedness is definitional and exactly known, so \(P_c\) is exact (no
fitted-\(V_c\) fragility, no reference-panel approximation). That direction is
proposed, not yet validated. This section records the honest scope so the theory
(correct and clean) is not overstated as a practical LDSC replacement.

## 14. The local/non-local decomposition (organizing framework)

Reframed 2026-07-18. The results above are organized by a single decomposition of
the LD score into a linkage (local) and a structure (non-local) component:

\[
\ell_j=\sum_k r_{jk}^2\;\approx\;
\underbrace{\ell_j^{\text{local}}}_{\sum_k (r_{jk}^{(P_c)})^2}
+\underbrace{\ell_j^{\text{struct}}}_{\text{structure}},
\]

where \(r_{jk}^{(P_c)}\) is the structure-whitened correlation (Sections 2, 7).
The non-local component has an **exact** form, not an approximation. With
column-standardized genotypes \(X\) (\(N\times M\)), \(r_{jk}=x_j^\top x_k/N\) and
GRM \(K=XX^\top/M\) with eigenpairs \((\lambda_p,u_p)\); then the genome-wide LD
score is exactly the **spectral decomposition of the GRM**:

\[
\sum_k r_{jk}^2=\frac{M}{N^2}\,x_j^\top K x_j
=\frac{M}{N}\sum_p \lambda_p\,\operatorname{corr}(x_j,u_p)^2 .
\]

This is textbook linear algebra, not an empirical discovery. The 2017 **"PC
score"** \(s_j=\sum_{p\le k}(\lambda_p/\lambda_t)\operatorname{corr}(g_j,u_p)^2\)
is precisely the **leading-\(k\)-eigenvector (renormalized) truncation** of this
identity — an estimate of \(\ell^{\text{struct}}\) from the top structure axes.

**What the experiment actually shows** (`experiments/ld_decomposition.py`;
re-derived and re-checked independently, EXTERNAL_REVIEW §3.3): the top-10-PC
truncation correlates \(0.99\) with the non-local remainder, **but captures only
\(\approx37\%\) of its magnitude** (correlation is scale-invariant; the truncation
misses the eigenvalue tail). The additive split that script prints
(\(\ell^{\text{ord}}=\ell^{\text{local}}+\ell^{\text{struct}}\)) is **by
construction** — \(\ell^{\text{struct}}\) is *defined* as the difference — so it is
not evidence of additivity; the only empirical content is that \(0.99\)
correlation. **Honest headline:** within the *windowed* LD score LDSC actually
uses, structure is only \(\approx5\%\) (genome-wide, cross-chromosome-inclusive,
it is \(\approx50\%\), but LDSC never regresses on that). The small windowed share
is exactly *why* cov-LDSC suffices and why the estimator angle died (Section 13).

**Unification.** Every structure-correction is an operation on this
decomposition: cov-LDSC (\(V_c=\sigma^2 I\)) and \(P_c\) compute
\(\ell^{\text{local}}\); the legacy WhiLD score is the special case of Section 8;
the PC score is the leading-PC truncation of \(\ell^{\text{struct}}\). Standard
LDSC's structure bias is
the conflation of \(\ell^{\text{local}}\) (heritability) with
\(\ell^{\text{struct}}\) (confounding); because \(\ell^{\text{struct}}\) is a small
fraction of a *local* window it surfaces mainly in the intercept and long-range
LD, not the windowed slope (consistent with Sections 10-11 and
`map_ldsc_miscalibration.py`).

**Applications (diagnostics).** Decomposing \(\ell^{\text{struct}}=\sum_p(\cdot)\)
gives a per-PC inflation spectrum (`detect_pc_inflation.py`); adding
\(\ell^{\text{struct}}\)/squared PC loadings as an S-LDSC annotation gives a
structure-robustness test for functional enrichment (`sldsc_pc_annotation.py`).

**Remark (structure aliasing — the identifiability limit).** Take marginal,
structure-*uncorrected* per-SNP statistics \(z_j=x_j^\top y/\sqrt N\) under
\(y=Xb+\sum_{p\in S}u_p\gamma_p+e\): an isotropic polygenic background
(\(E[bb^\top]=\tau^2 I\)); structure-aligned coefficients
\(\gamma_p=\gamma_p^{\mathrm{gen}}+\gamma_p^{\mathrm{env}}\) — a genuine
ancestry-differentiated *directional* genetic effect plus environmental
stratification along axis \(p\), variances \(V_p^{\mathrm{gen}},V_p^{\mathrm{env}}\);
mutually independent \(b,\gamma,e\); \(\sigma_e^2=1\). Using
\(x_j^\top u_p=\sqrt N\,\operatorname{corr}(x_j,u_p)\) (which needs \(x_j\)
standardized and \(\bar u_p=0\) — both guaranteed by column-centering, since then
\(K\mathbf 1=0\)) and the spectral identity above,
\[
E[z_j^2]=1+\underbrace{\tau^2 N\,\ell_j}_{\text{ordinary-LDSC signal}}
+\sum_{p\in S}\bigl(V_p^{\mathrm{gen}}+V_p^{\mathrm{env}}\bigr)\operatorname{corr}(x_j,u_p)^2,
\]
and the ordinary term *already* deposits genuine heritability on every structure
axis: \(\tau^2 N\ell_j=\tau^2 M\sum_p\lambda_p\operatorname{corr}(x_j,u_p)^2\), so
the full coefficient on \(\operatorname{corr}(x_j,u_p)^2\) is
\(\tau^2 M\lambda_p+V_p^{\mathrm{gen}}+V_p^{\mathrm{env}}\) (numerically confirmed:
an isotropic-only background inflates \(E[z_j^2]\) along the top PC at exactly this
rate). The confounding \(V_p^{\mathrm{env}}\) and the *directional* genetic
variance \(V_p^{\mathrm{gen}}\) enter every \(z_j\) through the identical realized
scalar \(\gamma_p\); they are therefore **exactly aliased**, not merely collinear
— no moment, joint distribution, sample size \(N\), or SNP count \(M\) of
structure-only statistics separates them.

**Scope.** The aliasing binds only the *rank-1, axis-directional* component. Two
genuine-signal channels stay identifiable and must not be claimed as lost:
(i) the isotropic footprint \(\tau^2 M\lambda_p\) is pinned by the local LDSC
slope \(\tau^2\); (ii) a *heteroscedastic-along-axis* architecture (per-SNP
\(h^2\) elevated for high-loading SNPs) rides on local LD and leaves a
\(\operatorname{corr}(x_j,u_p)^2\times\ell_j^{\mathrm{local}}\) interaction that a
fixed environmental axis shift does not. So the honest statement is
*non-identifiability from the second moments of structure-only statistics*,
confined to the directional component.

**Corollary.** Every structure-only correction (cov-LDSC, \(P_c\), PC-score
subtraction, an S-LDSC PC annotation, Bercovich residualization) removes the whole
\(\operatorname{corr}(x_j,u_p)^2\) contribution — so it removes the *directional*
genuine heritability \(V_p^{\mathrm{gen}}\) (and the aligned share of
\(\tau^2 M\lambda_p\)) together with \(V_p^{\mathrm{env}}\). It is conservative:
it de-confounds only by co-attenuating real directional ancestry-aligned \(h^2\);
the LD-borne heritability survives. Diagnostics therefore *flag* axes with
\((V_p^{\mathrm{gen}}+V_p^{\mathrm{env}})>0\); they cannot *certify* confounding.

**What breaks it.** Only external information. Within-family designs do **not**
recover \(V_p^{\mathrm{gen}}\): the axis is shared between sibs (\(u_p\) lies in the
family span), so family-mean-centering annihilates it — exactly the
\(M_BRM_B=0.5\,M_B\) of `within_family_killtest.py` — nulling the aligned genetic
*and* environmental parts together and estimating a strictly smaller,
confound-robust *within-family* heritability that omits \(V_p^{\mathrm{gen}}\)
(this is why within-family \(h^2<\) population \(h^2\)). Genuine separation of
direct genetic from stratification needs cross-population replication or parental
genotypes (transmitted/non-transmitted, [Kong et al. 2018]; RDR,
[Young et al. 2018]).

**Positioning.** This is a *formalization*, in the LD-score metric, of the
long-known polygenic-signal-vs-stratification confounding — the senior author's
own "genetic background correlated with structure" point
([Vilhjálmsson & Nordborg 2013](https://doi.org/10.1038/nrg3382)) and the
polygenic-selection debate (Berg & Coop 2014; Berg 2019; Sohail 2019; Zaidi &
Mathieson 2020, adjacent); Bhatia 2016 is the *first-moment* (signed) analogue.
It is **not** a new impossibility result — and, per the prior-art sweep
(PRIOR_ART.md), **not even the corollary is new**: [Huang et al.
2025](https://doi.org/10.1093/genetics/iyaf100) already state in print that
structure/PC correction and cov-LDSC (by name) co-attenuate genuine
ancestry-aligned genetic variance, and [Sotoudeh et al.
2026](https://doi.org/10.1073/pnas.2525790123) give a first-moment version (PC
adjustment diminishes direct genetic effects). The only thin genuinely-new tile
is the *exact marginal single-regressor* form — writing the second-moment
coefficient on \(\operatorname{corr}(x_j,u_p)^2\) as one aliased sum — and even
that is exposition; Huang 2025 must be the primary predecessor, not a footnote.

**Positioning (revised after the prior-art sweep, PRIOR_ART.md).** The qualitative
concept (LD has a structure component) is established (cov-LDSC; Bhatia et al.
2016; Price et al. 2008; Bercovich et al. 2025); the exact spectral identity is
textbook algebra (a GRM quadratic form; [de Vlaming et al.
2017](https://doi.org/10.1101/211821) already put LDSC and the GRM eigen-spectrum
in one frame, and [LDER, Song et al. 2022](https://doi.org/10.1016/j.ajhg.2022.03.013)
is a spectral LDSC reformulation). **Correction the sweep forced:** the published
structure-corrected-LD operators are *not* all hard projections. Bercovich 2025
and cov-LDSC are hard rank-K projections, but **[Mangin et al.
2012](https://doi.org/10.1038/hdy.2011.73) is already an eigenvalue-weighted
GLS / mixed-model-whitened operator** (\(V^{-1/2}\) with a \(K+Q\) covariance,
pseudoinverse-projected) matched to LMM association-test power — operator-identical
to \(P_c\) at the pairwise level. So \(P_c\) is *not* the first GLS generalization;
the genuinely-thin new tile is only **carrying that GLS operator from pairwise
\(r^2\) to LD *scores* (row sums) with a LOCO covariance matched to a LOCO-LMM
score statistic** — and even that collapses to \(\approx\)cov-LDSC in the windowed
regime (Section 13). The unification-of-estimators move itself also has a precedent
([MQS/GEMMA, Zhou 2017](https://doi.org/10.1214/17-AOAS1052), fusing HE regression
and LDSC). **Honest contribution:** an *organizing map* — one LD-score spectral
decomposition under which cov-LDSC, Mangin's \(r^2_{VS}\), the WhiLD/Mahalanobis
whitened GRM, the PC score, and a LOCO-LMM-matched \(P_c\) appear as choices of
*metric* (identity → hard projection → eigenvalue-weighted GLS) × *target*
(pairwise \(r^2\) → LD-score row sum) — plus the honest negative results. This is a
clarifying/unifying exposition, not a set of novel contributions. See
[PRIOR_ART.md](PRIOR_ART.md).

## References

1. Bulik-Sullivan et al. (2015), [LD Score regression distinguishes confounding
   from polygenicity](https://doi.org/10.1038/ng.3211).
2. Luo et al. (2021), [Estimating heritability and its enrichment in
   tissue-specific gene sets in admixed populations](https://doi.org/10.1093/hmg/ddab130).
3. Huang et al. (2025), [Interpreting SNP heritability in admixed
   populations](https://doi.org/10.1093/genetics/iyaf100).
4. Vilhjálmsson & Nordborg (2013), [The nature of confounding in genome-wide
   association studies](https://doi.org/10.1038/nrg3382) (the qualitative origin
   of the structure-aliasing remark: the confounder is environment *and* the
   genetic background correlated with structure).
5. Bhatia et al. (2016), [Subtle stratification confounds estimates of
   heritability from polygenic models](https://doi.org/10.1101/048181) (the
   first-moment / signed PC-loading analogue of the §14 remark).
6. Young et al. (2018), [Relatedness disequilibrium regression estimates
   heritability without environmental bias](https://doi.org/10.1038/s41588-018-0178-9)
   (RDR — an external anchor that breaks the aliasing).
7. Kong et al. (2018), [The nature of nurture: effects of parental genotypes](https://doi.org/10.1126/science.aan6877)
   (transmitted/non-transmitted — the other external anchor).
