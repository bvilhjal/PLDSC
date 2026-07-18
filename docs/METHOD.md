# PLDSC method: the mixed-model-matched LD score

Date: 2026-07-17

This is the M0 "freeze the theorem" document (see
[the development plan](DEVELOPMENT_PLAN.md)). It states the model, derives the
oracle conditional moment on which PLDSC rests, records the exact assumptions,
and shows the reductions to covariate-adjusted and ordinary LD Score regression.
The identity in Section 4 has been derived independently three times and
confirmed numerically by the M-spike (`experiments/m_spike_kill_or_continue.py`;
Section 11).

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
\(E[u_j^2]=q_j\) and is a required simulation axis.

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

## References

1. Bulik-Sullivan et al. (2015), [LD Score regression distinguishes confounding
   from polygenicity](https://doi.org/10.1038/ng.3211).
2. Luo et al. (2021), [Estimating heritability and its enrichment in
   tissue-specific gene sets in admixed populations](https://doi.org/10.1093/hmg/ddab130).
3. Huang et al. (2025), [Interpreting SNP heritability in admixed
   populations](https://doi.org/10.1093/genetics/iyaf100).
