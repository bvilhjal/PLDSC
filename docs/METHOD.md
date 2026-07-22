# PLDSC method: the mixed-model-matched LD score

Date: 2026-07-17; M0 contract completed and re-audited 2026-07-22;
theory addendum added 2026-07-22

This is the M0 "freeze the theorem" document (see
[the development plan](DEVELOPMENT_PLAN.md)). It states the model, derives the
oracle conditional moment on which PLDSC rests, records the exact assumptions,
and states the precise conditions under which its geometry or full regression
moment reduces to covariate-adjusted and ordinary LD Score regression.
The identity in Section 4 has been derived independently three times and
confirmed numerically by the M-spike (`experiments/m_spike_kill_or_continue.py`;
Section 11).

**Organizing framework (reframed 2026-07-18, qualified 2026-07-22).** These
results are best read through the operator-dependent local/non-local
decomposition in **Section 14**. The total LD score has an exact GRM spectral
identity; calling part of that identity "structure" additionally requires a
declared window, operator, or externally justified PC subspace. Cov-LDSC,
WhiLD, the PC score, and \(P_c\) are related operations on that geometry, not a
unique algebraic separation of linkage from structure.

**Theory addendum (added 2026-07-22).** Section 15 asks what the same
moments mean under assortative mating and parental genetic nurture. Assortment
creates signed, effect-aligned long-range gametic-phase disequilibrium (GPD),
whereas pure nurture can preserve the ordinary LDSC line while making its slope
entirely non-direct. Neither mate choice nor parental nurture is causally
identified from offspring GWAS plus offspring LD alone.

## 1. Frozen scales, model, and variant sets

Fix a target chromosome \(c\). For an imputed diploid dosage \(g_k\), estimate
the allele frequency in the analysis sample and use the same HWE scaling as the
reference implementation:

\[
\widehat p_k=\frac{\mathbf 1^Tg_k}{2N},\qquad
x_k=\frac{g_k-2\widehat p_k\mathbf 1}
          {\sqrt{2\widehat p_k(1-\widehat p_k)}},\qquad
d_k=\frac{x_k^Tx_k}{N}. \tag{S1}
\]

Missing dosages are mean-imputed only after the declared genotype QC;
monomorphic variants are excluded. Thus \(\mathbf 1^Tx_k=0\), but \(d_k\) is
not forced to one. HWE scaling and sample-RMS scaling are not interchangeable
in an admixed sample. All exact formulas below retain \(d_k\); sample-RMS
standardization is the special convention \(d_k=1\).

Let \(Z_{-c}\) contain the similarly scaled variants in the LOCO-GRM set
\(B_{-c}\), and freeze the GRM and phenotype scale as

\[
K_{-c}=\frac{Z_{-c}Z_{-c}^T}{M_{-c}},\qquad
\frac{\operatorname{tr}(K_{-c})}{N}
 =\frac{1}{M_{-c}}\sum_{m\in B_{-c}}d_m,\qquad
y=\frac{y_{\rm raw}}{s_y},\quad
s_y^2=\frac{y_{\rm raw}^TM_Cy_{\rm raw}}{N}, \tag{S2}
\]

where \(M_C=I-C(C^TC)^+C^T\) and \(C\) includes an intercept. No secondary
trace normalization is applied to \(K_{-c}\). Analyses on the raw phenotype
scale are algebraically valid, but must say so; then every variance and effect
below is in raw squared-phenotype units. Under (S2), \(\beta_k\) is measured in
conditional-phenotype SD per HWE-standardized genotype unit. If background
coefficients have covariance \(\sigma_g^2I/M_{-c}\), their genetic-value
covariance is \(\sigma_g^2K_{-c}\); their average per-individual variance is
\(\sigma_g^2\operatorname{tr}(K_{-c})/N\), not necessarily \(\sigma_g^2\).

**Table 1. Variant-set contract**

| Symbol | Role | May differ from |
|---|---|---|
| \(B_{-c}\) | off-chromosome variants used only to build \(K_{-c}\) | every target-chromosome set |
| \(S_c\) | modeled target-chromosome effect/tagging universe in \(X_{S_c}\) | \(T_c,R_c\) |
| \(T_c\) | variants whose score statistics \(z_j\) are defined | \(S_c,R_c\) |
| \(R_c\subseteq T_c\) | predeclared QC-passing regression observations | \(S_c,T_c\) |
| \(W_j\subseteq S_c\) | effect variants retained by a score approximation for test variant \(j\) | the complete \(S_c\) |

Write \(M_{-c}=|B_{-c}|\) and \(M_c=|S_c|\). The dense oracle and current
small-array API specialize to
\(S_c=T_c=R_c\); that is an implementation choice, not a theorem. Regression
QC may filter \(R_c\), but it does not silently remove the corresponding terms
from \(S_c\). Let \(C\in\mathbb R^{N\times K}\) be the complete GWAS
fixed-effect design, with rank \(K'\leq K\). The oracle null covariance and
data-generating model are

\[
V_c=\sigma_g^2K_{-c}+\sigma_e^2I,\qquad \sigma_e^2>0. \tag{M1}
\]

\[
y=C\alpha+X_{S_c}\beta_{S_c}+\varepsilon. \tag{M2}
\]

Conditional on the fixed genotype/design data and the declared \(V_c\), assume

\[
\begin{aligned}
E[\beta_{S_c}\mid X,V_c]&=0,&
\operatorname{Cov}(\beta_{S_c}\mid X,V_c)&=
   \operatorname{diag}(\tau_k^2:k\in S_c),\\
E[\varepsilon\mid X,V_c]&=0,&
\operatorname{Cov}(\varepsilon\mid X,V_c)&=V_c,\\
\operatorname{Cov}(\beta_{S_c},\varepsilon\mid X,V_c)&=0.
\end{aligned} \tag{M3}
\]

Only these conditional first and second moments are needed for the mean
identity. Gaussianity is needed later for the exact covariance of squared
scores, not for Equation (M6). Genotype-derived covariates, PCs, and GRMs are
compatible with this conditioning; a phenotype-derived \(C\) or
\(\widehat V_c(y)\) is not fixed conditional on \(X\).

## 2. The generalized residual precision \(P_c\)

Write \(W=V_c^{-1}\). Define

\[
P_c=W-WC(C^{\top}WC)^{+}C^{\top}W, \tag{M4}
\]

where \((\cdot)^{+}\) is the Moore–Penrose pseudoinverse (so a rank-deficient
\(C\) is handled). Equivalently \(P_c=W Q\) with the oblique projector
\(Q=I-C(C^{\top}WC)^{+}C^{\top}W\).

The four identities used below are

\[
P_cC=0,\qquad P_cV_cP_c=P_c,\qquad P_c^T=P_c,\qquad
P_c\succeq0,\quad
x^TP_cx=0\Longleftrightarrow x\in\operatorname{span}(C). \tag{P1--P4}
\]

For completeness, put \(B=V_c^{-1/2}C\). The middle term in (M4) is the
orthogonal projector onto \(\operatorname{col}(B)\), conjugated by
\(V_c^{-1/2}\). This proves symmetry, positive semidefiniteness, and P1 even
when \(C\) is rank-deficient. The projector is idempotent, and conjugating that
identity back gives P2. Hence the self-weight below is nonnegative and vanishes
only for a tested genotype already spanned by \(C\).

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

**Theorem.** Under (M1)–(M3), for every \(j\in T_c\), conditional on \(X\) and
the fixed oracle \(V_c\),

\[
\boxed{\;
E\!\left[z_j^2\mid X,V_c\right]
=1+\sum_{k\in S_c}\tau_k^2\,\frac{(x_j^{\top}P_cx_k)^2}{q_j}
=1+\sum_{k\in S_c}\tau_k^2\,q_k\,\bigl\{r_{jk}^{(P_c)}\bigr\}^2 .
\;} \tag{M6}
\]

*Proof.* Insert (M2) into the numerator of (M5) and use P1
(\(x_j^{\top}P_cC=0\)):

\[
x_j^{\top}P_cy=\underbrace{x_j^{\top}P_cC\alpha}_{0}
+\sum_{k\in S_c}\beta_k\,(x_j^{\top}P_cx_k)+x_j^{\top}P_c\varepsilon
=\sum_{k\in S_c} a_{jk}\beta_k+u_j. \tag{M6a}
\]

with \(a_{jk}=x_j^{\top}P_cx_k\) and \(u_j=x_j^{\top}P_c\varepsilon\). Taking the
conditional second moment and using the conditional zero means and zero
effect-noise covariance in (M3):

- signal term \(E\bigl[(\sum_k a_{jk}\beta_k)^2\bigr]
  =\sum_{k,l}a_{jk}a_{jl}E[\beta_k\beta_l]=\sum_k a_{jk}^2\tau_k^2\)
  by (M3);
- cross term \(2\,E[(\sum_k a_{jk}\beta_k)\,u_j]=0\);
- noise term \(E[u_j^2]=x_j^{\top}P_c\,E[\varepsilon\varepsilon^{\top}]\,P_cx_j
  =x_j^{\top}P_cV_cP_cx_j=x_j^{\top}P_cx_j=q_j\) by P2.

Hence \(E[(x_j^{\top}P_cy)^2\mid X,V_c]=q_j+\sum_{k\in S_c} a_{jk}^2\tau_k^2\). Dividing
by \(q_j\) and substituting \(a_{jk}^2=q_jq_k\{r_{jk}^{(P_c)}\}^2\) gives (M6).
\(\square\)

Two structural facts follow immediately: the intercept is **exactly 1** (from
P2, i.e. the noise contributes exactly its null variance), and, when
\(j\in S_c\), the \(j=k\) term is present and contributes \(\tau_j^2q_j\)
(since \(r_{jj}^{(P_c)}=1\)).

## 5. The candidate LD score and the role of \(q_k\)

Under a homoscedastic effect model \(\tau_k^2\equiv\tau^2\), the complete
effect-set score gives the exact single-regressor line

\[
E[z_j^2\mid X,V_c]=1+\tau^2\,L_j^{(P_c)},\qquad
L_j^{(P_c)}=\sum_{k\in S_c}q_k\bigl\{r_{jk}^{(P_c)}\bigr\}^2. \tag{M7}
\]

A windowed version replaces \(S_c\) by \(W_j\) and is an approximation
(Section 10). The \(q_k\) terms do not cancel. If \(q_k=q_0\) for all \(k\), the
common factor is absorbed into the slope. Variantwise normalization
\(x_k^*=x_k/\sqrt{q_k}\) instead transforms the coefficient to
\(\beta_k^*=\sqrt{q_k}\beta_k\), so
\(\operatorname{Var}(\beta_k^*)=q_k\tau^2\). An unweighted row sum therefore
represents the same homoscedastic architecture only when the original \(q_k\)'s
are constant; otherwise it silently changes the prior.

The theorem also supplies the annotation extension. If
\(\tau_k^2=\sum_{a=1}^A\theta_a a_{ka}\), with the resulting variance
nonnegative for every \(k\), then

\[
E[z_j^2\mid X,V_c]=1+\sum_{a=1}^A\theta_aL_{ja}^{(P_c)},\qquad
L_{ja}^{(P_c)}=\sum_{k\in S_c}a_{ka}q_k
  \{r_{jk}^{(P_c)}\}^2. \tag{M7a}
\]

## 6. Exactness assumptions

(M6) is exact only if all of the following hold; each is a precondition, not
bookkeeping:

1. **Fixed oracle \(V_c\).** The covariance generating \(\varepsilon\) equals
   the \(V_c\) inside \(P_c\). An independent plug-in \(\widehat V_c\) may be
   conditioned on as fixed, but oracle exactness still requires it to equal the
   true residual covariance; otherwise the baseline is (M15).
2. **Conditional effect moments.** Equation (M3), including zero conditional
   means, diagonal conditional covariance, and zero effect-noise covariance,
   holds given \(X,V_c\). Independence and Gaussianity are not needed for the
   mean in (M6).
3. **Same \(C\).** The complete fixed-effect design used by the GWAS appears in
   (M4); phenotype-selected covariates require separate selection inference.
4. **Complete modeled-effect set.** The sum is over \(S_c\), whether or not its
   variants are tested. Omitted target signal leaves an omitted nonnegative
   term under the diagonal architecture. Filtering \(R_c\) does not alter that
   sum.
5. **Defined matched score.** \(j\in T_c\) and \(q_j>0\). A scale-aware
   threshold may filter \(T_c\) or \(R_c\), but not silently redefine \(S_c\).
6. **Matched statistic.** Equation (M5) is a score statistic. Wald,
   likelihood-ratio, or engine-calibrated statistics need their own moment or a
   demonstrated equivalence.

A violation of (1) or the infinitesimal-\(V_c\) structure by large,
sparse, or directional effects on the **other** chromosomes (which enter
\(\varepsilon\) and can be in structure-LD with \(x_j\)) breaks the noise term
\(E[u_j^2]=q_j\). This is now tested (`experiments/misspecified_background.py`):
the null intercept is **robust** to sparse (down to 0.5% causal → 0.993) and
directional (top-axis variance doubled → 1.004) off-chromosome backgrounds,
because both live mostly in the low-rank structure subspace \(P_c\) already
whitens. It remains sensitive to variance-component/scale misspecification (M4a).
The harder case — true-covariance structure *absent* from the GRM used to build
\(V_c\) — is also now tested (`experiments/misspecified_orthogonal.py`): an axis
orthogonal to the genotype structure axes is harmless (intercept 1.000–1.003 up to
\(v=2\sigma_e^2\); genotypes carry no loading on it), so the realistic failure
route is an under-covering GRM — a GRM built from every 10th SNP (standard
LD-thinning) triples the intercept (3.19), and family-shared environment under a
pruned GRM inflates it 4.4–7.7×. The intercept's calibration therefore depends on
the GRM *spanning* the sample's actual covariance, which folds into the known
fitted-\(V_c\) sensitivity (M4a) with GRM density as a first-order requirement.

## 7. Reductions

Let \(M_C=I-C(C^TC)^+C^T\) and
\(s_{C,k}^2=N^{-1}x_k^TM_Cx_k\).

**Covariate-adjusted geometry.** If \(V_c=\sigma^2I\), then
\(P_c=\sigma^{-2}M_C\), \(q_k=Ns_{C,k}^2/\sigma^2\), and
\(r_{jk}^{(P_c)}\) is the correlation of the residual genotypes. Hence

\[
E[z_j^2\mid X,V_c]
=1+\frac{N}{\sigma^2}\sum_{k\in S_c}
  \tau_k^2s_{C,k}^2\{r_{jk}^{(C)}\}^2. \tag{M8}
\]

Thus the **correlation geometry** is exactly cov-LDSC's fixed-covariate
geometry. The full one-slope row-sum moment is conventional cov-LDSC only if
\(\tau_k^2s_{C,k}^2\) is constant: for example, residual norms are common and
effects are homoscedastic, or an inverse-information effect prior is declared.
Arbitrary covariates generally make the residual norms differ. Luo et al.
(2021) motivate the residual-genotype geometry; Equation (M8) records the
additional scale condition needed for equality of regression moments.
Literal equivalence to a conventional cov-LDSC slope also requires the same
phenotype/null-variance normalization for the association statistic.

**Ordinary LDSC.** With only an intercept, (S1) gives
\(M_Cx_k=x_k\), \(s_{C,k}^2=d_k\), and

\[
E[z_j^2\mid X,V_c]
=1+\frac{N\tau^2}{\sigma^2}
  \sum_{k\in S_c}d_k r_{jk}^2
\longrightarrow 1+\frac{N\rho_c}{M_c}\ell_j,\qquad
\ell_j=\sum_{k\in S_c}r_{jk}^2,\quad
\rho_c=\frac{M_cd_0\tau^2}{\sigma^2}. \tag{M9}
\]

The arrow is exact under sample-RMS scaling \(d_k=1\), or under any common
\(d_k=d_0\). Under sample-RMS scaling \(d_0=1\). Under the frozen HWE scale
it is a population/asymptotic limit when \(d_k\to1\), not a finite-sample
identity. Equation (M9) is the ordinary-LDSC-shaped moment of Bulik-Sullivan et
al. (2015) when \(\sigma^2\) has the conventional association-statistic
normalization. Otherwise \(\rho_c\) is a signal-to-residual-variance ratio, not
automatically conventional SNP heritability.
The earlier \(N-1\) claim mixed a different standardization into (S1) and is
withdrawn. In short: both limiting geometries hold; the conventional regression
moments require the stated scale/prior conditions.

## 8. Relation to the legacy WhiLD score

The 2017 WhiLD construction computed ordinary LDSC after whitening genotypes by
a covariance-derived metric. This is a useful **metric analogue** of (M7), but
not literally the same homoscedastic-effect model with \(q_k\equiv1\). Choose a
square-root factor \(A_c\) with \(A_c^TA_c=P_c\), and let
\(\widetilde x_k=A_cx_k\). Re-standardizing that whitened column changes its
coefficient scale:

\[
x_k^*=\widetilde x_k/\sqrt{q_k}
     =A_cx_k/\sqrt{x_k^TP_cx_k},\qquad
\beta_k^*=\sqrt{q_k}\beta_k,\qquad
\operatorname{Var}(\beta_k^*\mid X,V_c)=q_k\tau^2. \tag{M10}
\]

A q-free WhiLD row sum therefore corresponds to the same prior only when the
original \(q_k\)'s are common; otherwise it corresponds approximately to an
inverse-information prior on the original scale. The qualitative ablation
ladder—identity metric, covariance metric, fixed-effect projection, and explicit
\(q_k\)—remains useful. The unstable pseudoinverse-plus-Cholesky estimator is
not implemented; WhiLD remains a paper-only comparator.

## 9. Estimand, calibration, and inference

Under the complete homoscedastic model, a correctly specified regression of
\(z_j^2\) on (M7) has slope \(\tau^2\): coefficient variance per
HWE-standardized genotype, in squared phenotype units. It is not a
"per-\(P_c\)-metric" variance; \(P_c\) controls association geometry, not the
input effect scale. The corresponding model-based target-chromosome genic
variance on the (S2) phenotype scale is

\[
v_{\mathrm{genic},c}
=E[\beta_{S_c}^T(X_{S_c}^TX_{S_c}/N)\beta_{S_c}\mid X,V_c]
=\sum_{k\in S_c}\tau_k^2d_k, \tag{M11}
\]

which becomes \(M_c\tau^2\) only when \(d_k=1\). A realized effect vector also
contains directional-LD covariance,

\[
\frac{\|X_{S_c}\beta_{S_c}\|^2}{N}
=\sum_k d_k\beta_k^2
+2\sum_{k<l}\frac{x_k^Tx_l}{N}\beta_k\beta_l. \tag{M12}
\]

Consequently neither \(\tau^2\) nor (M11) is automatically total realized
heritability in an admixed population (Huang et al. 2025). For a general
correlated-effect covariance \(\Omega_\beta\), with
\(a_j=X_{S_c}^TP_cx_j\), the exact signal is

\[
E[z_j^2\mid X,V_c]
=1+q_j^{-1}a_j^T\Omega_\beta a_j, \tag{M13}
\]

which contains cross-effect terms and generally has no scalar LD-score slope.
Section 15 gives two mechanistically distinct cautions: assortative mating makes
long-range genotype covariance align with the generating trait, whereas
parental nurture can produce direct-like marginal associations even when
\(\beta_{\rm direct}=0\).

If the working score, effect architecture, or window is misspecified, a
free-intercept WLS fit over \(R_c\) targets the design-dependent projection

\[
b_{R,w}^*
=\frac{\sum_{j\in R_c}w_j(L_j-\bar L_w)(\mu_j-\bar\mu_w)}
       {\sum_{j\in R_c}w_j(L_j-\bar L_w)^2},\qquad
\mu_j=E[z_j^2\mid X,V_c]. \tag{M14}
\]

It depends on the regression set and weights; it is \(\tau^2\) only when the
mean model is correct. Fixing the intercept at one changes this pseudo-target.

Under the oracle, the intercept is exactly one. For a fixed misspecified
operator \(P_*\), however, the noise baseline is variant-specific:

\[
b_j(P_*)=
\frac{x_j^TP_*V_{\mathrm{true}}P_*x_j}{x_j^TP_*x_j}. \tag{M15}
\]

It can correlate with \(L_j\), so a single fitted intercept is not uniquely
"residual environmental confounding"; it can also reflect covariance or scale
misspecification. If \(\widehat P_c\) is fitted from the same phenotype, even
(M15) is not a fixed-operator conditional identity. The fit must be handled
jointly or repeated inside an outer resampling scheme. Section 13 retains the
empirical verdict for the intended application, but the algebra does not grant
the intercept a stronger interpretation.

### 9.1 Regression variance, weights, and uncertainty

Let \(Q_T=\operatorname{diag}(q_j:j\in T_c)\),
\(A=X_{T_c}^TP_cX_{S_c}\), \(B=Q_T^{-1/2}A\), and define

\[
R_T=Q_T^{-1/2}X_{T_c}^TP_cX_{T_c}Q_T^{-1/2},\qquad
\Gamma=\operatorname{Cov}(z_{T_c}\mid X,V_c)
=R_T+B\Omega_\beta B^T. \tag{M16}
\]

For the diagonal architecture, its off-diagonal entries are

\[
\Gamma_{jl}=r_{jl}^{(P_c)}+
\sum_{k\in S_c}\tau_k^2q_k
r_{jk}^{(P_c)}r_{lk}^{(P_c)},\qquad
\Gamma_{jj}=1+\sum_{k\in S_c}\tau_k^2q_k
\{r_{jk}^{(P_c)}\}^2. \tag{M17}
\]

Under jointly Gaussian effects and residual noise, Isserlis' identity gives

\[
\operatorname{Cov}(z_j^2,z_l^2\mid X,V_c)=2\Gamma_{jl}^2,\qquad
\operatorname{Var}(z_j^2\mid X,V_c)=2\Gamma_{jj}^2. \tag{M18}
\]

Non-Gaussian fourth cumulants add to (M18); they do not alter (M6). With
\(s_R=(z_j^2:j\in R_c)\), design \(D_R=(\mathbf1,L_R)\), and
\(H_R=2\Gamma_R\circ\Gamma_R\), oracle optimal linear weighting on the
estimable subspace is

\[
\widehat\theta_{\mathrm{GLS}}
=(D_R^TH_R^{+}D_R)^{-1}D_R^TH_R^{+}s_R. \tag{M19}
\]

For singular \(H_R\), (M19) additionally requires
\(\operatorname{col}(D_R)\subseteq\operatorname{col}(H_R)\) and
\(D_R^TH_R^+D_R\) nonsingular. Otherwise the zero-variance directions impose
exact constraints that must be included explicitly (or redundant outcomes must
be removed). When \(H_R\) is positive definite, \(H_R^+\) is the ordinary
inverse. Because \(H_R\) depends on unknown parameters and \(z^2\) is not
Gaussian, (M19) is an oracle linear-weighting result, not a claim of globally
efficient likelihood inference.

There is no extra standalone \(q_j^{-1}\) regression weight: standardizing the
score in (M5) already makes its null variance one. The \(q_k\)'s enter the mean
and \(\Gamma\). A feasible diagonal Gaussian working weight is

\[
w_j^{\mathrm{diag}}=\{2\widehat\mu_j^2\}^{-1},\qquad
\widehat\mu_j=\widehat\alpha+\widehat\tau^2L_j. \tag{M20}
\]

This is what `two_step=True` in `pldsc/regression.py` approximates. It ignores
the off-diagonal overlap in (M18). Conventional LDSC's additional LD-overlap
weight is not automatically valid after changing both the score and the
metric; exact GLS needs the full \(\Gamma\), and any cheaper overlap
approximation needs separate validation. This follows the distinction between
heteroskedastic and LD-dependent weighting in Bulik-Sullivan et al. (2015); a
full squared-score covariance is also used by g-LDSC (Xiong et al. 2024).

**Table 2. M0 decision on regression and inference**

| Item | M0 conclusion |
|---|---|
| Mean model | Equation (M6) is exact for fixed oracle \(P_c\) and complete \(S_c\) |
| Exact Gaussian weights | Use (M18)–(M19); they require the full score covariance, not only \(L_j\) |
| Cheap working weights | Equation (M20) is permissible for point fitting, but is not efficient GLS |
| Genomic uncertainty | Use chromosome/position-sorted, LD-aware delete blocks and recompute every preliminary fit and data-dependent weight within each replicate |
| Fitted \(V_c\) | An ordinary variant-block jackknife conditions on \(\widehat P_c\); use an outer refit/bootstrap or a joint sandwich to propagate operator uncertainty |
| Reference LD | Raw in-sample products define the theorem; finite-reference corrections and errors-in-variables remain deployment research |

The existing regression helper is therefore a fixed-operator prototype, not a
validated genomic inference engine. A block method needs blocks spanning the
dependence range, ordered genomic input, and a block-size sensitivity analysis;
"200 blocks" is a convention, not a theorem (Künsch 1989). Cross-chromosome
scores can remain dependent because they share a phenotype and use different
LOCO operators, so pooling chromosomes also needs modeling or validation.

## 10. Window approximation

The theorem itself has no locality assumption: the exact score sums over all
of \(S_c\). For a predeclared physical/genetic window \(W_j\), the omitted
signal is

\[
\Delta_j(W_j)=\sum_{k\in S_c\setminus W_j}
\tau_k^2q_k\{r_{jk}^{(P_c)}\}^2\geq0. \tag{M21}
\]

A threshold on each \(|r_{jk}^{(P_c)}|\) does not control the cumulative tail
in (M21), and because it is data-dependent it also changes the estimand.
Window convergence must therefore be checked against the full dense oracle or
a valid bound. A large omitted remainder shows truncation error; by itself it
does not diagnose misspecification of \(V_c\).

Equation (M21) contains only omitted variants in \(S_c\), the modeled target-
chromosome effect set. The off-chromosome effect partners in cross-autosome GPD
lie outside \(S_c\); they enter through the LOCO background/model assumptions
and must not be mislabeled as the (M21) window tail. Section 15.1 treats that
separate mechanism.

**Empirical locality finding.** Within a physical window,
\(r_{jk}^{(P_c)}\) was close to raw local correlation in the simulations:
tight recombination-driven LD was nearly orthogonal to the few genome-wide
axes down-weighted by \(V_c^{-1}\). This is an empirical result, not a
consequence of (M6). In real multi-ancestry 1000 Genomes data, the LOCO-clean
structure share of a ±250 kb windowed score was about 17% in the full panel and
about 0% in EUR (`experiments/realdata_1000g_benchmarks.py`, B2); a few fixed
PCs absorbed most of that difference (B3). The map experiment (Section 11)
found no measurable relatedness effect on windowed regression after LMM
association, while population structure remained the regime cov-LDSC already
corrected.

### 10.1 Relatedness in the LD panel: the residual-df correction

Section 13 concludes that relatedness "does not propagate into the local,
windowed LD-score regression." That is correct about the *slope*, but it is not
the whole story about the LD score itself, and the remaining channel is
mis-specified in the current implementation.

The finite-sample adjustment used by `pldsc.ld.local_ld_scores`,
\(\widetilde r^2=r^2-(1-r^2)/(\nu-1)\) with \(\nu=N-\operatorname{rank}(C)\)
(Equation (7), [IDEA_REVIEW.md](IDEA_REVIEW.md) §5), subtracts the noise floor
that two *unlinked* variants generate by chance. That floor is \(1/(N-1)\) only
when the rows are independent. For a panel with relationship matrix \(A\) and LD
covariates \(C\), taking the expectation over one variant at a time gives
\(E[(x_j^\top x_k)^2]=\operatorname{tr}(\widetilde A^2)\) and
\(E\|x_j\|^2=\operatorname{tr}(\widetilde A)\) with
\(\widetilde A=P_CAP_C\), so

\[
E[\widehat r_{jk}^2]\;\approx\;
\frac{\operatorname{tr}(\widetilde A^2)}{\operatorname{tr}(\widetilde A)^2}
\;=\;\frac1{N_{\mathrm{eff}}},\qquad
N_{\mathrm{eff}}:=\frac{\operatorname{tr}(\widetilde A)^2}
{\operatorname{tr}(\widetilde A^2)}.
\tag{M22}
\]

\(N_{\mathrm{eff}}\) is the participation ratio of the projected relationship
spectrum. For \(A=I\) with an intercept it is exactly \(N-1\), recovering the
usual convention up to the customary \(O(1)\) offset; with relatives it is
strictly smaller, so \(\nu-1\) **under-subtracts**. It is exposed as
`pldsc.effective_pair_df` and accepted by `local_ld_scores(..., pair_df=...)`.

**Verified** (`experiments/relatedness_effective_df.py`, real Mendelian
transmission, unlinked SNPs, \(N=400\) throughout): (M22) is accurate to
\(\le1\%\) across full-sib designs spanning a 4.7× inflation of the floor
(\(N_{\mathrm{eff}}=399\to228\to145\to84\)), while \(\nu-1\) is wrong by that
same 4.7×. The unrelated control reproduces \(1/(N-1)\).

**Why this is a correctness fix and not a heritability result.** Accumulated
over a window of \(S\) variants the omitted floor is
\(S(1/N_{\mathrm{eff}}-1/(\nu-1))\) — up to \(\approx19\) LD-score units at
\(S=2000\) in the 16-sib design, which is large relative to a real local score.
But it is *near-uniform across variants at fixed \(S\)*, so to first order it
shifts every \(\ell_j\) by a constant and therefore moves the LDSC **intercept**,
not the slope. \(h^2\) is protected; an intercept read as a confounding
barometer on a related in-sample panel is not. The residual non-uniformity
enters only through variation in \(S\) (variant density per window), which is a
density artifact that must be controlled for on its own account.

**Scope.** This is a statement about panels that contain relatives — in-sample
LD scores in family or founder cohorts. The standard fix remains to compute LD
from unrelated individuals, which sets \(N_{\mathrm{eff}}=N-1\) and makes the
correction moot; (M22) matters when that is not possible. It does not revive a
relatedness-based case for \(P_c\) as an estimator (Section 13): the channel it
describes is an intercept effect with a one-line remedy, whereas the biases
\(P_c\) was built for are slope effects.

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
  one differentiated niche; **kill-tested, NEGATIVE** (2026-07-18,
  `experiments/within_family_killtest.py`): the surviving within-family
  genetic covariance is a homoscedastic scale (\(M_BRM_B=0.5\,M_B\)) that the
  standard cluster-robust / effective-\(N\) standard error already removes.

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
advantage there. PLDSC's genuine differentiated value was then sought *upstream*
— at the association step, and specifically in **within-family / sib-GWAS
designs** where relatedness is definitional and exactly known, so \(P_c\) is
exact (no fitted-\(V_c\) fragility, no reference-panel approximation). That
direction was kill-tested (2026-07-18, `experiments/within_family_killtest.py`)
and came back **negative**: for full sibs the surviving within-family genetic
covariance is a homoscedastic scale (\(M_BRM_B=0.5\,M_B\)) that the standard
cluster-robust / effective-\(N\) standard error (Howe et al. 2022) already
removes, so the standard sib-GWAS pipeline is already calibrated and \(P_c\)
only reproduces it. This section records the honest scope so the theory
(correct and clean) is not overstated as a practical LDSC replacement.

That negative PLDSC-versus-cov-LDSC verdict is for the direct/isotropic model
tested here. The mechanisms in Section 15 are violations of causal or
architecture assumptions, not gaps repaired merely by swapping one offspring
LD metric for another.

## 14. The local/non-local decomposition (organizing framework)

Reframed 2026-07-18 and qualified 2026-07-22. The results above use a declared,
operator-dependent bookkeeping split:

\[
\ell_j=\ell_j^{(P_c)}+\Delta_j^{(P_c)},\qquad
\ell_j^{(P_c)}=\sum_k\{r_{jk}^{(P_c)}\}^2,\quad
\Delta_j^{(P_c)}=\ell_j-\ell_j^{(P_c)}.
\tag{D1}
\]

Here \(r_{jk}^{(P_c)}\) is the operator-adjusted correlation (Sections 2, 7).
Equation (D1) is exact only because the remainder is defined as the difference;
the two squared row sums do not in general form an orthogonal projection
decomposition, and \(\Delta_j^{(P_c)}\) need not be nonnegative.
Calling that remainder "structure" requires an externally justified structure
subspace or an empirical validation; it does not follow from linear algebra.

**The split is operator-dependent** (real-data finding,
`experiments/realdata_1000g_benchmarks.py` B2): an in-chromosome GRM whitens
local LD along with structure (for EUR it whitens the genome-wide score down to
the \(M/N\) noise floor), while a LOCO GRM acts differently on local LD. These
are two operational definitions, not proven upper and lower bounds on a unique
truth.

What *is* exact is a spectral identity for the **total** genome-wide correlation
score. Let \(\widetilde x_k=x_k/\sqrt{d_k}\), so
\(\widetilde x_k^T\widetilde x_k=N\), let
\(\widetilde X=[\widetilde x_1,\ldots,\widetilde x_M]\), and define
\(K_{\mathrm{corr}}=\widetilde X\widetilde X^T/M\) with eigenpairs
\((\lambda_p,u_p)\). Then

\[
\sum_k r_{jk}^2
=\frac{M}{N^2}\,\widetilde x_j^T K_{\mathrm{corr}}\widetilde x_j
=\frac{M}{N}\sum_p \lambda_p\,
  \operatorname{corr}(\widetilde x_j,u_p)^2. \tag{D2}
\]

This is textbook linear algebra, not an empirical discovery, and
\(K_{\mathrm{corr}}\) need not equal the HWE-scaled LOCO GRM in (S2). The 2017
**"PC score"** is a renormalized leading-eigenvector truncation of (D2). It is
an estimate of a structure component only when those selected eigenvectors are
independently justified as structure axes; eigenvalue rank alone does not make
that classification.

Assortative-mating GPD can appear as a diffuse low-rank contribution to (D2),
but a genotype-only eigenspace cannot assign that contribution to mate choice:
the same covariance geometry can arise from geography, ancestry, or selection
(Section 15.1).

**What the experiment actually shows** (`experiments/ld_decomposition.py`;
re-derived and re-checked independently, EXTERNAL_REVIEW §3.3): the top-10-PC
truncation correlates \(0.99\) with the non-local remainder, **but captures only
\(\approx37\%\) of its magnitude** (correlation is scale-invariant; the truncation
misses the eigenvalue tail). The additive split that script prints
(\(\ell^{\text{ord}}=\ell^{\text{local}}+\ell^{\text{struct}}\)) is **by
construction** — \(\ell^{\text{struct}}\) is *defined* as the difference — so it is
not evidence of additivity; the only empirical content is that \(0.99\)
correlation. **Honest headline (simulation, fst=0.1):** within the *windowed* LD
score LDSC actually uses, structure is only \(\approx5\%\) (genome-wide,
cross-chromosome-inclusive, it is \(\approx50\%\), but LDSC never regresses on
that). **On real 1000G data** (`experiments/realdata_1000g_benchmarks.py`) the
LOCO-clean windowed share is \(\approx0\%\) (EUR) to \(\approx17\%\)
(multi-ancestry panel) and the genome-wide share \(\approx90\%\) in the panel;
the in-chromosome operator overstates both (B2), and the share rises with MAF
(B5). The small single-ancestry windowed share
is exactly *why* cov-LDSC suffices and why the estimator angle died (Section 13).

**Unification.** Cov-LDSC, \(P_c\), WhiLD, and the PC score can all be placed on
the same metric/eigenspace map, but they do not all compute the same uniquely
defined \(\ell^{\text{local}}\) or \(\ell^{\text{struct}}\). The experiments show
that the declared remainder is small in the local windows studied, which is
consistent with Sections 10–11 and `map_ldsc_miscalibration.py`; the theorem in
Section 4 does not supply that empirical conclusion.

**Applications (diagnostics).** Decomposing \(\ell^{\text{struct}}=\sum_p(\cdot)\)
gives a per-PC inflation spectrum (`detect_pc_inflation.py`); adding
\(\ell^{\text{struct}}\)/squared PC loadings as an S-LDSC annotation gives a
structure-robustness test for functional enrichment (`sldsc_pc_annotation.py`).

**Remark (structure aliasing — the identifiability limit).** For this remark,
sample-RMS standardize the columns in (D2) and write them as \(X\), with
\(K=XX^T/M\). Take marginal, structure-*uncorrected* per-SNP statistics
\(z_j=x_j^\top y/\sqrt N\) under
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
+\sum_{p\in S}\bigl(V_p^{\mathrm{gen}}+V_p^{\mathrm{env}}\bigr)\operatorname{corr}(x_j,u_p)^2.
\tag{D3}
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
fixed environmental axis shift does not *asymptotically* — in finite samples the
environmental shift leaks a spurious interaction (median t \(\approx3.8\) at
\(N=600\)), so channel (ii) requires a calibrated null and has limited power at
realistic architectures (`experiments/local_adaptation_channel.py`, verdict
MARGINAL). So the honest statement is
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
confound-robust *within-family* heritability that omits \(V_p^{\mathrm{gen}}\).
This is one reason within-family \(h^2\) can be lower than population \(h^2\),
not a complete explanation. Direct-versus-stratification separation needs
cross-population replication or other external structure information;
direct-versus-nurture separation needs parental genotypes
(transmitted/nontransmitted, Kong et al. 2018) or another family anchor.
RDR (Young et al. 2018) can estimate a direct variance component robust to
environmental bias, but it does not by itself identify a specific parental-
nurture pathway (Section 15.2).

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

## 15. Theory addendum: assortative mating and genetic nurture

The M0 theorem conditions on a genotype matrix and a declared effect model. It
does not determine *why* offspring genotype predicts phenotype. This section
specializes established assortative-mating and indirect-genetic-effect results
to the LDSC moments above. It is an explanatory synthesis, not a new
identification theorem or a claim that the archived \(P_c\) estimator solves
either mechanism. Neither mate choice nor parental nurture is causally
identified from offspring GWAS plus offspring LD alone.

### 15.1 Assortative mating and directional GPD

Let \(G\) be a componentwise-standardized diploid genotype vector over \(M\)
loci, with \(\Sigma=\operatorname{Cov}(G)\), and let mates sort on the additive
phenotype

\[
S=G^Tb+\varepsilon_S,\qquad
a=\operatorname{Cov}(G,S)=\Sigma b,\qquad
\varepsilon_S\perp G. \tag{AM1}
\]

Under primary phenotypic assortment, equal genotype distributions in the two
sexes, linear-Gaussian conditional means, and conditional independence of mate
genotypes given their phenotypes, define
\(\rho_m=\operatorname{Corr}(S^{(m)},S^{(f)})\). Then

\[
\Sigma_{mf}:=\operatorname{Cov}(G^{(m)},G^{(f)})
=\frac{\rho_m}{\operatorname{Var}(S)}aa^T. \tag{AM2}
\]

The covariance is proportional to the **signed product** of the loci's loadings
on the mating phenotype. Here \(a_j=\operatorname{Cov}(G_j,S)\), not generally
the direct effect \(b_j\). Under alternative mating models (AM2) need not hold;
in some weak-assortment models it is only a first-order approximation.

Mendelian transmission converts cross-mate covariance into within-person
gametic-phase disequilibrium (GPD). Let \(\Xi_t\) and \(\Xi_{mf,t}\) be,
respectively, within-person and cross-mate covariances of centered allele counts
on one fixed, generation-invariant per-locus scale. For different chromosomes,

\[
(\Xi_{t+1})_{jk}
=\frac12(\Xi_t)_{jk}
+\frac14\{(\Xi_{mf,t})_{jk}+(\Xi_{mf,t})_{kj}\},\qquad
\operatorname{chr}(j)\ne\operatorname{chr}(k). \tag{AM3}
\]

If \(\Xi_{mf,t}=\Xi_{mf}\) is fixed and symmetric, \(\Xi_t\) approaches
\(\Xi_{mf}\) geometrically; at a stationary equilibrium with symmetric
\(\Xi_{mf,\infty}\), (AM3) implies
\(\Xi_{\infty,jk}=\Xi_{mf,\infty,jk}\). In the more specific
Fisher--Nagylaki model used by Border et al. (2022a), the exact equilibrium
causal-locus covariance is

\[
\Upsilon_\infty=D+2\phi\phi^T,\qquad
\phi_j=O(M^{-1/2}), \tag{AM4}
\]

where \(D\) is diagonal in that model and \(\phi\) depends on substitution
effects and mating strength. Real-genome local linkage and structure cannot be
inserted into (AM4) without a separate model.

For an analyzed trait with causal effects \(\beta\), equal to \(b\) only when it
is the mating trait, the cross-autosome contribution to genetic variance is
signed:

\[
V_{\rm GPD}^{\times}
=2\sum_{\substack{j<k\\
                  \operatorname{chr}(j)\ne\operatorname{chr}(k)}}
  \beta_j\beta_k r_{jk}. \tag{AM5}
\]

Positive assortment increases this variance for the mating trait, and can do so
for another trait aligned with its loading direction; there is no such claim
for an arbitrary \(\beta\). Squared LD scores discard the sign. If, on a
declared long-range set \(\mathcal F_j\),
\(r_{jk}=r^0_{jk}+\kappa a_ja_k\), then

\[
\Delta\ell_{j,\mathcal F}^{\rm AM}
=2\kappa a_j\sum_{k\in\mathcal F_j}r^0_{jk}a_k
+\kappa^2a_j^2\sum_{k\in\mathcal F_j}a_k^2. \tag{AM6}
\]

Only when background long-range covariance is negligible does the first term
vanish. This is a diagnostic identity, not automatically a corrected LDSC
regressor; its cross-term is another instance of the non-orthogonality in (D1).

For marginal signal \(\eta=\Sigma\beta\), the realized square is

\[
\eta_j^2
=\sum_k \sigma_{jk}^2\beta_k^2
+2\sum_{k<l}\sigma_{jk}\sigma_{jl}\beta_k\beta_l. \tag{AM7}
\]

Let \(E[\beta\mid\Sigma]=0\) and
\(\Omega_{\beta\mid\Sigma}=E[\beta\beta^T\mid\Sigma]\). On the same standardized
phenotype and association-statistic scale as (M6),

\[
\begin{aligned}
E[\chi_j^2\mid\Sigma]
&\approx1+N\mathbf e_j^T\Sigma\Omega_{\beta\mid\Sigma}\Sigma\mathbf e_j,\\
\Omega_{\beta\mid\Sigma}
&=\operatorname{diag}(\tau_1^2,\ldots,\tau_M^2)
\ \Longrightarrow\
1+N\sum_k\tau_k^2\sigma_{jk}^2,\\
\Omega_{\beta\mid\Sigma}
&=\tau^2I
\ \Longrightarrow\
1+N\tau^2\ell_j .
\end{aligned} \tag{AM8}
\]

A diagonal effect ensemble sets the expected cross terms in (AM7) to zero;
isotropy is the additional condition that yields an ordinary one-slope LDSC
line. Because assortment makes \(\Sigma\) depend on the generating trait,
conditional isotropy is not justified as a causal averaging model for that
trait, even though imposing it recovers the algebra in (M6). Equations
(M12)--(M14) describe the remaining projection problem. A local-window LD score
also excludes all cross-autosome and distant same-chromosome terms.

The genomic signature is established prior art: Yengo et al. (2018) detected
odd/even-chromosome GPD, Border et al. (2022a) derived (AM4), and Zhang et al.
(2025) introduced DGREML for directional-GPD variance. GPD alone does not prove
mate choice; geography/ancestry, some forms of selection, social homogamy, and
selection into the analyzed sample can also generate long-range covariance.
Pure vertical environmental transmission instead creates genotype--environment
covariance, which is the distinct problem below.

### 15.2 Pure parental nurture can preserve the LDSC line

Let \(G_o,G_m,G_f\) denote componentwise-standardized offspring, maternal, and
paternal genotype vectors on one common, stable allele-frequency scale, and write
\(\Sigma_o=\operatorname{Cov}(G_o)\) and
\(\Sigma_p=\operatorname{Cov}(G_m)=\operatorname{Cov}(G_f)\). Under Mendelian
transmission and random mating,

\[
E[G_o\mid G_m,G_f]=\frac{G_m+G_f}{2},\qquad
\operatorname{Cov}(G_o,G_m)
=\operatorname{Cov}(G_o,G_f)=\frac12\Sigma_p. \tag{GN1}
\]

Consider a phenotype caused entirely by an additive parental environmental
pathway:

\[
Y=G_m^T\gamma_m+G_f^T\gamma_f+\varepsilon,\qquad
\beta_{\rm direct}=0,\qquad
\operatorname{Cov}(G_o,\varepsilon)=0. \tag{GN2}
\]

The ordinary offspring-GWAS marginal association vector is

\[
\eta:=\operatorname{Cov}(G_o,Y)
=\Sigma_p\theta,\qquad
\theta=\frac{\gamma_m+\gamma_f}{2}. \tag{GN3}
\]

If LD is stable across generations and the reference is matched,
\(\Sigma_p=\Sigma_o=\Sigma\), this is identical at the marginal-moment level to
a direct model with effect \(\theta\). Suppose additionally that the phenotype
has unit variance, association statistics are correctly calibrated,
\(E[\theta\theta^T\mid\Sigma]=\tau_\theta^2I\), and the orthogonal parental
component is homoscedastic among unrelated offspring. Then

\[
E[\chi_j^2\mid\Sigma]
\approx1+N\tau_\theta^2\ell_j
=1+\frac{Nh_{\rm assoc}^2}{M}\ell_j,\qquad
h_{\rm assoc}^2=M\tau_\theta^2. \tag{GN4}
\]

Thus the textbook LDSC line and intercept can hold while direct SNP
heritability is zero. Its slope is offspring-genotype-tagged association
variance, not direct heritability or total nurture variance.

For equal maternal and paternal nurture effects
\(\gamma_m=\gamma_f=\gamma\), independent random mates, and no other pathway,

\[
V_{\rm nurture}=2\gamma^T\Sigma\gamma,\qquad
V_{\rm tagged}=\gamma^T\Sigma\gamma
=\frac12V_{\rm nurture}. \tag{GN5}
\]

Thus a unit-variance phenotype consisting entirely of symmetric biparental
nurture has \(h_{\rm assoc}^2=1/2\) in this toy model despite zero direct
heritability. With one nurturing parent, the tagged fraction is \(1/4\).
These are thought-experiment scalings, not empirical expectations.

With direct effects \(d\) plus symmetric parental nurture, offspring GWAS
identifies only their sum and its tagged variance:

\[
\theta_{\rm assoc}=d+\gamma,\qquad
V_{\rm assoc}
=(d+\gamma)^T\Sigma(d+\gamma)
=d^T\Sigma d+\gamma^T\Sigma\gamma+2d^T\Sigma\gamma. \tag{GN6}
\]

The variance contributed by the direct and parent-genotype-mediated pathways is
different: under random mating, the coefficient on
\(\gamma^T\Sigma\gamma\) is two, not one. Neither the LDSC slope nor intercept
separates the terms in (GN6). A nontransmitted parental allele supplies the
needed external nurture anchor (Kong et al. 2018); its reported effect ratio is
not a fraction of heritability.

For a shared parental pathway, correctly centered within-sibship association
removes (GN2), provided the centered residual is genotype-uncorrelated:

\[
\begin{aligned}
\operatorname{Cov}
\{G_{o,i}-\bar G_{o,\rm family},\,
  Y_i-\bar Y_{\rm family}\}&=0,\\
\operatorname{Cov}
\{G_{o,i}-\bar G_{o,\rm family},\,
  \varepsilon_i-\bar\varepsilon_{\rm family}\}&=0 .
\end{aligned} \tag{GN7}
\]

This does not cover sibling effects, child-genotype-evoked or differential
parenting, or other non-shared indirect effects. Population-minus-within-family
attenuation is therefore not itself a nurture estimator. Lee et al. (2018)
analyzed LDSC under parental environmental effects; family and
multigenerational methods such as DONUTS and JODIE already address direct versus
indirect decomposition.

### 15.3 Assortment plus nurture requires cross-generational geometry

On one common fixed per-locus scale, let
\(\Xi_o=\operatorname{Cov}(G_o)\),
\(\Xi_p=\operatorname{Cov}(G_m)=\operatorname{Cov}(G_f)\),
\(\Xi_{mf}=\operatorname{Cov}(G_m,G_f)\), and
\(\Xi_{oP}=\operatorname{Cov}(G_o,G_m+G_f)\). Mendelian transmission gives

\[
\Xi_{oP}
=\Xi_p+\frac12(\Xi_{mf}+\Xi_{mf}^T). \tag{GN8}
\]

For equal parental nurture effects \(\gamma\), \(\eta=\Xi_{oP}\gamma\).
Under isotropic nurture effects and the calibration assumptions of (GN4), the
matched offspring--parent row score is

\[
\ell_j^{\rm OP}
:=\frac{(\Xi_{oP}\Xi_{oP}^T)_{jj}}{(\Xi_o)_{jj}}
=\frac{\sum_k(\Xi_{oP})_{jk}^2}{\operatorname{Var}(G_{o,j})},
\qquad
E[\chi_j^2\mid\Xi_{oP},\Xi_o]
\approx1+N\tau_\gamma^2\ell_j^{\rm OP}. \tag{GN9}
\]

Random mating, stable LD, and componentwise standardization reduce
\(\Xi_{oP}\) to offspring LD and (GN9) to
ordinary LDSC. Assortment adds the long-range component in (AM2), so the
cross-generational operator need not equal offspring LD. Distinct maternal and
paternal architectures require separate offspring--mother and
offspring--father operators and effect covariances. Offspring genotypes do not
identify \(\Xi_{oP}\) without mating and transmission assumptions; a
directly matched score requires parent--offspring reference data.

For compactness, Table 3 writes \(\Sigma_{\rm AM}\) for the analyzed offspring
LD matrix under assortment.

**Table 3. What offspring GWAS plus offspring LD can identify**

| Generative case | Marginal association object | Ordinary scalar LDSC? | Causal mechanism identified? |
|---|---|---|---|
| Direct, isotropic effects | \(\Sigma d\) | Yes, under the usual assumptions | Only with the direct-effect model imposed |
| Pure parental nurture, random mating | \(\Sigma\theta\) from (GN3) | Yes, under the same effect assumptions | No; aliased with a direct coefficient at the marginal-moment level |
| Trait-aligned assortative mating | \(\Sigma_{\rm AM}\beta\), with LD aligned to \(b\) | Only if conditional isotropy and a full matched LD score are imposed; otherwise M14 is a projection | No; signed GPD does not identify its cause |
| Assortment plus nurture | \(\Xi_{oP}\gamma\) | Only under isotropic nurture effects with a matched score such as (GN9) | No; requires family and mate anchors |

Changing ordinary LD to \(P_c\)-weighted LD changes association geometry; it
does not manufacture the family or historical information needed to distinguish
direct biology, mate choice, and parental environment. RDR or within-family
designs can estimate a direct component under their assumptions;
transmitted/nontransmitted parental or trio data identify parental nurture; and
spouse data provide an assortment anchor. Any summary-statistic directional-GPD
or cross-generational extension would be a new branch beside PLDSC and should be
benchmarked against DGREML, DONUTS/JODIE, and family GWAS.

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
8. Isserlis (1918), [On a formula for the product-moment coefficient of any
   order of a normal frequency distribution in any number of variables](https://doi.org/10.1093/biomet/12.1-2.134).
9. Künsch (1989), [The jackknife and the bootstrap for general stationary
   observations](https://doi.org/10.1214/aos/1176347265).
10. Xiong et al. (2024), [Improved estimation of functional enrichment in SNP
    heritability using feasible generalized least squares](https://doi.org/10.1016/j.xhgg.2024.100272).
11. Yengo et al. (2018), [Imprint of assortative mating on the human
    genome](https://doi.org/10.1038/s41562-018-0476-3).
12. Border et al. (2022a), [Assortative mating biases marker-based heritability
    estimators](https://doi.org/10.1038/s41467-022-28294-9).
13. Border et al. (2022b), [Cross-trait assortative mating is widespread and
    inflates genetic correlation estimates](https://doi.org/10.1126/science.abo2059).
14. Zhang et al. (2025), [The contribution of gametic phase disequilibrium to
    the heritability of complex traits](https://doi.org/10.1038/s41588-025-02192-4).
15. Lee et al. (2018), [The accuracy of LD Score regression as an estimator of
    confounding and genetic correlations in genome-wide association
    studies](https://doi.org/10.1002/gepi.22161).
16. Howe et al. (2022), [Within-sibship genome-wide association analyses
    decrease bias in estimates of direct genetic
    effects](https://doi.org/10.1038/s41588-022-01062-7).
17. Wu et al. (2021), [Estimating genetic nurture with summary statistics of
    multigenerational genome-wide association
    studies](https://doi.org/10.1073/pnas.2023184118) (DONUTS).
18. Krätschmer et al. (2026), [Separating direct, indirect, and
    parent-of-origin genetic effects in the human
    population](https://doi.org/10.1016/j.xgen.2026.101277) (JODIE).
