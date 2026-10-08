# LD scores for multi-ancestry GWAS

Date: 2026-09-25. The derivations behind `pldsc/multiancestry.py`, checked
three ways:

- exact in-sample identities (`tests/test_multiancestry.py`);
- a Monte Carlo with out-of-sample reference panels
  ([multiancestry_meta_sim.py](../experiments/multiancestry_meta_sim.py));
- a first measurement on real 1000 Genomes panels
  ([multiancestry_1000g_meta.py](../experiments/multiancestry_1000g_meta.py)).

Equation labels MA1–MA12 are the ones the code cites.

## 1. Setting

There are ancestries a = 1, …, A. In ancestry a, variant j has genotype
variance σ²_aj and population LD correlations r^a_jk with every other variant
k. Every panel counts the same allele, so the sign of r^a_jk is comparable
across ancestries.

**Assumptions.**

- **(A1) Infinitesimal effects, standardized within each ancestry.** The
  effects β^a_k = σ_ak b^a_k have mean 0, where b^a_k is the per-allele effect.
  Their covariance is Cov(β^a_k, β^b_l) = δ_kl τ_ab, with τ_aa = h_a²/M and
  τ_ab = ρ_ab h_a h_b/M. Here ρ_ab is the standardized-scale trans-ancestry
  genetic correlation, the quantity S-LDXR estimates. A per-allele-scale model
  would make τ_ab depend on the SNP. The formulas below carry through with
  τ_ab,j inside the sums.
- **(A2) Independent samples.** Each ancestry has its own GWAS on N_a
  unrelated individuals, with no sample overlap between ancestries.
- **(A3) Structure is controlled within each GWAS.** Phenotypes have unit
  variance in each ancestry.

The **cross-ancestry LD score** is

  (MA1) ℓ^ab_j = Σ_k r^a_jk r^b_jk,

so ℓ^aa_j is the ordinary LD score. The self term contributes 1 to every
ℓ^ab_j.

**Lemma 1** (per-ancestry moments; Bulik-Sullivan 2015 on the diagonal, and
Popcorn / S-LDXR off it). Write z_aj = √N_a · r̂_a(g_j, y). Then

  E[z_aj z_bj] = δ_ab + √(N_a N_b) τ_ab ℓ^ab_j.

*Proof.* z_aj = √N_a (Σ_k r^a_jk β^a_k + ε_aj), where the sampling errors ε
have variance 1/N_a and are independent across ancestries by A2. Take the
product and use A1. ∎

## 2. Fixed-effect meta-analysis

A meta-analysis forms z_j = Σ_a c_aj z_aj with weights satisfying
Σ_a c_aj² = 1. Two weight schemes are in common use:

- **Sample-size (Stouffer; METAL SAMPLESIZE):** c_aj² = N_a/N_j, where
  N_j = Σ_a N_a over the ancestries that report j.
- **Inverse-variance on the per-allele scale:** β̂_aj = z_aj/(√N_a σ_aj) has
  variance 1/(N_a σ²_aj). Hence c_aj² = N_a σ²_aj / Σ_b N_b σ²_bj.

**Theorem 1** (the meta-analysis LD score). Under A1–A3, with
η_a = h_a²/h² the heritability of ancestry a relative to a reference value h²,

  (MA2) E[χ²_j] = 1 + (N_j h²/M) · L_j,  L_j = Σ_a Σ_b ω_ab,j ρ_ab √(η_a η_b) ℓ^ab_j,

  (MA3) ω_ab,j = c_aj c_bj √(N_a N_b) / N_j.

The two weight schemes give:

- sample-size weights: ω_ab = f_a f_b, where f_a = N_a/N_j;
- inverse-variance weights: ω_ab,j = f_a f_b σ_aj σ_bj / σ̄²_j, where
  σ̄²_j = Σ_c f_c σ²_cj.

*Proof.* Expand E[z_j²] = Σ_ab c_aj c_bj E[z_aj z_bj] and apply Lemma 1 and
Σ c² = 1. ∎

Two bounds follow.

- **(MA4)** Σ_ab ω_ab,j equals 1 under sample-size weights and is at most 1
  under inverse-variance weights. For the latter this is Cauchy–Schwarz:
  (Σ f σ)² ≤ Σ f σ².
- **(MA5)** Since |ℓ^ab| ≤ √(ℓ^aa ℓ^bb) (Cauchy–Schwarz over k) and
  |ρ| ≤ 1, with η = 1,
  L_j ≤ (Σ_a f_a √ℓ^aa_j)² ≤ Σ_a f_a ℓ^aa_j under either scheme. **The
  correct score never exceeds the composition-weighted mean of the
  within-ancestry scores.** Equality requires identical LD rows in every
  ancestry and ρ = 1.

## 3. What a single-ancestry score costs

LDSC regresses χ²_j on (N_j/N̄) ℓ^R_j for a reference score ℓ^R, e.g.
1000G EUR. By MA2 the (weighted) population regression gives

  (MA6) slope = (N̄ h²/M) κ, κ = Cov_w(x_L, x_R)/Var_w(x_R);
      intercept = 1 + (N̄ h²/M) · (Ē_w x_L − κ Ē_w x_R),

where x_L = (N_j/N̄) L_j and x_R = (N_j/N̄) ℓ^R_j. So ĥ²/h² → κ, which is
computable from the panels, the reported ancestry mix, and an assumed ρ. No
GWAS is needed. For a single-ancestry GWAS with R = a, κ = 1 exactly. This is
`single_reference_bias(L, ell_R, sample_sizes=N_j)`.

## 4. Estimation from reference panels

The estimators are

  (MA7) ℓ̂^aa_j = Σ_{k∈W_j} [r̂²_ajk − (1 − r̂²_ajk)/ν_a],  ℓ̂^ab_j = Σ_{k∈W_j} r̂^a_jk r̂^b_jk  (a ≠ b),

with ν_a = n_a − 2 for unrelated panels, or `effective_pair_df` for related
ones (METHOD §10.1).

**(MA8) The noise floor is asymmetric.** A squared within-panel correlation
has E[r̂²] = r² + (1 − r²)/n + O(n⁻²), which is why the diagonal needs a
correction. The cross terms multiply correlations estimated in *independent*
panels, so E[r̂^a r̂^b] = E[r̂^a] E[r̂^b] = r^a r^b (1 + O(1/n)). There is
nothing to subtract. Both obvious mistakes are material at 1000G sizes:
subtracting a floor from the cross terms, and forgetting it on the diagonal.
Each misstates L by about ω·|W|/n, i.e. roughly 2 LD-score units per term for
~1000 window SNPs and n ≈ 500.

**(MA9) Within and cross parts.** Split MA2 at ρ:
L = D + ρ X for a common ρ, where D uses ρ = I and X is the cross part. In
code, `meta_ld_scores(..., rho=np.eye(A))` returns D and the difference from
`rho=ones` returns X. In principle a two-regressor LDSC of the meta χ² on
(D, X) identifies h² and h²ρ. §6 Table 5 shows that on real panels D and X
are too collinear for this to be useful.

## 5. Pooled (mega-) analysis

A pooled analysis regresses on all N = Σ N_a individuals with covariates C.
Let x̃_j be g_j residualized on C and scaled to unit mean square over the
pooled sample. Then, exactly given genotypes,

  (MA10) m^a_jk = (1/(N_a s_ak)) Σ_{i∈a} x̃_ij g_ik,

  (MA11) E[χ²_j | G] = (noise) + N (h²/M) L^pool_j,  L^pool_j = Σ_ab f_a f_b ρ_ab √(η_a η_b) Σ_k m^a_jk m^b_jk.

**Proposition 2** (MA12: a pooled analysis ≡ an inverse-variance meta for
discrete groups). If C contains the ancestry-group indicators, or PCs whose
span includes them, then x̃_ij = (g_ij − ḡ_aj)/σ̃_j for i ∈ a, and
m^a_jk = s_aj r̂^a_jk/σ̃_j. Therefore L^pool_j is exactly the
inverse-variance L_j of MA2. This is the LD-score form of the known
equivalence between fixed-effect meta- and mega-analysis (Lin & Zeng 2010).

A pooled analysis differs from a meta only through what the covariates miss:

- group mean differences left unadjusted, which give Σ_{i∈a} x̃_ij ≠ 0. This
  is the structure channel PLDSC measured, 29% of the windowed score at
  ±1 cM in the 1000G FULL panel;
- admixed individuals, who do not fall into discrete groups.

For those cases, `pooled_ld_scores` evaluates MA11 on a composition-matched
pooled panel. It applies no finite-sample correction.

## 6. Results

**Table 1. Monte Carlo check** ([multiancestry_meta_sim.py](../experiments/multiancestry_meta_sim.py)).

Setup: two simulated ancestries with different LD lengths; GWAS
N = 3000 + 1500; M = 2000; h² = 0.5; 3 × 100 replicates. Entries are mean
ĥ²/h², and "ratio" is the slope on the used score divided by the slope on the
correct score, which is what MA6 predicts.

| E share | ρ | weights | correct score | E-only: ratio vs predicted κ | ρ=1 score: ratio vs predicted |
|---|---|---|---|---|---|
| 1.00 | 1.0 | sample size | 1.041 | 1.000 vs 1.000 | 1.000 vs 1.000 |
| 0.67 | 1.0 | sample size | 0.966 | 0.608 vs 0.610 | 1.000 vs 1.000 |
| 0.67 | 0.6 | sample size | 1.047 | 0.547 vs 0.548 | 0.890 vs 0.889 |
| 0.67 | 1.0 | IVW | 0.994 | 0.738 vs 0.742 | 1.000 vs 1.000 |
| 0.67 | 0.6 | IVW | 0.962 | 0.679 vs 0.681 | 0.905 vs 0.906 |

These rows are the **post-failure rerun** with 6,000-person reference panels.
The rule-bearing default run used 400-person panels and **FAILED** its
absolute-calibration checks. Even its single-ancestry arm returned 0.80. The
cause is errors-in-variables from small panels on this toy's very narrow
LD-score spread (SD 0.32 around a mean of 1.9). On one realization
(`--baseline-check`), single-ancestry ĥ²/h² was 0.85 / 0.97 / 1.02 at panel
sizes 400 / 1,500 / 6,000, and 1.04 with in-sample LD. The ratio check (iii) matched within 0.015
in that run too, but it was observed afterwards, not pre-stated. The rerun's
rules were written before it ran.

**Table 2. Mean LD scores on 1000G chr22** (every 10th biallelic SNP, ±1 cM,
9,654 regression SNPs with full-panel MAF ≥ 0.05). The diagonal is corrected
by MA7; off-diagonal entries are cross scores.

| | EUR | AFR | EAS | SAS | AMR |
|---|---|---|---|---|---|
| EUR | 9.76 | 5.08 | 6.68 | 7.92 | 8.65 |
| AFR | | 7.32 | 4.32 | 4.70 | 5.95 |
| EAS | | | 8.64 | 7.04 | 6.83 |
| SAS | | | | 8.72 | 7.66 |
| AMR | | | | | 12.33 |

AMR's large within-ancestry score is its admixture LD. A 1-in-3 denser
cache scales every entry by about 3 (EUR 30.4) and leaves κ essentially
unchanged (Table 3).

**Table 3. κ = ĥ²/h² for a meta regressed on the 1000G EUR score.** These are
the corrected numbers (Table 4): N_j-scaled regressors, regression SNPs with
EUR MAF ≥ 0.05, and sample-size weights. Each cell gives OLS / LDSC-style
weights, 1-in-10 cache.

| mix | ρ = 1.0 | ρ = 0.8 | ρ = 0.6 |
|---|---|---|---|
| EUR 0.90 + AFR | 0.908 / 0.906 | 0.890 / 0.889 | 0.872 / 0.873 |
| EUR 0.75 + AFR | 0.788 / 0.786 | **0.750 / 0.750** | 0.713 / 0.714 |
| EUR 0.50 + AFR | 0.643 / 0.644 | 0.592 / 0.595 | 0.541 / 0.546 |
| EUR 0.90 + EAS | 0.955 / 0.955 | 0.930 / 0.930 | 0.904 / 0.906 |
| EUR 0.75 + EAS | 0.896 / 0.898 | **0.841 / 0.845** | 0.787 / 0.792 |
| EUR 0.50 + EAS | 0.828 / 0.840 | 0.753 / 0.766 | 0.677 / 0.693 |
| EUR 0.90 + SAS | 0.971 / 0.970 | 0.940 / 0.940 | 0.910 / 0.910 |
| EUR 0.75 + SAS | 0.932 / 0.933 | **0.868 / 0.870** | 0.804 / 0.807 |
| EUR 0.50 + SAS | 0.879 / 0.889 | 0.794 / 0.805 | 0.709 / 0.720 |
| EUR 0.90 + AMR | 0.987 / 0.981 | 0.954 / 0.949 | 0.920 / 0.917 |
| EUR 0.75 + AMR | 0.977 / 0.963 | **0.907 / 0.896** | 0.838 / 0.829 |
| EUR 0.50 + AMR | 0.985 / 0.957 | 0.892 / 0.867 | 0.800 / 0.778 |

The 1-in-3 cache (32,446 regression SNPs) gives, at the bold cells:
AFR 0.743 / 0.748, EAS 0.830 / 0.845, SAS 0.870 / 0.870, AMR 0.909 / 0.898.
The corrected intercept shift is small: +0.00 to +0.03 of the LD-driven
excess χ² at EUR 0.75, and at most +0.08 at EUR 0.50 + AFR.

**Reading Table 3.**

- κ < 1 even at ρ = 1. The dominant cost is LD mismatch, not effect
  heterogeneity: a EUR panel overstates the tagging a mixed sample actually
  experiences, which is MA5 in numbers.
- Lower ρ adds a further, roughly linear, loss.
- AFR matters most, because its LD is shortest. AMR matters least, because
  its LD is EUR-like plus admixture LD.

**Table 4. The kill rule** ([multiancestry_1000g_meta.py](../experiments/multiancestry_1000g_meta.py);
stated before running). The rule cell is EUR 0.75 of N + one other
superpopulation, ρ = 0.8, sample-size weights. KEEP if |κ − 1| ≥ 0.10 for any
X.

| computation | AFR | EAS | SAS | AMR | verdict |
|---|---|---|---|---|---|
| as pre-stated (constant N, full-panel regression set) | 0.716 | 0.811 | 0.821 | 0.850 | **KEEP** |
| corrected, OLS | 0.750 | 0.841 | 0.868 | 0.907 | KEEP |
| corrected, LDSC-style weights | 0.750 | 0.845 | 0.870 | 0.896 | KEEP |

The correction came from a robustness check run after the rule. Under
LDSC-style weights the pre-stated computation gave κ ≈ 0.1–0.4, which was
implausible. The cause was two defects in how I first computed κ, not in MA6:

- it ignored per-SNP N_j, although LDSC scales each LD score by N_j;
- it kept regression SNPs that are rare or absent in EUR. No EUR LD-score file
  contains them, and their EUR score of ~0 made them dominate a weighted fit.

The verdict is the same either way. The pre-stated computation overstated
both the slope loss and the intercept inflation (+0.08 to +0.11 there,
+0.00 to +0.03 corrected).

**Table 5. Is ρ identifiable from the meta alone?** This is the correlation of
the N_j-scaled within part D and cross part X across SNPs (MA9), at EUR 0.75.

| X | AFR | EAS | SAS | AMR |
|---|---|---|---|---|
| corr(D, X), 1-in-10 | 0.847 | 0.813 | 0.980 | 0.988 |
| corr(D, X), 1-in-3 | 0.871 | 0.798 | 0.983 | 0.988 |

For SAS and AMR the two parts are collinear, so ρ is not identified. For AFR
and EAS it is weakly identified at best (variance inflation of about 3–4).
**A meta-analysis alone cannot tell LD mismatch from effect heterogeneity.**
ρ has to come from ancestry-specific summary statistics (Popcorn, S-LDXR) or
be bracketed, as in Table 3.

**Table 6. Genome-wide κ at the rule cell** (EUR 0.75 + X, ρ = 0.8;
[multiancestry_genome_kappa.py](../experiments/multiancestry_genome_kappa.py),
2026-10-08). The data are bigsnpr's 1000G phase-3 set (2,490 samples,
1.66M variants), with 1,190,723 regression SNPs on 22 autosomes and the MHC
excluded. Entries are OLS / LDSC-style weights. The SE is a 200-block
jackknife over the genome.

| X | chr22 (Table 3) | genome-wide | SE (OLS) | per chromosome (OLS) | MHC kept | IVW weights |
|---|---|---|---|---|---|---|
| AFR | 0.750 / 0.750 | **0.737 / 0.726** | 0.005 | 0.709–0.819 | 0.728 / 0.726 | 0.746 / 0.735 |
| EAS | 0.841 / 0.845 | **0.846 / 0.837** | 0.008 | 0.753–0.936 | 0.820 / 0.836 | 0.849 / 0.841 |
| SAS | 0.868 / 0.870 | **0.877 / 0.866** | 0.005 | 0.839–0.938 | 0.848 / 0.865 | 0.876 / 0.866 |
| AMR | 0.907 / 0.896 | **0.896 / 0.891** | 0.004 | 0.849–0.939 | 0.875 / 0.890 | 0.893 / 0.890 |

- **The step-1 gate passes** (AFR κ ≤ 0.90). Genome-wide κ is within 0.025
  of chr22 in every cell.
- **The bridge to the step-1 protocol missed its tolerance.** That protocol
  is a 1-in-10 thinning of the full VCFs, and only chr10/21/22 are local. On
  those three chromosomes the same code on the two variant sets differs by
  up to 0.039 (chr10 AFR), against a pre-stated 0.03. So this table is
  specific to its variant set, and the headline stays chr22's until the
  thinned-VCF run.
- Restricting regression SNPs to HapMap3 moves κ by ≤ 0.002.
  η_X = 0.8 / 1.25 gives AFR 0.716 / 0.761 (OLS). Genome-wide corr(D, X) is
  0.87 / 0.85 / 0.98 / 0.98, so Table 5's conclusion holds.
- The per-chromosome outliers (chr12 EAS 0.75, chr14 0.82–0.94) come from a
  few extreme-LD regions that dominate OLS. On chr14, 0.1% of SNPs carry 10%
  of the EUR-score variance.

**Table 7. The GBMI asthma casualty**
([gbmi_asthma_casualty.py](../experiments/gbmi_asthma_casualty.py),
2026-10-08). Data: GBMI release 052021 asthma, i.e. the all-biobank meta
plus the five ancestry metas, joined to Table 6's scores by rsID. That gives
953,596 regression SNPs. N is the effective N implied by each file's SEs,
and at GBMI's K = 0.085 the liability conversion factor is 1.00. Rules were
stated before the files were read. Composition by N: EUR 0.80, EAS 0.14,
AFR, SAS and AMR 0.03 or less.

| all-biobank meta regressed on | h² (SE) |
|---|---|
| the 1000G EUR score (GBMI's method) | 0.133 (0.012) |
| the meta score MA2, ρ̂ EUR–EAS, other pairs ρ = 0.8 | **0.156 (0.014)** |
| … other pairs ρ = 0.6 / 1.0 | 0.159 / 0.153 |
| post hoc: N-weighted mean of the per-stratum h², each on its own score | 0.154 (0.012) |
| GBMI reported (UKB LD scores, liability) | 0.085 |

- **ρ is estimable only for EUR–EAS: 0.755 (SE 0.046).** The other strata
  are too small (h² z = 2.1 / 3.4 / 0.4 for AFR / SAS / AMR), so their pairs
  are bracketed. Their point estimates (EUR–AFR 0.73, EUR–SAS 0.83) lie
  inside the bracket.
- **The step-3 gate is not met robustly.** h² moves by 1.9 SE of the
  EUR-score estimate at the central bracket, and by 1.6–2.2 SE across it. By
  STRATEGY §7, asthma alone makes a practice note. The shift itself is
  precise: +0.023 with a paired jackknife SE of 0.002. The gate's SE is
  mostly sampling noise that both estimates share.
- **Does the data agree with MA2?** The observed ratio of the two estimates
  is 0.852 against a predicted κ of 0.850. That agreement is nearly
  automatic, because the two regressors correlate 0.9965. The real check is
  the post-hoc row: per-stratum estimates on matched scores average 0.154,
  which agrees with the meta score and not with the EUR score.
- **0.133 is not GBMI's 0.085.** The LD scores differ (UKB, denser) and
  GBMI's N convention is not stated. If the ratio carries over, GBMI's value
  would be about 0.10 on the meta score. That is an extrapolation, not a
  measurement.

## 7. Scope and caveats

- **Variant sets.** Tables 2–5 are chr22, a 1-in-10 thinning of the phase-3
  VCF. Table 6 is genome-wide on bigsnpr's set. 1000G superpopulations stand
  in for ancestries (n = 347–661): AMR is admixed and SAS is heterogeneous.
- **Single map.** One population-averaged genetic map (SHAPEIT4 GRCh37) is
  used for every ancestry. The phase-3 chr22 file is the S3 mirror's `v5a`
  release.
- **The effect model A1 is one choice.** A MAF-dependent architecture
  (α ≠ −1) or per-ancestry heritabilities (η ≠ 1) change κ. Both enter MA2
  explicitly (`h2_ratio`, and τ_ab,j), but they are not explored here.
- **A random 1-in-10 thinning is not HapMap3.** The density check (1-in-3)
  moves κ by ≤ 0.015.
- **The LDSC-style weights are a simplified version** (overcounting ×
  heteroscedasticity at mean χ² 1.5), not ldsc's two-step fit.
- **κ is a population-regression quantity.** Finite-GWAS noise adds variance,
  not bias.

## 8. Prior art (not yet swept)

- **Occupied.** Lemma 1 is Popcorn's (Brown et al. 2016) and S-LDXR's (Shi et
  al. 2021) cross-population moment. Proposition 2 is the LD-score analogue of
  Lin & Zeng 2010. MAMA (Turley et al. 2021) uses cross-ancestry LD scores to
  combine ancestries.
- **Plausibly new, pending a sweep** (STRATEGY §4):
  - the explicit meta-analysis LD score MA2, with its two weight schemes;
  - the panel-only bias factor κ (MA6) and its bound MA5;
  - the noise-floor asymmetry MA8, stated as a warning;
  - the measured κ on 1000G.

  The one recorded empirical warning that LDSC underestimates h² on
  meta-analysed summary statistics — a 2026 bioRxiv preprint, now citable
  (doi 10.64898/2026.07.05.736573, PMID 42465239; reports the phenomenon
  and one magnitude, AD 0.265 → 0.160 with imputed variants, but not MA2
  or a panel-only κ) — points the same way as Table 3. Momin et al. 2023
  (*Nat Commun* 14:1152) and cov-LDSC (Yao et al. 2019) are the other two
  the step-2 sweep must position against (STRATEGY §5).
