# Prior-art positioning

Date: 2026-07-18

**Status: Bercovich 2025 now read in full** (Oxford HTML + PDF, cross-checked
term census); the novelty-decisive gate flagged by the external review is closed.
The *broader* landscape verdict is still from ~6 targeted web searches plus the
applications-panel findings — the six-cluster subagent workflow
(`workflows/scripts/pldsc-prior-art-*.js`) hit a session usage limit and should be
re-run for a fully source-checked sweep. Treat the Bercovich-specific conclusions
below as verified and the rest as provisional.

## Headline

The reframed project (local/non-local LD-score decomposition; non-local part
approximated by the PC score; unification of the structure-corrections; structure
diagnostics) is a **framework/synthesis + theory-transfer + diagnostics**
contribution, **not a new estimator**. The core idea "PC-residualization makes an
LD quantity structure-invariant" has a clear published lineage that we must credit
and not claim:

- **Mangin et al. 2012** — originator: regress genotypes on inferred ancestry
  proportions, correlate the residuals → structure-corrected LD.
- **Bercovich et al. 2025** — the modern, proved version for *pairwise r²* used in
  **pruning / clumping / FST / PCA QC** (Albrechtsen group).
- **Luo et al. 2021 (cov-LDSC)** — the same residualization applied to *LD scores*
  for **heritability** — i.e. the direct competitor for our application.

What is left after removing those is genuinely narrower than the first pass
suggested, but real: (1) the **operator** — an eigenvalue-weighted GLS/mixed-model
residual-precision `P_c`, not a hard rank-K projection (see below); (2) the
**exact LD-score spectral decomposition** `Σ_k r_jk² = (M/N)Σ_p λ_p corr(x_j,u_p)²`
used as the organizing identity, with the 2017 PC score identified as its
leading-eigenvector *truncation* (top-10 PCs: corr ≈ 0.99 but only ~37% of the
non-local magnitude — the split is by construction, not an empirical discovery);
(3) the **unification** of cov-LDSC / WhiLD / PC-score / `P_c` under that identity;
and (4) **LDSC-specific diagnostics + the identifiability limit** + efficiency.

## Bercovich et al. 2025 — full read (verified)

DOI [10.1093/genetics/iyaf009](https://doi.org/10.1093/genetics/iyaf009), *Genetics*
229(3), iyaf009 (bioRxiv 2024.05.02.592187). Authors: Bercovich, Rasmussen, Li,
Wiuf, Albrechtsen (Copenhagen).

- **Operator.** Residual matrix `R = G − 2Π̂` with `Π̂ = Q F̂` (predicted individual
  allele frequencies). For known/estimated ancestry design `Q`, the adjusted
  covariance is `D̂_adj = 1/(2(n−k))·Gᵀ(I−P)G` with `P = Q(QᵀQ)⁻¹Qᵀ`. This is a
  **hard rank-k orthogonal projection** off the admixture / top-K-PC subspace
  (ordinary OLS residualization). **No eigenvalue weighting.** Adjusted r² is the
  squared correlation of the residuals.
- **Theory.** *Theorem 1:* `(D_adj)_st = Σ_ℓ w_st^ℓ (D^dℓ)_st` — the adjusted
  covariance is a weighted average of *within-ancestry* LD, hence invariant to
  structure/admixture. *Theorem 3:* the estimator is consistent for that
  population quantity. Assumptions: admixture `Π = QF`, `G_is ~ Bin(2, Π_is)`,
  independence across individuals and between the two haplotypes within an
  individual (HWE within ancestry).
- **Scope — what it does NOT touch.** Term census, two independent fetches in
  agreement: **heritability 0, LD score 0, LDSC 0, genetic correlation 0, variance
  component 0, polygenic 0, partitioning 0.** One passing "mixed models" reference
  (Zhou & Stephens) and four "GWAS" mentions, all as clumping motivation. The
  paper is a **QC-for-population-genetics** method (unbiased pruning/clumping,
  unbiased FST and PCA after pruning). It has no heritability-estimation content.
- **Their own prior art:** they credit **Mangin et al. 2012** for the
  residual-correlation LD idea and differentiate by avoiding explicit confounder
  estimation and linearity assumptions.

**Precise delta vs PLDSC.** Two independent axes separate us from Bercovich:
1. *Different operator.* Bercovich = **hard projection** `I − P` off the top-K
   subspace (idempotent, unweighted). PLDSC `P_c = V_c⁻¹ − V_c⁻¹C(CᵀV_c⁻¹C)⁺CᵀV_c⁻¹`
   with `V_c = σ_g²K_{-c} + σ_e²I` is an **eigenvalue-weighted GLS whitening** over
   the *whole* GRM spectrum (down-weights direction `p` by `1/(σ_g²λ_p+σ_e²)`), of
   which the hard projection is the `σ_g²→∞`-on-top-eigenvalues limit. cov-LDSC
   also uses the hard-projection form. So the GLS-weighted LD-score operator is
   genuinely un-published — though we have shown that in the *windowed* regime it
   collapses to ≈cov-LDSC, so the practical gap is small.
2. *Different target.* Bercovich corrects **pairwise r²** for pruning/clumping;
   PLDSC uses row-sums (**LD scores**) as the regressor in **heritability LD-score
   regression**. Non-overlapping downstream objects.

**Bottom line:** Bercovich does *not* collide with any heritability-facing claim
(C4/C5/C6) and only partially overlaps C1/C2 (same "PC-residualization ⇒
structure-invariance" intuition, but a different operator and a different target).
It *does* pre-empt any claim that "structure-corrected LD from PC residuals" is
new — that belongs to Mangin/Bercovich (r²) and cov-LDSC (LD scores).

## Per-claim assessment

| Claim | Verdict | Nearest prior work |
|---|---|---|
| C1 LD = local (linkage) + non-local (structure) | **published (concept)** | Bercovich 2025 (pairwise r², hard projection); cov-LDSC (LD-score, "long-range admixture LD"); Mangin 2012; Price 2008. Ours differs only as the *exact spectral form* of the LD **score**. |
| C2 non-local LD *score* ≈ PC score (leading-eigenvector truncation of an exact identity) | **apparently novel** (thin, framing) | Bercovich/cov-LDSC correct LD via PCs but never state the LD-*score* = local + eigenvalue-weighted-PC-loadings identity, nor tie the 2017 PC score to it as a truncation. |
| C3 unification (cov-LDSC / WhiLD / PC-score / `P_c`) | **novel as synthesis** | ingredients all exist; no unifying decomposition found. |
| C4 squared PC loadings as S-LDSC annotations | novel-ish, incremental; **no Bercovich collision** | S-LDSC/baseline-LD (Finucane 2015, Gazal 2017), graphREML 2024, FGLS-enrichment 2024 occupy the space; PC-loading annotation specifically not found. |
| C5 per-PC inflation spectrum (χ² scale) | **adjacent**; **no Bercovich collision** | Bhatia 2016 (signed, on β̂) already maps stratification to PCs; χ²-scale is the thin new slice. |
| C6 identifiability limit as a theorem | **adjacent** (known); **no Bercovich collision** | Huang 2025 (directional-LD estimand in admixture); differential-selection debates; novelty = clean formalization. |

## Nearest prior work (cite and position against)

- **Mangin et al. 2012** — the origin of the idea Bercovich proves: LD measures
  (`r²_S`, `r²_V`, `r²_VS`) that correct `r²` for **structure *and* relatedness**
  via a kinship matrix. The relatedness-corrected `r²_VS` is the closest published
  ancestor to `P_c`'s `V_c = σ_g²K + σ_e²I` whitening philosophy — but on pairwise
  `r²`, not on LD scores, and via a moment correction rather than a GLS-precision
  operator ([10.1038/hdy.2011.73](https://doi.org/10.1038/hdy.2011.73)).
- **Bercovich et al. 2025** — structure-adjusted *pairwise r²* via a hard
  projection off the ancestry subspace, for pruning/clumping/FST/PCA; the closest
  competitor to C1 but a different operator and no heritability content
  ([10.1093/genetics/iyaf009](https://doi.org/10.1093/genetics/iyaf009);
  bioRxiv 2024.05.02.592187).
- **Luo et al. 2021** — cov-LDSC; PC-residualized **LD scores** for heritability —
  the direct application competitor
  ([10.1093/hmg/ddab130](https://doi.org/10.1093/hmg/ddab130)).
- **Vilhjálmsson & Nordborg 2013** — "The nature of confounding in GWAS"; the
  qualitative origin of the identifiability limit (the confounder is environment
  *and* the genetic background correlated with structure) — and, notably, by the
  project's own senior author, so the natural lineage to own
  ([10.1038/nrg3382](https://doi.org/10.1038/nrg3382)). The §14 structure-aliasing
  remark is its LD-score-second-moment formalization.
- **Bhatia et al. 2016** — PC-loading regression; per-PC stratification detection
  — but on the *first-moment* (signed) scale, framing the aligned component as
  confounding to be subtracted (the opposite stance to the aliasing remark's "it
  also contains genuine V_p^gen") ([bioRxiv 10.1101/076133](https://doi.org/10.1101/076133)).
- **Young et al. 2018 (RDR)** and **Kong et al. 2018 (transmitted/non-transmitted)**
  — external anchors (relatedness disequilibrium; parental genotypes) that *break*
  the structure-only aliasing; polygenic-selection debate (Berg & Coop 2014;
  Berg 2019; Sohail 2019; Zaidi & Mathieson 2020, adjacent) is the empirical
  backdrop.
- Bulik-Sullivan et al. 2015 (LDSC, [10.1038/ng.3211](https://doi.org/10.1038/ng.3211));
  Finucane et al. 2015 (S-LDSC, [10.1038/ng.3404](https://doi.org/10.1038/ng.3404));
  Gazal et al. 2017 (baseline-LD).
- Price et al. 2008 — long-range LD confounds genome scans in admixed populations
  (AJHG; PMID 18606306). (The "structure/relatedness-corrected LD is not new"
  point is the Mangin 2012 entry above — same paper.)
- Active enrichment/heritability space to acknowledge: **graphREML**
  (medRxiv 2024.11.04.24316716), **FGLS functional enrichment** (Cell HGG Adv
  2024, S2666-2477(24)00011-3), **LDSC++** (medRxiv 2025.04.07.25324946),
  **Huang 2025** (admixed-h² estimand, [10.1093/genetics/iyaf100](https://doi.org/10.1093/genetics/iyaf100)).
- Efficiency: a Python-3 LDSC *platform* preprint (Dec 2025) — fast tooling alone
  is not novel. "Susceptibility of PGS to ancestry stratification" (Dec 2025) —
  adjacent to the diagnostics.

## Positioning recommendation

Frame as a **unifying-perspective methods paper**: the exact LD-score decomposition
identity (`ℓ = ℓ_local + ℓ_struct`, with the PC score as its leading-eigenvector
truncation), the unification of the existing structure-corrections under it, the
GLS-weighted `P_c` operator as the natural generalization of the hard-projection
corrections, and the LDSC-specific diagnostics + identifiability limit.
**Credit the structure-corrected-LD idea to Mangin 2012 / Bercovich 2025 (pairwise
r²) and cov-LDSC (LD scores); do not claim it, and do not claim a new estimator.**
Per-PC detection goes to Bhatia. Given the narrowed slice and the crowded active
space, the honest bar is: is *identity + unification + operator generalization +
diagnostics* enough of a clarifying contribution? The Bercovich gate is now
closed; the remaining open item is the full verified landscape sweep (re-run the
prior-art workflow after the session reset).
