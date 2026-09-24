# PLDSC strategy: LD Score regression on multi-ancestry data

Date: 2026-09-25. This replaces the three-track version of 2026-09-24, which
remains in git history. The GC-versus-stratification probe that version ran
is recorded in [experiments/README](../experiments/README.md).

## 1. The direction

**PLDSC asks one question: what goes wrong when LDSC-based tools run on
multi-ancestry data, how large is the error, and what is the fix?**

Four reasons for choosing this question:

1. **It is the only regime where this project's own measurements are
   material.** The windowed structure share is about 0% in EUR and 29% in the
   1000G multi-ancestry panel (README). The new result in §2 is larger still.
2. **The field is there.** Multi-ancestry meta-analyses routinely run LDSC
   with European LD scores. GBMI did so, with a heuristic justification
   (landscape check, 2026-09-24).
3. **The sibling family is going there.** LDpred3-atw is the cross-ancestry
   successor. SMARTpred's ready-made references are UKB-EUR, and it gates on
   the LDSC intercept. LDpred3 initializes h² from LDSC.
4. **The claims are measurements and corrections.** Those are the only kind
   of claim that has survived this project's prior-art sweeps.

The estimator line (`P_c`, the theory note) stays concluded.

## 2. Where it stands

**Table 1. Built and measured (2026-09-25).**

| piece | where | status |
|---|---|---|
| Math: the LD score a meta-analysis follows (MA2), the bias factor κ of a single-ancestry score (MA6), the noise-floor asymmetry (MA8), pooled ≡ inverse-variance meta (MA12) | [MULTI_ANCESTRY.md](MULTI_ANCESTRY.md) | derived |
| Code: `cross_ancestry_ld_scores`, `meta_weights`, `meta_ld_scores`, `single_reference_bias`, `pooled_ld_scores` | `pldsc/multiancestry.py` | 15 unit tests (81 in the suite), including two exact in-sample identities |
| Simulation with out-of-sample panels | [multiancestry_meta_sim.py](../experiments/multiancestry_meta_sim.py) | default run FAILED its pre-stated calibration rule (small-panel errors-in-variables). The post-failure rerun passed; predicted κ matched within 0.004 |
| First real-data step, 1000G chr22 | [multiancestry_1000g_meta.py](../experiments/multiancestry_1000g_meta.py) | pre-stated kill rule: **KEEP** |

**The headline so far.** Take a meta-analysis that is 75% EUR and 25% from one
other 1000G superpopulation, with trans-ancestry genetic correlation ρ = 0.8.
If LDSC regresses it on the EUR LD score, it underestimates h² by:

| other 25% | h² underestimated by |
|---|---|
| AFR | 25% |
| EAS | 16% |
| SAS | 13% |
| AMR | 9–10% |

These are κ = 0.75 / 0.84 / 0.87 / 0.90–0.91. They are stable to SNP
density (±0.015) and to LDSC-style weights. Most of the loss is **LD
mismatch**: even at ρ = 1, κ = 0.79 for EUR + AFR. The intercept moves little
(+0.00 to +0.03 of the excess χ² at 75% EUR).

Two negatives from the same run:

- A meta-analysis alone cannot separate LD mismatch from effect heterogeneity.
  The within and cross parts correlate 0.80–0.99 across SNPs
  (MULTI_ANCESTRY Table 5).
- The first κ computation overstated the loss. It had two defects, which a
  post-rule robustness check found. The pre-stated verdict and the corrected
  numbers are both recorded (MULTI_ANCESTRY Table 4).

## 3. The project

**Table 2. One question, three parts.**

| part | question | status | existing PLDSC work |
|---|---|---|---|
| **1. LD scores (core)** | Which LD score does a multi-ancestry GWAS follow, and what does the EUR score cost in h² and r_g? | math, code and chr22 done | spectral engine, cM windows, the 29% measurement |
| 2. Enrichment | Do ancestry-informative annotations (LLD_AFR, allele age) show false S-LDSC enrichment in multi-ancestry data? | pending: the D1(b) kill-test, rule already written | [D1_ANNOTATION_GATE](D1_ANNOTATION_GATE.md), raw-loading control |
| 3. Intercept | How much do relatives in the LD panel shift the intercept? | small effect (D2: ~0.005 typical, ~0.03 at biobank slope); one section | `effective_pair_df` |

## 4. Next steps, each with its gate

**Table 3. Order of work.**

| # | task | cost | gate |
|---|---|---|---|
| 1 | κ genome-wide: all 22 autosomes from the S3 mirror (~15 GB; the EBI FTP is blocked here) | 1–2 days | the rule-cell AFR κ stays ≤ 0.90; otherwise chr22 was unrepresentative and the claim shrinks |
| 2 | Prior-art sweep scoped to MA2, MA6, MA8 and the measured κ (Popcorn, S-LDXR, MAMA, the 2026 meta-LDSC preprint); verify the snippet-level citations | 1 day | if MA2/κ are in print, the paper becomes the measurement plus the casualty |
| 3 | **The casualty:** re-run one published multi-ancestry meta-analysis (GBMI or similar, reported ancestry mix) with the meta score instead of the EUR score | 1–2 weeks | h² moves by more than 2 SE, or the paper becomes a practice note |
| 4 | Genetic correlation: the cross-trait version of MA2, and what EUR scores do to r_g between two meta-analyses | 1 week | only if step 3 survives |
| 5 | Part 2 (the D1(b) kill-test) | 1–2 weeks | as written in D1_ANNOTATION_GATE |
| 6 | Family integration: SMARTpred "reference fit" (κ from a GWAS's reported ancestry mix) and LDpred3-atw Phase 0 | 1 week | after steps 1–3 |

## 5. Landscape (single-agent check, 2026-09-24; snippet-level unless noted)

- **LD mismatch in meta-analyses: partial.** GBMI used UKB-EUR LD scores on
  its multi-ancestry meta. A 2026 bioRxiv preprint reports LDSC underestimating
  h² on meta-analysed summary statistics without a correction; this matches the
  direction of §2. Popcorn and S-LDXR use cross-population LD scores for
  ancestry-specific summary statistics, not for a pooled meta. No
  composition-matched meta score and no panel-only bias predictor were found.
- **Relatives in LD panels: open.** Pan-UKB's LD code uses raw n in its r²
  correction (code read). The NYGC 3,202-sample 1000G set includes 602 trios.
- **False S-LDSC enrichment on ancestry-correlated annotations: mostly open.**
  Bhatia 2016 is prior art for the raw-loading control.
- **Admixture LD under global PCs: partial.** cov-LDSC recommends 20 cM
  windows for admixed samples. HAMSTA's LAD score is the nearest
  local-ancestry analogue.

## 6. What we stop

- New estimators or constructions, including the local-ancestry LD score.
  Four sweeps have closed every construction claim this project has made.
- `P_c` as a contribution.
- Prior-art sweeps on theory claims.
- The standalone WHEN_STRUCTURE_MATTERS note. Its content folds into Part 1's
  pooled case (MA11–MA12), where the 29% structure share is the residual after
  imperfect ancestry adjustment.

## 7. Risks

**Table 4.**

| risk | consequence | mitigation |
|---|---|---|
| chr22 is unrepresentative | κ moves | step 1 |
| MA2 or κ is already in print (MAMA, S-LDXR supplements) | novelty shrinks to the measurement | step 2 before writing; the casualty is the differentiator |
| No published meta moves by > 2 SE | no casualty | publish as a practice note with the κ table and a `pldsc` function |
| ρ is unknown for most traits, and a meta cannot identify it (MULTI_ANCESTRY Table 5) | κ is a bracket, not a point | report κ over ρ ∈ [0.6, 1]; use S-LDXR/Popcorn ρ where published |
| The model assumptions (A1 effect scale, η = 1, one genetic map) | κ shifts | MA2 carries η and τ_ab,j explicitly; add sensitivity runs in step 1 |
