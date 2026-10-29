# When does population structure actually break LDSC? — note outline

Date: 2026-07-23. Status: **stopped as a standalone note (2026-09-25)** — its
content folds into the live project as MULTI_ANCESTRY.md §5's pooled case (see
STRATEGY.md §6, "what we stop"); §4's two-regime intercept numbers remain
source material. The line below records how it started. Direction chosen after the
relatedness-pivot idea (taking `P_c` to livestock/plant/founder panels) was
argued down — see "Why this and not the estimator" below.

This is a **different note** from [NOTE_OUTLINE.md](NOTE_OUTLINE.md). That one
leads with theory (decomposition / unification / identifiability) and the
prior-art sweep graded all four of its claims *partial* or *occupied*. This one
leads with **empirical magnitudes and a decision rule**, and its central content
— how large the windowed structure share actually is, and which moment each
nuisance lands on — is not in print anywhere. It reuses the same machinery and
most of the same text; it changes what is foregrounded.

## The question the field actually has

Practitioners run cov-LDSC, or don't, largely by superstition. The literature
establishes *that* LD has a structure component (Price 2008; Mangin 2012;
cov-LDSC; Bercovich 2025) and *that* structure correction co-attenuates genuine
ancestry-aligned h² (Huang 2025). Nobody has published **how big the effect is
in the regime LDSC operates in**, or a rule for when it matters. That gap is
this note's contribution, and it is an empirical gap, not a theoretical one —
which is precisely why the prior-art wall that stopped the theory note does not
stand in front of it.

## The one-paragraph answer

Two nuisances, two different moments, two different remedies:

| nuisance | magnitude | lands on | remedy |
|---|---|---|---|
| **Relatedness** in the LD panel | `O(1/N_eff)` per pair; up to ~19 LD-score units accumulated at S=2000 in a 16-sib design | **intercept** (near-uniform across variants) | compute LD from unrelateds; if impossible, the `N_eff` correction of METHOD §10.1 |
| **Population structure** | `O(1)`, does not vanish with N — but only ~0% (EUR) to ~29% (multi-ancestry) of the *windowed* score at LDSC's ±1 cM window | **slope** | cov-LDSC, and a handful of PCs suffices (B3) |

The headline practitioner result: **in a single-ancestry cohort the windowed
structure share is ~0%, so cov-LDSC is unnecessary; in a multi-ancestry panel it
is ~29% at LDSC's default ±1 cM window, and 2–10 hard-projected PCs capture
essentially all of what a full GLS whitening would remove.** The genome-wide
share (~90% in a multi-ancestry panel) is the number people intuit, and it is
irrelevant, because LDSC never regresses on the genome-wide score.

**Quote the cM number, not the physical one.** The windowed share depends on both
the width *and the unit* of the window, because a fixed physical window spans
variable genetic distance: on the FULL panel it is **29.3% at ±1 cM**, 17% at
±250 kb, and ~11% at ±0.25 cM
([realdata_genome_loco.py](../experiments/realdata_genome_loco.py);
[realdata_1000g_benchmarks.py](../experiments/realdata_1000g_benchmarks.py) B2).
LDSC's default window is ±1 cM, so **29%** is the practitioner-facing figure and
the ±250 kb result belongs in the window-unit remark. Earlier drafts of this note
led with 17%, understating the effect by ~1.7×.

## Section structure

1. **The two-moment framing.** Which nuisance hits the intercept, which hits the
   slope, and why that determines whether you should care. This is the note's
   organizing idea and it is cheap to state.
2. **Structure: how much is actually in the window?** The B1–B6 benchmarks
   ([realdata_1000g_benchmarks.py](../experiments/realdata_1000g_benchmarks.py)):
   density-robustness (B1), the LOCO-vs-in-chromosome operator bracket (B2), the
   ancestry gradient (B6), the MAF gradient (B5), plus the genome-wide LOCO
   driver's cM windows. Headline: genome-wide ~90% vs windowed ~29% (±1 cM) in
   the multi-ancestry panel, ~0% in EUR — with the window-unit dependence
   (29% at ±1 cM, 17% at ±250 kb, ~11% at ±0.25 cM) stated as its own result.
3. **How many PCs is enough?** B3: top-2 eigenvectors capture 92% of the
   genome-wide score; hard-projecting 2 PCs leaves 8% of what `P_c` removes.
   Real-data support for "cov-LDSC ≈ the GLS operator," stated as practical
   guidance rather than as a negative result about `P_c`.
4. **Relatedness is an intercept problem.** METHOD §10.1 + Eq. (M22) +
   [relatedness_effective_df.py](../experiments/relatedness_effective_df.py).
   The `ν−1` correction is wrong by up to 4.7× on a related panel; the resulting
   LD-score bias is large but near-uniform, so h² survives and the intercept
   does not. Practical consequence: **an LDSC intercept computed against an
   in-sample related LD panel is not a clean confounding barometer.** The
   2026-08-11 gate round quantifies the practical size on real panels
   ([relatedness_d2_killtest.py](../experiments/relatedness_d2_killtest.py),
   [d2_intercept_propagation.py](../experiments/d2_intercept_propagation.py)):
   the naive denominator under-subtracts on every standard panel (N_eff ratios
   EUR 0.80 / AFR 0.86 / AMR 0.88 / FULL 0.55), which moves the intercept by
   ~0.005 at typical chi² slope (immaterial, SE ~0.01) but ~0.03 at biobank
   slope with the FULL multi-ancestry panel — so the barometer caveat is
   material exactly at biobank scale.
5. **The decision rule.** A short table: cohort type → is the windowed share
   material → what to run. Plus the honest boundary: the operator-dependence of
   the split (B2), and the fact that the correction is conservative (it removes
   real directional h² too — Huang 2025, METHOD §14 corollary).
6. **What we could not settle.** The identifiability limit as a *scope
   statement* rather than a headline claim: diagnostics flag axes, they do not
   certify confounding. Brief, credited to Huang 2025 / Sotoudeh 2026 /
   Vilhjálmsson & Nordborg 2013.

## Why this and not the estimator

The relatedness-pivot idea — take `P_c`, whose [m3](../experiments/RESULTS.md)
advantage over cov-LDSC is real and specific to relatedness, to populations where
relatedness dominates (cattle, maize, *A. thaliana*, founder isolates) — fails
for three reasons, in increasing order of severity:

1. **Wrong moment.** The relatedness contamination of the LD score is
   near-uniform, so it lands on the intercept, not the slope. The mechanism that
   makes it large is the mechanism that makes it harmless for h².
2. **The competitor is a filtering flag.** Relatedness in an LD panel is fixed by
   dropping relatives or using an external panel, not by a GLS operator.
3. **Motivational.** The rebuttal to (2) — that in a breed or a selfing panel you
   *cannot* drop relatives — is correct, and fatal to the pivot rather than to
   the objection: those are exactly the settings where you own the individual-
   level genotypes and should run GREML/GBLUP instead. LDSC exists to work from
   summary statistics when individual data is unavailable; that constraint does
   not bind in breeding programs.

Untested residual: human founder isolates working from consortium sumstats
(FinnGen, deCODE, Sardinia) sit between the regimes — high relatedness *and* a
genuine summary-statistics constraint. Probably still dead (FinnGen runs SAIGE,
so the association step is already LMM-corrected, and at N≈500k the `1/N_eff`
floor is negligible again), but it has not been kill-tested.

Argument (1) is not merely a reason to drop the pivot — it is **Section 1 of this
note**, and the sharpest thing to come out of the post-archive work.

## Open tasks

- [ ] **Prior-art check on the two-moment framing.** The individual facts are
      known; the framing and the empirical magnitudes are the claim. Given this
      project's history — four theory claims graded partial-or-occupied — do this
      *before* writing. Specifically: does any cov-LDSC/LDSC-methods paper already
      state "relatedness → intercept, structure → slope" explicitly?
- [ ] **Genome-wide LOCO run.** B1–B6 are three chromosomes with physical-kb
      windows and an intercept-only `P_c` metric. The note needs all 22
      autosomes, a recombination map (cM windows), and LOCO throughout. This is
      the single largest piece of work and it is shared with the
      [#4 per-axis study](REAL_DATA_SURVEY.md).
- [ ] **One real trait end-to-end**, to make the decision rule concrete rather
      than panel-only: same sumstats, LD scores from an EUR panel vs a
      multi-ancestry panel, show the h² and intercept difference.
- [ ] Venue/length — a short methods/practice note; *Bioinformatics* application
      note, a Genetics methods note, or bioRxiv. User call.

## Honest self-assessment

The theory content is all pre-existing; the novelty is *quantitative and
practical*. That is a weaker claim than the archived note attempted and a much
more defensible one — the failure mode of the archived note was claiming
theoretical novelty that a sweep could refute, and empirical magnitudes cannot be
refuted that way, only superseded. The risk is different: a reviewer asks "so
this is just a cov-LDSC benchmarking paper?" The answer has to be the two-moment
framing in Section 1, which is why the prior-art check on it is the gating task.
