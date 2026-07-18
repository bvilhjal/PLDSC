# Theory/methods note — outline (working draft)

Date: 2026-07-18. Direction chosen by the user after the within-family kill-test
came back negative: **(b) a clean theory/methods note** (no differentiated
application). Working title candidates:

- "The local/non-local decomposition of the LD score, and what structure-only
  corrections can and cannot do."
- "One decomposition behind cov-LDSC, WhiLD, the PC score, and the mixed-model
  LD score — and its identifiability limit."

## What the note claims (and does not)

**Claims (the standing contributions):**
1. **Decomposition.** The genome-wide LD score is exactly the GRM spectral
   decomposition, `sum_k r_jk^2 = (M/N) sum_p lambda_p corr(x_j,u_p)^2`
   (METHOD §14). Local (linkage) vs non-local (structure) is a partition of this
   sum; the 2017 PC score is the leading-eigenvector truncation of the non-local
   part.
2. **Unification.** cov-LDSC, the legacy WhiLD score, the PC score, and the
   mixed-model operator `P_c` are all operations on this decomposition (removing
   or measuring `ell^struct`); the published structure-corrected-LD operators
   (Mangin 2012 / Bercovich 2025 on r^2, cov-LDSC on LD scores) are *hard rank-K
   projections*, and `P_c` is their eigenvalue-weighted GLS generalization.
3. **Identifiability limit (the intellectual core) — VERIFIED + corrected.**
   Formalized and adversarially verified (3 lenses; see METHOD §14 "structure
   aliasing" remark). Corrected statement: in the marginal `E[z_j^2]` the
   coefficient on `corr(x_j,u_p)^2` is `τ²Mλ_p + V_p^gen + V_p^env`; the
   *directional* genuine variance `V_p^gen` and the confounding `V_p^env` enter
   through the identical realized scalar, so they are **exactly aliased** (not
   merely collinear) — unidentifiable from the second moments of any structure-only
   statistic, at any N or M. **Scope:** only the rank-1 *directional* component is
   aliased; LD-borne heritability (isotropic footprint `τ²Mλ_p`, pinned by the
   local slope; and a heteroscedastic-along-axis architecture, which leaves a
   `corr²×ℓ^local` interaction env lacks) stays identifiable. **Corollary:** every
   structure-only correction (cov-LDSC/P_c/PC-score/Bercovich) conservatively
   co-attenuates the *directional* real ancestry-aligned h². **Breaks only** with
   external anchors — cross-population replication or parental genotypes (RDR
   Young 2018; transmitted/non-transmitted Kong 2018). Within-family does NOT
   separate them — it *sidesteps* the shared axis (why within-family h² < pop h²;
   the repo's own `M_B R M_B = 0.5 M_B`). **Lineage to own:** Vilhjálmsson &
   Nordborg 2013 (the senior author's own paper) is the qualitative origin; state
   this as a named *remark*, a formalization, not a new impossibility.

**Explicitly does NOT claim:** a new estimator; a novel structure-corrected LD
(that is Mangin/Bercovich/cov-LDSC); a differentiated LDSC application (the
summary-stat, S-LDSC-annotation, and within-family angles were all tested and do
not differentiate PLDSC — Sections 13, `map_ldsc_miscalibration.py`,
`sldsc_pc_annotation.py`, `within_family_killtest.py`); a new impossibility
result (the identifiability limit formalizes a known confounding).

## Proposed section structure

1. **Introduction** — LD score has a structure component; the qualitative fact is
   old (Price 2008, cov-LDSC, Bercovich 2025); what is useful is the exact
   decomposition as an organizing device and its consequences.
2. **The decomposition** — the spectral identity; local vs non-local; PC score as
   truncation; the honest magnitudes (windowed ~5%, genome-wide ~50%).
3. **Unification** — the four methods as one ladder; hard projection vs GLS `P_c`;
   the exactness/assumptions table.
4. **The identifiability limit** — the proposition + proof; the corollary
   (conservatism of every structure-only correction); what breaks it
   (within-family), tying to why within-family h^2 is the ancestry-robust
   standard and why PLDSC adds nothing there.
5. **Diagnostics (illustrative, not headline)** — per-PC inflation spectrum;
   structure-robust S-LDSC enrichment as a *sensitivity analysis* (conservative,
   per the limit), with honest caveats.
6. **Scope / negative results** — the honest record: the estimator angle reduces
   to cov-LDSC; within-family kill-test negative. This is a feature of the note
   (it says precisely where the framework does and does not buy anything).
7. **Prior art** — Mangin 2012, Bercovich 2025, cov-LDSC, Bhatia 2016, the
   selection-vs-stratification literature, Huang 2025 (PRIOR_ART.md).

## Open tasks for the note

- [x] Finalize the identifiability proposition into METHOD §14 — DONE (verified +
      corrected; the "structure aliasing" remark). The most defensible genuinely-new
      slice is the **corollary** (the whole LD-score correction family inherits the
      conservatism) — foreground it; frame the rest as formalizing V&N 2013.
- [x] A clean numerical illustration of the *scoped* limit — DONE
      (`experiments/identifiability_aliasing.py`): directional-gen and env are
      identical to 2 dp (exact aliasing); heteroscedastic-along-axis leaves a
      larger `corr²×ℓ^local` interaction (identifiable via the local channel).
      Corrected moment coefficient `τ²Mλ_p` separately confirmed (corr 0.998).
- [ ] Decide venue/length (short methods note / preprint vs a Genetics-style
      perspective) — user call.
- [ ] Full prior-art landscape sweep to confirm the unification niche is
      unoccupied (deferred; the Bercovich gate is closed).

## Honest self-assessment (carry into the write-up)

The theory is real but textbook-adjacent (the decomposition is standard linear
algebra; the limit formalizes a known confounding). The note's value is
*clarity/unification + an unusually honest map of what does not work*, not a new
method. Pitch and venue should match that — a concise, honest unifying note, not
an over-sold methods paper.
