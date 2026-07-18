# Prior-art positioning

Date: 2026-07-18

**Status: inline scan (fuller than the first pass, still not the full verified
workflow).** The six-cluster subagent workflow hit a session usage limit; this
verdict is from ~6 targeted web searches plus the applications-panel findings.
Re-run `workflows/scripts/pldsc-prior-art-*.js` after the reset for a fully
source-checked version. The inline scan already changed the verdict once (see
Bercovich, below), so treat the novelty estimate as provisional.

## Headline

The reframed project (local/non-local LD-score decomposition; non-local part =
PC score; unification of the structure-corrections; structure diagnostics) is a
**framework/synthesis + diagnostics** contribution. The novel slice is **thinner
than the first pass suggested**: the qualitative framework and even
"structure-corrected LD from PC residuals" are already published (esp. **Bercovich
et al. 2025**), and the heritability/enrichment space is actively crowded
(cov-LDSC, graphREML, FGLS-enrichment, LDSC++). The genuinely distinct pieces are
the **explicit LD-*score* decomposition identity tying `ell_struct` to the PC
score**, the **unification** of cov-LDSC / WhiLD / PC-score / `P_c`, and the
**LDSC-specific diagnostics + identifiability limit**. This can still be a
unifying-perspective paper, but it must position hard against Bercovich 2025 and
cov-LDSC, and it is not a new estimator.

## Per-claim assessment

| Claim | Verdict | Nearest prior work |
|---|---|---|
| C1 LD = local (linkage) + non-local (structure) | **largely published** | **Bercovich 2025** (structure-corrected LD via PC-residualized genotypes; two-locus Wahlund effect); cov-LDSC ("long-range admixture LD"); Price 2008; "Novel LD measures correcting structure/relatedness", Heredity 2011 |
| C2 non-local LD *score* = PC score (explicit identity) | **apparently novel** (thin) | closest is Bercovich (structure LD via PCs) + the 2017 manuscripts' empirical correlation; but the explicit LD-score = local + PC-score identity, tied to LDSC, not found stated |
| C3 unification (cov-LDSC / WhiLD / PC-score / P_c) | **novel as synthesis** | ingredients all exist; no unifying decomposition found |
| C4 squared PC loadings as S-LDSC annotations | novel-ish, incremental | S-LDSC/baseline-LD (Finucane 2015, Gazal 2017), graphREML 2024, FGLS-enrichment 2024 occupy the enrichment space; PC-loading annotation specifically not found |
| C5 per-PC inflation spectrum (chi^2 scale) | **adjacent** | Bhatia 2016 (signed, on beta-hat) already maps stratification to PCs; chi^2-scale is the thin new slice |
| C6 identifiability limit as a theorem | **adjacent** (known) | Huang 2025 (directional-LD estimand in admixture); differential-selection debates; novelty = clean formalization |

## Nearest prior work (cite and position against)

- **Bercovich et al. 2025** — LD from PC-residualized genotypes in structured
  populations; the closest competitor to C1/whitening
  ([10.1093/genetics/iyaf009](https://doi.org/10.1093/genetics/iyaf009);
  bioRxiv 2024.05.02.592187).
- **Luo et al. 2021** — cov-LDSC; structure component of LD removed via PCs
  ([10.1093/hmg/ddab130](https://doi.org/10.1093/hmg/ddab130)).
- **Bhatia et al. 2016** — PC-loading regression; per-PC stratification detection
  ([bioRxiv 10.1101/076133](https://doi.org/10.1101/076133)).
- Bulik-Sullivan et al. 2015 (LDSC, [10.1038/ng.3211](https://doi.org/10.1038/ng.3211));
  Finucane et al. 2015 (S-LDSC, [10.1038/ng.3404](https://doi.org/10.1038/ng.3404));
  Gazal et al. 2017 (baseline-LD).
- Price et al. 2008 — long-range LD confounds genome scans in admixed populations
  (AJHG; PMID 18606306). "Novel LD measures correcting structure/relatedness",
  Heredity 2011 (structure/relatedness-corrected LD is not new).
- Active enrichment/heritability space to acknowledge: **graphREML**
  (medRxiv 2024.11.04.24316716), **FGLS functional enrichment** (Cell HGG Adv
  2024, S2666-2477(24)00011-3), **LDSC++** (medRxiv 2025.04.07.25324946),
  **Huang 2025** (admixed-h² estimand, [10.1093/genetics/iyaf100](https://doi.org/10.1093/genetics/iyaf100)).
- Efficiency: a Python-3 LDSC *platform* preprint (Dec 2025) — fast tooling alone
  is not novel. "Susceptibility of PGS to ancestry stratification" (Dec 2025) —
  adjacent to the diagnostics.

## Positioning recommendation

Frame as a **unifying-perspective methods paper**: the LD-score decomposition
identity (`ell_struct = PC score`, validated r=0.99), the unification of the
existing structure-corrections under it, and the LDSC-specific diagnostics + the
identifiability limit. **Credit the structure-corrected-LD idea to Bercovich 2025
and cov-LDSC** (do not claim it), and per-PC detection to Bhatia. Do not claim a
new estimator or novel structure-corrected LD. Given the narrowed slice and the
crowded active space, the honest bar is: is the *identity + unification + the
diagnostics* enough of a clarifying contribution? Re-run the full verified scan
(esp. reading Bercovich 2025 in full) before committing.
