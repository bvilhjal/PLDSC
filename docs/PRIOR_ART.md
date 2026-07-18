# Prior-art positioning (preliminary)

Date: 2026-07-18

**Status: preliminary.** A light inline literature scan (targeted web searches +
the web-verified prior-art work from the applications panel). The full six-cluster
verified scan was launched as a workflow but hit a session usage limit before
running; re-run it (`workflows/scripts/pldsc-prior-art-*.js`) after the reset for
a complete, source-checked verdict.

## Headline

The reframed project (a local/non-local decomposition of the LD score, with the
non-local part identified as the PC score, unifying existing structure-corrections,
plus structure diagnostics) is a **framework/synthesis + diagnostics** contribution
— **not already published as a whole, and not a new estimator/tool.** The
qualitative framework is established; the novel slices are real but thin, and the
paper lives or dies on careful positioning against cov-LDSC and Bhatia et al.

## Per-claim assessment

| Claim | Verdict | Nearest prior work |
|---|---|---|
| C1 LD score = local (linkage) + non-local (structure) | adjacent (concept known) | cov-LDSC frames "long-range admixture LD" vs local; Price 2008; Bercovich 2025 |
| C2 non-local LD score = PC score (explicit identity) | apparently novel (as a stated identity) | 2017 manuscripts noted the empirical correlation only |
| C3 unification of cov-LDSC / WhiLD / PC-score / P_c | novel (as synthesis) | each exists separately; no unifying decomposition found |
| C4 squared PC loadings as S-LDSC annotations | novel-ish, incremental | Gazal 2017 baseline-LD uses LD/MAF annotations, not PC-loading |
| C5 per-PC inflation spectrum (chi^2 scale) | adjacent | Bhatia 2016 maps stratification to PCs (signed, on beta-hat); our chi^2-scale is the new slice |
| C6 identifiability limit as a theorem | adjacent (known folklore) | differential-selection / ancestry-gradient debates; novelty = clean formalization |

## Nearest prior work (cite and position against)

- Bulik-Sullivan et al. 2015 — LDSC ([10.1038/ng.3211](https://doi.org/10.1038/ng.3211)).
- Finucane et al. 2015 — S-LDSC / partitioned h² ([10.1038/ng.3404](https://doi.org/10.1038/ng.3404)).
- Gazal et al. 2017 — baseline-LD model (LD/MAF annotations) [DOI to verify].
- Luo et al. 2021 — cov-LDSC; **owns** "structure component of LD, removed via PCs / long-range admixture LD" ([10.1093/hmg/ddab130](https://doi.org/10.1093/hmg/ddab130)).
- Bhatia et al. 2016 — PC-loading regression; **owns** per-PC stratification detection/correction (signed, on effect estimates) ([bioRxiv 10.1101/076133](https://doi.org/10.1101/076133)).
- Price et al. 2008 — long-range LD and population structure [AJHG; DOI to verify].
- Bercovich et al. 2025 — LD / pruning-clumping in structured populations ([10.1093/genetics/iyaf009](https://doi.org/10.1093/genetics/iyaf009)).
- Recent preprints to check: a Python-3 LDSC *platform* (Dec 2025) — efficiency tooling exists, so a fast implementation alone is not novel; "Quantifying the susceptibility of polygenic scores to ancestry stratification" (Dec 2025) — adjacent to the diagnostics.

## Positioning recommendation

Frame as a **unifying-perspective methods paper**: the decomposition identity
(`ℓ_struct = PC score`, validated r=0.99), the unification of cov-LDSC / WhiLD /
PC-score / P_c under it, and the diagnostics (per-PC spectrum, structure-robust
S-LDSC enrichment) with the identifiability limit as a stated result. Credit the
qualitative concept to cov-LDSC (structure LD component), Bhatia (PC-loading
regression), and Price/Bercovich (long-range structure LD). Do **not** claim a
novel estimator or a killer tool; the value is clarity, unification, and honest
diagnostics. Re-run the full verified scan before submission.
