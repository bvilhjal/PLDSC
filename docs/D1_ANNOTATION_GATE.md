# D1(a) gate: baseline-LD annotations vs ancestry PC loadings

**Date:** 2026-08-11. **Script:** [experiments/d1_annotation_loading_gate.py](../experiments/d1_annotation_loading_gate.py).
**Question (PROJECT_DIRECTIONS §4, D1, kill-test a):** does a real
baseline-LD-style annotation set correlate with 1000G ancestry PC loadings at
all? "If essentially none, D1 collapses to a footnote."

## What was run

- **Loadings.** Multi-ancestry 1000G FULL panel (phase 3, N = 2504). LOCO
  discipline: the GRM eigenspace is built from chr21 only (4,057 SNPs, MAF ≥
  0.05), and chr22 SNP loadings `a_j = corr(z_j, u_p)`, p = 1..10, are evaluated
  against it — no SNP's loading sits inside the GRM that defined the axes.
- **Annotations.** The real baseline-LD v2.2 set (Gazal et al. 2017;
  `baselineLD.22.annot.gz`, Zenodo record 10515792), matched by position to the
  cached chr22 SNP set: 8,781 of 9,654 SNPs matched, 96 non-constant
  annotations. Nothing recomputed or approximated.
- **Statistics per annotation:** max |corr| over the 10 loadings, and the joint
  R² of the annotation on all 10 (the quantity S-LDSC's design matrix sees).

## Results

| annotation | max \|corr\| | R² (10 PCs) |
|---|---|---|
| MAF_Adj_Predicted_Allele_Age | 0.124 | 0.018 |
| MAF_Adj_LLD_AFR | 0.123 | 0.025 |
| MAFbin1, MAFbin6–10 | 0.06–0.07 | ≈0.01 |
| Repressed_Hoffman, Intron_UCSC | ≈0.055 | <0.01 |
| **median over 96** | **0.025** | — |

Counts: |corr| > 0.05: 10 annotations; > 0.1: 2; > 0.2: 0. No annotation has
joint R² > 0.025.

## Reading

The signal is weak and concentrated exactly where it should be. The two
annotations above 0.1 are ancestry instruments *by construction* — `LLD_AFR`
is African low-LD, and MAF-adjusted predicted allele age is ancestry-linked.
Mainstream functional annotations top out at ≈0.055, and the MAF bins sit at
0.06–0.07 (expected: allele-frequency spectra differ by ancestry). In a
single-ancestry (EUR) analysis the correlations would be smaller still — these
are FULL-panel continental axes.

Against the gate's binary wording this is close to "essentially none" for the
mainstream baseline-LD set. What keeps it from being a clean kill is the
mechanism D1 is about: in
[raw_loading_control.py](../experiments/raw_loading_control.py) the S-LDSC
self-term pathology appeared when the *underlying* SNP-level correlation
between a PC loading and its own self-dropped LD-smeared score was only
0.10–0.14 — the self term inflates that to 0.82 collinearity inside the
LD-smeared annotation. That 0.10–0.14 is a *different quantity* from what this
gate measures (corr(annotation, PC loading)); the two are related only by
analogy, and the analogy is the weakest link in this verdict: it equates
"an annotation correlated 0.12 with an ancestry axis" with "a loading whose
smearing self-correlates at 0.10–0.14", without measuring the annotation's
own LD-smearing behaviour. `LLD_AFR` and predicted allele age sit at 0.12,
inside that window under the analogy. The pathology threshold is not "high
correlation"; it is "correlation in the range the self term can weaponize" —
and whether an annotation-loading at 0.12 is in that range is precisely what
kill-test (b) must establish rather than assume.

## Caveats on this measurement

- **chr22 only, chr21-only GRM (~4k SNPs).** Loading estimation noise
  attenuates all correlations downward. The top continental axes have enormous
  eigenvalues, so attenuation should be modest for the leading PCs, but the
  numbers here are lower bounds, not exact values.
- 91% position match between the annotation file and the cached SNP set; no
  strand/allele reconciliation was needed (annotations are positional), but the
  unmatched 9% is not characterized.

## Verdict

**D1 survives in a narrowed form.** The claim shrinks from "a general
correction to S-LDSC annotation practice" to: *ancestry-informative annotations
(LLD_AFR, predicted allele age, MAF bins) used in multi-ancestry or admixed
S-LDSC carry enough loading signal that the self-term collinearity can
produce confounding-driven τ — and the raw-loading control is the fix.* Whether
0.12 is enough to produce τ false positives of practical size is an empirical
question the gate cannot answer; it is exactly what kill-test (b) measures
(weighted S-LDSC on simulated multi-ancestry data with the *real* LLD_AFR /
allele-age annotations, τ_CONF with vs without the raw-loading control,
M ≥ 10⁵, block-jackknife SEs).

If (b) shows τ_CONF ≈ 0 even without the control, D1 is a footnote after all.
If it shows a false positive the control removes, the narrowed claim is the
paper, and the scoped prior-art sweep (Finucane 2015, Gazal 2017, Bhatia 2016,
the repo's own PC-annotation experiment) gates the writing.
