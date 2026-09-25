# Project review: why it stalled, and how to make it interesting

Date: 2026-07-25. Independent review of the full record (README, [METHOD](METHOD.md),
[DEVELOPMENT_PLAN](DEVELOPMENT_PLAN.md), [PRIOR_ART](PRIOR_ART.md),
[NOTE_OUTLINE](NOTE_OUTLINE.md), [WHEN_STRUCTURE_MATTERS](WHEN_STRUCTURE_MATTERS.md),
[REAL_DATA_SURVEY](REAL_DATA_SURVEY.md), `experiments/`), plus one new experiment
run for this review ([raw_loading_control.py](../experiments/raw_loading_control.py)).

## 1. Verdict

**The archive decision was right for the project as posed — and the project was
posed the same way four times.** Every framing that died (`P_c` as an estimator,
the LMM-matched LD score, the within-family angle, the theory note) was a variant
of *"a better estimator of `h²` under structure."* That framing has a computable
ceiling, and this repo is the thing that computed it: the structure share of the
**windowed** LD score — the only part LDSC ever regresses on — is ~0%
single-ancestry and 17–29% multi-ancestry. An estimator improvement is bounded by
that number, and in the regime where the number is large, 2–10 hard-projected PCs
already capture it (B3). Each negative kill-test was a rediscovery of the same
ceiling from a new angle.

What the record does **not** establish is a ceiling on the project's other output.
Note the asymmetry in what survived the sweeps:

| kind of claim | examples | sweep outcome |
|---|---|---|
| **construction** (an operator, an identity, a unification, a theorem) | `P_c`, DECOMP, UNIF, ALIAS | all four graded *partial* or *occupied* |
| **measurement / correction** (a number, an error) | windowed vs genome-wide shares; window-unit dependence; `ν−1` → `tr(Ã)²/tr(Ã²)` | none refuted; two are new facts |

Construction claims are the most refutable kind — one adversarial sweep kills
them, and this project ran four. Measurement claims can only be *superseded*, and
this project is unusually well equipped to make them: an exact `O(M N²)`
genome-wide spectral score with no `M × M` matrix, LOCO, cM windows, 59 passing
tests (81 as of 2026-09-25 — see [README.md](../README.md)), and a demonstrated
willingness to publish its own negative numbers.

**The one-line recommendation: stop asking "can I estimate `h²` better?" and start
asking "what is the field computing wrong, and by how much?"** The instrument is
built. Nothing in the archive record argues against pointing it at that question.

## 2. Three structural reasons it stalled

1. **The estimand never changed.** Every angle competed with cov-LDSC on `h²`.
   Because the windowed structure share is small in the single-ancestry regime and
   PC-capturable in the multi-ancestry one, *any* operator lands in the same place.
   §13 reads as a discovery; it was a consequence of the estimand choice.
2. **Novelty was sought where the field is 15 years old.** LD-with-structure has
   Price 2008, Mangin 2012, Bhatia 2016, cov-LDSC, Bercovich 2025, Huang 2025.
   A new *construction* there is a coin flip against a sweep. A new *number* is not.
3. **The kill-tests were run at the scale of the sim, not the scale of the claim.**
   [local_adaptation_channel.py](../experiments/local_adaptation_channel.py) closed
   the one surviving channel at `N = 600–1200`, `M ≈ 2000`, with a hand-rolled
   3-term regression whose null the experiment itself documents as badly
   anti-conservative (ENV median `t ≈ 3.8`). That is a fair verdict *for
   wild-population studies*, which is what it was scoped to. It is not a verdict on
   the same channel at `M ≈ 10⁶`, and §4.1 below shows the anti-conservatism was an
   artifact of the parameterization, not of the channel.

## 3. What the repo owns

Worth being explicit, because the directions below all spend these:

- `pldsc/genome.py` — exact spectral genome-wide score (`O(M N²)`, no `M × M`),
  cM windows, LOCO, validated to `3e-13` against brute force on real chr22.
- `pldsc.effective_pair_df` — the corrected LD-score noise floor under relatedness.
- A calibrated simulator with structure, admixture, and family designs.
- An honest, dated, adversarially-reviewed negative-results record. This is worth
  more than it looks: it makes any positive claim from this repo credible.

## 4. Four directions, ranked

### D1 (recommended). Population structure inherits into S-LDSC enrichment — and the fix is un-smearing, not projection

**The question.** S-LDSC's `τ_C` is estimated by regressing `χ²` on annotation LD
scores `ℓ_j(C) = Σ_{k∈W_j} r_jk² C_k`. If an annotation is correlated with a
structure axis, does stratification leak into `τ_C`, and how much?

**Why it is interesting.** It is a *false-positive* question about a method in
very heavy use, not an efficiency question about an estimator nobody runs. It has
a concrete failure mode, a cheap fix, and a testable magnitude. And it sits
exactly where this project's identifiability work has something to say.

**The observation this review adds.** Environmental confounding along axis `u_p`
enters `E[z_j²]` through the SNP's **own** loading `a_j = corr(x_j,u_p)²` and is
*not propagated by LD*. An ancestry-differentiated genetic architecture enters
through the **LD-smeared** score `ℓ_j(a) = Σ_{k∈W_j} r_jk² a_k`. Those are
separable regressors — but only once the self term is dropped, because `r_jj² = 1`
puts a whole copy of `a_j` inside `ℓ_j(a)`. **S-LDSC's annotation LD scores include
that self term.** That is the mechanical reason a smeared PC annotation cannot
separate the two channels, and it is fixable by adding the raw loading as its own
regressor. This is LDSC's own intercept-versus-slope logic (uniform = confounding,
LD-scaling = polygenic), applied per axis.

**Evidence** ([raw_loading_control.py](../experiments/raw_loading_control.py),
run for this review; `N=600`, `M≈2000`, admixed, 1 cM windows, seed 20260725):

| | measured |
|---|---|
| `corr(a, ℓ(a))` **without** the self term | +0.10 to +0.14 |
| `corr(a, ℓ(a))` **with** the self term (S-LDSC convention) | **+0.82** |
| median `t` under a *purely environmental* axis shift, S-LDSC-convention statistic | **+7.6 to +13.4** |
| median `t`, same data, raw/smeared split | **−0.31 to −0.08** |
| median `t`, same data, the naive interaction test | +0.65 → +1.34 (worsens with `N`) |

And on the practical case — `τ` on a functional annotation's LD score, under
`CONF` (no true enrichment, environmental axis shift) versus `ENRICH` (true 4×
enrichment, no confounding):

| design | annotation ancestry-correlated: `τ` CONF (want 0) | `τ` ENRICH (want > 0) | control annotation: `τ` ENRICH retained |
|---|---|---|---|
| A standard S-LDSC | **2.54e+01** | 2.77e+00 | 1.00× |
| B + smeared PC annot (the repo's fix) | −2.61e+00 | 1.15e−01 (0.04×) | 0.97× |
| C + **raw** loading | **8.5e−03 (0.00×)** | 3.23e−01 (0.12×) | **1.02×** |
| D + both | 5.4e−02 | 3.88e−01 (0.14×) | 1.01× |

Three readings, in decreasing order of confidence:

1. **Standard S-LDSC on an ancestry-correlated annotation reports a coefficient
   ~9× larger under pure environmental confounding (25.4) than under a genuine 4×
   enrichment (2.77).** This is the same failure
   [sldsc_pc_annotation.py](../experiments/sldsc_pc_annotation.py) found at ×3.5;
   at a stronger gradient it is much worse. It is the most publishable number in
   the repo right now.
2. **The raw loading dominates the smeared PC annotation on both axes.** It removes
   the false enrichment completely (0.00× vs B's −0.10×, i.e. B over-corrects into
   negative territory) and is essentially free on an ancestry-independent
   annotation (1.02× of genuine enrichment retained, vs B's 0.97×).
3. **Neither fix rescues genuine ancestry-correlated enrichment** (0.12× and
   0.04×). That is the §14 identifiability limit, unchanged and not claimed away.
   The honest framing is a *sensitivity analysis with a much better null*, not an
   unbiased estimator.

**Honest negative from the same run.** The split does **not** rescue the
local-adaptation channel's *power*: at a calibrated threshold, power on the
heteroscedastic arm is 0.06–0.19, no better than the naive test. It fixes the
null, not the sensitivity. Whether the channel has power at `M ≈ 10⁶` (where `t`
scales roughly as `√M`) is untested and is the real open question.

**Kill-tests, in order.** (a) How many real annotations are materially correlated
with 1000G PC loadings? Compute `corr(C, a)` for a baseline-LD-style annotation
set — one afternoon with `pldsc/genome.py`, and if the answer is "essentially
none," D1 collapses to a footnote. (b) Redo the above with proper weighted S-LDSC,
block-jackknife SEs, and `M ≥ 10⁵`. (c) Re-run one published enrichment with the
raw-loading control and see whether anything moves.

**Prior-art risk: MEDIUM.** Nearest known: S-LDSC (Finucane 2015), baseline-LD
(Gazal 2017, which already models LD-dependent architecture), Bhatia 2016
(first-moment PC-loading regression), and the repo's own PC-annotation experiment.
The specific move — *raw annotation alongside its LD score, self term removed, as
the confounding control* — I did not find, but this review ran no sweep and this
project's history says do one before writing a line.

**Cost:** kill-test (a) ≈ 1 day; the full study ≈ 3–4 weeks.

### D2. The LDSC intercept is contaminated by the LD panel, not just by the GWAS

**The question.** [§10.1](METHOD.md#101-relatedness-in-the-ld-panel-the-residual-df-correction)
found the `ν−1` noise-floor subtraction under-corrects by up to 4.7× on a related
panel, and that the resulting bias is near-uniform — so it lands on the
**intercept**. The intercept and attenuation ratio are used across the literature
as *the* confounding barometer and as a QC gate. So: for the LD panels people
actually use, how far off is the intercept?

**Why it is interesting.** It is a correctness claim about a number that gates
published conclusions, and it is already derived and implemented here. It is the
cheapest high-information thing in the whole repo.

**Kill-test (do this first, before anything else in this document).** Compute
`tr(Ã)²/tr(Ã²)` versus `ν` for the panels in real use: 1000G EUR, 1000G ALL,
UKB-unrelated, and a founder/sib design. If the ratio is ≈1 for the standard
panels, the finding is confined to founder cohorts and breeding panels and should
be written up as a one-paragraph caveat, not a project. **One day of work, and it
decides whether D2 exists.** The repo's own argument in
[WHEN_STRUCTURE_MATTERS](WHEN_STRUCTURE_MATTERS.md) ("Why this and not the
estimator") already suspects the answer is unfavourable at biobank `N`; that
suspicion is worth converting into a number rather than carrying as a prior.

**Prior-art risk: LOW–MEDIUM** on the specific correction; the "compute LD from
unrelateds" folklore is universal, but the quantitative intercept statement is not
something I found in print.

**Cost:** 1 day to decide; ≈2 weeks to write up if it survives.

### D3. Turn WHEN_STRUCTURE_MATTERS from a benchmark into an audit

[WHEN_STRUCTURE_MATTERS](WHEN_STRUCTURE_MATTERS.md) has the right frame — two
nuisances, two moments, two remedies — and it correctly anticipates the reviewer
question ("so this is a cov-LDSC benchmarking paper?"). The answer to that
question is not a better framing paragraph; it is **an audit with names attached**.
The interesting version asks: *the field applies EUR-panel LD scores to
increasingly multi-ancestry meta-analyses — what does that cost, and can it be
predicted from the panel alone?* The decomposition answers exactly this: the
applied LD score carries the **panel's** non-local component, not the cohort's,
and that mismatch is a slope error. A panel-computable predictor of the bias would
need no cohort genotypes — only the published ancestry composition.

Two changes to the existing plan:

- **The headline number should be ~29%, not 17%.** The 17% figure is a ±250 kb
  *physical* window. [realdata_genome_loco.py](../experiments/realdata_genome_loco.py)
  gives **29.3%** at 1.0 cM — and 1 cM is the window LDSC actually uses. The doc's
  central empirical claim is currently understated by a factor of ~1.7.
- **Add one real trait end-to-end** (already on the open-task list) and at least
  one *published* `h²`/`r_g` that moves. A benchmark with no casualty is a
  benchmark; a benchmark with a casualty is a result.

**Prior-art risk: MEDIUM** (LD reference-panel mismatch is an active area; the
decomposition-based *predictor* is the differentiated part). **Cost:** the
22-autosome run is the gating expense (~15 GB, already scripted via `--plan`).
Do it *after* D1(a) and D2's kill-test, since both are cheaper and either could
redirect the effort.

### D4 (high variance, listed for completeness). The minimal-anchor question

[§15](METHOD.md#15-theory-addendum-assortative-mating-and-genetic-nurture) ends on
a genuinely unresolved point: offspring GWAS plus offspring LD cannot identify
mate choice or parental nurture. The interesting follow-up is not another
impossibility statement but its converse — **what is the minimal external moment
that restores identification?** That has a clean answer shape (a small table:
anchor → what becomes identified), and it is the kind of question the §15 algebra
is already set up to answer.

**Prior-art risk: HIGH.** Border 2022, Kong 2018, Young 2018, Wu 2021 (DONUTS),
Zhang 2025 (DGREML) occupy most of the neighbourhood, and §15 is already graded
*occupied* line by line. Listed because the user asked what would be *interesting*,
and this is the most intellectually interesting thing left in the repo — not
because it is the best bet.

## 5. What to stop doing

- **Stop sweeping prior art on construction claims.** Any remaining one is dead on
  arrival; the sweeps are now confirming a settled result at ~700k tokens each.
- **Stop treating `P_c` as a candidate contribution.** Keep it as a comparator and
  as the thing that generates the operator bracket in B2. It is good infrastructure
  and a bad headline.
- **Don't spend the 22-autosome download until a decision depends on it.** Right
  now it is infrastructure with no attached decision. D3 attaches one; D1(a) and
  D2 are cheaper and might redirect it. *(Superseded, 2026-09-25: the live
  direction did attach a decision — the genome-wide κ is STRATEGY §4 step 1,
  so the download is now the first thing to spend. See
  [scripts/fetch_1000g.py](../scripts/fetch_1000g.py).)*

## 6. Suggested order

| # | task | cost | decides |
|---|---|---|---|
| 1 | `tr(Ã)²/tr(Ã²)` vs `ν` for real LD panels | 1 day | D2 exists or is a caveat |
| 2 | `corr(annotation, PC loading)` for a baseline-LD-style set | 1 day | D1 exists or is a footnote |
| 3 | D1 at scale: weighted S-LDSC, jackknife SEs, `M ≥ 10⁵` | 1–2 weeks | the headline number |
| 4 | prior-art sweep — **only** on whichever of D1/D2 survived | 1 day | write or don't |
| 5 | 22-autosome LOCO run + one real trait end-to-end | 2–3 weeks | D3's audit |

Steps 1 and 2 are each a day and each can close a direction. Both were run on
2026-08-11; outcomes and narrowed verdicts in §8 below.

## 7. Corrections to the current record — APPLIED (2026-07-25)

- [x] **README, "29 tests"** — the suite is **59** (`test_genome` 21, `test_ld` 19,
  `test_lmm` 13, `test_regression` 6). Fixed.
- [x] **The "~17%" windowed share quoted as a headline** — that is the ±250 kb
  *physical*-window figure. LDSC's default window is ±1 cM, where the FULL-panel
  share is **29.3%** ([realdata_genome_loco.py](../experiments/realdata_genome_loco.py)),
  so the practitioner-facing claim was understated by ~1.7×. The share depends on
  the window's **unit as well as its width**: 29.3% at ±1 cM, 17% at ±250 kb,
  ~11% at ±0.25 cM. Every place the number was quoted as a general headline now
  leads with the cM figure and keeps the physical one as the window-unit remark:
  [WHEN_STRUCTURE_MATTERS](WHEN_STRUCTURE_MATTERS.md) (two-moment table, headline
  paragraph, §2 outline), [README](../README.md),
  [DEVELOPMENT_PLAN](DEVELOPMENT_PLAN.md) §1, [NOTE_OUTLINE](NOTE_OUTLINE.md) §2,
  [METHOD](METHOD.md) §10, [REAL_DATA_SURVEY](REAL_DATA_SURVEY.md), and the
  `realdata_1000g_structure.py` correction block and B1–B6 verdict in
  [experiments/README](../experiments/README.md).
  **Not changed:** the B1–B6 result tables themselves and
  [EXTERNAL_REVIEW](EXTERNAL_REVIEW.md), which are dated records of specific
  ±250 kb runs and of a review round — correct as written, and rewriting them
  would falsify a history.
- [x] **`identifiability_aliasing.py`'s "finite-sample orthogonalization
  residual"** — localized: it is the self-term collinearity between `a_j` and
  `ℓ_j(a)` (`corr = +0.82` with the self term, `+0.14` without), and dropping the
  self term centres that null at zero. Recorded as a named correction in
  [METHOD §14](METHOD.md#14-the-localnon-local-decomposition-organizing-framework)'s
  scope paragraph, which no longer claims channel (ii) needs an empirically
  calibrated threshold. The MARGINAL verdict is unaffected — the fix is to the
  null, not the power.

## 8. Gate outcomes (2026-08-11)

The two one-day gates of §6 have now both been run, on real 1000G phase-3 data.

- **D2 — CLOSED at typical slope; LIVE at biobank-scale slope.** The tr-ratio
  gate ([relatedness_d2_killtest.py](../experiments/relatedness_d2_killtest.py))
  found the naive `nu - 1` under-subtraction is real on every standard panel
  (PC-adjusted `N_eff` ratios: EUR 0.80, AFR 0.86, AMR 0.88, FULL 0.55), but
  propagating it to the LDSC intercept
  ([d2_intercept_propagation.py](../experiments/d2_intercept_propagation.py))
  gives a shift of ~0.005 at chi2 slope 0.05 and ~0.03 at slope 0.3, against a
  typical real intercept SE of ~0.01. Against the propagation script's own
  pre-stated rule (a shift well under ~0.01 kills; at or above keeps the
  claim alive), the slope-0.05 arm kills but the slope-0.3 arm — ~3x the
  yardstick, and FULL (ratio 0.55) is the standard multi-ancestry reference —
  meets the keep-alive condition. The earlier "CLOSED" verdict (and the
  2026-08-11 commit title) overreached for that arm; corrected here. So: the
  intercept shift is immaterial at typical slopes, material at biobank-scale
  slope with a multi-ancestry panel — a quantified two-regime caveat for the
  note, not a paper. (The "corrected scores recover 1.00" line is by
  construction — the simulation defines the corrected score as truth — so the
  naive-vs-corrected delta, not that half, is the evidence.)
- **D1(a) — SURVIVES, narrowed.** The annotation-loading gate
  ([d1_annotation_loading_gate.py](../experiments/d1_annotation_loading_gate.py),
  full write-up [D1_ANNOTATION_GATE](D1_ANNOTATION_GATE.md)): on the real
  baseline-LD v2.2 set, max |corr| with FULL-panel PC loadings is 0.12, carried
  by the ancestry-informative annotations (LLD_AFR, predicted allele age; MAF
  bins 0.06-0.07; mainstream functional <= 0.055). The claim narrows to
  ancestry-informative annotations in multi-ancestry S-LDSC. Next gate: (b)
  at-scale S-LDSC with the real LLD_AFR / allele-age annotations, tau_CONF with
  vs without the raw-loading control — then the scoped prior-art sweep on the
  survivor, per §6.
