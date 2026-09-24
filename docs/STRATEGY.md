# Making PLDSC relevant: a strategy

Date: 2026-09-24. A review of the whole record (README,
[PROJECT_DIRECTIONS](PROJECT_DIRECTIONS.md), [D1_ANNOTATION_GATE](D1_ANNOTATION_GATE.md),
[WHEN_STRUCTURE_MATTERS](WHEN_STRUCTURE_MATTERS.md), [PRIOR_ART](PRIOR_ART.md),
`experiments/`) and of the sibling packages that now use LDSC in production.
It adds one small probe, [gc_vs_stratification.py](../experiments/gc_vs_stratification.py) (§3.1),
and a landscape check (§4). The probe **failed** its pre-stated rule, and §3.1
reports that result as it stands.

## 1. The verdict in one paragraph

The archive decision is right for the question it asked: can we estimate h²
better under structure? I do not recommend reopening it. The estimator's
ceiling is measured. On real 1000G data the windowed structure share is about
0% in EUR and 29% in the multi-ancestry panel at LDSC's ±1 cM window, and 2–10
PCs already capture it. What goes unused is the **instrument**. It has five
parts: an exact genome-wide spectral LD-score engine with LOCO and cM windows;
a corrected noise floor under relatedness; a mechanism (the `r_jj² = 1` self
term) that explains a false-positive mode in S-LDSC; a real-data gate showing
where that mechanism can bite (baseline-LD's ancestry-informative annotations);
and a credible negative-results record. Meanwhile the sibling family now
depends in production on the LDSC quantities this instrument audits:

- SMARTpred has an intercept gate.
- LDpred3 initializes h² from LDSC.
- gwfm and LDpred3-funct learn functional priors from all 96 baseline-LD v2.2
  columns, including LLD_AFR and allele age.
- The default reference is UKB-EUR.

The strategy is therefore **to make PLDSC the family's LD-reference auditor,
and to publish one measurement paper for the multi-ancestry era from it.**
Relevance comes from having consumers. What makes it cool is one sharp result
with named casualties.

## 2. Why the project reads as irrelevant today

**Table 1. Symptoms, causes, remedies.**

| symptom | cause | remedy |
|---|---|---|
| The first screen of the README says "CONCLUDED / archived … no differentiated contribution" | The record is honest, but it is placed where a visitor decides whether to read on | Lead with what the instrument measures and move the archive record below it. Keep every word of it (§7) |
| Nothing in the family imports PLDSC | `ldpred3/docs/family.md` lists it only as "concluded: PLDSC (archived)" | Give it three named consumers (Table 2) |
| Four framings died the same death | Every one competed with cov-LDSC on h², which is bounded by the windowed structure share | Change the question from *estimate* to *attribute and audit* ([PROJECT_DIRECTIONS](PROJECT_DIRECTIONS.md) §1 made this diagnosis; this document adds the consumers) |
| Every number comes from a sim or from 1000G chr21/22 | The 22-autosome run and one real trait were never attached to a decision | Attach them to Track C (§5), whose gate is a published estimate that moves |
| The API is research-grade (NumPy arrays in, no CLI, no file formats) | The package was never meant to be used by anyone else | Add a CLI and a sidecar format that SMARTpred can register (Track A) |

## 3. Where the family already depends on what PLDSC measures

**Table 2. Consumers in the sibling family.** Every row cites a file that
exists today.

| # | consumer | where | what it uses | what can go wrong that PLDSC measures | PLDSC module |
|---|---|---|---|---|---|
| C1 | SMARTpred trait QC | `smartpred/_runner_trait_qc.py`; `templates/_ldsc_caveat.html` | Free-intercept LDSC. It warns if the intercept falls outside [0.8, 1.2] or the ratio exceeds 0.5 (only with ≥ 10,000 variants) | The caveat says a free intercept "can reflect confounding, ancestry/LD mismatch, or model misspecification", but the service cannot say which | A2 fingerprint (which axis?), A1 panel card (is it the panel?) |
| C2 | LDpred3 `ldsc_init` | `ldpred3/pipeline.py`; `ldpred3/ldsc.py` | The LDSC slope initializes `h2_init` for LDpred3-auto | A slope error from panel/cohort composition mismatch (the D3 predictor) | A1, Track C |
| C3 | ppb z calibration | `ppb/docs/CALIBRATION.md` C2–C3 | `z_true = z_obs / √intercept`. The document says the intercept "cannot separate" GC from stratification | A raw-loading regressor might separate them. **Probe failed its rule, so not recommended (§3.1)** | none for now |
| C4 | gwfm / LDpred3-funct priors | `gwfm/gwfm/annotations.py`; `smartpred/_runner_funct.py` | Priors learned from all 96 baseline-LD v2.2 columns | D1: MAF_Adj_LLD_AFR and predicted allele age correlate up to 0.12 with ancestry loadings, so a stratified GWAS may teach the prior ancestry | A3 flags, Track B |
| C5 | LDpred3-atw | `README.md`, roadmap phases 0–1 | "Gap decomposed into LD-mismatch and effect-heterogeneity components"; "reproducible reference builds for at least three ancestry groups" | Every new reference build needs its structure share, N_eff/ν and λ spectrum recorded | A1 panel card |
| C6 | SMARTpred reference registration | `docs/setup.md` | Default references are UKB-EUR HM3/HM3+. The docs ask for "an external manifest with ancestry, genome build, sample selection" | The manifest is free text today; nothing in it is measured | A1 as a machine-readable manifest sidecar |

**Figure 1. Where the instrument sits.**

```text
  LD reference panel (PLINK + genetic map)              GWAS summary statistics
                 │                                                │
      ┌──────────┴──────────┐                                     │
      │  A1  panel card     │ structure share (±1 cM, LOCO),      │
      │  pldsc audit-panel  │ N_eff/ν, λ₁..λ₁₀, intercept-shift   │
      └──────────┬──────────┘ forecast, loadings sidecar a_jp     │
     ┌───────────┼────────────────┬───────────────────┐           │
     ▼           ▼                ▼                   ▼           ▼
 SMARTpred   LDpred3-atw     gwfm / LDpred3-funct   A2 fingerprint: χ²_j ~ 1 + ℓ_j + Σ_p a_jp
 manifest    phase-1 builds  A3 annotation flags         │
 (C6)        (C5)            (C4) ◄── Track B            ▼
                                                     SMARTpred trait-QC card (C1)
```

### 3.1 The one probe run for this review: GC versus stratification (C3)

The claim tested: ppb's two-term intercept confounds genomic-control (GC)
deflation with stratification. The PLDSC mechanism predicts that
stratification enters through the raw loading `a_j`, not uniformly. If so,
`χ²_j ~ 1 + ℓ_j + a_j` should put stratification into the loading coefficient
and leave an intercept of 1/λ. The rule was stated before the run: PASS if,
with GC (λ = 1.2) and stratification together, the three-term intercept lands
within 0.02 of 0.833 while the two-term intercept misses by more than 0.05.

The toy is an admixed cohort with N = 1200, M ≈ 2000, 1 cM windows, h² = 0.5,
no PC adjustment, and an in-sample (oracle) loading. The stratification shift
adds 0.2 to mean χ². Both magnitudes are illustrative, not calibrated to a
real GWAS.

**Table 3. Probe outcome.** Values are means with Monte Carlo SE in
parentheses. The "default run" column is the rule-bearing run (4 realizations
× 100 phenotypes). The "diagnostic rerun" column is the higher-precision rerun
(12 × 500), labelled as diagnostic because it was run after the failure. Its
last row comes from `--homogeneous-check` (8 × 300).

| quantity | default run | diagnostic rerun | reading |
|---|---|---|---|
| GC+STRAT three-term intercept − 1/λ | −0.044 (0.023) | −0.036 (0.007) | **FAIL**: misses the 0.02 tolerance |
| GC+STRAT two-term intercept − 1/λ | +1.311 | +1.445 | the unadjusted admixed GWAS puts structure-mediated polygenic signal into the intercept |
| paired shift from stratification, two-term intercept | +0.266 (0.070) | +0.515 (0.021) | stratification lands in the intercept |
| paired shift from stratification, three-term intercept | −0.002 (0.011) | −0.008 (0.003) | stratification is absorbed by the loading term |
| structure-free cohort (fst → 0), plain two-term intercept | — | 0.956 (0.013) | the toy LD model's own LDSC intercept bias |

What this says:

1. **The mechanism holds.** In paired comparisons, stratification moves the
   two-term intercept by +0.27 to +0.52, and the three-term intercept by
   less than 0.01.
2. **The rule fails, and the verdict stands.** The −4% miss is not produced by
   the raw-loading term. Plain LDSC shows the same bias in a cohort with no
   structure at all (0.956), so it is baseline LDSC bias in this toy's LD
   model. The design choice was poor: the tolerance should have been set
   against the toy's structure-free intercept, not against 1, and at 3 MC SE.
3. **Consequence for the strategy.** I am not recommending the ppb use. To come
   back in, it needs a re-stated rule (reference = the structure-free
   intercept) or, better, a real-data check: a GC-corrected GIANT height GWAS
   against a non-GC GWAS of the same trait, using 1000G loadings. A toy cannot
   settle ppb's question, because ppb's problem is about external loadings on
   real sumstats.

## 4. Landscape check (2026-09-24)

A single-agent search pass (about 33 queries) was run on the five open
questions the tracks depend on. Most full-text hosts were blocked, so the
**checked** column matters. "Code read" means the agent opened the source.
"Snippet" means it saw only search-result text: verify before citing.

**Table 4. What the check found.**

| question | nearest work found | checked | status |
|---|---|---|---|
| Q1: does admixture LD survive global-PC adjustment inside ±1 cM? | cov-LDSC (Luo 2021) recommends a **20 cM** window for admixed samples because of residual admixture LD (README read). HAMSTA (Chan 2023, AJHG) defines a local-ancestry "LAD score" for admixture-mapping h². Huang 2025 (Genetics) reports LDSC bias in admixed samples from directional LD | code read / snippet | **partial**. Nobody measured the residual inside ±1 cM. Local-ancestry-aware LD scores for LDSC are **open** (but that is a construction, §6) |
| Q2: do relatives in the LD panel bias LD scores or the intercept? | Mangin 2012 (kinship-corrected r², not LDSC). Pan-UKB's LD code uses `r2_adj = ((n−1)/(n−2))·r2 − 1/(n−2)` with raw n and no effective-N or covariate-df term (code read). The NYGC 3202-sample 1000G set adds 698 relatives, completing 602 trios, to the 2504 unrelateds (snippet) | code read / snippet | **open**: the most open of the five |
| Q3: false S-LDSC τ from ancestry-correlated annotations? | Gazal 2017 notes stratification risk only for signed analyses. Bhatia 2016 (bioRxiv) "PC loading regression" corrects z by reference loadings, which is prior art for the raw-loading *control*, but applied to z rather than as an S-LDSC regressor. Alexander & Curtis 2020: LD scores correlate with between-population frequency differences | snippet | **mostly open** for the false-τ demonstration; **partial** for the control |
| Q4: LDSC on multi-ancestry meta-analyses | GBMI (Zhou 2022, Cell Genomics) ran LDSC on its multi-ancestry meta with UKB-EUR LD scores, justified by an LDSC ratio not significantly above 0.2. A 2026 bioRxiv reports LDSC h² underestimation on meta-analysed Alzheimer's sumstats without proposing a correction | snippet | **partial**: no composition-matched LD scores, and no panel-only bias predictor found |
| Q5: per-axis stratification from summary statistics alone | Blanc & Berg 2025 (Genetics): one axis per polygenic-score test, needs GWAS-panel genotypes. Bhatia 2016 remains the only summary-statistics-only per-axis method found | snippet | **partial**. The fingerprint (A2) is a QC feature, not a paper claim |

This was one pass, not the adversarial multi-agent sweep this project
normally runs. It ranks the tracks; it does not clear any of them. Step 4 of
Table 6 is that check, and it is scoped to whichever claims survive their
gates.

## 5. The strategy: three tracks and a wildcard

**Table 5. Tracks.**

| track | deliverable | primary reader | why relevant | why cool | cost | gate | prior-art risk |
|---|---|---|---|---|---|---|---|
| **A. Instrument** | `pldsc audit-panel` (A1), `pldsc fingerprint` (A2), annotation flags (A3); JSON card + loadings sidecar with an identity hash | the family: SMARTpred, LDpred3-atw, gwfm | Consumers C1, C2, C4–C6 exist today | "Why is my intercept 1.15?" gets a measured answer on the results page | ~2 weeks | none; it is infrastructure the other tracks use | low (engineering on known ideas; no novelty claimed) |
| **B. Sharp result** | "Do functional PGS priors learn ancestry?" D1(b) at scale, then propagation into LDpred3-funct / gwfm | statistical-genetics methods community | Functional priors from baseline-LD are standard (SBayesRC, LDpred-funct, and this family) | A two-line mechanism (the self term) reaches a production PGS pipeline. The toy number exists: standard S-LDSC gives τ = 25.4 under pure confounding versus 2.77 under a genuine 4× enrichment | 1–2 weeks for the gate; ~4 weeks to write | B1 as already pre-stated in [D1_ANNOTATION_GATE](D1_ANNOTATION_GATE.md) | medium (Q3: the false-τ demonstration is mostly open; the control has Bhatia 2016 precedent) |
| **C. Measurement paper** | "What an LD reference costs you in the multi-ancestry era": the D3 audit with casualties | GWAS consortia and applied users (GBMI, MVP, All of Us, Pan-UKB) | Multi-ancestry meta-analyses are the growth area, and they run LDSC with EUR panels (Q4) | A published h² or r_g that moves, with the panel-only predictor that anticipated it | 3–5 weeks | C3: at least one published estimate moves by > 2 SE, otherwise it becomes a practice note | medium (Q4 partial; Q2 open for the relatives section) |
| W. Wildcard | A short perspective: pre-registered kill-tests and adversarial prior-art sweeps as a research method, with this repo as the dataset | meta-science / methods readers | Agent-assisted research needs examples of *stopping* | The dated record has four sweeps, about 15 vetted ideas, pre-registered rules, and at least three self-corrected overreaches (in-chromosome GRM, 17% → 29%, D2 "CLOSED") | ~1 week | only if it says something quantitative (verdict flips, cost per verdict) | low novelty risk; high audience risk |

### 5.1 Track A — the instrument

- **A1 panel card.** Input: PLINK genotypes, a genetic map, optional kinship.
  Output: n, ν, N_eff = tr(Ã)²/tr(Ã²) and N_eff/ν; λ₁…λ₁₀; the structure
  share at ±1 cM (LOCO) and genome-wide, with the window-unit remark; a
  forecast of the intercept shift at χ² slopes 0.05 and 0.3 (the D2
  propagation formula, slope × mean bias); and a loadings sidecar `a_jp`,
  p = 1…K. Almost all of this already exists in `pldsc/genome.py`,
  `pldsc.effective_pair_df`, `realdata_genome_loco.py` and
  `d2_intercept_propagation.py`. What is new is the CLI, the I/O and the JSON
  schema.
- **A2 stratification fingerprint.** `χ²_j ~ 1 + ℓ_j + Σ_p a_jp` with block
  jackknife, reporting a z-score for each axis. This is Bhatia 2016's idea
  deployed as QC, and the doc and UI must say so. The limit (METHOD §14)
  belongs on the card: signal on an axis can be confounding, genuine
  ancestry-aligned architecture, or, in an unadjusted GWAS, structure-mediated
  polygenic LD. §3.1 shows that last one dominating. The fingerprint flags;
  it does not certify.
- **A3 annotation flags.** The D1(a) computation, generalized: for every
  column of a functional panel, max |corr| with the loadings and R² on K
  loadings, with LLD_AFR and allele age as the known positives. Whether to add
  raw loadings as nuisance columns in funct fits waits on Track B.
- **A4 (optional, outreach).** An interactive "anatomy of an LD score" page.
  Sliders for panel composition and window width would show the linkage share
  against the structure share, built from A1 cards. It is cheap once A1
  exists, and it is the most shareable thing in this list.

**An honest alternative to Track A.** If maintaining another package is
unattractive, A1–A3 belong in `ldpred3`, which already owns `ldsc.py` and the
I/O, and PLDSC stays the research record. The most relevant form of PLDSC may
be three functions in ldpred3 plus one paper.

### 5.2 Track B — "do functional priors learn ancestry?"

- **B1 = D1 kill-test (b), unchanged.** Weighted S-LDSC with block-jackknife
  SEs and M ≥ 10⁵, on real 1000G genotypes with the real LLD_AFR and
  allele-age annotations. Compare τ_CONF with and without the raw-loading
  control. The rule is the one already written in D1_ANNOTATION_GATE: if
  τ_CONF ≈ 0 without the control, D1 is a footnote.
- **B2, only if B1 survives.** Run the same phenotypes through LDpred3-funct
  and gwfm. Does the learned coefficient on LLD_AFR or allele age move under
  confounding, and does the held-out PGS correlate more with PC1? The rule, to
  be stated before running: kill if |Δ corr(PGS, PC1)| < 0.02 at a
  stratification that adds 0.2 to mean χ² after the GWAS's own PC adjustment.
- **Why B before C.** It is the cheapest gate with the largest upside, its
  rule is already written, and a positive result is also a bug report against
  the family's own functional mode, which is useful whatever it does for a
  paper.

### 5.3 Track C — "what an LD reference costs you"

- **C1.** The 22-autosome LOCO structure shares at ±1 cM for the 1000G
  superpopulations and the admixed subgroups (ASW, ACB, MXL, PUR, CLM, PEL).
  This is already scripted (`realdata_genome_loco.py --plan`).
- **C2.** A panel-only predictor of the slope error from a meta-analysis's
  ancestry composition. The applied LD score carries the *panel's* non-local
  component, not the cohort's (PROJECT_DIRECTIONS D3).
- **C3, the casualty.** Re-run one GBMI trait, or a Pan-UKB multi-ancestry
  meta, with composition-matched LD scores instead of UKB-EUR, and report what
  moves and whether C2 predicted it.
- **C4, relatives in the panel (Q2, the most open).** Compare 2504 with 3202
  1000G. Check Pan-UKB's raw-n correction for the small ancestry groups (AMR,
  MID, EAS), after confirming what its "hq samples" list excludes. D2 already
  says the size of this effect: ~0.005 of intercept at typical slope and ~0.03
  at biobank slope. It is one section, not the headline.
- **C5, optional.** Residual admixture LD inside ±1 cM after global PCs (Q1).
  This is a measurement, which is the kind of claim that survives. The
  local-ancestry LD-score *construction* that it suggests is not recommended
  (§6).
- The decision-rule table from WHEN_STRUCTURE_MATTERS §5 closes the paper. The
  gate is C3: a benchmark with no casualty is a practice note, and that is a
  fine fallback.

## 6. What not to do

- **Do not look for "cool" in a new estimator or construction.** Four sweeps
  graded every construction partial or occupied. The local-ancestry-aware LD
  score (Q1) is open today, but it is exactly the kind of claim this project
  has watched die, and Huang 2025 and HAMSTA are close.
- **Do not reopen `P_c` as a contribution.** Keep it as the comparator that
  generates the B2 operator bracket.
- **Do not build a UI inside PLDSC.** SMARTpred is the UI. PLDSC ships cards
  and sidecars.
- **Do not claim that the fingerprint or the flags separate confounding from
  real ancestry-aligned signal.** §14 says they cannot, and §3.1 shows the
  structure-mediated polygenic term dominating an unadjusted GWAS.
- **Do not publish the WHEN_STRUCTURE_MATTERS note on its own.** Without
  consumers and a casualty it invites the question "so this is a cov-LDSC
  benchmark?", and that note says so itself.

## 7. Repositioning the repository (half a day, no science)

1. **README, first screen.** State what the instrument measures and for whom
   (Table 2). Status line: "estimator: concluded (2026-07-19); instrument:
   active." Move the archive paragraph directly beneath it, unedited.
2. **`ldpred3/docs/family.md`.** Once A1 exists, move PLDSC from "concluded"
   to a leaf that SMARTpred and LDpred3-atw consume. That is a change to a
   sibling repository, proposed here and not made.
3. **Keep the name.** A rename breaks links for little gain. The subtitle
   carries the new role.

## 8. Order of work

**Table 6. The first six to eight weeks.**

| # | task | cost | decides |
|---|---|---|---|
| 1 | README repositioning (§7.1) | ½ day | — |
| 2 | Track B1 (D1 kill-test b) | 1–2 weeks | whether Track B exists |
| 3 | A1 panel card on 1000G phase 3: superpopulations, FULL, 2504 vs 3202 (in parallel with 2) | 1 week | nothing; every track uses it |
| 4 | Scoped prior-art sweep on whichever of B and C survived, including verification of Table 4's snippet-level citations | 1 day | write or don't |
| 5 | A2 + A3 into SMARTpred QC behind a flag (a SMARTpred PR) | 1 week | whether you want the service to show it |
| 6 | Track C1–C3 with one real multi-ancestry meta-analysis | 2–3 weeks | paper or practice note |
| 7 | W and A4, only if wanted | ~1 week each | — |

**Table 7. Risks.**

| risk | my estimate | consequence | mitigation |
|---|---|---|---|
| B1 comes back small (the real annotation–loading correlations are only ≤ 0.12) | about even odds | Track B becomes a footnote | A3 still ships as a QC flag; C is independent of B |
| The family does not want another dependency | moderate | A has no users | Fold A1–A3 into ldpred3 (§5.1) |
| C's differentiator is already published (LD-reference mismatch is an active area) | moderate | C shrinks to a practice note | Step 4 before writing; the panel-only predictor plus a casualty is the differentiator |
| The toy-to-real gap: every mechanism number is from N ≤ 2400, M ≈ 2000 sims or chr21/22 | certain today | numbers move at scale | B1 and C1 are exactly the at-scale versions, and both have rules |
| The identifiability limit is read as a caveat rather than a design constraint | moderate | overclaiming | Put "flags, not certificates" on the card itself (A2) |

## 9. What I would tell you plainly

- PLDSC will not become cool by becoming a better estimator. It can become
  relevant by being *used*, and the family already has the uses.
- The strongest single bet is Track B's gate, because its rule is already
  written and a positive answer would matter to your own functional-PGS mode.
  The strongest paper is Track C, and only with a casualty.
- The one idea I tested for this review (C3, helping ppb) failed its
  pre-stated rule. The mechanism held, but the rule was badly posed, and I am
  not recommending the idea on that evidence.
- It would be reasonable to decide that the most relevant form of this work is
  three functions in ldpred3 and one measurement paper, with this repository
  kept as the honest record behind them.
