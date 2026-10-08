# PLDSC as a student project — the milestone ladder

Date: 2026-09-25. What PLDSC looks like taken on as an **MSc thesis
(14–17 weeks / 30 ECTS)** or a **PhD rotation / 10-week project**, with a
ladder whose every rung has a deliverable that survives even if its gate
fails — the property that makes a research project safe to grade in a
repo whose own history teaches that many directions die. The supervisor
prep for this plan (README restructure, fetch script, the committed
example bundle, CONTRIBUTING, reading order) is done; if any of it turns
out broken, **stop and fix the repo, not the student** (see M0's gate).

## The ladder

| M | weeks | task | deliverable | gate | if the gate fails |
|---|---|---|---|---|---|
| **M0** | 1 | Environment (`CONTRIBUTING.md`), run the tests, run `experiments/multiancestry_meta_sim.py --n-ref 6000` (and once without, to *see* the recorded failure), reproduce the headline from `examples/mini`; read STRATEGY + MULTI_ANCESTRY §1–§6 | a one-page "what I reproduced" note + **their first figure** (`pldsc.plotting.kappa_heatmap` on the bundle) | sim reproduces Table 1; bundle reproduces Table 3 to printed precision | the repo is broken — fix it (that is M0's finding, not the student's failure) |
| **M1** | 2–4 | **κ genome-wide** (STRATEGY step 1; the gate passed on 2026-10-08 on bigsnpr's 1000G set, MULTI_ANCESTRY Table 6, so what is left here is the thinned-VCF bridge and 1b). 1a: the fast, progress-reporting, resumable VCF → cache parser (the current one is a nested pure-Python loop, chr21/22/10 only). 1b: the two unexplored assumption axes of MULTI_ANCESTRY §7 (η ≠ 1, α ≠ −1) as sensitivity runs | κ table v2 + the per-chromosome profile figure (`kappa_profile`) + a sensitivity section | the rule-cell AFR κ stays ≤ 0.90 | chr22 was unrepresentative → the claim narrows; the *measurement* is still the contribution |
| **M2** | 3–4 | **Prior-art sweep** (step 2): generalize `workflows/scripts/pldsc-prior-art-sweep.js` from the four archived claims to MA2/MA6/MA8/κ; position against PMID 42465239, Momin 2023, cov-LDSC; verify STRATEGY §5's snippet-level citations | a new PRIOR_ART section with per-claim occupied / partial / unoccupied verdicts | are MA2/κ unoccupied? | if occupied: the paper becomes the measurement + the casualty, which is M3 anyway |
| **M3** | 4–8 | **The casualty** (step 3, upgraded): GBMI's ancestry-stratified sumstats are open — estimate ρ instead of bracketing it (first: verify per-trait availability in the manifest), then recompute h² with the EUR score vs the meta score vs the published value | the before/after forest plot (`forest_plot`) over GBMI traits + the point-estimate κ | h² moves by more than 2 SE | publish as a practice note with the κ table — explicitly STRATEGY §7's fallback |
| **M4** | 9–11 | **Ship it**: `pldsc kappa --mix EUR=0.75,AFR=0.25 --rho 0.8` (κ needs no GWAS — MA6; this creates the package's first `[project.scripts]` entry point), then the SMARTpred "reference fit" hook and the `ldpred3.ldsc.ldsc_h2` correction path (step 6) | a released function + docs + a one-line change in a sibling that consumes it | does the family adopt it? | the CLI stands alone as a tool |
| **M5** | 11–15 | **Choose one:** Part 2's S-LDSC false-enrichment kill-test (rule already written in `D1_ANNOTATION_GATE.md`; `raw_loading_control.py` is compute-heavy, budget it), *or* r_g for two meta-analyses (step 4), *or* the local-adaptation power question at real M | a second, independent result | as pre-written | a negative result with a pre-stated rule is a **pass** in this repo's culture |
| **M6** | 15–17 | Write-up: the note (fold WHEN_STRUCTURE_MATTERS §4's two-regime intercept numbers in as Part 3) or a thesis chapter | draft | — | — |

**Sizing.** MSc thesis: M0–M4 plus one M5 option — a complete arc
(theory done, measurement, external validation on open data, a shipped
tool). Rotation / summer: M0–M2 (reproduce, genome-wide κ, novelty
verdict — ends with a figure, a table, and a written assessment).
Undergraduate final-year: weaker fit — the theory load in MULTI_ANCESTRY
§1–§5 is real (fixed-effect meta weights, Cauchy–Schwarz bounds,
errors-in-variables) and needs a supervisor walk-through of MA2/MA6.

**What a pass looks like, stated up front:** a reproduced genome-wide κ
table, a prior-art verdict, one figure set, and a working `pldsc kappa`.
**What a distinction looks like:** the GBMI casualty moving a published
h² by > 2 SE with an estimated ρ, adopted by a sibling package.

## Do-not-do list (protect yourself from this repo's own history)

1. **No new estimators or constructions** — STRATEGY §6 stops them; four
   prior sweeps closed every construction claim this project made. This
   includes the local-ancestry LD score.
2. **Do not revive `P_c`** as a contribution, and do not re-open the
   within-family/sib angle (`within_family_killtest.py` closed it).
3. **Do not re-litigate the archived theory note** (NOTE_OUTLINE,
   IDEA_REVIEW, METHOD §1–§9/§11–§13/§15). It is a record.
4. **Do not start from `DEVELOPMENT_PLAN.md`** — it is the archived plan
   and its §1.3 scope excludes the current project.
5. **Do not read `experiments/README.md` top-down** — its opening
   "Verdict: CONTINUE" belongs to the dead line; use its dated index.

The archive is the training, not the backlog: read PRIOR_ART and
MULTI_ANCESTRY §6 to learn how this repo runs a kill rule, then run yours
the same way.

## If you want a side question, the best-scoped unclaimed ones

- **The local-adaptation power question at real M.**
  `local_adaptation_channel.py` came back MARGINAL by its own pre-stated
  rule (power 0.38 at N = 600, κ = 1); whether it has power at M ~ 10⁶
  "is untested, and is the open question". Compute-bound, no new theory.
- **The per-PC inflation spectrum on real sumstats.**
  `detect_pc_inflation.py` works but was graded "largely Bhatia +
  incremental"; `REAL_DATA_SURVEY.md` Q2 (LD Hub / GWAS Catalog) is unrun.
  Sumstats-only — a clean one-semester audit.
- **An audit of shipped LD references for the noise-floor denominator.**
  STRATEGY §5 records two concrete targets: Pan-UKB's LD code using raw n
  in its r² correction, and the 602 trios inside the NYGC 3,202-sample
  1000G set. Small, sharp, and lands on tools thousands of people use.
