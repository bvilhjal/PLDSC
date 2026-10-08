# PLDSC

PLDSC studies LD Score regression when GWAS samples contain population
structure. The live project (since 2026-09-25) asks one question: **what
goes wrong when LDSC-based tools run on multi-ancestry data, how big is
the error, and what is the fix?** — plan and gates in
[docs/STRATEGY.md](docs/STRATEGY.md), math in
[docs/MULTI_ANCESTRY.md](docs/MULTI_ANCESTRY.md), how to read this repo in
[docs/READING_ORDER.md](docs/READING_ORDER.md).

## The headline so far (1000G chr22)

Take a fixed-effect meta-analysis that is 75% EUR and 25% one other 1000G
superpopulation, with trans-ancestry genetic correlation ρ = 0.8. If LDSC
regresses it on the standard 1000G **EUR** LD score, it understates h² by:

| other 25% | h² understated by | κ = ĥ²/h² |
|---|---|---|
| AFR | 25% | 0.750 |
| EAS | 16% | 0.841 |
| SAS | 13% | 0.868 |
| AMR | 9–10% | 0.907 / 0.896 |

**Genome-wide (2026-10-08):** the step-1 gate passes. AFR κ is 0.737 on
1.19M regression SNPs over 22 autosomes, and every cell is within 0.025 of
the table above (MULTI_ANCESTRY Table 6). That run used a different variant
set, and its bridge to the chr22 protocol missed the pre-stated tolerance,
so the table above stays the chr22 one for now. The first real casualty
test was GBMI asthma (Table 7). There, h² moves 0.133 → 0.156 when the
EUR score is replaced by the meta score. That is 1.9 SE, just short of the
pre-stated 2 SE.

Most of the loss is **LD mismatch, not effect heterogeneity**: even at
ρ = 1, κ = 0.79 for EUR + AFR. The intercept moves little (+0.00 to +0.03
of the LD-driven excess χ² at 75% EUR). Two honest negatives from the same
run: a meta cannot separate LD mismatch from effect heterogeneity
(MULTI_ANCESTRY Table 5), and the first κ computation overstated the loss
— the pre-stated verdict and the corrected numbers are both recorded
(Table 4). κ needs no GWAS: it is a property of the panels, the ancestry
mix and ρ (MA6).

## Run it in 60 seconds

No download, no genotypes — reproduce the headline tables from the
committed 2 MB LD-score bundle ([examples/mini](examples/mini)):

```bash
python examples/reproduce_chr22_tables.py     # Tables 3-5 in ~0.6 s
```

The Monte Carlo validation of the machinery (no data, seeded, ~20-35 s) —
**pass `--n-ref 6000`**:

```bash
python experiments/multiancestry_meta_sim.py --n-ref 6000   # 3x PASS
```

The no-flag default is the 400-person-panel run that *records the
pre-registered failure* in MULTI_ANCESTRY §6: it prints (i) FAIL (ii) FAIL
by design. Running it once to see a recorded failure is instructive; it
is not a bug.

The real-data chr22 experiment itself (needs `.scratch_1000g/`, see
below) then runs in seconds:

```bash
python experiments/multiancestry_1000g_meta.py    # Tables 2-5, ~4 s
```

## Install and environment

- Python ≥ 3.11, numpy ≥ 2.0. Editable install with the test extras:
  `python -m pip install -e ".[test]"` (also `.[plot]` for
  `pldsc.plotting`).
- **There is no `pldsc` conda environment.** Older experiment docstrings
  say `conda run -n pldsc ...` — that env never existed on this machine;
  use any 3.11+ interpreter. On the maintainer's machine the package is
  installed editable in the shared `ldpred3-accelerate` conda env
  (Python 3.14). The ambient base anaconda (3.10) is below the floor.
- Checks (see [CONTRIBUTING.md](CONTRIBUTING.md)):
  `ruff check .` and, with BLAS threads pinned,
  `pytest -q` or `python -m unittest discover -s tests -v` (81 tests,
  <1 s).
- Real-data experiments read `.scratch_1000g/` (gitignored, 100%
  data-free fresh clone). Build it with
  `scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22`:
  S3 mirror primary (the EBI FTP has been unreachable from here and serves
  a different re-release, v5b), sha256-verified against the copies that
  produced the recorded results, both genetic-map formats included.

## First numerical slice

Local LD scores with fixed covariates projected out:

```python
from pldsc import local_ld_scores

scores = local_ld_scores(
    genotypes,
    positions_cm,
    covariates=principal_components,
    window_cm=20.0,  # illustrative; select by an LD-score convergence check
)
```

It operates on one chromosome at a time, includes the LD Score self-term,
and retains negative finite-sample-adjusted terms. Covariates are treated
as fixed; the appropriate window is population- and data-dependent (20 cM
is not a universal default).

Multi-ancestry LD scores for a meta-analysis, and what a single-ancestry
score costs it (docs/MULTI_ANCESTRY.md):

```python
from pldsc import (cross_ancestry_ld_scores, meta_weights, meta_ld_scores,
                   single_reference_bias)

s = cross_ancestry_ld_scores([eur, afr], positions_cm, window_cm=1.0)
omega, n_j = meta_weights([n_eur, n_afr], s.available, scheme="sample_size")
L = meta_ld_scores(s, omega, rho=[[1, 0.8], [0.8, 1]])
kappa = single_reference_bias(L, s.scores[0, 0], sample_sizes=n_j).kappa
```

The panels must share variants and the counted allele. Within-ancestry
scores get the usual noise-floor correction; cross-ancestry scores need
none (MA8).

## History and the archived record

PLDSC is a clean restart of a 2017 research prototype. Its first line —
an estimator built on the mixed-model residual-precision operator `P_c` —
was **concluded and archived on 2026-07-19**: the machinery was correct
and well-tested, but four adversarial prior-art sweeps and an exhaustive
contribution search found no differentiated contribution, so the repo
keeps an unusually candid negative-results record instead of a paper.

That record is deliberate and worth reading — it is the rigour culture
the live project runs on (pre-registered kill rules that were honoured
when they failed, post-rule robustness checks that broke the author's own
headline, corrections appended above verdicts). Start from
[docs/READING_ORDER.md](docs/READING_ORDER.md), which classifies every
document as live or archival. The short map:

- **Archive** (do not treat as backlog or current scope):
  [docs/METHOD.md](docs/METHOD.md) (§10.1 and §14 are live),
  [docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md),
  [docs/PRIOR_ART.md](docs/PRIOR_ART.md),
  [docs/EXTERNAL_REVIEW.md](docs/EXTERNAL_REVIEW.md),
  [docs/IDEA_REVIEW.md](docs/IDEA_REVIEW.md),
  [docs/NOTE_OUTLINE.md](docs/NOTE_OUTLINE.md),
  [docs/GEMMA_PARITY.md](docs/GEMMA_PARITY.md),
  [docs/WHEN_STRUCTURE_MATTERS.md](docs/WHEN_STRUCTURE_MATTERS.md)
  (source material for MULTI_ANCESTRY §5).
- **Live, part of the archive-to-now arc:**
  [docs/PROJECT_DIRECTIONS.md](docs/PROJECT_DIRECTIONS.md) (the
  construction/measurement asymmetry),
  [docs/REAL_DATA_SURVEY.md](docs/REAL_DATA_SURVEY.md),
  [docs/D1_ANNOTATION_GATE.md](docs/D1_ANNOTATION_GATE.md)
  (Part 2's pre-written kill rule), [experiments/README.md](experiments/README.md)
  (dated index; the two `multiancestry_*` entries are current).

Durable value of the archived line, kept and tested: the `P_c` operator
and LD-score decomposition (`pldsc/`), the verified theory write-up
(METHOD §1–§15, including the adversarially checked structure-aliasing
remark §14), and the negative-results record itself.
