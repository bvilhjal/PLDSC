# Reading order: what is live and what is archive

The standalone LDSC implementation and the multi-ancestry research are live.
The documents from the concluded estimator line are retained as a historical
record, rather than a backlog or plan. The line counts below are a historical
snapshot from 2026-09-25.

**For ordinary and bivariate LDSC**, start with [LDSC.md](LDSC.md).

**For multi-ancestry research**, start with [STRATEGY.md](STRATEGY.md) (15 min), then
[MULTI_ANCESTRY.md](MULTI_ANCESTRY.md) §1–§6 (the math and every current
table), then `experiments/README.md`'s live entries. The archived line is
worth reading once you know it is an archive — it is the rigour culture
the live project runs on.

## Live (current project)

**Table 1. Live documentation.**

| doc | lines | role |
|---|---|---|
| [LDSC.md](LDSC.md) | — | Installation, univariate and bivariate LDSC, QC, reference scores and examples |
| [Methods report](../report/methods.pdf) | — | LDSC equations, uncertainty, implementation and limitations |
| [STRATEGY.md](STRATEGY.md) | 123 | **Read first.** The one question, the headline, six next steps each with a pre-stated gate, landscape, risks |
| [MULTI_ANCESTRY.md](MULTI_ANCESTRY.md) | 299 | **Read second, fully.** MA1–MA12, Tables 1–5, §7 caveats (each unexplored item is a bounded task) |
| [D1_ANNOTATION_GATE.md](D1_ANNOTATION_GATE.md) | 88 | Part 2's pre-written kill rule (the S-LDSC false-enrichment test) |
| [PROJECT_DIRECTIONS.md](PROJECT_DIRECTIONS.md) | 313 | ~55% live: §1 the construction/measurement asymmetry, §3 what the repo owns, §4 D1, §8 gate outcomes; §4 D2–D4, §5, §6 superseded by STRATEGY |
| [REAL_DATA_SURVEY.md](REAL_DATA_SURVEY.md) | 118 | Mixed: the spectral engine and status blocks are live; Q2–Q4 unclaimed |
| [STUDENT_PLAN.md](STUDENT_PLAN.md) | — | The milestone ladder (M0–M6) if this is being taken on as a thesis/rotation project |
| `experiments/README.md` | 936 | ~30% live: read the two `multiancestry_*` entries and `realdata_genome_loco` / `raw_loading_control` / the D2 trio. **Do not read top-down** — it opens with the dead line's verdicts; use its dated index |

`examples/ordinary_ldsc.py` runs ordinary and bivariate LDSC on synthetic
files. For the multi-ancestry research, `examples/reproduce_chr22_tables.py`
reproduces the headline tables, and the Monte Carlo validation runs with
`python experiments/multiancestry_meta_sim.py --n-ref 6000`.

## Archive (historical record of the concluded estimator line)

Each of these carries a banner at the top of the file. §10.1 and §14 of
METHOD.md are the two live sections inside an otherwise-archived note
(the effective-pair-df correction, used by `pldsc.effective_pair_df`, and
the spectral identity `pldsc.genome` implements).

**Table 2. Archived documentation.**

| doc | lines | what it is |
|---|---|---|
| [METHOD.md](METHOD.md) | 1,203 | The archived theory note; §10.1 and §14 are live and consumed by `pldsc/` |
| [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) | 605 | The archived plan; its §1.3 scope excludes the current project — do not start from it |
| [PRIOR_ART.md](PRIOR_ART.md) | 203 | The do-not-re-open register (six adversarial sweeps); the model for STRATEGY step 2 |
| [EXTERNAL_REVIEW.md](EXTERNAL_REVIEW.md) | 661 | Independent re-derivation of the archived claims; §5.1/§6.1 reusable as method |
| [IDEA_REVIEW.md](IDEA_REVIEW.md) | 220 | Origin of the `P_c` line; its Eq. (7) noise floor is superseded by METHOD §10.1 |
| [NOTE_OUTLINE.md](NOTE_OUTLINE.md) | 122 | The clarifying note "was not pursued" |
| [GEMMA_PARITY.md](GEMMA_PARITY.md) | 114 | Parity check for the archived estimator |
| [WHEN_STRUCTURE_MATTERS.md](WHEN_STRUCTURE_MATTERS.md) | 145 | Stopped as a standalone note (STRATEGY §6); source material for MULTI_ANCESTRY §5's pooled case |

## One-paragraph orientation for newcomers

PLDSC provides standalone univariate and bivariate LDSC. Its current research
measures what an EUR LD score costs a multi-ancestry meta-analysis in h² (κ)
on 1000G panels, with a test on a published meta-analysis (STRATEGY §4 step
3). Everything dated 2026-07-17 … 2026-09-24 that talks about `P_c`, a
new estimator, a theory note, or a "PLDSC becomes a validation project"
decision rule is the concluded line: read it to learn why constructions
kept dying, not for things to do. The do-not-do list is STRATEGY §6; the
best-scoped unclaimed side questions are listed in
[STUDENT_PLAN.md](STUDENT_PLAN.md).
