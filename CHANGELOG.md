# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions are
single-sourced from `pldsc.__version__`; the package is pre-release
(`0.1.0.dev1`) and the API may still change.

## [Unreleased]

### Added

- Standalone univariate and bivariate LDSC with optional Numba acceleration,
  dense and low-rank block scores, and allelign-backed GWAS QC/harmonization.
  The regression kernels are adapted from LDpred3 and bipred.
- A known-truth finite-reference score-noise sensitivity experiment.
- A methods report explaining univariate and bivariate LDSC, the file
  workflow, computation and the relationship to the multi-ancestry research.

- `experiments/multiancestry_genome_kappa.py`: STRATEGY step 1, κ on all 22
  autosomes (bigsnpr's 1000G phase-3 set). It runs the unchanged
  `cross_ancestry_ld_scores` on segments that contain every window in full,
  so it matches one whole-chromosome call to 1e-12 under 1 GB of memory.
  Rules were pre-stated. G1, the step-1 gate, passes: AFR κ 0.737 / 0.726.
  G2, the bridge to the thinned-VCF protocol, missed its 0.03 tolerance
  (0.039), so the README headline stays chr22's. Results are in
  `docs/MULTI_ANCESTRY.md` Table 6, with STRATEGY, STUDENT_PLAN, README and
  the experiments index updated to match.
- `experiments/gbmi_asthma_casualty.py`: STRATEGY step 3 on GBMI asthma.
  ρ comes from the ancestry metas (EUR–EAS 0.755; the smaller strata are
  bracketed), and h² of the all-biobank meta is fitted on the EUR score and
  on the meta score. Rules were pre-stated. The gate is not met robustly:
  h² moves 0.133 → 0.156, which is 1.9 SE at the central bracket. Results
  are in `docs/MULTI_ANCESTRY.md` Table 7. `.scratch_gbmi/` is gitignored.

- `examples/mini/chr22_ldscores.npz` (2.0 MB): the chr22 cross-ancestry
  LD-score tensor (5×5×9,654) plus masks, variances and per-ancestry MAFs,
  committed so MULTI_ANCESTRY Tables 3–5 reproduce in under a second with
  no download (`examples/reproduce_chr22_tables.py`; regenerate with
  `examples/make_mini.py`). Provenance and data terms in
  `examples/mini/README.md`.
- `scripts/fetch_1000g.py`: the documented path to a working
  `.scratch_1000g/` — S3 mirror primary (EBI fallback, with the v5a/v5b
  release decision made explicit), sha256 verification against the copies
  that produced the recorded results, both genetic-map formats (SHAPEIT4
  gmap + plink-format maps converted from it, because the Beagle host is
  gone), baselineLD, and the thinned-cache build.
- `pldsc.plotting` (optional `.[plot]` extra): `kappa_heatmap`,
  `kappa_profile` (with the step-1 gate drawn), and the before/after
  `forest_plot` for the casualty. matplotlib is imported lazily; the
  package core stays numpy-only.
- `CONTRIBUTING.md` (interpreter floor, editable install, the three checks
  CI runs, the "there is no `pldsc` conda env" note, the data section, and
  how contributions work in a repo with pre-registered gates);
  `docs/READING_ORDER.md` (the live/archive classification of all docs);
  `docs/STUDENT_PLAN.md` (the M0–M6 milestone ladder, sizing, do-not-do
  list, and the best-scoped unclaimed side questions); `CITATION.cff`.
- Packaging: `[project.optional-dependencies]` (`test`, `plot`), license,
  authors and repository URLs in the metadata; CI now installs
  `.[test]`, runs `ruff check .`, tests on Python 3.11/3.13/3.14 via
  pytest (was: unittest on 3.11/3.13, no lint); pytest config with
  `testpaths`.
- This changelog. (The 2026-09-02 review had recorded its absence as
  "deliberately omitted, archived"; the project now has a live direction
  and a second contributor coming.)

### Changed

- Put fast univariate and bivariate LDSC in the README quick start; added a
  practical LDSC guide and a runnable synthetic GWAS example.
- Updated citation metadata and included documentation, examples, scripts,
  research notes and report files in the source distribution.
- README reorganized the multi-ancestry question, headline κ table,
  reproduction commands, environment and API slice; the ~180-line
  archived-record section is collapsed into a history paragraph with links
  (every word survives in the linked docs).
- `docs/STRATEGY.md`: step 1 now prices the VCF→cache parser sub-task (the
  1–2 days covered only the download); step 3 upgraded from "GBMI or
  similar" to a point estimate — ρ estimated from GBMI's openly
  downloadable ancestry-stratified sumstats instead of bracketed, with
  per-trait manifest verification as the first sub-task; §5 landscape
  adds the now-citable 2026 meta-LDSC preprint (PMID 42465239), Momin et
  al. 2023 and cov-LDSC positioning; the ρ risk row's mitigation updated
  accordingly. `docs/MULTI_ANCESTRY.md` §8 mirrors the citation fix.
- `experiments/README.md`: dated index at the top with the live entries
  first (it previously opened with the dead line's "Verdict: CONTINUE"),
  plus the no-`pldsc`-env and data-fetch notes;
  `realdata_genome_loco.download_plan()` now points at
  `scripts/fetch_1000g.py` instead of EBI `v5b` URLs and a self-reference.
- Archive banners at the top of DEVELOPMENT_PLAN, PRIOR_ART,
  EXTERNAL_REVIEW, IDEA_REVIEW, NOTE_OUTLINE and GEMMA_PARITY
  ("historical record of the concluded estimator line"); stale claims
  annotated in place: DEVELOPMENT_PLAN §1.3/§1.4 (scope/decision rule
  superseded) and §5 ("no .github/workflows"), WHEN_STRUCTURE_MATTERS's
  status ("proposed" → stopped, folded into MULTI_ANCESTRY §5),
  IDEA_REVIEW Eq. (7) (superseded by METHOD §10.1), PROJECT_DIRECTIONS §5
  (the 22-autosome work order, superseded by STRATEGY step 1), and the
  test-count mentions now point at the canonical count.

### Fixed

- Allow CI to check out the public allelign reference with the default GitHub
  token when no repository-specific read token is configured.
- Stable centered WLS, rank-aware effective degrees of freedom, finite-score
  validation, and positive-semidefinite genetic-correlation matrices.
- Score caches bind map, data, and implementation hashes; old caches require
  an explicit legacy opt-in for the GBMI experiment. Refit estimated rho in
  its paired jackknife and replace the indefinite upper correlation bracket.
- Correct claims about finite-panel κ, correlated regressors and the scope of
  the pooled/meta identity.

- `experiments/m_spike_kill_or_continue.py`: the two simulated populations
  now share one ancestral spectrum per chromosome. Each had been drawn from
  its own, so the realized Hudson F_ST was 0.26 at nominal 0.08, and 0.21 in
  the "homogeneous" control, which is now one population. Readout B's MCSE
  now uses Var(z²), not Var(z), and the script prints the realized F_ST.
  The re-run keeps the CONTINUE verdict, with ordinary-OLS null inflation
  under structure of 2.25 (was 12.3) and ~1 for all three models in the
  control (ordinary had read 9.07). The experiments README, METHOD §11 and
  the DEVELOPMENT_PLAN status are updated.
- `pldsc.genome.read_plink_genetic_map` skips `#` comment lines; the
  plink-format maps written by `scripts/fetch_1000g.py` carry a provenance
  header and previously crashed the parser (regression test added).
- `experiments/README.md`: "maps and panel already local" was false for a
  fresh clone; corrected to the fetch-script command.
- Verified end-to-end after the data fixes (2026-09-25): all three
  previously data-blocked experiments run — `multiancestry_1000g_meta.py`
  reproduces Tables 2–5 cell-for-cell (~4 s), `realdata_1000g_benchmarks.py`
  completes B1–B6 (97 s), `realdata_genome_loco.py` completes and keeps the
  29 %/90 % headline on the converted map (38 s).

### Tests

- 81 → 86: `tests/test_plotting.py` (matplotlib-gated smoke tests for the
  three figures) and the map-comment regression test in
  `tests/test_genome.py`.
