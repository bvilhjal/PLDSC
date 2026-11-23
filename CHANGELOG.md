# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions are
single-sourced from `pldsc.__version__`; the API may still change before
version 1.0.

## [Unreleased]

## [0.1.3] - 2026-11-23

### Added

- Report section 11 and Figure 4 on Cochran's Q: a pooled statistic with its
  Q recovers cross-ancestry genetic covariance, and within-ancestry Q gives
  a ceiling for a cross-ancestry correlation.
- `experiments/q_heterogeneity_sim.py`: known-truth simulation on 1000
  Genomes block LD with three stated rules, and its recorded results.
- `experiments/gbmi_asthma_heterogeneity.py`: the exploratory GBMI asthma
  analysis of the same two questions, and its recorded estimates.
- `ldsc_compat=True` for `ldsc_h2` and `ldsc_rg`, `ldsc_h2_compat`,
  `ldsc_rg_compat` and `pldsc h2/rg --ldsc-compat`: an opt-in estimator that
  follows the reference `ldsc` as `ldsc --h2` and `ldsc --rg` run by default
  (h2 two-step at chi-square 30, r_g one-step; fixed-weight jackknife, 200
  blocks, separate `w_ld` via `--w-ld-scores`). It agrees
  with `bulik/ldsc` commit 2fdeeb3 to about 1e-15 on the committed fixture
  `tests/data/ldsc_parity.npz` (`scripts/make_ldsc_parity_fixture.py`).
  Defaults are unchanged.
- `experiments/q_model.py`: the Z-and-Q design with three stratum intercepts
  (old shared-intercept form kept), a simulation arm with unequal intercepts
  and a test showing that the shared form leaks -9.4% into the covariance
  under inverse-variance weights on the chr22 bundle.
- Command line: reads ldsc `munge_sumstats` output (SNP, A1, A2, N, Z);
  refuses fits on fewer than `--min-snps` (default 200,000) regression SNPs
  unless `--allow-low-coverage`; checks declared `--sumstats-build` and
  `--reference-build` against each other and against the GRCh37 region
  filter; checks `--out` before reading data. JSON output adds `coverage`,
  `m_snps` (value and convention) and `build_check`; r_g output adds
  `gcov_se` and `h2_se`.
- A `UserWarning` when the array functions default `m_snps` to the row count.

### Changed

- Report derivations, regression details and the three Q propositions and
  proofs now appear in Supplementary Methods S1–S7. The main text retains
  ten key equations with intuitive explanations; the PDF and LaTeX source
  include the supplement.
- Report section 11 now states what Q discards, with Propositions 1-3 and
  proofs; the remaining excess-loss statement is marked as a conjecture.
  The earlier wording that Q "loses only the sign" and "identifies g" is
  replaced.
- The asthma inverse-variance Z-and-Q fit was rerun with three stratum
  intercepts: 0.1008 (0.0187), still 19% below the two-strata estimate
  (1.6 SE). The intercept leak moved g by only +0.003, and no cause is
  attributed; the allele-frequency explanation is withdrawn.
  `report/asthma_heterogeneity.json` adds the three-intercept fits; every
  earlier value is unchanged.
- The Q simulation was rerun with the three-intercept fit in every cell and
  the arm with stratum intercepts 1.08 and 1.0 under a new rule R4, stated
  in the driver before the run. R4 passed: the three-intercept fit averaged
  0.2216 against g = 0.222 (−0.2%), while one intercept per equation gave
  0.2098 (−5.5%). At 85% EUR the extra intercepts raise the SD of g from 1.7
  to 2.0 times that of the two strata. Every earlier value is unchanged.
- README: the multi-ancestry section leads with the null prediction result
  and the Q result, defines the estimand, and moves the 0.133 to 0.155 shift
  to supporting results. CLI help and docs carry an ancestry-mix caveat.
- Wording: removed claims that the default estimator matches `ldsc` or is
  faster than it ("Fast" dropped from the package description); documented
  the M convention (bundled 9,997,231 versus `M_5_50` 5,961,159, a factor of
  1.677 in h2) in `docs/LDSC.md`, with parity and timing tables.

### Fixed

- `ldsc_rg_compat` and `--ldsc-compat` r_g now default to one-step fits.
  Reference `ldsc --rg` never applies its two-step default, because
  `estimate_rg` tests `intercept_h2 is None` after turning it into a list;
  the earlier default reproduced `ldsc --rg --two-step 30` instead (r_g
  0.3284 against 0.3924 on the parity fixture). A two-step cutoff with a
  fixed cross-trait intercept now raises an error, as in ldsc; the docs no
  longer say that ldsc raises one by default.
- `tests/test_ldsc_port.py::test_provider_parity` and
  `experiments/ldsc_provider_parity_bench.py` failed with
  `KeyError: 'gcov_se'` when bipred is installed; they now compare the fields
  both results have.
- The report PDF and LaTeX source are rebuilt from the edited Markdown, and
  `report/reproduction.txt` records the current Section 11 drivers and
  results.
- The README and `docs/LDSC.md` said that PLDSC never detects builds or
  lifts over; with allelign 0.3 or later, alignment does both. A declared
  build mismatch now also names that route. The parse log of an ldsc
  munged file names its N as the total sample size, and the docs give the
  consequence for case-control traits.
- Citations: Popcorn authors in `CITATION.cff`; PRIOR_ART DOIs (Lin, Xing &
  Xing 2012; Mathew, Léon & Sillanpää 2018); the Bhatia et al. title; Xue &
  Zhao (2025, Annals of Statistics) cited as direct prior art in the report,
  `docs/MULTI_ANCESTRY.md`, `PRIOR_ART.md` and `CITATION.cff`.

## [0.1.2] - 2026-11-19

### Added

- Unified [methods and evidence report](report/methods.md) with the
  simulation results and three reproducible figures on chr22 tagging,
  held-out chromosome prediction for GBMI asthma and simulation power.
- Outcome test of the matched LD-score pattern, with aggregate chromosome
  cross-products for reanalysis and a two-population known-truth experiment.
- `cross_ancestry_ld_scores(..., effect_scales=...)`: per-panel scales of
  the standardized causal effects, giving weighted scores
  `sum_k v_ak v_bk r^a_jk r^b_jk` for effect models other than A1.
  The default (`None`) leaves the scores unchanged.
- `experiments/multiancestry_alpha_kappa.py`: kappa under shared per-allele
  effects with frequency-dependent variance, with pre-stated rules.

## [0.1.1] - 2026-10-29

### Changed

- First release with a non-development version number; the code is that of
  0.1.0.dev5.

## [0.1.0.dev5] - 2026-10-29

### Added

- Bundled `hm3-eur` LD reference: 1,190,268 HapMap3 regression SNPs with
  1000 Genomes Phase 3 EUR LD scores (GRCh37), loadable with
  `load_ld_reference`, describable with `pldsc ld-reference`, and usable
  in `pldsc h2`/`pldsc rg` via `--ld-reference` instead of explicit
  reference, score and `m_snps` inputs. The bundle is built from the
  public LDSC reference deposit by `scripts/build_ld_reference.py` and
  carries its provenance, including both published variant counts; the
  loader checks the bundle's recorded checksum. With `--ld-reference`,
  `--m-snps` replaces the bundled variant count (for example with the
  MAF 5–50% count).

### Changed

- `pldsc h2` and `pldsc rg` leave out reference variants in allelign's
  GRCh37 long-range LD regions, including the MHC, by default;
  `--include-long-range-ld-grch37` keeps them, replacing the opt-in
  `--exclude-long-range-ld-grch37`. The regions are found from reference
  coordinates, so GWAS files matched by ID no longer need positions, and
  the JSON `long_range_ld` entry is a single record for the reference.
  `load_ld_reference` leaves the same regions out unless
  `include_long_range_ld=True`.

## [0.1.0.dev4] - 2026-10-29

### Changed

- Simplify installation and the LDSC and research documentation for the public
  repository; preserve the detailed experiment results in a separate record.
- Remove obsolete environment-specific commands and the private-dependency
  token fallback from CI. Scientific computations and results are unchanged.
- Use current SPDX license metadata for package builds.

## [0.1.0.dev3] - 2026-10-29

### Added

- Command-line h² and bivariate LDSC with allelign processing, optional
  GRCh37 region exclusions, JSON diagnostics, and checked conversion of
  standard LDSC text scores to reference order.
- An executable notebook and a command-line walkthrough with synthetic data.

## [0.1.0.dev2] - 2026-10-29

### Changed

- Route the historical GBMI asthma GWAS reader and allele alignment through
  allelign, preserving its pre-specified selection and sample-size rules.
- Expand file-workflow tests and explain allelign column overrides and QC logs.

## [0.1.0.dev1] - 2026-10-29

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
