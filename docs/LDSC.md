# LDSC and bivariate LDSC

PLDSC estimates SNP heritability with univariate LD Score regression (LDSC)
and genetic correlation with bivariate LDSC. The default implementations
adapt the regression kernels of LDpred3 and bipred; they are **not** the
estimator of the reference `ldsc` software, and their defaults differ from
it ([LDSC equivalence](#ldsc-equivalence)). Optional Numba compilation
accelerates the delete-block jackknife; bivariate fits can also run its blocks
in parallel. Allelign reads, filters and harmonizes GWAS summary statistics.
The [unified report](../report/methods.md) ([PDF](../report/methods.pdf))
explains the regressions and their uncertainty, then tests how their
interpretation changes with study design and reference population.

## Install and try it

With Python 3.11 or newer, clone the repository and install PLDSC. Pip installs
the public allelign dependency automatically:

```bash
git clone https://github.com/bvilhjal/PLDSC.git
cd PLDSC
python -m pip install ".[fast]"
python examples/ordinary_ldsc.py
```

The example creates temporary reference and GWAS files, estimates h² and
r_g, and prints QC and allele-alignment counts. Its synthetic z scores are
independent across variants: it illustrates the API, not the calibration of
an analysis. `.[fast]` installs Numba; without it the NumPy implementation
remains available. Compilation takes extra time on the first call.

## Bundled HapMap3 reference (EUR, GRCh37)

For European-ancestry GWAS on GRCh37 coordinates, PLDSC bundles the
HapMap3 LD scores published for LDSC, so no score files are needed:

```python
from pldsc import load_ld_reference, ldsc_rg_sumstats

ref = load_ld_reference("hm3-eur")
rg = ldsc_rg_sumstats("trait1.tsv.gz", "trait2.tsv.gz",
                      ref.variants, ref.scores, m_snps=ref.m_snps, ncores=4)
```

```bash
pldsc rg --sumstats1 trait1.tsv.gz --sumstats2 trait2.tsv.gz \
  --ld-reference hm3-eur --ncores 4 --out rg.json
pldsc h2 --sumstats trait1.tsv.gz --ld-reference hm3-eur --m-snps 5961159
pldsc ld-reference   # provenance of the bundle, as JSON
```

The bundle holds 1,190,268 HapMap3 regression SNPs with LD scores from
the 1000 Genomes Phase 3 EUR QC panel (489 people, 1 cM windows), as
published for LD Score regression (Zenodo record 7768714; the build is
reproducible with `scripts/build_ld_reference.py`). Its `m_snps` is
9,997,231, the number of panel variants summed into each score. The
original LDSC software defaults instead to the MAF 5–50% count
(5,961,159 here, available as `ref.m_snps_5_50`); the choice rescales
h² and genetic covariance but leaves rg unchanged when the same count
is used in all three regressions, so state which count a reported h²
used. On the command line, `--m-snps` replaces the bundled count, as in
the second example. Fifty-three printed SNPs whose published score is
not positive are excluded from the regression set and listed in the
provenance. By default `load_ld_reference` and the commands also leave
out the 30,437 variants in allelign's GRCh37 long-range LD regions,
including the MHC, and record the count in `ref.long_range_ld` (the
`--m-snps` count is unchanged). LDSC runs on this deposit leave out the MHC
through their `weights_hm3_no_MHC` regression weights. To keep these
regions, pass `include_long_range_ld=True` or
`--include-long-range-ld-grch37`; `include_apoe=True` or `--include-apoe`
also leaves out the APOE region.
Do not use this reference for non-European or GRCh38 summary statistics;
supply matched scores instead, as below.

## Command-line walkthrough

Create a small set of synthetic files and keep them in a **new** directory:

```bash
python examples/ordinary_ldsc.py --out-dir demo-data
```

The directory contains `trait1.tsv.gz`, `trait2.tsv.gz`, `reference.bim`, and
`reference.ldscores.npy`. Here `M = 4000` because the synthetic scores were
summed over 4000 reference variants. Fit h² and then r_g:

```bash
pldsc h2 --sumstats demo-data/trait1.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --allow-low-coverage --out demo-data/h2.json

pldsc rg --sumstats1 demo-data/trait1.tsv.gz --sumstats2 demo-data/trait2.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --ncores 2 --allow-low-coverage --out demo-data/rg.json
```

Both commands call allelign to load and QC the GWAS, then harmonize alleles
against the reference. The output JSON records estimates, package versions,
input paths, settings, per-trait QC and alignment logs, the number of common
variants before the chi-square cap, and the number retained. It also records
`coverage` (regression SNPs, reference variants and the `--min-snps`
threshold), `m_snps` (the value used and its convention, e.g. whether it is
the bundled count or `M_5_50`) and `build_check`; r_g output includes the
genetic covariance SE (`gcov_se`) and both marginal h² SEs (`h2_se`). An
undefined estimate or SE is `null`. Read the full logs, for example:

```bash
python -m json.tool demo-data/rg.json
```

`--out` creates a new file and refuses to overwrite an existing one; the
path is checked before any data are read, so a bad `--out` fails at once.
Omit it to write JSON to standard output; warnings, such as an undeclared
reference build, go to standard error.

The commands fail closed (exit status 2) when fewer than 200,000 SNPs remain
for the regression, the count below which ldsc warns that estimates are
unreliable. Change the threshold with `--min-snps`, or add
`--allow-low-coverage` to report the fit anyway, as the 4,000-variant demo
above does. To see all
options, run `pldsc h2 --help` or `pldsc rg --help`. `ncores` parallelizes
bivariate jackknife fits when Numba is installed; it does not accelerate every
part of file loading.

### Use your own files

Provide a `.bim` or `.pvar` reference variant table, sorted by genomic
position within contiguous chromosome blocks; a one-dimensional `.npy` LD
score vector with **one score in exactly the same row order**; and one or two
GWAS files. Scores must be suitable for the GWAS populations and variant set.
For bivariate LDSC, the same reference scores must be suitable for both
traits. Make sure GWAS and reference use the same genome build. Declare them
with `--sumstats-build` (e.g. `GRCh38`, `hg19`) and, for an explicit
reference, `--reference-build` (the bundled reference is GRCh37); a declared
mismatch is an error. Undeclared builds are left to allelign: from version
0.3 its alignment detects a file on the other build when it matches less than
half of the reference, lifts it over (downloading a UCSC chain file once) and
records this under `liftover` in the alignment log. Older allelign versions
do neither, and such a file then fails the `--min-snps` check. PLDSC itself
does not infer reference suitability.

**Ancestry.** The scores must describe the population that produced the
GWAS. For a GWAS or meta-analysis that pools ancestries, a single-population
reference such as `hm3-eur` mis-scales h² and the genetic covariance, and the
quantity being estimated is a composition-weighted combination of the groups'
heritabilities and cross-ancestry covariances, not any one population's
heritability. See [MULTI_ANCESTRY.md](MULTI_ANCESTRY.md).

Set `--m-snps` to the number of variants summed into each LD score, **not**
the number of regression rows remaining after QC. For example, scores summed
over one million variants use `--m-snps 1000000` even if only 800,000 GWAS
rows survive. The command does not read the original LDSC text score files or
their `.M` files directly during fitting. To convert standard LDSC score
files containing `CHR`, `SNP`, `BP` and `L2`, use:

```bash
pldsc prepare-scores --reference reference.bim \
  --score-files ldscores/*.l2.ldscore.gz --out reference.ldscores.npy
```

The converter reads the reference with allelign, joins by unique variant ID,
checks chromosome and position, and writes scores in reference order. It
rejects missing or duplicated reference matches and non-positive scores;
score-file rows outside the reference are counted and ignored. For different
column names, see `pldsc prepare-scores --help`. Select a reference table
whose rows are the desired regression variants; filtering its rows does not
change `M`. Obtain `M` from the score construction or the matching LDSC `.M`
or `.M_5_50` files, making sure it counts the variants actually summed into
the scores. The converter does not infer `M` or verify genome builds.

Allelign recognizes many GWAS layouts. Inspect its detected columns before a
fit with `allelign inspect trait1.tsv.gz`. Output of ldsc's `munge_sumstats.py`
(`SNP`, `A1`, `A2`, `N`, `Z`, tab- or whitespace-delimited, optionally gzip
compressed) is also accepted: it is read as z scores with unit SE, rows with
a missing Z or N are dropped and counted, and the parse log reports the
format as `ldsc-munged`. Such files carry no positions, MAF or INFO, so the
MAF and INFO filters do not apply to them. Their `N` is ldsc's total sample
size and is used as the effective N. For a case-control GWAS h² is then on
the observed scale at the sample's case fraction φ, which is 4φ(1 − φ) times
the balanced-scale value an effective N gives; pass a constant effective N
with `--n-eff` (`--n-eff1`, `--n-eff2` for r_g) instead. For nonstandard column names,
pass mappings to PLDSC; repeat the option for each field:

```bash
pldsc h2 --sumstats trait1.tsv.gz --reference reference.bim \
  --ld-scores reference.ldscores.npy --m-snps 1000000 \
  --column beta=EFFECT --column se=STDERR --column n_eff=NEFF --out h2.json

pldsc rg --sumstats1 trait1.tsv.gz --sumstats2 trait2.tsv.gz \
  --reference reference.bim --ld-scores reference.ldscores.npy --m-snps 1000000 \
  --column1 n_eff=NEFF --column2 n_cases=CASES --column2 n_controls=CONTROLS \
  --ncores 4 --out rg.json
```

Use `--column` for h² and `--column1`/`--column2` independently for r_g.
Fields include `id`, `chrom`, `pos`, `ea`, `oa`, `beta`, `or`, `se`, `pval`,
`n_eff`, `n_cases`, `n_controls`, `eaf` and `info`. To supply a justified
constant effective sample size instead, use `--n-eff 100000` for h² or
`--n-eff1`/`--n-eff2` for r_g. A constant replaces **every** per-variant N.
The commands reject conflicting explicit N sources. Case/control-derived N
needs special care in meta-analyses, as explained below.

By default the commands use allelign's file QC, drop strand-ambiguous
variants, and apply the LDSC chi-square cap. They also leave out reference
variants in allelign's GRCh37 long-range LD regions, including the MHC; add
`--include-apoe` to leave out the APOE locus as well. The regions are found
from the reference coordinates, so GWAS files matched by ID need no
positions, and the JSON records the number of reference variants left out.
The default assumes a GRCh37 reference. A declared `--reference-build`
other than GRCh37 is an error unless the regions are kept, and an undeclared
build draws a warning. Add
`--include-long-range-ld-grch37` to keep the regions, for example for another
build; then filter with a build-appropriate region list using allelign in
Python if needed.

`--min-n-ratio`, `--min-maf`, and `--min-info` adjust both traits' file QC.
`--allow-ambiguous` and `--no-chi2-cap` change variant selection and should
only be used with an appropriate justification. `--constrain-intercept` fixes
the univariate intercept for h² and the **cross-trait** intercept for r_g.

## Notebook walkthrough

Install Jupyter in the same Python environment, then launch the executable
[notebook](../examples/ldsc_notebook.ipynb) from the repository root:

```bash
python -m pip install jupyterlab ipykernel
python -m ipykernel install --user --name pldsc --display-name "Python (PLDSC)"
python -m jupyterlab examples/ldsc_notebook.ipynb
```

Choose **Python (PLDSC)** as the kernel. Run all cells to create synthetic
files, load them with allelign, fit both models, inspect QC and alignment,
and demonstrate optional long-range LD exclusion. The notebook explains how
to replace the paths and `M` with real inputs.

## Fit GWAS files in Python

Suppose `reference.bim` is a sorted reference variant table and
`reference.ldscores.npy` is a NumPy vector with one LD score for each row of
that table, in the same order. The scores must be based on an appropriate
reference population and compatible variants. The GWAS files must use the
same genome build as the reference. These functions do not check builds;
allelign 0.3 or later lifts a file on the other build over during alignment
(see above). In this example the scores were summed over **exactly the
variants in the reference table**:

```python
import allelign
import numpy as np
from pldsc import ldsc_h2_sumstats, ldsc_rg_sumstats

reference = allelign.read_variants("reference.bim")
scores = np.load("reference.ldscores.npy", allow_pickle=False)
M = len(reference)  # valid here because this is the full score-sum set

h2 = ldsc_h2_sumstats("trait1.tsv.gz", reference, scores, m_snps=M)
rg = ldsc_rg_sumstats("trait1.tsv.gz", "trait2.tsv.gz", reference, scores,
                      m_snps=M, ncores=4)

print("h²:", h2.estimate.h2, "SE:", h2.estimate.h2_se)
print("r_g:", rg.estimate.rg, "SE:", rg.estimate.rg_se)
print("genetic covariance:", rg.estimate.gcov)
print("trait 1 QC:", rg.qc[0])
print("trait 2 alignment:", rg.alignment[1])
print("common before chi-square filter:", rg.n_common,
      "retained:", len(rg.reference_index),
      "removed by chi-square filter:", rg.n_chi2_dropped)
```

The `.bim` can instead be passed as a path, or replaced by a `.pvar` file.
The reference must have chromosome blocks in contiguous order, with positions
sorted within each chromosome. The scores need not be in `.npy` format, but
you must load and align them yourself: PLDSC does not read LDSC's text score
files directly. Check variant identities and alleles, not just vector length.

Allelign accepts tab-, comma- or whitespace-delimited files, including gzip,
and recognizes common GWAS column names. Supply effect and other alleles,
a signed effect (beta or odds ratio) and its standard error or a p-value,
and an appropriate sample size. See `allelign.read_sumstats` for column
aliases and overrides. Allelign can use a per-variant N column or derive
effective N from case/control columns; inspect the `sumstats_parse` entry in
the QC log for the selected source. Check that it is appropriate for the
study. If a file uses nonstandard names, pass them explicitly, for example
`read_options={"beta": "EFFECT", "se": "STDERR", "n_eff": "NEFF"}`.
When a `TEST` column is present, allelign retains `ADD` and `ADD-CONDTL`
rows and reports the others in the parse log. For a file with no N column
and a justified constant sample size:

```python
h2 = ldsc_h2_sumstats("trait1.tsv.gz", reference, scores, m_snps=M,
                       read_options={"n_eff": 100_000})
```

A numeric `n_eff` **overrides every row**, even when the file contains a
per-variant N column. To select a particular N column without flattening it,
use `read_options={"n_eff": "column_name"}`. The bivariate wrapper accepts
`read_options1` and `read_options2` separately. For case-control GWAS,
choose an effective sample size appropriate to the analysis; the usual
per-study value is `4 / (1 / N_cases + 1 / N_controls)`. For a compatible
meta-analysis of non-overlapping studies with different case fractions, sum
the studies' effective sample sizes rather than applying this formula to
pooled counts.

## QC and variant selection

The file wrappers use allelign to remove non-finite effects, non-positive
standard errors or sample sizes, and duplicated identifiers; filter N below
70% of the median; and, when available, filter MAF below 0.01 and INFO below
0.7. They align alleles to the reference, flip effects as needed, and drop
ambiguous A/T and C/G sites by default. Overrides are passed as `qc_options`
for one trait or `qc_options1` and `qc_options2` for two; see
`allelign.qc_sumstats`. The `qc` and `alignment` tuples contain a log per trait.
Inspect the alignment logs for unmatched variants and allele mismatches before
interpreting the estimates. When both files have coordinates, allelign checks
them even for matching rsIDs; this catches some build mismatches but does not
establish that the builds agree. Use `drop_ambiguous=False` only when the
strand and orientation of palindromic alleles are known independently.

After alignment, the wrappers retain variants with chi-square at most
`max(80, 0.001 * N)` for each trait, using the individual variant's effective
`N` when it varies. Set `chi2_cap=False` only when a
separately justified selection has already been made. The wrappers do **not**
exclude long-range LD regions themselves; a bundled reference from
`load_ld_reference` already leaves them out. For your own reference, if the
input and reference are GRCh37 and those exclusions are appropriate:

```python
ss = allelign.read_sumstats("trait1.tsv.gz")
in_region = allelign.in_long_range_ld(ss.chrom, ss.pos,
                                      genome_build="GRCh37", include_apoe=True)
h2 = ldsc_h2_sumstats(ss.subset(~in_region), reference, scores, m_snps=M)
```

The built-in region list is for GRCh37. For another build, provide an
appropriate `regions` list to allelign. An in-memory `allelign.Sumstats` object
can be passed to either wrapper; `read_options` apply only to file paths.

## LD scores and the variant count

`m_snps` is the number of variants **summed into each LD score**, not the
number of GWAS rows surviving QC or the number of rows used in the regression.
Keep that number fixed when removing regression rows. For example, if scores
sum over one million variants but 800,000 regression variants are retained,
use `m_snps=1_000_000`. The score vector must still contain one value per
reference table row, in its order. Scores of zero or less cannot be used in
these regressions, including scores made non-positive by a finite-reference
correction; investigate the reference or variant selection rather than silently
clipping them.

PLDSC can calculate scores from an LD correlation matrix:

```python
import numpy as np
from pldsc import ld_scores

R = np.array([[1.0, 0.3, 0.1], [0.3, 1.0, 0.2], [0.1, 0.2, 1.0]])
scores = ld_scores(R, n_ref=1_000)  # n_ref is optional
```

For large data, pass disjoint `(matrix, indices)` blocks covering the full
variant map. Dense float and scaled int8 matrices are supported, as are
`LowRankLD` factors, which avoid constructing a dense matrix. Block scores
omit between-block LD. The optional `n_ref` applies a finite-reference
correction. A small or arbitrarily partitioned example is not a validated
LDSC reference; for real analyses use scores appropriate to the population,
variant set and LD window. `local_ld_scores` is also available to calculate
windowed scores from genotypes on one chromosome at a time.

## Use arrays directly

For already QC'd and allele-aligned arrays in genomic order, the low-level
functions avoid file processing. Here `z1`, `z2`, and `ell` are aligned
one-dimensional arrays and `M` is the count used to construct `ell`:

```python
import numpy as np
from pldsc import ldsc_chi2_keep, ldsc_h2, ldsc_rg_z

N1 = np.full(len(z1), 100_000.0)  # replace with the actual per-variant N
N2 = np.full(len(z2), 80_000.0)
keep1, _ = ldsc_chi2_keep(z1**2, N1)
keep2, _ = ldsc_chi2_keep(z2**2, N2)
keep = keep1 & keep2

h2 = ldsc_h2(z1[keep]**2, ell[keep], N1[keep], m_snps=M)
rg = ldsc_rg_z(z1[keep], z2[keep], ell[keep], N1[keep], N2[keep],
               m_snps=M, ncores=4)
print(h2.h2, h2.h2_se, rg.rg, rg.rg_se)
```

These functions do not filter, harmonize or sort inputs. They accept scalar
or per-variant sample sizes. `ldsc_rg_z` takes signed z scores; the separate
`ldsc_rg` takes **LDpred3-standardized effects**, not ordinary GWAS beta
coefficients. Array functions default `m_snps` to the number of supplied
rows, which is inappropriate when the LD scores sum over a larger set, and
issue a `UserWarning` when they do; the compat functions require it.

## LDSC equivalence

Two estimators are available. Neither changes the other's defaults.

| | Default (`ldsc_h2`, `ldsc_rg`, `pldsc h2/rg`) | `ldsc_compat=True`, `--ldsc-compat` | Reference `ldsc` |
|:--|:--|:--|:--|
| Free-intercept h² | one-step IRWLS on all rows | two-step: intercept from χ² < 30, slope on all rows | `--h2`: two-step at 30 |
| r_g and its marginal h² | one-step IRWLS | one-step | `--rg`: one-step unless `--two-step` is given |
| Over-counting weights | regression scores | separate `w_ld` (default: regression scores) | `--w-ld` |
| Fitted mean in the weights | floored at 1 | h² clipped to [0, 1], no floor | same as compat |
| Jackknife | weights refit in each deletion | final weights fixed | fixed |
| Blocks | 100 (h²), 200 (r_g) | 200 | 200 |
| χ² cap max(80, 0.001 N) | on in the sumstats functions and CLI | off | off for one-annotation h² and r_g |
| M | bundled count 9,997,231 unless `--m-snps` | must be given | `M_5_50` in published runs |

The compat mode agrees with the estimator code of `bulik/ldsc` (commit
2fdeeb3, run under Python 3 with import changes only) to about 1e-15 for h²,
intercept, r_g, genetic covariance and every SE, on the committed fixture
`tests/data/ldsc_parity.npz` (3,000 SNPs with per-SNP N, `w_ld` different
from the regression scores and outliers that trigger the two-step filter),
and on 1,159,831 bundled HapMap3 scores with simulated z. It does not
reproduce ldsc's file handling, allele merging or `.sumstats` munging, and the
bundle has no separate `weights_hm3_no_MHC` scores, so with `--ld-reference`
the regression scores also serve as `w_ld` (the region filter below removes
the MHC).

ldsc applies its default two-step cutoff to `--h2` only. In `--rg` it tests
`intercept_h2 is None` after turning that option into a list, so the test
always fails and every r_g regression, including the marginal h² printed in
an r_g log, is one-step unless `--two-step` is given. Compat mode follows the
program, not the intent: on the parity fixture `ldsc --rg` gives r_g 0.3924
(SE 0.0900), as does `ldsc_rg_compat`, and `ldsc --rg --two-step 30` gives
0.3284 (0.0895), as does `ldsc_rg_compat(..., twostep=30)`. As in ldsc, a
two-step cutoff cannot be combined with a fixed cross-trait intercept.

On equal inputs (`w_ld` equal to the regression scores, the same N and M, no
χ² cap) the default estimator and ldsc's one-step fit agree to about 1e-4 in
h², but their standard errors differ more: with 200 blocks each, 0.0592
against 0.0512 on the parity fixture, because the default jackknife refits
the weights in every deletion.

**The M convention is the largest difference.** At a given slope h² and the
genetic covariance scale with M. The bundled count 9,997,231 is 1.677 times
`M_5_50` = 5,961,159, so default h² from `--ld-reference hm3-eur` is 1.677
times what the same fit gives with ldsc's convention; r_g is unchanged. In a
review run on 2026-11-23 (not part of the test suite; regression scores used
as `w_ld`), the 2010 GIANT height GWAS gave h² 0.325 with the default
settings and 0.193 with ldsc's estimator and `M_5_50`, a factor of 1.68; for
the 2010 HDL GWAS the factor was 1.19, because the two-step estimator and the
χ² cap also matter for traits with large loci. Report which M and which
estimator an h² used.

**Speed.** The default jackknife refits the weights in every deletion; the
compat jackknife solves one small system per block with fixed weights. Warm
timings on 1,159,831 rows (the bundled scores with simulated z), 8 cores,
Numba 0.62+, 2026-11-23:

| | h², 100 blocks | h², 200 blocks | r_g, 200 blocks |
|:--|--:|--:|--:|
| Default (one core) | 1.2 s | 2.3 s | 5.4 s |
| `ldsc_compat=True` | — | 0.23 s | 0.39 s |
| `bulik/ldsc` estimator code, same inputs | — | 0.32 s | 0.46 s |

The r_g timings for compat mode and ldsc are for one-step fits, as
`ldsc --rg` runs by default; with a two-step cutoff of 30 they were 0.69 s
and 0.82 s, measured separately the same day. The default is therefore not
faster than ldsc's estimator: its cost grows linearly with the number of
blocks, because each deletion is a full refit.

## Interpretation and performance

By default both fits use one-step iteratively reweighted regression and contiguous
block-jackknife standard errors; weights are refitted within each deletion.
The defaults are 100 blocks for h², 200 for r_g and two weight updates.
`n_blocks` and `n_iter` can be set explicitly. Intercepts are free by default.
For h², `constrain_intercept` fixes the univariate intercept. For bivariate
LDSC it fixes the **cross-trait** intercept, while the two marginal
intercepts remain free. A cross-trait intercept can reflect sample overlap,
correlated confounding or model/reference mismatch and does not by itself
identify their contributions.

The estimates are not clipped: a noisy h² can be negative and r_g can exceed
the interval [-1, 1]. r_g is NaN if a marginal h² is non-positive, and its
SE is NaN if a required jackknife replicate is invalid or fewer than two
blocks are available. Binary-trait h² is not automatically transformed to
the liability scale. Neither estimator implements partitioned h², liability
conversion or ldsc's file munging; see [LDSC equivalence](#ldsc-equivalence)
for what the compat mode does reproduce.

The Numba kernels avoid materializing a separate dataset for each deletion;
`ncores` parallelizes bivariate deletions. Compilation and cache setup affect
the first call, and speed depends on input size, blocks, hardware and thread
settings. A paired warm timing and parity driver against LDpred3 and bipred
is available when those sibling packages are installed. On macOS it requires
AC power and Low Power Mode off:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python experiments/ldsc_provider_parity_bench.py
```
