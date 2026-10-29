# Fast LDSC and bivariate LDSC

PLDSC estimates SNP heritability with univariate LD Score regression (LDSC)
and genetic correlation with bivariate LDSC. The implementations adapt the
fast regression kernels from LDpred3 and bipred. Optional Numba compilation
accelerates the delete-block jackknife; bivariate fits can also run its blocks
in parallel. Allelign reads, filters and harmonizes GWAS summary statistics.
The [methods report](../report/methods.pdf) explains the regression equations,
weights, jackknife and interpretation.

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
  --n-blocks 40 --out demo-data/h2.json

pldsc rg --sumstats1 demo-data/trait1.tsv.gz --sumstats2 demo-data/trait2.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --ncores 2 --out demo-data/rg.json
```

Both commands call allelign to load and QC the GWAS, then harmonize alleles
against the reference. The output JSON records estimates, package versions,
input paths, settings, per-trait QC and alignment logs, the number of common
variants before the chi-square cap, and the number retained. An undefined
estimate or SE is `null`. Read the full logs, for example:

```bash
python -m json.tool demo-data/rg.json
```

`--out` creates a new file and refuses to overwrite an existing one. Omit it
to write JSON to standard output; warnings go to standard error. To see all
options, run `pldsc h2 --help` or `pldsc rg --help`. `ncores` parallelizes
bivariate jackknife fits when Numba is installed; it does not accelerate every
part of file loading.

### Use your own files

Provide a `.bim` or `.pvar` reference variant table, sorted by genomic
position within contiguous chromosome blocks; a one-dimensional `.npy` LD
score vector with **one score in exactly the same row order**; and one or two
GWAS files. Scores must be suitable for the GWAS populations and variant set.
For bivariate LDSC, the same reference scores must be suitable for both
traits. Make sure GWAS and reference use the same genome build. PLDSC cannot
infer reference suitability or perform liftover.

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
fit with `allelign inspect trait1.tsv.gz`. For nonstandard column names,
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
variants, and apply the LDSC chi-square cap. They do not remove long-range LD
regions unless requested. For GWAS and reference coordinates known to be
GRCh37, add `--exclude-long-range-ld-grch37`; add `--include-apoe` if that
locus should also be excluded. Each GWAS row must have a valid chromosome and
position for this option. The JSON records the number removed per trait. This
flag declares the build; it does not verify it. For a different build, filter
with a build-appropriate region list using allelign in Python.

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
same genome build as the reference; PLDSC does not check builds or lift over
coordinates. In this example the scores were summed over **exactly the
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
`max(80, 0.001 * N)` for each trait. Set `chi2_cap=False` only when a
separately justified selection has already been made. The wrappers do **not**
automatically exclude long-range LD regions. For example, if the input and
reference are GRCh37 and those exclusions are appropriate:

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
rows, which is inappropriate when the LD scores sum over a larger set.

## Interpretation and performance

Both fits use one-step iteratively reweighted regression and contiguous
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
the liability scale. These one-step fits do not implement every feature of
the original standalone LDSC software.

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
