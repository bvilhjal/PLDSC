# PLDSC

PLDSC provides LD Score regression (LDSC) for SNP heritability and
bivariate LDSC for genetic correlation. Its default regression kernels are
adapted from LDpred3 and bipred and are not the reference `ldsc` estimator;
an opt-in `--ldsc-compat` mode reproduces that estimator and is tested
against `ldsc`'s own code ([what differs](docs/LDSC.md#ldsc-equivalence)).
Optional Numba compilation accelerates block jackknifing. [Allelign](https://github.com/bvilhjal/allelign) loads, filters
and aligns GWAS summary statistics to the reference alleles.

The package also includes tools and research on LD scores in structured and
multi-ancestry populations. Its API may still change before version 1.0.

**[Which LD score matches a GWAS?](report/methods.pdf)**
([read online](report/methods.md)). The unified report connects ordinary
LDSC to multi-ancestry meta-analysis, the simulation results, the GBMI
asthma application, and the earlier investigations into mixed models,
relatedness, enrichment, assortative mating and parental effects.

## Install and run

With Python 3.11 or newer, clone the repository and install it:

```bash
git clone https://github.com/bvilhjal/PLDSC.git
cd PLDSC
python -m pip install ".[fast]"
```

Allelign and NumPy are installed as dependencies. The `fast` extra installs
Numba; without it, `python -m pip install .` uses the NumPy implementation.
The first accelerated call includes compilation overhead.

Create synthetic inputs in a new directory and run both analyses:

```bash
python examples/ordinary_ldsc.py --out-dir demo-data

pldsc h2 --sumstats demo-data/trait1.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --allow-low-coverage --out demo-data/h2.json

pldsc rg --sumstats1 demo-data/trait1.tsv.gz --sumstats2 demo-data/trait2.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --ncores 2 --allow-low-coverage --out demo-data/rg.json
```

The example also prints the estimates and QC logs. Its independent synthetic
z scores illustrate the workflow; they do not establish calibration. The
JSON files include processing and alignment counts as well as the estimates.

For real data, use reference LD scores appropriate to the populations and
genome build, in exactly the order of the reference variants. `--m-snps` is
the number of variants summed into each LD score, not the number of GWAS
rows retained after filtering. h² scales with it: the bundled
`--ld-reference hm3-eur` count (9,997,231) gives h² 1.677 times that of the
MAF 5–50% count (5,961,159) that published `ldsc` results use, so pass
`--m-snps 5961159` to compare with them, and add `--ldsc-compat` to use
ldsc's estimator as well. Check the QC and alignment logs before
interpreting estimates. The command line refuses runs with fewer than
200,000 regression SNPs (`--min-snps`, `--allow-low-coverage`) and, when
`--sumstats-build` is given, a build that differs from the reference's.
Without that option, allelign 0.3 or later detects a GWAS on the other build
during alignment and lifts it over, recording this in the alignment log. It
reads allelign-supported files and ldsc `munge_sumstats` output (SNP, A1, A2,
N, Z); the latter's N is a total sample size, so for case-control traits see
[the LDSC guide](docs/LDSC.md#use-your-own-files).

For a GWAS that pools ancestries, a single-population reference LD score
mis-scales h² and the genetic covariance; see
[Multi-ancestry research](#multi-ancestry-research) below.

## Guides

- [LDSC guide](docs/LDSC.md): real-data inputs, command-line options,
  allelign QC, Python APIs, interpretation and performance.
- [Runnable notebook](examples/ldsc_notebook.ipynb): create example GWAS,
  preprocess them with allelign and fit both models.
- [Contributing](CONTRIBUTING.md): development setup and checks.

## Multi-ancestry research

**What the data show so far.** On 953,596 GBMI asthma SNPs, replacing the
European LD score with a composition-matched score did not improve
prediction of held-out chromosomes: the pooled prediction error fell by
0.0007% (95% interval −0.020% to +0.022%), and the matched score did better
on 13 of 22 chromosomes. A fixed-effect meta-analysis also reports
Cochran's Q, the other half of a rotation of the two ancestry statistics. In
simulations on 1000 Genomes LD the pooled statistic and Q together recovered
the cross-ancestry genetic covariance, at 1.7 times the standard deviation of
the two-strata estimate for an 85% EUR study (2.0 times with the three
stratum intercepts that inverse-variance weights need); part of that loss is
provably unavoidable ([Proposition 3](report/methods.md#11-what-does-cochrans-q-add)).
In asthma the fit matched the two-strata estimate under sample-size weights
and fell 22% below it under inverse-variance weights. Refitting with three
stratum intercepts ([experiments/q_model.py](experiments/q_model.py)), which
stops unequal intercepts from leaking into the covariance, left it 19% below
(1.6 SE), so the leak explains little of the gap.

**What is being estimated.** The matched score targets the heritability of
the meta-analysed GWAS itself: a composition-weighted combination of the
ancestry groups' heritabilities and cross-ancestry covariances, with weights
set by each group's share of the sample (or inverse-variance weight). It is
not the heritability of any one population. The asthma fits assume equal
heritabilities across groups and fixed or bracketed cross-ancestry
correlations. Their
estimates are on the balanced observed scale of a case-control study, so
they depend on prevalence and need an explicit conversion before comparison
with liability-scale values.

**Supporting results.** On 1000 Genomes chromosome 22, a 75% EUR and 25%
AFR composition with effect correlation 0.8 gives a panel plug-in slope ratio
of 0.750. Only 1.98% of the matched predictor's variation remains after
rescaling the EUR predictor. In simulations with known LD, the power to
detect the matched pattern rose from 12% at a 90:10 mix to 99% at 50:50,
holding total sample size fixed. In asthma the fixed-correlation matched
score moved the balanced-scale heritability estimate from 0.133 to 0.155.
Only 0.734% of the matched predictor's weighted variation remained after
rescaling the European one, so the change in scale did not measurably
change prediction.
[Asthma prediction](report/methods.md#61-can-the-scores-predict-an-unseen-chromosome),
[power study](report/methods.md#10-when-can-the-two-scores-be-distinguished)
and [full multi-ancestry derivation](docs/MULTI_ANCESTRY.md).

![Estimated asthma heritability rises from 0.133 to 0.155. Prediction-error changes vary in sign across held-out chromosomes.](report/asthma_prediction.png)

**Figure 2 from the report.** Each point on the right holds out a whole
chromosome. The diamond pools the errors on all 953,596 SNPs.

Reproduce the chromosome 22 tables from the committed aggregate LD-score
bundle, without downloading genotypes:

```bash
python examples/reproduce_chr22_tables.py
```

See the [experiment index](experiments/README.md) for simulations, data
requirements and historical results, and the [reading guide](docs/READING_ORDER.md)
for the current research and the archived estimator line.
