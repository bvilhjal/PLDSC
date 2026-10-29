# PLDSC

PLDSC provides fast LD Score regression (LDSC) for SNP heritability and
bivariate LDSC for genetic correlation. Its regression kernels are adapted
from LDpred3 and bipred; optional Numba compilation accelerates block
jackknifing. [Allelign](https://github.com/bvilhjal/allelign) loads, filters
and aligns GWAS summary statistics to the reference alleles.

The package also includes tools and research on LD scores in structured and
multi-ancestry populations. Its API may still change before version 1.0.

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
  --n-blocks 40 --out demo-data/h2.json

pldsc rg --sumstats1 demo-data/trait1.tsv.gz --sumstats2 demo-data/trait2.tsv.gz \
  --reference demo-data/reference.bim \
  --ld-scores demo-data/reference.ldscores.npy --m-snps 4000 \
  --n-blocks 40 --ncores 2 --out demo-data/rg.json
```

The example also prints the estimates and QC logs. Its independent synthetic
z scores illustrate the workflow; they do not establish calibration. The
JSON files include processing and alignment counts as well as the estimates.

For real data, use reference LD scores appropriate to the populations and
genome build, in exactly the order of the reference variants. `--m-snps` is
the number of variants summed into each LD score, not the number of GWAS
rows retained after filtering. Check the QC and alignment logs before
interpreting estimates. PLDSC does not verify genome builds or perform
liftover.

## Guides

- [LDSC guide](docs/LDSC.md): real-data inputs, command-line options,
  allelign QC, Python APIs, interpretation and performance.
- [Runnable notebook](examples/ldsc_notebook.ipynb): create example GWAS,
  preprocess them with allelign and fit both models.
- [Methods report](report/methods.pdf): regression equations, weights,
  uncertainty and assumptions.
- [Contributing](CONTRIBUTING.md): development setup and checks.

## Multi-ancestry research

[Multi-ancestry methods and results](docs/MULTI_ANCESTRY.md) describes the
model, the 1000 Genomes panel analyses and the GBMI asthma example, including
their limitations. The panel plug-in slope comparison is a sensitivity
analysis, not a calibrated estimate of absolute heritability bias. A
genome-wide bridge between two reference protocols missed its pre-stated
tolerance, and the asthma comparison did not robustly meet its 2-SE gate.

Reproduce the chromosome 22 tables from the committed aggregate LD-score
bundle, without downloading genotypes:

```bash
python examples/reproduce_chr22_tables.py
```

See the [experiment index](experiments/README.md) for simulations, data
requirements and historical results, and the [reading guide](docs/READING_ORDER.md)
for the current research and the archived estimator line.
