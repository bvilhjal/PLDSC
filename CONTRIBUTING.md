# Contributing

Use Python 3.11 or newer. From a PLDSC checkout, install the package in
editable mode with the test and optional acceleration dependencies:

```bash
python -m pip install -e ".[test,fast]"
python -m pip install ruff
```

This installs the public `allelign` package automatically. To develop
against a sibling checkout instead, install `../allelign` in editable mode
first. Add the `plot` extra when working on plots.

## Checks

CI tests Python 3.11, 3.13 and 3.14, with an explicitly pinned allelign
revision. Run its lint and test checks from the repository root:

```bash
ruff check .
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 pytest -q
```

Ruff checks pyflakes rules only; long lines and multi-statement lines are
intentional. Pin BLAS threads when comparing numerical results. When changing
an accelerated path, also check its fallback:

```bash
NUMBA_DISABLE_JIT=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  pytest -q tests/test_ldsc_port.py tests/test_review_fixes.py
```

## Experiments and results

The [experiment index](experiments/README.md) describes the current
simulations, real-data inputs and historical results. The committed
[aggregate chromosome 22 bundle](examples/mini) reproduces selected tables
without downloads. Individual-level genotypes and GWAS data belong in
gitignored scratch directories, not in the repository.

For the chromosome 22 genotype-based experiments, fetch the exact 1000
Genomes phase-3 sources used in the recorded results:

```bash
python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22
```

The script uses the S3 `v5a` mirror and checks the recorded SHA-256 digests;
other experiments may require additional data. Record the source, inputs,
configuration, seed and outcome for new experiments, including failed
pre-stated gates. Preserve results and distinguish observed performance from
untested expectations. The earlier estimator research is a historical record,
as explained in the [reading guide](docs/READING_ORDER.md).

Versions are single-sourced from `pldsc.__version__`; update `CITATION.cff`
and `CHANGELOG.md` with a release. Keep benchmarks separate from tests.
