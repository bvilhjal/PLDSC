# Development environment

PLDSC requires Python ≥ 3.11 (`pyproject.toml`), NumPy ≥ 2.0 and the private
`allelign` package, and is
developed on Python 3.14. A default Anaconda base environment is 3.10 and
is below both floors.

**There is no dedicated `pldsc` conda environment.** Some archived
experiment docstrings say `conda run -n pldsc python ...`; that
environment never existed on this machine — plain `python` is correct. On
the maintainer's machine the package is installed editable in the shared
`ldpred3-accelerate` conda environment (Python 3.14, alongside the
sibling packages). Any 3.11+ interpreter works:

```bash
python -m pip install -e ../allelign
python -m pip install -e ".[test,fast]"  # add plot for the plotting extra
```

## Checks that CI runs

`.github/workflows/ci.yml` runs on Python 3.11, 3.13 and 3.14:

```bash
ruff check .                          # pyflakes only (select = ["F"]), by design
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 pytest -q
NUMBA_DISABLE_JIT=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  pytest -q tests/test_ldsc_port.py tests/test_review_fixes.py
```

Pin BLAS threads as
above when comparing numbers across machines: BLAS reductions are the one
place threading can move the last digits. Long lines and multi-statement
lines are intentional; do not reformat or add style rules.

## Real-data experiments and `.scratch_1000g/`

A fresh clone is 100% data-free; `.scratch_1000g/` is gitignored and holds
the 1000 Genomes phase-3 files. Build only what you need (chr22 for the
headline; chr21/22/10 for the structure-share suite):

```bash
python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22
```

The S3 mirror is primary (the EBI FTP has been unreachable from this
environment and serves the `v5b` re-release; all recorded numbers used
`v5a`, whose sha256 the script verifies). The script also derives the
plink-format genetic maps the LOCO driver needs, because the historical
Beagle host for them is gone. No-download reproduction of the headline
tables does not need any of this: `examples/reproduce_chr22_tables.py`
runs off the committed bundle in [examples/mini](examples/mini) (data
provenance and terms are documented there).

## How contributions work in this repo

- The live plan and its pre-stated gates live in
  [docs/STRATEGY.md](docs/STRATEGY.md) §4; the doc classification (what is
  live, what is archive) in [docs/READING_ORDER.md](docs/READING_ORDER.md).
  **The archive is a record, not a backlog** — see STRATEGY §6 ("what we
  stop") before proposing a new estimator or construction.
- Experiments print their own PASS/FAIL against rules that are written
  *before* they run. If a rule fails, record the failure and the verdict
  honestly (see `docs/MULTI_ANCESTRY.md` §6 Table 1 for the model); a
  pre-stated gate closing is a result.
- New experiments get an entry in `experiments/README.md` (dated, at the
  top with the other live entries). Record runtimes — several archived
  experiments have none.
- `CHANGELOG.md` is Keep-a-Changelog style with one `[Unreleased]`
  section; versions are single-sourced from `pldsc.__version__`.
- Benchmarks and comparisons that measure performance belong under
  `benchmarks/` (not yet present) with archived results, not in `tests/`.
