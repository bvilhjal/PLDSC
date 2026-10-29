#!/usr/bin/env python3
"""Run univariate and bivariate LDSC on temporary synthetic GWAS files.

The simulated z scores have LDSC-like marginal moments; they are independent
across variants. This is an API demonstration, not a realistic GWAS simulation
or a calibration study. Run from the repository root::

    python examples/ordinary_ldsc.py

Use ``--out-dir demo-data`` to keep the generated files for command-line use.
"""
import argparse
from contextlib import nullcontext
from pathlib import Path
from tempfile import TemporaryDirectory

import allelign
import numpy as np

from pldsc import ld_scores, ldsc_h2_sumstats, ldsc_rg_sumstats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", help="create a new directory and keep the synthetic inputs")
    args = parser.parse_args()
    rng = np.random.default_rng(723)
    m = 4000
    positions = 1000 + 1000 * np.arange(m)
    ids = np.array([f"rs{i + 1}" for i in range(m)])
    chrom = np.full(m, "1")
    a1, a2 = np.full(m, "A"), np.full(m, "C")
    n1, n2 = 80_000.0, 70_000.0

    # The scores sum r² over the full variant map, in independent 2-SNP
    # blocks. For an actual study use a validated, ancestry-matched reference.
    correlations = rng.uniform(0.05, 0.85, m // 2)
    blocks = [(np.array([[1.0, r], [r, 1.0]]), np.array([2*i, 2*i + 1]))
              for i, r in enumerate(correlations)]
    scores = ld_scores(blocks)

    mu1 = 1.0 + n1 * 0.20 * scores / m
    mu2 = 1.0 + n2 * 0.20 * scores / m
    cross = np.sqrt(n1 * n2) * 0.10 * scores / m
    a, b = rng.standard_normal((2, m))
    z1 = np.sqrt(mu1) * a
    z2 = cross / np.sqrt(mu1) * a + np.sqrt(mu2 - cross**2 / mu1) * b
    se = np.full(m, 0.01)

    def summary(z, n, swap):
        ea = np.where(swap, a2, a1)
        oa = np.where(swap, a1, a2)
        beta = np.where(swap, -z, z) * se
        return allelign.Sumstats(ids, chrom, positions, ea, oa, beta, se,
                                 np.full(m, n), np.full(m, 0.3), np.ones(m))

    # Half of trait 2's effects count the opposite allele. Allelign aligns
    # them to the reference, including the sign changes, before fitting.
    trait1 = summary(z1, n1, np.zeros(m, dtype=bool))
    trait2 = summary(z2, n2, np.arange(m) % 2 == 0)
    if args.out_dir:
        destination = Path(args.out_dir)
        destination.mkdir(parents=True, exist_ok=False)
        directory = nullcontext(destination)
    else:
        directory = TemporaryDirectory(prefix="pldsc-example-")
    with directory as tmp:
        tmp = Path(tmp)
        bim = tmp / "reference.bim"
        with bim.open("w") as out:
            for snp, pos in zip(ids, positions):
                out.write(f"1\t{snp}\t0\t{pos}\tA\tC\n")
        file1, file2 = tmp / "trait1.tsv.gz", tmp / "trait2.tsv.gz"
        allelign.write_sumstats(file1, trait1)
        allelign.write_sumstats(file2, trait2)
        np.save(tmp / "reference.ldscores.npy", scores)

        h2 = ldsc_h2_sumstats(file1, bim, scores, m_snps=m, n_blocks=40)
        rg = ldsc_rg_sumstats(file1, file2, bim, scores, m_snps=m,
                              n_blocks=40, ncores=2)
        print("Illustrative independent synthetic z scores; not a calibration study")
        print("Trait 1:", h2.estimate)
        print("Traits 1 and 2:", rg.estimate)
        print("Trait 1 QC:", h2.qc[0])
        print("Bivariate retained:", len(rg.reference_index),
              "of", rg.n_common, "common variants;",
              rg.n_chi2_dropped, "dropped by chi-square cap")
        print("Trait 2 alignment:", rg.alignment[1])
        if args.out_dir:
            print("Synthetic inputs written to:", tmp)


if __name__ == "__main__":
    main()
