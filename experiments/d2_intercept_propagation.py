"""D2 follow-up: does the noise-floor under-subtraction move an LDSC intercept?

The D2 kill-test (experiments/relatedness_d2_killtest.py) found the naive
``nu - 1`` finite-sample correction under-subtracts the unlinked-pair r^2 noise
floor on every real 1000G panel (PC-adjusted N_eff ratios 0.55-0.88), but left
the consequence undecided: the record's verdict rule needs the bias propagated
into the quantity D2 is about — the LDSC intercept read as a confounding
barometer.

This script measures that propagation directly. On real 1000G chr22 genotypes
(EUR and FULL panels, PC-residualized scores so only the relatedness floor
differs):

1. Compute windowed LD scores twice: with the package-default denominator
   ``residual_df - 1`` (naive) and with ``pair_df = effective_pair_df(K, PCs)``
   (corrected, METHOD.md §10.1 Eq. M22).
2. Calibrate against truth: simulate chi2_j = z_j^2 with
   E[chi2_j] = 1 + slope * ell_j, taking the *corrected* score as ell (the
   noise-floor-free population score it estimates), at two realistic chi2
   slopes: 0.05 (N=100k, h2=0.5, M=1M) and 0.3 (biobank scale).
3. Regress chi2 on the naive and the corrected score (pldsc.ld_score_regression,
   block jackknife) over Monte Carlo replicates. The mean intercept shift of
   the naive regression vs 1 is the D2 bias; the replicate SD (independent-noise
   replicates, so an *underestimate* of real LDSC SEs) and the canonical real
   intercept SE ~0.01 are the yardsticks.

Run:  python experiments/d2_intercept_propagation.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import realdata_1000g_benchmarks as bench  # noqa: E402
from pldsc import effective_pair_df, local_ld_scores, ld_score_regression  # noqa: E402

MAF = 0.05
K_PC = 10
WINDOW_BP = 250_000.0   # the ±250 kb window of benchmarks B1-B6 (positions in bp)
N_REPS = 200
SLOPES = (0.05, 0.3)
N_BLOCKS = 50
SEED = 20260811


def propagate(name: str, G: np.ndarray, pos: np.ndarray, rng) -> None:
    p = G.mean(0) / 2.0
    keep = (np.minimum(p, 1 - p) >= MAF) & (G.var(0) > 0)
    G, pos = G[:, keep], pos[keep]
    Z = bench.std(G).astype(np.float32)
    K = ((Z @ Z.T) / Z.shape[1]).astype(np.float64)
    pcs = np.linalg.eigh(K)[1][:, -K_PC:]
    n_eff = effective_pair_df(K, pcs)

    # PC-residualized scores (cov-LDSC style): what differs between the two
    # denominators is then purely the relatedness noise floor.
    kw = dict(covariates=pcs, window_cm=WINDOW_BP)
    ell_corr = local_ld_scores(G, pos, pair_df=n_eff, **kw)
    ell_naive = local_ld_scores(G, pos, **kw)  # default pair_df = residual_df - 1
    bias = ell_naive - ell_corr
    print(f"  {name}: N={G.shape[0]}, M={G.shape[1]}, N_eff={n_eff:.0f}, "
          f"mean score bias {bias.mean():+.4f} (sd {bias.std():.4f})")

    ell_true = np.clip(ell_corr, 0.0, None)
    for slope in SLOPES:
        ints_naive = np.empty(N_REPS)
        ints_corr = np.empty(N_REPS)
        amp = np.sqrt(slope * ell_true)
        for r in range(N_REPS):
            z = amp * rng.standard_normal(ell_true.size) + rng.standard_normal(ell_true.size)
            chi2 = z * z
            ints_naive[r] = ld_score_regression(chi2, ell_naive, n_blocks=N_BLOCKS).intercept
            ints_corr[r] = ld_score_regression(chi2, ell_corr, n_blocks=N_BLOCKS).intercept
        print(f"    slope={slope:<5} intercept naive  {ints_naive.mean():.4f} "
              f"+- {ints_naive.std(ddof=1):.4f}   corrected {ints_corr.mean():.4f} "
              f"+- {ints_corr.std(ddof=1):.4f}   "
              f"shift {ints_naive.mean() - 1.0:+.4f} "
              f"(analytic {slope * bias.mean():+.4f})")


def main():
    t0 = time.time()
    print("D2 intercept propagation: naive vs corrected noise-floor subtraction")
    samples, G22, pos22 = bench.load_chrom(22)  # t25 cache from the kill-test
    s21, G21, _ = bench.load_chrom(21)
    assert s21 == samples
    pop, superpop = bench.load_panel()
    masks = bench.cohort_masks(samples, pop, superpop)
    rng = np.random.default_rng(SEED)
    for cname in ("EUR", "FULL"):
        propagate(cname, G22[masks[cname]], pos22.astype(np.float64), rng)
        print(f"    ({time.time() - t0:.0f}s)")
    print("\nRead: |shift| well under ~0.01 (a typical real intercept SE) kills "
          "D2; a shift at or above it keeps the barometer claim alive.")


if __name__ == "__main__":
    main()
