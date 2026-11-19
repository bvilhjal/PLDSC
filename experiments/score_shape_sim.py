"""Known-truth calibration of the residual LD-score coefficient.

Design fixed before running the pilot. The two group proportions are 90:10,
75:25 and 50:50, rho=0.8, Nh²/M=1. For each composition use true LD, 400 or
6,000 reference people per group. The negative control generates a EUR-only
GWAS while supplying a distinct matched candidate. The positive control
generates the specified meta-GWAS. Every replicate draws new GWAS statistics
and new independent reference samples on the same population LD landscape.

With 1,000 replicates per cell, the existing design targets are 92–98%
coverage of nominal 95% intervals and at least 80% two-sided power. Also
record how often held-out mean squared error favors the matched predictor.
Use one observation per independent two-SNP block, averaged over its two
correlated z² values. Thirty-two regions each contain 512 independent blocks.
The Gaussian design is a test of the method, not a model of the GBMI cohort.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_shape import centered_products, chromosome_losses, gain_percent, residual_slope  # noqa: E402


GROUPS = 32
BLOCKS = 16_384
SHARES = (0.10, 0.25, 0.50)
PANEL_SIZES = (400, 6000)
RHO = 0.8
T_CRIT = 2.0395134464  # two-sided t, 31 df (32 delete-region estimates)


def reference_correlations(rng, true_r, n):
    """Exact distribution of Pearson r from n independent bivariate Gaussian rows."""
    orthogonal = np.sqrt(1 - true_r**2)
    along = true_r * np.sqrt(rng.chisquare(n - 1, size=true_r.size)) + orthogonal * rng.normal(size=true_r.size)
    across = orthogonal * np.sqrt(rng.chisquare(n - 2, size=true_r.size))
    return along / np.hypot(along, across)


def scores(r_e, r_f, panel_n, f):
    """Squared and cross-LD row sums for two-locus blocks; f is the F share."""
    square_e, square_f = r_e**2, r_f**2
    if panel_n:
        square_e -= (1 - square_e) / (panel_n - 2)
        square_f -= (1 - square_f) / (panel_n - 2)
    eur = 1 + square_e
    matched = (1 - f)**2 * eur + f**2 * (1 + square_f) + 2 * f * (1 - f) * RHO * (1 + r_e * r_f)
    return eur, matched


def block_chi_square(r_e, r_f, f, draws):
    """Mean of correlated z², using the two exact eigenmodes of the meta-GWAS."""
    e, g = 1 - f, f
    mean_square = np.zeros_like(r_e)
    for sign, chi_square in zip((1, -1), draws):
        a, b = 1 + sign * r_e, 1 + sign * r_f
        sampling = e * a + g * b
        genetic = (e * a)**2 + (g * b)**2 + 2 * e * g * RHO * a * b
        mean_square += (sampling + genetic) * chi_square / 2
    return mean_square


def simulate(replicates, seed):
    population, references, outcomes = [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(3)]
    r_e = population.uniform(0.1, 0.9, BLOCKS)
    r_f = 0.25 * r_e + 0.75 * population.uniform(-0.35, 0.8, BLOCKS)
    methods = (0, *PANEL_SIZES)  # zero denotes exact population LD
    rows = np.zeros((replicates, len(SHARES), len(methods), 2, 3))
    # Fields per replicate: residual coefficient, delete-region SE, held-out gain.
    design = np.ones((BLOCKS, 5))
    by_region = design.reshape(GROUPS, -1, 5)
    for i in range(replicates):
        estimates = [(r_e, r_f)] + [
            (reference_correlations(references, r_e, n), reference_correlations(references, r_f, n))
            for n in PANEL_SIZES]
        draws = outcomes.normal(size=(2, BLOCKS))**2
        design[:, 3] = block_chi_square(r_e, r_f, 0, draws)
        for j, share in enumerate(SHARES):
            design[:, 4] = block_chi_square(r_e, r_f, share, draws)
            for k, (n, (ref_e, ref_f)) in enumerate(zip(methods, estimates)):
                design[:, 1], design[:, 2] = scores(ref_e, ref_f, n, share)
                all_moments = by_region.swapaxes(1, 2) @ by_region
                for truth in (0, 1):
                    cols = [0, 1, 2, 3 + truth]
                    moments = all_moments[:, cols][:, :, cols]
                    total = moments.sum(axis=0)
                    theta = residual_slope(centered_products(total))
                    deleted = residual_slope(centered_products(total - moments))
                    se = np.std(deleted) * np.sqrt(GROUPS - 1)
                    gain = gain_percent(chromosome_losses(moments))
                    rows[i, j, k, truth] = (theta, se, gain)
    cells = []
    for j, share in enumerate(SHARES):
        for k, n in enumerate(methods):
            for truth, model in enumerate(("EUR", "matched")):
                theta, se, gain = rows[:, j, k, truth].T
                cells.append(dict(other_share=share, panel_n=n or None, generating_model=model,
                                  theta_truth=truth, mean_theta=float(theta.mean()), mean_se=float(se.mean()),
                                  sd_theta=float(theta.std(ddof=1)),
                                  covered=int(np.sum(np.abs(theta - truth) <= T_CRIT * se)),
                                  rejected_zero=int(np.sum(np.abs(theta) > T_CRIT * se)),
                                  matched_lower_loss=int(np.sum(gain > 0)),
                                  mean_gain_percent=float(gain.mean())))
    return dict(seed=seed, replicates=replicates, blocks=BLOCKS, regions=GROUPS, rho=RHO,
                expected_slope=1.0, panel_n_unit="people per group; null denotes population LD",
                interval="two-sided 95% delete-region t(31)", numpy_version=np.__version__, cells=cells)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261119)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("at least two independent replicates are required")
    result = simulate(args.replicates, args.seed)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"{args.replicates} replicates per cell; seed {args.seed}")
    for row in result["cells"]:
        print(f"F {row['other_share']:.0%} / panel {row['panel_n'] or 'exact'} / {row['generating_model']}: "
              f"mean slope {row['mean_theta']:.3f}, covered {row['covered']}, "
              f"rejected zero {row['rejected_zero']}, matched won {row['matched_lower_loss']}")


if __name__ == "__main__":
    main()
