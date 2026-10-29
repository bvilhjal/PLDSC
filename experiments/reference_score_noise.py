"""Known-truth sensitivity to reference-panel score error.

Two populations have the same block LD and perfectly correlated effects, so
there is no population LD mismatch. Compare the panel plug-in kappa with slopes
against the known expected chi-square score. This is a calibration demonstration,
not an estimate of the 1000 Genomes results' bias.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pldsc import cross_ancestry_ld_scores, meta_ld_scores, meta_weights, single_reference_bias  # noqa: E402


def run(n_ref=500, replicates=30, seed=73):
    rng = np.random.default_rng(seed)
    groups, width = 80, 5
    corr = np.linspace(.05, .5, groups)
    m = groups * width
    truth = np.repeat(1 + (width - 1) * corr**2, width)
    rows = []
    for _ in range(replicates):
        def draw():
            common = rng.normal(size=(n_ref, groups, 1))
            noise = rng.normal(size=(n_ref, groups, width))
            return (np.sqrt(corr)[None, :, None] * common
                    + np.sqrt(1 - corr)[None, :, None] * noise).reshape(n_ref, m)
        panel = cross_ancestry_ld_scores([draw(), draw()], np.zeros(m))
        omega, _ = meta_weights([7500, 2500], panel.available)
        meta = meta_ld_scores(panel, omega)
        eur = panel.scores[0, 0]
        plug_in = single_reference_bias(meta, eur).kappa
        eur_slope = single_reference_bias(truth, eur).kappa
        meta_slope = single_reference_bias(truth, meta).kappa
        rows.append((plug_in, eur_slope, meta_slope, eur_slope / meta_slope))
    return np.mean(rows, axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-ref", type=int, default=500)
    ap.add_argument("--replicates", type=int, default=30)
    ap.add_argument("--seed", type=int, default=73)
    args = ap.parse_args()
    if args.n_ref < 4 or args.replicates < 1:
        ap.error("n-ref must be at least 4 and replicates must be positive")
    result = run(args.n_ref, args.replicates, args.seed)
    print(f"identical population LD, rho=1; n_ref={args.n_ref}, replicates={args.replicates}, seed={args.seed}")
    for label, value in zip(("panel plug-in kappa", "true expected EUR h2/h2",
                             "true expected meta h2/h2", "true expected EUR/meta"), result):
        print(f"  {label}: {value:.5f}")


if __name__ == "__main__":
    main()
