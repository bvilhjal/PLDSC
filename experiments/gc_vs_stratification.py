"""Does a raw-loading regressor separate genomic control from stratification?

Strategy-review probe (2026-09-24, docs/STRATEGY.md §3, consumer C3). The
sibling package ppb rescales z by the LDSC intercept (`z_true = z_obs /
sqrt(intercept)`, ppb docs/CALIBRATION.md C2-C3) and records that the
intercept "is one number standing for two mechanisms": genomic-control (GC)
deflation divides every chi2 by lambda, while stratification ADDS to chi2.
When both are present they can cancel in the intercept, and the deflation
goes undetected.

The PLDSC mechanism (raw_loading_control.py) says stratification along axis
u_p does not add a constant: it enters E[chi2_j] through the SNP's own raw
loading a_j = corr(x_j, u_p)^2. So a three-term regression

    chi2_j ~ 1 + ell_j + a_j

should return intercept 1/lambda under GC (which scales all three terms) and
push stratification into the a_j coefficient, where two-term LDSC
(chi2_j ~ 1 + ell_j) folds it into the intercept.

Design: admixed cohort (pldsc.simulate), N=1200, 5 x 400 variants, 1 cM
windows, h2=0.5, GWAS without PC adjustment, the cohort's own top axis as the
loading (an ORACLE: real use would take loadings from an external panel,
which attenuates them). 2 x 2 factorial: GC lambda in {1, 1.2} x
stratification in {none, an environmental shift along u1 tuned so it adds
~0.2 to mean chi2}. Both magnitudes are illustrative, not calibrated to a
real GWAS.

Pre-stated rule: PASS if, in the GC+STRAT arm, the three-term intercept is
within 0.02 of 1/lambda = 0.833 while the two-term intercept is off by more
than 0.05. Otherwise the ppb use of the raw loading is dropped from the
strategy.

Outcome (recorded in docs/STRATEGY.md §3.1): the default run (4 x 100) FAILS
that rule, -0.044 against a 0.02 tolerance -- but the tolerance was smaller
than the design's own Monte Carlo SE (~0.023), and the GC arms are exact
rescalings of the no-GC arms (dividing chi2 by lambda divides every OLS
coefficient by lambda), so the miss is the no-GC arm's deviation from 1
divided by lambda. The paired rows below (same phenotype draw with and
without the shift) are the mechanistic evidence. `--realizations 12 --reps
500` is the higher-precision DIAGNOSTIC rerun, run after the failure and
labelled as such; it does not replace the pre-stated verdict.

`--homogeneous-check` runs the post-failure diagnostic that localizes the
miss: the same LD model with fst -> 0 (no structure at all), where ordinary
two-term LDSC is the right model. Its intercept is also ~0.95, so the -4% is
the toy LD model's baseline LDSC intercept bias (AR1 LD, rho = 0.9 per
0.15 cM, extends well past the 1 cM window), not something the raw-loading
term introduces.

Run:  python experiments/gc_vs_stratification.py [--realizations R --reps K]
      python experiments/gc_vs_stratification.py --homogeneous-check
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from raw_loading_control import CHROMS, WINDOW, _tune, annot_score, cohort  # noqa: E402

from pldsc.genome import grm_eigh, unit_columns  # noqa: E402
from pldsc.simulate import simulate_admixture, standardize  # noqa: E402

SEED = 20260924
N_INDIV = 1200
REALIZATIONS = 4
REPS = 100
H2 = 0.5
LAMBDA_GC = 1.2
STRAT_EXCESS = 0.2          # mean-chi2 excess the stratification adds (illustrative)


def fit(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def one_realization(rng, reps):
    Z, Xn, chrom, pos, u1, a = cohort(rng, N_INDIV)
    N, M = Z.shape
    ell = annot_score(Xn, chrom, pos, WINDOW)
    two = np.column_stack([np.ones(M), ell])
    three = np.column_stack([np.ones(M), ell, a])
    ax = u1 - u1.mean()
    ax /= ax.std()
    za = Z.T @ ax / np.sqrt(N)

    def draw():
        g = Z @ rng.normal(0.0, np.sqrt(H2 / M), M)
        g -= g.mean()
        g *= np.sqrt(H2) / g.std()
        e = rng.normal(size=N)
        e -= e.mean()
        e *= np.sqrt(1 - H2) / e.std()
        return g + e

    # Tune the shift once per realization on fixed draws (deterministic in s).
    base = [draw() for _ in range(20)]
    clean_mean = np.mean([(((Z.T @ y / np.sqrt(N)) / y.std()) ** 2).mean() for y in base])

    def excess(s):
        vals = []
        for y in base:
            ys = y + s * ax
            vals.append((((Z.T @ ys / np.sqrt(N)) / ys.std()) ** 2).mean())
        return np.mean(vals) - clean_mean

    s_star = _tune(excess, STRAT_EXCESS, hi=5.0)

    arms = {(gc, st): [] for gc in (False, True) for st in (False, True)}
    for _ in range(reps):
        y0 = draw()
        for (gc, st) in arms:
            y = y0 + s_star * ax if st else y0
            chi2 = ((Z.T @ y / np.sqrt(N)) / y.std()) ** 2
            if gc:
                chi2 = chi2 / LAMBDA_GC
            b2, b3 = fit(two, chi2), fit(three, chi2)
            arms[(gc, st)].append((b2[0], b3[0], b3[2] * a.mean(), chi2.mean()))
    return arms, M, s_star, float(np.corrcoef(a, ell)[0, 1]), float(za.std())


def homogeneous_check(rng, realizations=8, reps=300):
    """Clean phenotypes only; fst ~ 0 vs 0.15; two- vs three-term intercept."""
    for fst in (1e-6, 0.15):
        i2, i3 = [], []
        for _ in range(realizations):
            G, chrom, pos = simulate_admixture(rng, CHROMS, N_INDIV, fst=fst, rho=0.9)
            Z = standardize(G)
            Xn, dead = unit_columns(Z)
            keep = ~dead
            Z, Xn, chrom, pos = Z[:, keep], Xn[:, keep], chrom[keep], pos[keep]
            _, U = grm_eigh(Z)
            N, M = Z.shape
            a = (Z.T @ U[:, -1] / np.sqrt(N)) ** 2
            ell = annot_score(Xn, chrom, pos, WINDOW)
            two = np.column_stack([np.ones(M), ell])
            three = np.column_stack([np.ones(M), ell, a])
            for _ in range(reps):
                g = Z @ rng.normal(0.0, np.sqrt(H2 / M), M)
                g -= g.mean()
                g *= np.sqrt(H2) / g.std()
                e = rng.normal(size=N)
                e -= e.mean()
                e *= np.sqrt(1 - H2) / e.std()
                y = g + e
                chi2 = ((Z.T @ y / np.sqrt(N)) / y.std()) ** 2
                i2.append(fit(two, chi2)[0])
                i3.append(fit(three, chi2)[0])
        i2, i3 = np.asarray(i2), np.asarray(i3)
        print(f"fst={fst:g}: two-term intercept {i2.mean():.3f} (se {i2.std() / np.sqrt(len(i2)):.3f}), "
              f"three-term intercept {i3.mean():.3f} (se {i3.std() / np.sqrt(len(i3)):.3f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--realizations", type=int, default=REALIZATIONS)
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--homogeneous-check", action="store_true")
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    if args.homogeneous_check:
        homogeneous_check(rng)
        return
    acc = {}
    Ms, shifts, cors = [], [], []
    for _ in range(args.realizations):
        arms, M, s_star, cor, _ = one_realization(rng, args.reps)
        Ms.append(M)
        shifts.append(s_star)
        cors.append(cor)
        for k, v in arms.items():
            acc.setdefault(k, []).extend(v)
    n_tot = args.realizations * args.reps
    print(f"N={N_INDIV}, M~{int(np.mean(Ms))}, {args.realizations} realizations x {args.reps} reps, "
          f"lambda_GC={LAMBDA_GC}, strat shift s*={np.mean(shifts):.3f} "
          f"(adds {STRAT_EXCESS} to mean chi2), corr(a, ell)={np.mean(cors):+.3f}")
    print(f"target intercept: 1.000 without GC, {1 / LAMBDA_GC:.3f} with GC\n")
    hdr = (f"{'arm':12s} {'mean chi2':>9s} {'2-term int':>11s} {'(se)':>7s} "
           f"{'3-term int':>11s} {'(se)':>7s} {'b_a*mean(a)':>12s}")
    print(hdr)
    print("-" * len(hdr))
    verdict = {}
    for (gc, st), rows in acc.items():
        r = np.asarray(rows)
        name = ("GC" if gc else "--") + "+" + ("STRAT" if st else "--")
        se = r.std(axis=0) / np.sqrt(n_tot)
        print(f"{name:12s} {r[:, 3].mean():9.3f} {r[:, 0].mean():11.3f} {se[0]:7.3f} "
              f"{r[:, 1].mean():11.3f} {se[1]:7.3f} {r[:, 2].mean():12.3f}")
        verdict[(gc, st)] = r[:, :2].mean(axis=0)
    print("\npaired effect of the stratification shift (same draws; want 0 for a clean readout)")
    for gc in (False, True):
        d = np.asarray(acc[(gc, True)]) - np.asarray(acc[(gc, False)])
        se = d.std(axis=0) / np.sqrt(n_tot)
        print(f"  {'GC' if gc else 'no GC'}: two-term intercept {d[:, 0].mean():+.3f} (se {se[0]:.3f}), "
              f"three-term intercept {d[:, 1].mean():+.4f} (se {se[1]:.4f}), "
              f"loading term {d[:, 2].mean():+.3f}")
    two_int, three_int = verdict[(True, True)]
    target = 1 / LAMBDA_GC
    ok = abs(three_int - target) <= 0.02 and abs(two_int - target) > 0.05
    print(f"\nGC+STRAT: two-term off by {two_int - target:+.3f}, "
          f"three-term off by {three_int - target:+.3f} -> {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
