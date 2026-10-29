"""Kill-test: does the local-LD channel have POWER to detect polygenic local
adaptation against a matched environmental axis shift, at wild-study scale?

Context. The structure-aliasing limit (METHOD.md sec.14) says the rank-1
DIRECTIONAL component of ancestry-aligned signal is exactly aliased with
environmental confounding, but a HETEROSCEDASTIC-along-axis architecture
(per-SNP h^2 elevated for high-loading SNPs -- what polygenic local adaptation
along an ancestry cline should look like) leaves a corr^2 x ell^local
interaction the environmental shift lacks (identifiability_aliasing.py). That
script showed the channel EXISTS for one seed at one setting. This kill-test
asks whether it has POWER across genotype realizations at wild-population-study
scale (N = 600 / 1200, discrete populations, no reference panel -- PCs and LD
scores are in-sample), because that is the only regime where a new test would
differentiate itself.

Design. Discrete multipop cohorts (Balding-Nichols freqs, AR(1) local LD),
target chr0. Three amplitude-matched phenotype arms per realization, all
explaining the same mean chi^2 (target_chi2):
  ENV : rank-1 environmental shift along the top axis u1 (the confounder; NULL)
  DIR : rank-1 directional genetic (projection of u1 into col(Xc); aliased
        control -- should sit on the ENV null)
  HET : heteroscedastic-along-axis genetic (per-SNP effect variance
        w_j = (c2_j/mean(c2))^kappa; the local-adaptation architecture; ALT)
Statistic: t of the interaction coefficient in
  z_j^2 ~ 1 + corr(x_j,u1)^2 + corr(x_j,u1)^2 * (ell^local - mean),
the same orthogonalized 3-regressor fit as identifiability_aliasing.py.
Thresholds are EMPIRICALLY CALIBRATED on the ENV arm per setting (95th
percentile, one-sided -- HET predicts a positive interaction), so type-I is
controlled by construction; power is the HET exceedance fraction. Power at the
fixed t>2 threshold is also reported.

PRE-REGISTERED VERDICT RULE (fixed before running):
  PASS     -- power >= 80% at the MODERATE setting (N=600, target_chi2=2.5,
              kappa=1) with calibrated type-I ~5%;
  KILL     -- power < 50% at the STRONGEST setting (target_chi2=4.0, kappa=2);
  MARGINAL -- otherwise (the channel is real but too weak to base a method on).

Run:  python experiments/local_adaptation_channel.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.simulate import ar1_haps, standardize  # noqa: E402

SEED = 20260719
CHROMS = [400, 600]         # chr0 = target; chr1 feeds the GRM
FST = 0.10
RHO = 0.6
REPS_PHENO = 2000           # phenotype replicates per arm per realization


def multipop_freqs(rng, m, fst, npops):
    anc = rng.uniform(0.1, 0.9, m)
    a, b = anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst
    return [rng.beta(a, b) for _ in range(npops)]


def simulate(rng, npops, n_pop):
    """Discrete multipop cohort (like wild-population studies)."""
    freqs = [multipop_freqs(rng, m, FST, npops) for m in CHROMS]
    blocks = []
    for g in range(npops):
        h = np.hstack([ar1_haps(rng, np.tile(freqs[c][g], (2 * n_pop, 1)), RHO)
                       for c in range(len(CHROMS))])
        blocks.append(h[0::2] + h[1::2])
    G = np.vstack(blocks)
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(CHROMS)])
    keep = G.var(0) > 0
    return G[:, keep], chrom[keep]


def windowed_ldscore(X, half=20):
    """ell_j^local = sum_{|k-j|<=half} r_jk^2 (in-sample; no reference panel)."""
    N, M = X.shape
    R = (X.T @ X) / N
    p = np.arange(M)
    mask = (np.abs(p[:, None] - p[None, :]) <= half).astype(float)
    return ((R ** 2) * mask).sum(1)


def fit3(c2, Lc, z2):
    """t of the interaction coefficient in the orthogonalized 3-regressor fit."""
    inter = c2 * Lc
    A = np.column_stack([np.ones_like(c2), c2, inter])
    coef, *_ = np.linalg.lstsq(A, z2, rcond=None)
    resid = z2 - A @ coef
    dof = max(1, len(z2) - 3)
    s2 = (resid @ resid) / dof
    cov = np.linalg.pinv(A.T @ A) * s2
    return coef[2] / np.sqrt(np.diag(cov))[2]


def run_arm(Xc, direction, c2, Lc, kappa, target_chi2, seed):
    """One amplitude-matched arm; returns (t_interaction, mean_chi2).

    direction=None -> HET arm; otherwise a fixed per-individual direction
    (u1 for ENV, its projection into col(Xc) for DIR)."""
    N, Mc = Xc.shape
    rloc = np.random.default_rng(np.random.PCG64(seed))
    if direction is None:
        w = (c2 / c2.mean()) ** kappa
        B = rloc.standard_normal((Mc, REPS_PHENO)) * np.sqrt(w)[:, None]
        Gsig = Xc @ B
        Gsig = Gsig / Gsig.std(0, keepdims=True)
    else:
        Gsig = (direction / direction.std())[:, None] * np.ones((1, REPS_PHENO))
    Ei = rloc.standard_normal((N, REPS_PHENO))
    Zsig = (Xc.T @ Gsig) / np.sqrt(N)
    Znoise = (Xc.T @ Ei) / np.sqrt(N)
    a = np.sqrt(max(target_chi2 - np.mean(Znoise ** 2), 0.0) / np.mean(Zsig ** 2))
    z2 = np.mean((a * Zsig + Znoise) ** 2, axis=1)
    return fit3(c2, Lc, z2), float(np.mean(z2))


def run_setting(label, npops, n_pop, kappa, target_chi2, n_real):
    t_env, t_dir, t_het, m_env, m_het = [], [], [], [], []
    for r in range(n_real):
        rng = np.random.default_rng(np.random.PCG64(SEED + 1000 + r))
        G, chrom = simulate(rng, npops, n_pop)
        N = G.shape[0]
        Xall = standardize(G)
        Xc = standardize(G[:, chrom == 0])
        K = (Xall @ Xall.T) / Xall.shape[1]
        _, U = np.linalg.eigh(K)
        u1 = U[:, -1]
        c2 = ((Xc.T @ u1) / np.sqrt(N)) ** 2
        Lc = windowed_ldscore(Xc)
        Lc = Lc - Lc.mean()
        g_dir = Xc @ np.linalg.lstsq(Xc, u1, rcond=None)[0]
        t, m = run_arm(Xc, u1, c2, Lc, kappa, target_chi2, SEED + 7)
        t_env.append(t); m_env.append(m)
        t, _ = run_arm(Xc, g_dir, c2, Lc, kappa, target_chi2, SEED + 7)
        t_dir.append(t)
        t, m = run_arm(Xc, None, c2, Lc, kappa, target_chi2, SEED + 7)
        t_het.append(t); m_het.append(m)
    t_env, t_dir, t_het = map(np.array, (t_env, t_dir, t_het))
    thr = np.quantile(t_env, 0.95)
    power = float(np.mean(t_het > thr))
    power_t2 = float(np.mean(t_het > 2.0))
    type1_t2 = float(np.mean(t_env > 2.0))
    print(f"  {label:26s} medENV {np.median(t_env):6.2f}  medDIR {np.median(t_dir):6.2f}"
          f"  medHET {np.median(t_het):6.2f}  thr95 {thr:6.2f}"
          f"  power {power:5.2f}  (t>2: env {type1_t2:4.2f}, het {power_t2:4.2f})"
          f"  meanChi2 {np.mean(m_env):.2f}/{np.mean(m_het):.2f}")
    return dict(label=label, thr=thr, power=power, med_het=float(np.median(t_het)))


def main():
    print("Local-adaptation channel kill-test (METHOD.md sec.14 scope boundary)")
    print(f"seed={SEED} bitgen=PCG64  chroms={CHROMS}  fst={FST}  rho={RHO}"
          f"  pheno reps={REPS_PHENO}")
    print("arms: ENV (null, rank-1 axis shift) / DIR (aliased genetic control) /"
          " HET (local-adaptation architecture)\n")

    results = []
    # gradient-strength sweep at moderate signal (N=600)
    for kappa in (0.5, 1.0, 2.0):
        results.append(run_setting(f"N600 chi2=2.5 kappa={kappa}", 3, 200,
                                   kappa, 2.5, 24))
    # signal-strength sweep at kappa=1 (N=600)
    for tc in (1.5, 4.0):
        results.append(run_setting(f"N600 chi2={tc} kappa=1", 3, 200, 1.0, tc, 24))
    # larger cohort at the moderate setting
    results.append(run_setting("N1200 chi2=2.5 kappa=1", 4, 300, 1.0, 2.5, 16))
    # strongest setting (the KILL criterion arm)
    results.append(run_setting("N600 chi2=4.0 kappa=2", 3, 200, 2.0, 4.0, 24))

    moderate = next(r for r in results if r["label"] == "N600 chi2=2.5 kappa=1.0")
    strongest = results[-1]
    print("\n" + "=" * 74)
    print("PRE-REGISTERED VERDICT RULE: PASS = power>=0.80 at N600/chi2=2.5/kappa=1;")
    print("KILL = power<0.50 at chi2=4.0/kappa=2; else MARGINAL.")
    if moderate["power"] >= 0.80:
        verdict = "PASS"
    elif strongest["power"] < 0.50:
        verdict = "KILL"
    else:
        verdict = "MARGINAL"
    print(f"moderate power = {moderate['power']:.2f}   strongest power = "
          f"{strongest['power']:.2f}   ->  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
