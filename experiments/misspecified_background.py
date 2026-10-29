"""Non-circular robustness test: where does the theory bend?

Every other experiment generates data under exactly the model Eq (3) assumes
(epsilon ~ (0, V_c) with V_c = sigma_g^2 K_all + sigma_e^2 I, i.e. an
INFINITESIMAL off-chromosome background). METHOD.md Section 6 flags the axis that
matters most and implements it nowhere: what if the off-chromosome architecture
is NOT infinitesimal -- sparse, or directional along the structure axis? Then the
true Cov(epsilon) != V_c, the P2 identity P_c V_c P_c = P_c no longer applies to
the true covariance, and the null intercept moves off 1.

This measures that bend EXACTLY (no Monte Carlo): the null-target-chromosome
intercept is mean_j x_j^T P_c V_true P_c x_j / (x_j^T P_c x_j) (= null_inflation),
where P_c is built from the ASSUMED V_c (all-SNP LOCO GRM) but V_true is the
covariance implied by a misspecified background. Infinitesimal (f=1.0, alpha=0) is
the control and must return 1.000.

Run:  python experiments/misspecified_background.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import null_inflation, residual_precision  # noqa: E402
from pldsc.simulate import ar1_haps, standardize  # noqa: E402

SEED = 20260718
CHROMS = [150, 600]          # target chr 0; chr 1 is the off-chromosome background
N_POP = 150
N_POPS = 3
FST = 0.10
RHO = 0.6
SG2, SE2 = 0.7, 0.3


def multipop_freqs(rng, m, fst, npops):
    anc = rng.uniform(0.1, 0.9, m)
    a, b = anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst
    return [rng.beta(a, b) for _ in range(npops)]


def simulate(rng):
    freqs = [multipop_freqs(rng, m, FST, N_POPS) for m in CHROMS]
    blocks = []
    for g in range(N_POPS):
        h = np.hstack([ar1_haps(rng, np.tile(freqs[c][g], (2 * N_POP, 1)), RHO)
                       for c in range(len(CHROMS))])
        blocks.append(h[0::2] + h[1::2])
    G = np.vstack(blocks)
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(CHROMS)])
    keep = G.var(0) > 0
    return G[:, keep], chrom[keep]


def main():
    print("Non-circular robustness: null intercept under a MISSPECIFIED off-chr background")
    print(f"seed={SEED} {N_POPS} pops x {N_POP}  fst={FST}  sigma_g2={SG2}  reps=6")
    print("intercept target = 1.000 exactly under the assumed (infinitesimal) model")
    reps = 6

    sparse_fracs = [1.0, 0.10, 0.02, 0.005]      # 1.0 = infinitesimal control
    dir_alphas = [0.0, 0.5, 1.0, 2.0]            # extra variance along top structure axis
    sp_rows = {f: [] for f in sparse_fracs}
    dir_rows = {a: [] for a in dir_alphas}

    for r in range(reps):
        rng = np.random.default_rng(np.random.PCG64(SEED + 500 + r))
        G, chrom = simulate(rng)
        off = chrom == 1
        tgt = chrom == 0
        n = G.shape[0]
        Xc = standardize(G[:, tgt])
        Zoff = standardize(G[:, off])
        m_off = Zoff.shape[1]
        K_all = (Zoff @ Zoff.T) / m_off
        Vc = SG2 * K_all + SE2 * np.eye(n)
        P_c = residual_precision(Vc, np.ones((n, 1)))

        # sparse background: true genetic covariance uses only a causal subset
        for f in sparse_fracs:
            n_c = max(2, int(round(f * m_off)))
            idx = rng.choice(m_off, size=n_c, replace=False)
            K_sparse = (Zoff[:, idx] @ Zoff[:, idx].T) / n_c
            V_true = SG2 * K_sparse + SE2 * np.eye(n)
            sp_rows[f].append(float(np.mean(null_inflation(P_c, Xc, V_true))))

        # directional background: extra variance along the top structure axis
        w, V = np.linalg.eigh(K_all)
        lam1, u1 = w[-1], V[:, -1]
        for a in dir_alphas:
            V_true = Vc + a * SG2 * lam1 * np.outer(u1, u1)
            dir_rows[a].append(float(np.mean(null_inflation(P_c, Xc, V_true))))

    print("\n  SPARSE off-chr background (fraction of off-chr SNPs causal):")
    print("    " + "  ".join(f"f={f:<6g}" for f in sparse_fracs))
    print("    " + "  ".join(f"{np.mean(sp_rows[f]):<8.3f}" for f in sparse_fracs)
          + "   <- mean null intercept")

    print("\n  DIRECTIONAL off-chr background (extra variance = alpha x polygenic")
    print("  variance along the top structure axis):")
    print("    " + "  ".join(f"a={a:<6g}" for a in dir_alphas))
    print("    " + "  ".join(f"{np.mean(dir_rows[a]):<8.3f}" for a in dir_alphas)
          + "   <- mean null intercept")

    print("\n" + "=" * 74)
    print("Read: control (f=1.0, alpha=0) must be 1.000 (P2 identity). If the")
    print("intercept moves as the background becomes sparse or directional, the")
    print("oracle calibration is FRAGILE to off-chromosome model misspecification")
    print("-- the failure axis METHOD Section 6 names. This is the one non-circular")
    print("result: it violates epsilon ~ (0, V_c) rather than assuming it.")


if __name__ == "__main__":
    main()
