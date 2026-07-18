"""The linchpin of the local/non-local LD-score decomposition.

Claim:  ell_j (ordinary LD score) = ell_j^local + ell_j^struct, where
  ell_j^local  = sum_k (r_jk^(P_c))^2   (structure-whitened; linkage/recombination)
  ell_j^struct = ell_j^ord - ell_j^local
and the theory's punchline is that the NON-LOCAL component equals the PC score:
  ell_j^struct  ~  s_j = sum_p (lambda_p/lambda_t) corr(g_j, PC_p)^2   (LOCO).

If ell^struct correlates tightly with the PC score, the decomposition unifies
the LD score, the P_c operator, cov-LDSC, and the 2017 "PC score" under one
identity. Test on a multi-population cohort, windowed and genome-wide.

Run:  conda run -n pldsc python experiments/ld_decomposition.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import residual_precision, score_pieces  # noqa: E402
from pldsc.simulate import ar1_haps, grm, standardize  # noqa: E402

SEED = 20260718
CHROMS = [400, 500]          # target chr 0; chr 1 builds the LOCO GRM / PCs
N_POP = 150
N_POPS = 4
FST = 0.10
RHO = 0.6
WINDOW_CM = 2.0
K_PC = 10                    # PCs used in the PC score
SG2, SE2 = 0.6, 0.4


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


def r2_ldscore(corr_like, mask):
    return (np.clip(corr_like, -1, 1) ** 2 * mask).sum(1)


def main():
    print("LD-score decomposition: is the non-local component the PC score?")
    print(f"seed={SEED} {N_POPS} pops x {N_POP}  fst={FST}  window={WINDOW_CM}cM  reps=10")
    reps = 10
    rows = []
    for r in range(reps):
        rng = np.random.default_rng(np.random.PCG64(SEED + 400 + r))
        G, chrom = simulate(rng)
        tgt = chrom == 0
        Xc = standardize(G[:, tgt])
        Xn = Xc / np.linalg.norm(Xc, axis=0)
        n, M = Xc.shape
        pos = np.arange(M) * 0.15
        winmask = (np.abs(pos[:, None] - pos[None, :]) <= WINDOW_CM).astype(float)
        allmask = np.ones((M, M))

        # ordinary correlations and structure-whitened (P_c, intercept-only) correlations
        Rord = Xn.T @ Xn
        K = grm(G[:, ~tgt])
        Vc = SG2 * K + SE2 * np.eye(n)
        P = residual_precision(Vc, np.ones((n, 1)))
        PXc, q = score_pieces(P, Xc)
        A = Xc.T @ PXc
        Rpc = A / np.sqrt(np.outer(q, q))          # r_jk^(P_c)

        # PC score (LOCO): eigenvalue-weighted squared SNP-PC correlations
        w, V = np.linalg.eigh(K)
        lam = w[::-1][:K_PC]
        pcs = V[:, ::-1][:, :K_PC]
        pcs = (pcs - pcs.mean(0)) / pcs.std(0)
        corr_pc = Xn.T @ (pcs / np.sqrt(n))         # (M x K)
        pc_score = (corr_pc ** 2 * (lam / lam.sum())).sum(1)

        for tag, mask in (("windowed", winmask), ("genome-wide", allmask)):
            ell_ord = r2_ldscore(Rord, mask)
            ell_local = r2_ldscore(Rpc, mask)
            ell_struct = ell_ord - ell_local
            # correlation of the non-local component with the PC score
            cc = np.corrcoef(ell_struct, pc_score)[0, 1]
            frac = ell_struct.mean() / ell_ord.mean()
            rows.append((tag, cc, frac,
                         ell_ord.mean(), ell_local.mean(), ell_struct.mean()))

    print(f"\n{'scope':13s} {'corr(ell_struct, PCscore)':>26s} {'struct/total':>13s}"
          f" {'ell_ord':>9s} {'ell_local':>10s} {'ell_struct':>11s}")
    for tag in ("windowed", "genome-wide"):
        sub = [x for x in rows if x[0] == tag]
        cc = np.mean([x[1] for x in sub])
        frac = np.mean([x[2] for x in sub])
        eo = np.mean([x[3] for x in sub])
        el = np.mean([x[4] for x in sub])
        es = np.mean([x[5] for x in sub])
        print(f"{tag:13s} {cc:>26.3f} {frac:>13.2f} {eo:>9.2f} {el:>10.2f} {es:>11.2f}")

    print("\n" + "=" * 74)
    print("Read: if corr(ell_struct, PC score) is high (esp. genome-wide), the")
    print("non-local part of the LD score IS the PC score -- the decomposition and")
    print("the unification hold. 'struct/total' shows how much of the LD score is")
    print("structure (large genome-wide, small within a local window).")


if __name__ == "__main__":
    main()
