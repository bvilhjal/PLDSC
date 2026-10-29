"""Orthogonal misspecification: a true-covariance axis ABSENT from the GRM.

misspecified_background.py tested sparse and directional off-chromosome
backgrounds and found the oracle intercept robust (<1% movement) -- but both
misspecifications live mostly in the low-rank structure subspace P_c already
whitens. The harder case it explicitly flagged as untested: a true-covariance
axis that is ABSENT from (orthogonal to) the GRM used to build V_c. That is the
case where V_c^{-1} does not down-weight the offending direction and the P2
identity has nothing to say.

Three arms (exact, no Monte Carlo; null-target-chromosome intercept =
mean_j null_inflation(P_assumed, X, V_true); 6 realizations):

  (a) BATCH axis orthogonal to ancestry: V_true = V_c + v u u^T with u a
      centered within-population batch indicator (shared microenvironment).
      By construction u has ~zero loading on K's structure axes -- verified
      numerically (projection onto the top-10 eigenvectors of K is printed).
  (b) PRUNED-GRM estimation error: P_c built from sigma_g2 K_pruned + sigma_e2 I
      with K_pruned from every 10th off-chromosome SNP (standard LD-thinning
      practice); truth uses K_full. This is not an exotic axis but the
      realistic way a covariance component ends up "absent from the GRM".
  (c) FAMILY-SHARED ENVIRONMENT + pruned GRM: on an admixed sibship cohort,
      V_true = sigma_g2 K_full + sigma_e2 I + v E_fam (block-diagonal sibship
      covariance); P_c built from the pruned GRM, which captures the family
      axes only weakly. Within-arm oracle (P_c from V_true) must return 1.000.

Run:  python experiments/misspecified_orthogonal.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import null_inflation, residual_precision  # noqa: E402
from pldsc.simulate import ar1_haps, simulate_admixture, standardize  # noqa: E402

SEED = 20260719
CHROMS = [150, 600]         # chr0 = target; chr1 = off-chromosome background
N_POP = 150
N_POPS = 3
FST = 0.10
RHO = 0.6
SG2, SE2 = 0.7, 0.3
V_SWING = [0.25, 0.5, 1.0, 2.0]     # v / sigma_e2
PRUNE = 10                          # keep every 10th off-chr SNP in K_pruned
REPS = 6


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


def intercept(P, Xc, V_true):
    return float(np.mean(null_inflation(P, Xc, V_true)))


def main():
    print("Orthogonal misspecification: axes absent from the GRM")
    print(f"seed={SEED} bitgen=PCG64  reps={REPS}  sigma_g2={SG2} sigma_e2={SE2}")

    # ---- arms (a) and (b): discrete multipop cohort -------------------------
    a_rows = {v: [] for v in V_SWING}
    b_rows = []
    proj_report = []
    for r in range(REPS):
        rng = np.random.default_rng(np.random.PCG64(SEED + 500 + r))
        G, chrom = simulate(rng)
        off, tgt = chrom == 1, chrom == 0
        n = G.shape[0]
        Xc = standardize(G[:, tgt])
        Zoff = standardize(G[:, off])
        K_full = (Zoff @ Zoff.T) / Zoff.shape[1]
        Vc = SG2 * K_full + SE2 * np.eye(n)
        one = np.ones((n, 1))

        # (a) batch axis: random half within each population, centered
        batch = np.zeros(n)
        for g in range(N_POPS):
            lo = g * N_POP
            batch[lo:lo + N_POP // 2] = 1.0
        u = batch - batch.mean()
        u = u / np.linalg.norm(u)
        w, U = np.linalg.eigh(K_full)
        proj = float(np.sum((U[:, -10:].T @ u) ** 2))   # ||proj of u onto top10||
        proj_report.append(proj)
        P = residual_precision(Vc, one)
        for v in V_SWING:
            V_true = Vc + v * SE2 * np.outer(u, u)
            a_rows[v].append(intercept(P, Xc, V_true))

        # (b) pruned GRM
        Zpr = Zoff[:, ::PRUNE]
        K_pr = (Zpr @ Zpr.T) / Zpr.shape[1]
        P_pr = residual_precision(SG2 * K_pr + SE2 * np.eye(n), one)
        b_rows.append(intercept(P_pr, Xc, Vc))

    print("\n(a) BATCH axis (shared microenvironment), projection of u onto the")
    print(f"    top-10 GRM eigenvectors = {np.mean(proj_report):.4f} "
          "(~0 => axis absent from the GRM)")
    print("    v/sigma_e2:      " + "  ".join(f"{v:<8g}" for v in V_SWING))
    print("    mean intercept:  " + "  ".join(f"{np.mean(a_rows[v]):<8.3f}"
                                               for v in V_SWING))

    print("\n(b) PRUNED-GRM estimation error (K from every "
          f"{PRUNE}th off-chr SNP; truth uses the full GRM):")
    print(f"    mean intercept:  {np.mean(b_rows):.3f}   (control with full GRM"
          " would be 1.000)")

    # ---- arm (c): admixed sibships + family-shared environment --------------
    c_oracle, c_full, c_pruned = [], {v: [] for v in V_SWING}, {v: [] for v in V_SWING}
    for r in range(REPS):
        rng = np.random.default_rng(np.random.PCG64(SEED + 9000 + r))
        G, chrom, _ = simulate_admixture(rng, CHROMS, n_indiv=400, fst=0.15,
                                         rho=RHO, families=True, sibs=4)
        keep = G.var(0) > 0
        G, chrom = G[:, keep], chrom[keep]
        off, tgt = chrom == 1, chrom == 0
        n = G.shape[0]
        Xc = standardize(G[:, tgt])
        Zoff = standardize(G[:, off])
        K_full = (Zoff @ Zoff.T) / Zoff.shape[1]
        Vc = SG2 * K_full + SE2 * np.eye(n)
        one = np.ones((n, 1))

        n_fam = n // 4
        F = np.zeros((n, n_fam))
        F[np.arange(n), np.arange(n) // 4] = 1.0
        E_fam = F @ F.T

        Zpr = Zoff[:, ::PRUNE]
        K_pr = (Zpr @ Zpr.T) / Zpr.shape[1]
        P_pr = residual_precision(SG2 * K_pr + SE2 * np.eye(n), one)
        P_full = residual_precision(Vc, one)
        for v in V_SWING:
            V_true = Vc + v * SE2 * E_fam
            c_oracle.append(intercept(residual_precision(V_true, one), Xc, V_true))
            c_full[v].append(intercept(P_full, Xc, V_true))
            c_pruned[v].append(intercept(P_pr, Xc, V_true))

    print("\n(c) FAMILY-SHARED environment on admixed sibships (100 sibships x4):")
    print(f"    oracle (P_c from V_true):            {np.mean(c_oracle):.3f}"
          "   <- must be 1.000")
    print("    v/sigma_e2:      " + "  ".join(f"{v:<8g}" for v in V_SWING))
    print("    full GRM:        " + "  ".join(f"{np.mean(c_full[v]):<8.3f}"
                                               for v in V_SWING))
    print("    pruned GRM:      " + "  ".join(f"{np.mean(c_pruned[v]):<8.3f}"
                                               for v in V_SWING))

    print("\n" + "=" * 74)
    print("Read: (a) tests an axis orthogonal to the ancestry subspace; (b),(c)")
    print("test the realistic route to 'absent from the GRM' -- thinning the SNPs")
    print("the GRM is built from. Compare with misspecified_background.py, where")
    print("in-subspace misspecification moved the intercept <1%.")


if __name__ == "__main__":
    main()
