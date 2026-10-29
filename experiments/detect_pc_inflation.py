"""Detecting structure inflation in GWAS and mapping it to specific PCs.

Regress the GWAS per-SNP statistic on per-SNP PC loadings, one coefficient per
PC, to get an inflation *spectrum* rather than an aggregate lambda_GC / LDSC
intercept. Signed version (beta-hat on loadings) = Bhatia et al. 2016; here we
also test the chi^2-scale (squared-loading) magnitude version. LOCO PCs to avoid
proximal contamination.

Three scenarios test both power and the identifiability limit:
  S1  confounding along PC1 and PC3, no PC-aligned true signal  -> should flag 1,3
  S2  genuine causal signal ALIGNED with PC2, NO confounding    -> false flag at 2?
  S3  null: no confounding, no PC-aligned signal                -> should flag none

Detection statistic per PC: joint regression coefficient / block-jackknife SE
(a per-PC z). |z| large => that PC is implicated.

Run:  python experiments/detect_pc_inflation.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.simulate import ar1_haps, grm, standardize  # noqa: E402

SEED = 20260718
CHROMS = [400, 500]          # target chr 0 (GWAS + loadings); chr 1 builds LOCO PCs
N_POP = 150
N_POPS = 4                   # 4 populations -> several structure PCs
FST = 0.10
RHO = 0.6
K_PC = 6                     # PCs to scan
H2_POLY = 0.3


def multipop_freqs(rng, m, fst, npops):
    anc = rng.uniform(0.1, 0.9, m)
    a, b = anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst
    return [rng.beta(a, b) for _ in range(npops)]


def simulate(rng):
    freqs = [multipop_freqs(rng, m, FST, N_POPS) for m in CHROMS]  # [chr][pop]
    blocks = []
    for g in range(N_POPS):
        h = np.hstack([ar1_haps(rng, np.tile(freqs[c][g], (2 * N_POP, 1)), RHO)
                       for c in range(len(CHROMS))])
        blocks.append(h[0::2] + h[1::2])
    G = np.vstack(blocks)
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(CHROMS)])
    keep = G.var(0) > 0
    return G[:, keep], chrom[keep]


def block_jackknife_se(x_cols, y, n_blocks=40):
    """SE of each joint-regression coefficient via delete-a-block jackknife."""
    A = np.column_stack([np.ones_like(y)] + x_cols)
    n = y.size
    full, *_ = np.linalg.lstsq(A, y, rcond=None)
    blocks = np.array_split(np.arange(n), min(n_blocks, n))
    ests = []
    for blk in blocks:
        keep = np.ones(n, bool); keep[blk] = False
        c, *_ = np.linalg.lstsq(A[keep], y[keep], rcond=None)
        ests.append(c)
    ests = np.array(ests)
    B = len(blocks)
    se = np.sqrt((B - 1) / B * np.sum((ests - ests.mean(0)) ** 2, axis=0))
    return full, se


def run(rng, scenario):
    G, chrom = simulate(rng)
    tgt = chrom == 0
    Xc = standardize(G[:, tgt])
    Xn = Xc / np.linalg.norm(Xc, axis=0)
    n, M = Xc.shape

    # LOCO sample PCs (from the other chromosome) and per-SNP loadings on them
    K = grm(G[:, ~tgt])
    w, V = np.linalg.eigh(K)
    pcs = V[:, ::-1][:, :K_PC]                       # sample PC values (n x K)
    pcs = (pcs - pcs.mean(0)) / pcs.std(0)
    loadings = Xn.T @ (pcs / np.sqrt(n))             # (M x K) signed loadings ~ corr

    # build the phenotype
    beta = rng.standard_normal(M)
    if scenario == "S2_aligned_signal":
        beta = beta + 1.5 * loadings[:, 1]           # true causal effects aligned w/ PC2
    y = Xc @ beta
    y = y / y.std() * np.sqrt(H2_POLY)
    if scenario == "S1_confound_1_3":
        y = y + 1.0 * pcs[:, 0] + 1.0 * pcs[:, 2]    # confounding along PC1 and PC3
    y = y + rng.standard_normal(n) * np.sqrt(1 - H2_POLY)

    # ordinary marginal GWAS
    yc = y - y.mean()
    bhat = Xn.T @ yc / (n - 1)                        # marginal effect estimate
    chi2 = (Xn.T @ yc / yc.std()) ** 2

    # per-PC detection: joint regression of the statistic on the loadings
    signed, signed_se = block_jackknife_se([loadings[:, p] for p in range(K_PC)], bhat)
    sq, sq_se = block_jackknife_se([loadings[:, p] ** 2 for p in range(K_PC)], chi2)
    z_signed = signed[1:] / np.where(signed_se[1:] > 0, signed_se[1:], np.inf)
    z_sq = sq[1:] / np.where(sq_se[1:] > 0, sq_se[1:], np.inf)
    return z_signed, z_sq


def main():
    print("Detecting structure inflation and mapping it to specific PCs")
    print(f"seed={SEED} {N_POPS} pops x {N_POP}  fst={FST}  K_PC={K_PC}  reps=15")
    print("per-PC detection z-scores (|z| large => PC implicated)")
    reps = 15
    scenarios = [
        ("S1: confound along PC1 & PC3", "S1_confound_1_3"),
        ("S2: TRUE signal aligned w/ PC2, no confound", "S2_aligned_signal"),
        ("S3: null (no confound, no aligned signal)", "S3_null"),
    ]
    for label, key in scenarios:
        zs = np.array([run(np.random.default_rng(np.random.PCG64(SEED + 300 + r)), key)
                       for r in range(reps)])
        z_signed = zs[:, 0, :].mean(0)
        z_sq = zs[:, 1, :].mean(0)
        print(f"\n{label}")
        print("  PC:            " + " ".join(f"{p+1:>6d}" for p in range(K_PC)))
        print("  signed |z|:    " + " ".join(f"{abs(v):>6.1f}" for v in z_signed))
        print("  chi2-sq |z|:   " + " ".join(f"{abs(v):>6.1f}" for v in z_sq))

    print("\n" + "=" * 74)
    print("Read: S1 should light up PC1 & PC3 (detection works). S2 should light up")
    print("PC2 too -- the identifiability limit: genuine ancestry-aligned signal is")
    print("indistinguishable from confounding. S3 should stay flat (specificity).")


if __name__ == "__main__":
    main()
