"""Within-family kill-test: does the exact within-family P_c beat the standard
sib-GWAS pipeline? Verdict: NO (for standard full-sib designs). This is the
honest, fair version -- a prior draft used a rigged metric (see NOTE below).

BACKGROUND. The one candidate niche left for PLDSC after the summary-statistics
angle collapsed (relatedness does not reach the windowed LD-score regression;
see map_ldsc_miscalibration.py and METHOD sec.13) was within-family / sib-GWAS,
where relatedness is exactly known by pedigree so P_c is exact with no fitted-V_c
fragility. The claim to test: the within-family genetic covariance that survives
family-mean-centering (M_B = I - B(B^T B)^{-1}B^T) inflates the standard sib-GWAS
LD-score regression, and P_c removes it.

WHY THE CLAIM FAILS (the crux, proved here in Part 1). For full sibs the expected
relatedness is R = 0.5 (I + B) with B the family block-of-ones. The FE projector
annihilates the family span (M_B B = 0), so

        M_B R M_B = 0.5 * M_B     (EXACTLY).

The within-family genetic covariance that survives mean-centering is therefore
PROPORTIONAL to M_B: a pure homoscedastic scale (0.5 sigma_g^2), not
LD-correlated confounding. A real within-sibship GWAS (Howe et al. 2022) reports
sibship-clustered / cluster-robust standard errors, which consistently estimate
the true score variance x_j^T M_B V_true M_B x_j and hence divide this scale out
-- the null chi^2 is 1 by construction, with NO modelling of V. There is nothing
left for P_c to correct.

NOTE ON A RIGGED METRIC. Standardizing by x_j^T M_B x_j (the OLS/V=I variance)
instead of the cluster-robust variance makes FE's null chi^2 read
sigma_e^2 + 0.5 sigma_g^2 (~0.55-0.65) and P_c read 1.0 -- a 'gap' that is pure
normalization, zero confounding. That is the analytic P2-identity metric
EXTERNAL_REVIEW.md 5.1 flagged as circular; do not use it as the headline.

Part 2 confirms, by Monte Carlo with each estimator given a consistent (best-case
cluster-robust-equivalent) SE, that FE and P_c BOTH recover a calibrated windowed
LD-score regression -- equivalent, so P_c adds nothing. Part 3 shows the tiny
structured remainder of the surviving covariance is orthogonal to the windowed
LD score (Wall-1), so it cannot bias the slope either.

Run:  conda run -n pldsc python experiments/within_family_killtest.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import residual_precision, score_pieces, ld_scores  # noqa: E402
from pldsc.simulate import (  # noqa: E402
    simulate_admixture,
    simulate_relatedness,
    standardize,
)

SEED = 20260718


def fam_block(n_fam, sibs):
    B = np.zeros((n_fam * sibs, n_fam))
    for f in range(n_fam):
        B[f * sibs:(f + 1) * sibs, f] = 1.0
    return B


def fe_projector(B):
    C = np.column_stack([np.ones(B.shape[0]), B])
    return np.eye(B.shape[0]) - C @ np.linalg.pinv(C.T @ C) @ C.T


def windowed_mask(m, half=15):
    p = np.arange(m)
    return (np.abs(p[:, None] - p[None, :]) <= half).astype(float)


def fit(x, y):
    A = np.column_stack([np.ones_like(x), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c[0], c[1]


def build(kind, rng, cs, n_fam, sibs, fst=0.15):
    if kind == "homog":
        G, chrom, _ = simulate_relatedness(rng, cs, n_fam=n_fam, sibs=sibs, rho=0.6)
    else:
        G, chrom, _ = simulate_admixture(rng, cs, n_indiv=n_fam * sibs, fst=fst,
                                         rho=0.6, families=True, sibs=sibs)
    keep = G.var(0) > 0
    return G[:, keep], chrom[keep]


# --------------------------------------------------------------------------- #
# Part 1: the surviving within-family covariance is a pure scale (M_B R M_B = c M_B)
# --------------------------------------------------------------------------- #
def part1_scale_not_confounding(rng, cs):
    print("=" * 78)
    print("PART 1  Surviving within-family covariance is a SCALE, not confounding")
    print("        M_B (sigma_g^2 K_-c) M_B  ~=  c * M_B   =>  FE's own SE divides")
    print("        it out (null chi^2 -> 1) with no modelling of V.")
    print("-" * 78)
    print(f"  {'scenario':32s} {'c_hat/sig_g2':>12s} {'||resid||/||M_B R M_B||':>24s}")
    for kind, nf, s in [("homog", 120, 2), ("homog", 60, 4),
                        ("admix", 120, 2), ("admix", 60, 4)]:
        G, chrom = build(kind, rng, cs, nf, s)
        N = G.shape[0]
        Zoff = standardize(G[:, chrom != 0])
        K = (Zoff @ Zoff.T) / Zoff.shape[1]
        B = fam_block(nf, s)
        M = fe_projector(B)
        MRM = M @ K @ M
        # best scale c minimising ||MRM - c M||_F  => c = <MRM,M>/<M,M>
        c = np.sum(MRM * M) / np.sum(M * M)
        resid = MRM - c * M
        frac = np.linalg.norm(resid) / (np.linalg.norm(MRM) + 1e-300)
        print(f"  {kind+f'  {nf}x{s} sibs':32s} {c:12.3f} {frac:24.3f}")
    print("  Reading: c_hat ~= 0.5 (the expected full-sib relatedness). The scale")
    print("  is what FE's residual-variance / cluster-robust SE absorbs. The")
    print("  Frobenius remainder is large ONLY because it is off-diagonal,")
    print("  cross-individual noise over many off-target markers (Part 3 shows it")
    print("  is orthogonal to the windowed LD score, so it never reaches a slope).")


# --------------------------------------------------------------------------- #
# Part 2: with a consistent SE, FE and P_c are BOTH calibrated  (MC)
# --------------------------------------------------------------------------- #
def consistent_z(Op, X, y, V_true):
    """z_j = (x_j^T Op y) / sqrt(x_j^T Op V_true Op x_j).

    The denominator is the exact score variance that a cluster-robust / effective-N
    SE consistently estimates -- so this is best-case (oracle-consistent-SE) for
    EVERY operator. Under the null E[z_j^2]=1 for all of them by construction;
    the question Part 2 asks is only whether the SLOPE is recovered unbiased.
    """
    OpX = Op @ X
    num = OpX.T @ y                                   # (M, reps)
    var = np.einsum("ij,ij->j", OpX, V_true @ OpX)    # (M,)
    return num / np.sqrt(var)[:, None]


def part2_slope_calibration(rng, cs):
    print("\n" + "=" * 78)
    print("PART 2  With a consistent SE (what cluster-robust/effective-N achieves),")
    print("        FE and PLDSC are BOTH calibrated. MC: regress z^2 on the operator's")
    print("        own windowed LD score -> (intercept, slope) vs true tau^2.")
    print("-" * 78)
    reps = 8000
    h2_target = 0.4
    print(f"  {'scenario':22s} {'method':16s} {'intercept':>9s} {'slope':>9s} {'tau^2':>10s}")
    for kind, nf, s, sg2, sc2, se2 in [("homog", 120, 2, 0.7, 0.0, 0.3),
                                       ("homog", 60, 4, 0.5, 0.2, 0.3),
                                       ("admix", 60, 4, 0.5, 0.2, 0.3)]:
        G, chrom = build(kind, rng, cs, nf, s)
        N = G.shape[0]
        tgt = chrom == 0
        Xc = standardize(G[:, tgt])
        Mc = Xc.shape[1]
        Zoff = standardize(G[:, ~tgt])
        K = (Zoff @ Zoff.T) / Zoff.shape[1]
        B = fam_block(nf, s)
        Bn = B @ B.T
        V_true = sg2 * K + sc2 * Bn + se2 * np.eye(N)
        wm = windowed_mask(Mc)

        # operators
        M_fe = fe_projector(B)
        Kped = np.zeros((N, N))
        for f in range(nf):
            Kped[f * s:(f + 1) * s, f * s:(f + 1) * s] = 0.5
        np.fill_diagonal(Kped, 1.0)
        one = np.ones((N, 1))
        P_ped = residual_precision(sg2 * Kped + sc2 * Bn + se2 * np.eye(N), one)
        P_or = residual_precision(V_true, one)

        # true per-SNP tau^2 from the oracle full LD score so h2 ~ target
        _, q_or = score_pieces(P_or, Xc)
        L_or_full = ld_scores(P_or, Xc)
        tau2 = h2_target / np.mean(L_or_full)

        Lc = np.linalg.cholesky(V_true)
        Beta = rng.standard_normal((Mc, reps)) * np.sqrt(tau2)
        Y = Xc @ Beta + Lc @ rng.standard_normal((N, reps))

        for name, Op in [("within-sib FE", M_fe), ("PLDSC pedigree", P_ped),
                         ("PLDSC oracle", P_or)]:
            z = consistent_z(Op, Xc, Y, V_true)
            z2 = np.mean(z ** 2, axis=1)
            L = ld_scores(Op, Xc, wm)
            a, b = fit(L, z2)
            print(f"  {kind+f' {nf}x{s}':22s} {name:16s} {a:9.3f} {b:9.4f} {tau2:10.4f}")
        print()
    print("  Reading: all three sit at intercept ~1 and slope ~tau^2. FE (the")
    print("  standard tool, given its normal cluster-robust SE) is already")
    print("  calibrated; PLDSC does not beat it. Pedigree-P_c matches only when its")
    print("  V approximates the truth; under admixture it is the most fragile.")


# --------------------------------------------------------------------------- #
# Part 3: the surviving structured remainder is orthogonal to the windowed score
# --------------------------------------------------------------------------- #
def part3_wall1(rng, cs):
    print("=" * 78)
    print("PART 3  Wall-1: even the structured part of the surviving covariance is")
    print("        orthogonal to the windowed LD score, so it cannot bias the slope.")
    print("-" * 78)
    print(f"  {'scenario':22s} {'corr(FE leverage, windowed L_FE)':>34s}")
    for kind, nf, s in [("homog", 120, 2), ("homog", 60, 4), ("admix", 60, 4)]:
        G, chrom = build(kind, rng, cs, nf, s)
        tgt = chrom == 0
        Xc = standardize(G[:, tgt])
        Zoff = standardize(G[:, ~tgt])
        K = (Zoff @ Zoff.T) / Zoff.shape[1]
        B = fam_block(nf, s)
        M = fe_projector(B)
        # per-SNP surviving off-target genetic leverage x_j^T M K M x_j
        MX = M @ Xc
        lev = np.einsum("ij,ij->j", MX, K @ MX) / np.einsum("ij,ij->j", Xc, MX)
        L = ld_scores(M, Xc, windowed_mask(Xc.shape[1]))
        r = np.corrcoef(lev, L)[0, 1]
        print(f"  {kind+f' {nf}x{s}':22s} {r:34.3f}")
    print("  Reading: ~0. The surviving covariance is off-target (LOCO); its")
    print("  realized fluctuation is independent of within-target-chromosome LD")
    print("  (independent assortment across chromosomes), so it moves at most a")
    print("  genome-wide intercept a consistent SE already absorbs -- never the")
    print("  windowed slope. Same structural reason the summary-stat angle died.")


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    cs = [80, 80, 80]
    print("PLDSC WITHIN-FAMILY KILL-TEST  (fair version)")
    print(f"seed={SEED} bitgen=PCG64\n")
    part1_scale_not_confounding(rng, cs)
    part2_slope_calibration(rng, cs)
    part3_wall1(rng, cs)
    print("\n" + "=" * 78)
    print("VERDICT: KILLED for standard full-sib designs. The within-family genetic")
    print("covariance surviving family-mean-centering is a homoscedastic SCALE the")
    print("standard sib-GWAS SE (cluster-robust / effective-N, Howe et al. 2022)")
    print("already removes; its structured remainder is orthogonal to the windowed")
    print("LD score. PLDSC's exact P_c reproduces the standard result, it does not")
    print("beat it, and pedigree-only P_c is fragile under admixture.")
    print("The ONLY unruled-out regime: HETEROGENEOUS pedigrees (mixed full/half")
    print("sibs, cousins) where M_B K M_B is NOT proportional to M_B, evaluated on a")
    print("genome-wide (not windowed) intercept -- but that estimand is already")
    print("owned by Sib-Regression (Visscher 2006) and RDR (Young 2018), and needs")
    print("a proper per-meiosis-recombination simulator this repo does not have.")


if __name__ == "__main__":
    main()
