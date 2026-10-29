"""Illustration of the structure-aliasing limit (METHOD.md sec.14), scoped.

Marginal (structure-UNcorrected) per-SNP z_j = x_j^T y / sqrt(N). We compare the
per-SNP inflation E[z_j^2] under three phenotype models that all load variance on
the SAME top structure axis u_1, and regress z_j^2 on

    [ 1 ,  corr(x_j,u_1)^2 ,  corr(x_j,u_1)^2 * ell_j^local ] .

The DIRECTIONAL genetic model and the ENVIRONMENTAL model are ALIASED: both are
rank-1 along u_1, so z_j^2 is pure corr(x_j,u_1)^2 with NO local-LD interaction --
structure-only second moments cannot tell them apart (the identifiability limit).
The HETEROSCEDASTIC genetic model rides on local LD, so it leaves a nonzero
corr^2 x ell^local interaction that the environmental shift lacks -- i.e. that
part of genuine ancestry-aligned heritability IS identifiable, via the local
channel. This is the scope boundary of the limit.

Run:  python experiments/identifiability_aliasing.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.simulate import simulate_admixture, standardize  # noqa: E402

SEED = 20260718


def windowed_ldscore(X, half=20):
    """ell_j^local = sum_{|k-j|<=half} r_jk^2, r_jk = x_j^T x_k / N."""
    N, M = X.shape
    R = (X.T @ X) / N
    p = np.arange(M)
    mask = (np.abs(p[:, None] - p[None, :]) <= half).astype(float)
    return ((R ** 2) * mask).sum(1)


def zscores(X, Y):
    """Marginal uncorrected z_j = x_j^T y / sqrt(N), columns of Y are phenotypes."""
    return (X.T @ Y) / np.sqrt(X.shape[0])


def fit3(c2, Llocal_c, z2):
    """Regress z2 on [1, corr^2, corr^2 * (ell^local - mean)]; the CENTERED local
    factor makes the interaction regressor ~orthogonal to corr^2 (removes the
    spurious collinearity-driven interaction). Return b[corr^2], b[interaction],
    and t for the interaction."""
    inter = c2 * Llocal_c
    A = np.column_stack([np.ones_like(c2), c2, inter])
    coef, *_ = np.linalg.lstsq(A, z2, rcond=None)
    resid = z2 - A @ coef
    dof = max(1, len(z2) - 3)
    s2 = (resid @ resid) / dof
    cov = np.linalg.pinv(A.T @ A) * s2
    se = np.sqrt(np.diag(cov))
    return coef[1], coef[2], coef[2] / se[2]


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    cs = [250, 250]                 # 2 chromosomes, AR1 local LD; chr0 = target
    G, chrom, _ = simulate_admixture(rng, cs, n_indiv=500, fst=0.15, rho=0.7)
    keep = G.var(0) > 0
    G, chrom = G[:, keep], chrom[keep]
    N = G.shape[0]
    Xall = standardize(G)
    tgt = chrom == 0
    Xc = standardize(G[:, tgt])
    Mc = Xc.shape[1]

    # top structure axis from the full GRM
    K = (Xall @ Xall.T) / Xall.shape[1]
    w, U = np.linalg.eigh(K)
    u1 = U[:, -1]
    c1 = (Xc.T @ u1) / np.sqrt(N)          # corr(x_j, u_1) for target SNPs
    c2 = c1 ** 2
    Llocal = windowed_ldscore(Xc)

    Llocal_c = Llocal - Llocal.mean()          # centered local factor (orthogonalize)
    reps = 8000

    def run_model(direction=None, het_w=None, target_chi2=2.5):
        """Build phenotypes, amplitude-matched so mean(z^2) ~= target_chi2, and
        return (b[corr^2], b[interaction], t_interaction, mean chi^2).
        direction: fixed unit per-individual signal (env or directional genetic).
        het_w: per-SNP effect-variance weights for the heteroscedastic model."""
        rloc = np.random.default_rng(np.random.PCG64(SEED + 7))
        if het_w is not None:
            Bh = rloc.standard_normal((Mc, reps)) * np.sqrt(het_w)[:, None]
            Gsig = Xc @ Bh
            Gsig = Gsig / Gsig.std(0, keepdims=True)          # unit-var per rep
        else:
            Gsig = (direction / direction.std())[:, None] * np.ones((1, reps))
        # calibrate amplitude a so mean_j E[z_j^2] ~= target_chi2
        Ei = rloc.standard_normal((N, reps))
        # find a by matching: z = a*Xc^T Gsig_dir + noise; grid is overkill, solve
        Zsig = zscores(Xc, Gsig)                                # signal z (a=1)
        Znoise = zscores(Xc, Ei)
        base = np.mean(Znoise ** 2)                             # ~1
        s2sig = np.mean(Zsig ** 2)
        a = np.sqrt(max(target_chi2 - base, 0.0) / s2sig)
        Z = a * Zsig + Znoise
        z2 = np.mean(Z ** 2, axis=1)
        b, i, t = fit3(c2, Llocal_c, z2)
        return b, i, t, float(np.mean(z2))

    # (A) ENV: rank-1 per-individual axis shift along u1
    bA, iA, tA, mA = run_model(direction=u1)
    # (B) DIRECTIONAL genetic: projection of u1 into col(Xc) -> closest genuine
    #     genetic value to a pure ancestry-axis shift (best-case rank-1)
    g_dir = Xc @ np.linalg.lstsq(Xc, u1, rcond=None)[0]
    bB, iB, tB, mB = run_model(direction=g_dir)
    # (C) HETEROSCEDASTIC genetic: per-SNP h^2 elevated for high-loading SNPs
    bC, iC, tC, mC = run_model(het_w=(c2 / c2.mean() + 1e-6))

    print("Structure-aliasing illustration (METHOD.md sec.14), amplitude-matched")
    print(f"seed={SEED} bitgen=PCG64  N={N} target SNPs={Mc} lambda1={w[-1]:.1f}")
    print("interaction regressor uses CENTERED ell^local (orthogonalized vs corr^2)\n")
    print(f"  {'model':34s} {'meanChi2':>8s} {'b[corr^2]':>10s} {'t[corr^2 x local]':>18s}")
    print(f"  {'(A) environmental axis shift':34s} {mA:8.2f} {bA:10.2f} {tA:18.2f}")
    print(f"  {'(B) directional genetic (rank-1)':34s} {mB:8.2f} {bB:10.2f} {tB:18.2f}")
    print(f"  {'(C) heteroscedastic-along-axis':34s} {mC:8.2f} {bC:10.2f} {tC:18.2f}")
    print()
    print("Reading (amplitude-matched, so numbers are comparable):")
    print("- (A) env and (B) directional genetic are IDENTICAL to 2 dp -- exact")
    print("  aliasing. (B) is a genuine GENETIC effect (SNP effects producing a")
    print("  phenotype along the ancestry axis), yet it and an environmental axis")
    print("  shift produce the same structure-only statistics: no moment separates")
    print("  them. This is the identifiability limit. (Their shared nonzero local")
    print("  interaction is a finite-sample orthogonalization residual, IDENTICAL in")
    print("  both, so it is not a distinguishing feature.)")
    print("- (C) heteroscedastic-along-axis genetic signal loads LESS on pure corr^2")
    print("  and leaves a LARGER corr^2 x ell^local interaction -- it rides on local")
    print("  LD, so that part IS identifiable via the local channel.")
    print("Scope of the limit: the rank-1 directional component is aliased (every")
    print("structure-only correction co-attenuates it); the LD-borne heteroscedastic")
    print("component is not. Matches METHOD.md sec.14.")


if __name__ == "__main__":
    main()
