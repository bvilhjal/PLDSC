"""M5: LD-score-regression inference, and the fitted-V_c intercept-SE caveat.

Builds on the pldsc.regression OLS path with a delete-one-SNP jackknife and the
pldsc._lmm oracle. Two arms per replicate on an admixture+relatedness cohort:

  ORACLE arm  z^2 and L from the true P_c   -> probe toy fixed-operator
              coverage of intercept (=1) and slope (=tau^2).
  FITTED arm  z^2 and L from a REML-fitted P_hat_c -> the jackknife SE (which
              holds V_hat_c fixed) MISSES the variance-component uncertainty, so
              the intercept may under-cover at scale. The displayed propagated
              term is illustrative, not a validated correction.

This illustrates the requirement M4 surfaced: a fitted-V_c intercept SE must
propagate variance-component uncertainty. It is not a validation of genomic
blocks, recomputed weights, or a deployable fitted-V_c correction.

Run:  conda run -n pldsc python experiments/m5_inference.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))  # repo root, for `pldsc`
sys.path.insert(0, _HERE)                    # experiments, for the REML helper

from pldsc._lmm import ld_scores, null_inflation, residual_precision  # noqa: E402
from pldsc.regression import ld_score_regression  # noqa: E402
from pldsc.simulate import grm, simulate_admixture, standardize  # noqa: E402
from m4_fitted_vc import reml_fit  # noqa: E402

SEED = 20260718


def coverage(est, truth, se, z=1.96):
    est, se = np.asarray(est), np.asarray(se)
    return float(np.mean(np.abs(est - truth) < z * se))


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    chrom_sizes = [80, 80, 80]
    sg2, se2 = 0.7, 0.3
    h2_true = sg2 / (sg2 + se2)
    R = 300

    G, chrom, pos = simulate_admixture(rng, chrom_sizes, n_indiv=400, fst=0.15,
                                       families=True, sibs=4)
    keep = G.var(0) > 0
    G, chrom = G[:, keep], chrom[keep]
    tgt = chrom == 0
    N = G.shape[0]
    Xc = standardize(G[:, tgt])
    K = grm(G[:, ~tgt])
    Vc = sg2 * K + se2 * np.eye(N)
    one = np.ones((N, 1))
    lam, U = np.linalg.eigh(K)
    Ct = U.T @ one
    Lv = np.linalg.cholesky(Vc)

    P_oracle = residual_precision(Vc, one)
    PX_o, q_o = P_oracle @ Xc, None
    q_o = np.einsum("ij,ij->j", Xc, PX_o)
    L_oracle = ld_scores(P_oracle, Xc)
    tau2 = 0.30 / float(np.mean(L_oracle))

    print("PLDSC M5: inference and the fitted-V_c intercept-SE caveat")
    print(f"seed={SEED} N={N} target-chr SNPs={Xc.shape[1]} R={R} "
          f"tau2={tau2:.4g} h2_true={h2_true:.2f}")

    o_int, o_slp, o_int_se, o_slp_se = ([] for _ in range(4))
    f_int, f_slp, f_int_se, f_slp_se, h2_hat = ([] for _ in range(5))

    for _ in range(R):
        beta = rng.standard_normal(Xc.shape[1]) * np.sqrt(tau2)
        y = Xc @ beta + Lv @ rng.standard_normal(N)

        # oracle arm
        z2 = (PX_o.T @ y / np.sqrt(q_o)) ** 2
        r = ld_score_regression(z2, L_oracle, n_blocks=Xc.shape[1])
        o_int.append(r.intercept); o_slp.append(r.slope)
        o_int_se.append(r.intercept_se); o_slp_se.append(r.slope_se)

        # fitted arm
        sg, se = reml_fit(U.T @ y, Ct, lam)
        h2_hat.append(sg / (sg + se))
        Ph = residual_precision(sg * K + se * np.eye(N), one)
        PXh = Ph @ Xc
        qh = np.einsum("ij,ij->j", Xc, PXh)
        Lh = ld_scores(Ph, Xc)
        z2h = (PXh.T @ y / np.sqrt(qh)) ** 2
        rf = ld_score_regression(z2h, Lh, n_blocks=Xc.shape[1])
        f_int.append(rf.intercept); f_slp.append(rf.slope)
        f_int_se.append(rf.intercept_se); f_slp_se.append(rf.slope_se)

    o_int, o_slp = np.array(o_int), np.array(o_slp)
    o_int_se, o_slp_se = np.array(o_int_se), np.array(o_slp_se)
    f_int, f_slp = np.array(f_int), np.array(f_slp)
    f_int_se, f_slp_se = np.array(f_int_se), np.array(f_slp_se)
    h2_hat = np.array(h2_hat)

    # data-driven propagation: dIntercept/dh^2 = empirical slope of intercept on
    # h^2_hat across replicates; Var contribution = (dg/dh2)^2 * Var(h2_hat).
    dg_dh2 = np.polyfit(h2_hat, f_int, 1)[0]
    var_vc = dg_dh2 ** 2 * np.var(h2_hat)
    f_int_se_total = np.sqrt(f_int_se ** 2 + var_vc)

    def line(label, est, truth, se):
        return (f"    {label:34s} coverage={coverage(est, truth, se)*100:5.1f}%  "
                f"empirical_sd={np.std(est):.4f}  mean_reported_se={np.mean(se):.4f}")

    print("\n  ORACLE V_c (toy delete-one-SNP coverage probe):")
    print(line("intercept (truth 1)", o_int, 1.0, o_int_se))
    print(line("slope (truth tau2)", o_slp, tau2, o_slp_se))

    print("\n  FITTED V_c (REML per replicate):")
    print(line("intercept, jackknife-only SE", f_int, 1.0, f_int_se))
    print(line("intercept, jackknife+V_c SE", f_int, 1.0, f_int_se_total))
    print(line("slope (truth tau2)", f_slp, tau2, f_slp_se))

    # variance decomposition of the fitted-arm intercept
    se_reg = float(np.mean(f_int_se))
    se_vc = float(np.sqrt(var_vc))
    M = Xc.shape[1]
    crossover_M = M * (se_reg / se_vc) ** 2  # where the two components are equal
    print(f"\n  intercept variance decomposition (fitted arm):")
    print(f"    regression SE (jackknife, shrinks ~1/sqrt(M)) = {se_reg:.3f}  at M={M}")
    print(f"    V_hat_c SE (shared systematic, ~const in M)   = {se_vc:.3f}  "
          f"[dIntercept/dh2={dg_dh2:.2f} x sd(h2_hat)={h2_hat.std():.3f}]")
    print(f"    -> V_hat_c term dominates once M > ~{crossover_M:,.0f} SNPs")
    print(f"    h2_hat = {h2_hat.mean():.3f} +/- {h2_hat.std():.3f} (true {h2_true:.2f})")

    print("\n" + "=" * 74)
    print("Reading: the delete-one-SNP jackknife covers in this toy fixed-operator")
    print("experiment (both arms ~95%); this is not a genomic-block validation.")
    print("At only M=80 target SNPs the regression intercept is")
    print("noise-dominated, so V_hat_c uncertainty is negligible and coverage is")
    print("fine either way. But that V_hat_c term is a SHARED systematic that does")
    print("NOT shrink with M, while the regression SE does (~1/sqrt(M)); past the")
    print("crossover M above it dominates and the jackknife-only intercept SE may")
    print("under-cover -- so at realistic M the M5 design must propagate variance-")
    print("component uncertainty. The extrapolation is illustrative, not a finished")
    print("correction; an outer refit/bootstrap or joint sandwich is still needed.")


if __name__ == "__main__":
    main()
