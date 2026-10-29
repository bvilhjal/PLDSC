"""M4: does PLDSC stay calibrated with a FITTED V_c (not the oracle)?

Everything so far (M-spike, M3) used an oracle V_c. Real use estimates the
variance components. This experiment quantifies the plug-in error:

  M4a  deliberate misspecification sweep: fix the true V_c, vary the ASSUMED
       h^2 used to build P_c, and read the null intercept. Exactly 1 at the
       truth; h^2=0 recovers cov-LDSC. This is the intercept-as-diagnostic
       curve (docs/METHOD.md Section M9).
  M4b  fitted plug-in error: over many null replicates, REML-fit (sigma_g^2,
       sigma_e^2) from the phenotype, build P_hat_c, and report the resulting
       intercept distribution. Tight around 1 => calibration survives fitting.

Note: PLDSC the package stays oracle-only (no REML). The REML here is local to
the experiment -- it stands in for the fitted null variance components that M4
would otherwise import from GEMMA. P_c depends on V_c only through its overall
scale AND the h^2 ratio, and REML fits both, so this is a faithful probe.

Run:  python experiments/m4_fitted_vc.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import null_inflation, residual_precision  # noqa: E402
from pldsc.simulate import grm, simulate_admixture, standardize  # noqa: E402

SEED = 20260718


def golden_min(f, lo, hi, tol=1e-6, maxit=200):
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(maxit):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
        if abs(b - a) < tol:
            break
    return 0.5 * (a + b)


def reml_fit(yt, Ct, lam):
    """REML for y = C beta + g + e, g~N(0, sg2 K), e~N(0, se2 I), profiled over
    delta = se2/sg2 in the eigenbasis of K (yt=U^T y, Ct=U^T C, lam=eig(K))."""
    n, q = Ct.shape

    def neg_ll(log_delta):
        w = lam + np.exp(log_delta)
        invw = 1.0 / w
        CtW = Ct * invw[:, None]
        CtWCt = Ct.T @ CtW
        CtWy = CtW.T @ yt
        beta = np.linalg.solve(CtWCt, CtWy)
        rss = float(np.sum(yt * yt * invw) - CtWy @ beta)
        _, logdet = np.linalg.slogdet(CtWCt)
        return (n - q) * np.log(rss / (n - q)) + np.sum(np.log(w)) + logdet

    log_delta = golden_min(neg_ll, np.log(1e-4), np.log(1e4))
    w = lam + np.exp(log_delta)
    invw = 1.0 / w
    CtW = Ct * invw[:, None]
    CtWCt = Ct.T @ CtW
    CtWy = CtW.T @ yt
    beta = np.linalg.solve(CtWCt, CtWy)
    rss = float(np.sum(yt * yt * invw) - CtWy @ beta)
    sg2 = rss / (n - q)
    se2 = np.exp(log_delta) * sg2
    return sg2, se2


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    chrom_sizes = [80, 80, 80]
    sg2_true, se2_true = 0.7, 0.3          # h^2_true = 0.7 (K diag ~ 1)
    h2_true = sg2_true / (sg2_true + se2_true)

    G, chrom, pos = simulate_admixture(rng, chrom_sizes, n_indiv=400, fst=0.15,
                                       families=True, sibs=4)
    keep = G.var(0) > 0
    G, chrom, pos = G[:, keep], chrom[keep], pos[keep]
    tgt = chrom == 0
    N = G.shape[0]
    Xc = standardize(G[:, tgt])
    K = grm(G[:, ~tgt])
    Vc_true = sg2_true * K + se2_true * np.eye(N)
    one = np.ones((N, 1))

    lam, U = np.linalg.eigh(K)
    Ct = U.T @ one
    Lv = np.linalg.cholesky(Vc_true)

    print("PLDSC M4: fitted-V_c plug-in error")
    print(f"seed={SEED} bitgen=PCG64  N={N}  target-chr SNPs={Xc.shape[1]}")
    print(f"true sigma_g2={sg2_true} sigma_e2={se2_true}  h2_true={h2_true:.2f}")

    def intercept(sg2, se2):
        P = residual_precision(sg2 * K + se2 * np.eye(N), one)
        return float(np.mean(null_inflation(P, Xc, Vc_true)))

    # ---- M4a: deliberate misspecification sweep over assumed h^2 ----
    print("\n[M4a] mean null chi^2 vs ASSUMED h^2 (true V_c fixed; want 1.000 at truth)")
    for h2 in [0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
        tag = "  <-- cov-LDSC (h2=0)" if h2 == 0.0 else (
            "  <-- truth" if abs(h2 - h2_true) < 1e-9 else "")
        print(f"    assumed h2={h2:4.2f}:  {intercept(h2, 1.0 - h2):.3f}{tag}")

    # ---- M4b: REML-fitted plug-in error over null replicates ----
    R = 300
    sg2_hat = np.empty(R)
    se2_hat = np.empty(R)
    icpt = np.empty(R)
    for r in range(R):
        y = Lv @ rng.standard_normal(N)          # null phenotype ~ N(0, Vc_true)
        yt = U.T @ y
        sg2, se2 = reml_fit(yt, Ct, lam)
        sg2_hat[r], se2_hat[r] = sg2, se2
        icpt[r] = intercept(sg2, se2)
    h2_hat = sg2_hat / (sg2_hat + se2_hat)

    print(f"\n[M4b] REML plug-in over {R} null replicates")
    print(f"    REML recovery: sigma_g2={sg2_hat.mean():.3f}+/-{sg2_hat.std():.3f} "
          f"(true {sg2_true})   h2={h2_hat.mean():.3f}+/-{h2_hat.std():.3f} "
          f"(true {h2_true:.2f})")
    print(f"    fitted-P_c intercept: mean={icpt.mean():.4f}  sd={icpt.std():.4f}  "
          f"[5%,95%]=[{np.quantile(icpt,0.05):.3f}, {np.quantile(icpt,0.95):.3f}]")
    print(f"    oracle intercept = 1.000 ; cov-LDSC (h2=0) = {intercept(0.0,1.0):.3f}")

    print("\n" + "=" * 74)
    print("Reading: M4a should pass through 1.000 at the true h^2 and inflate as")
    print("h^2 -> 0 (cov-LDSC). M4b tells the real story: if the fitted-P_c")
    print("intercept clusters tightly around 1, PLDSC's calibration survives")
    print("realistic variance-component estimation error.")


if __name__ == "__main__":
    main()
