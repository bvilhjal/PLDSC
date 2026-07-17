"""M-spike: kill-or-continue experiment for PLDSC.

Throwaway per docs/DEVELOPMENT_PLAN.md (Table 1a, milestone "M-spike"). Builds a
dense P_c oracle on ONE tiny structured cohort and runs:

  Readout B - correctness: does the empirical Monte-Carlo E[z_j^2] match the
              oracle moment Eq (3), and does the windowed L_j (Eq 4) converge to
              the full-chromosome sum?
  Readout A - does the problem matter: on the SAME simulated phenotypes, do
              ordinary and covariate-adjusted LDSC miscalibrate (intercept away
              from 1) where the mixed-model-matched (P_c) score calibrates?

The three association models share one set of simulated phenotypes and differ
only in the operator P applied to form the score statistic:
  * ordinary OLS         P = I - 1(1'1)^+1'                    -> ordinary LDSC
  * covariate OLS        P = I - C(C'C)^+C',  C=[1, LOCO PCs]  -> cov-LDSC
  * oracle LMM score     P = P_c (Eq 1) with V_c = sg2 K_-c + se2 I -> PLDSC

Intended to be RUN once as a go/no-go gate, not shipped. Reference-quality
float64; RNG bit generator pinned alongside the seed per the plan.

Run:  <miniconda-python> experiments/m_spike_kill_or_continue.py
"""
from __future__ import annotations

import numpy as np

SEED = 20260717


# --------------------------------------------------------------------------- #
# small numerical helpers
# --------------------------------------------------------------------------- #
def ndtri(p):
    """Inverse standard-normal CDF (Acklam's rational approximation, ~1e-9)."""
    p = np.asarray(p, dtype=np.float64)
    a = [-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00]
    b = [-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00]
    plow, phigh = 0.02425, 1 - 0.02425
    x = np.empty_like(p)
    lo, hi = p < plow, p > phigh
    mid = ~(lo | hi)
    q = np.sqrt(-2 * np.log(p[lo])) if lo.any() else np.array([])
    x[lo] = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = np.sqrt(-2 * np.log(1 - p[hi])) if hi.any() else np.array([])
    x[hi] = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p[mid] - 0.5
    r = q * q
    x[mid] = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
             (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    return x


def standardize(G):
    p = G.mean(0) / 2.0
    s = np.sqrt(2.0 * p * (1.0 - p))
    s = np.where(s == 0, 1.0, s)
    return (G - 2.0 * p) / s


def projector(C):
    """Symmetric OLS residual-forming projector I - C(C'C)^+C'."""
    n = C.shape[0]
    return np.eye(n) - C @ np.linalg.pinv(C.T @ C) @ C.T


def windowed_r2_score(Xstd, pos, window):
    """Sum of squared same-chromosome correlations within +/- window (incl self)."""
    m = Xstd.shape[1]
    Xn = Xstd / np.linalg.norm(Xstd, axis=0)
    L = np.empty(m)
    for j in range(m):
        sel = np.abs(pos - pos[j]) <= window
        r = Xn[:, sel].T @ Xn[:, j]
        L[j] = np.sum(np.clip(r, -1, 1) ** 2)
    return L


def fit_intercept_slope(x, y):
    A = np.column_stack((np.ones_like(x), x))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef[0], coef[1]


# --------------------------------------------------------------------------- #
# genotype simulation: two populations, AR(1) local LD, optional sib relatedness
# --------------------------------------------------------------------------- #
def bn_freqs(rng, m, fst):
    anc = rng.uniform(0.1, 0.9, size=m)
    a, b = anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst
    return rng.beta(a, b), rng.beta(a, b)


def ar1_haps(rng, freqs, rho, n_hap):
    """n_hap x m haplotypes: AR(1) latent field, exact marginal Bernoulli(freqs)."""
    m = freqs.size
    z = np.empty((n_hap, m))
    z[:, 0] = rng.standard_normal(n_hap)
    root = np.sqrt(1 - rho * rho)
    for s in range(1, m):
        z[:, s] = rho * z[:, s - 1] + root * rng.standard_normal(n_hap)
    thr = ndtri(freqs)[None, :]
    return (z < thr).astype(np.float64)  # allele=1 w.p. freqs


def simulate(rng, chrom_sizes, n_found_per_pop, n_sib_per_pop, fst, rho):
    """Return genotype matrix (N x M), chrom id per SNP, cM positions, pop labels."""
    npop = 2
    per_pop_freqs = []
    for _ in range(npop):
        per_pop_freqs.append([bn_freqs(rng, m, fst)[0] for m in chrom_sizes])
    # per-chromosome founder haplotypes, per population
    genos, pops = [], []
    for g in range(npop):
        n_found = n_found_per_pop
        # 2 haplotypes per founder, per chromosome
        hap_by_chr = [ar1_haps(rng, per_pop_freqs[g][c], rho, 2 * n_found)
                      for c in range(len(chrom_sizes))]
        found = [hap_by_chr[c][0::2] + hap_by_chr[c][1::2]
                 for c in range(len(chrom_sizes))]  # genotype per chr
        found_geno = np.hstack(found)
        # sib pairs: children of two random founders, independent chromosome assortment
        sib_rows = []
        for _ in range(n_sib_per_pop):
            pa, pb = rng.choice(n_found, size=2, replace=False)
            child = []
            for c in range(len(chrom_sizes)):
                ha = hap_by_chr[c][2 * pa + rng.integers(2)]
                hb = hap_by_chr[c][2 * pb + rng.integers(2)]
                child.append(ha + hb)
            sib_rows.append(np.concatenate(child))
        block = found_geno if not sib_rows else np.vstack([found_geno, np.array(sib_rows)])
        genos.append(block)
        pops.append(np.full(block.shape[0], g))
    G = np.vstack(genos)
    pop = np.concatenate(pops)
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(chrom_sizes)])
    pos = np.concatenate([np.arange(m) * 0.15 for m in chrom_sizes])  # cM
    # drop monomorphic SNPs
    keep = (G.var(0) > 0)
    return G[:, keep], chrom[keep], pos[keep], pop


# --------------------------------------------------------------------------- #
# the experiment for one scenario
# --------------------------------------------------------------------------- #
def run_scenario(name, rng, chrom_sizes, n_found, n_sib, fst, rho,
                 sigma_g2, sigma_e2, target_h2_infl, window, n_rep, n_pc=2):
    G, chrom, pos, pop = simulate(rng, chrom_sizes, n_found, n_sib, fst, rho)
    N, M = G.shape
    tgt = (chrom == 0)
    off = ~tgt
    Xc = standardize(G[:, tgt])          # target-chromosome standardized genotypes
    pos_c = pos[tgt]
    Mc = Xc.shape[1]

    # LOCO GRM and V_c from the OTHER chromosomes
    Zoff = standardize(G[:, off])
    Kloco = (Zoff @ Zoff.T) / Zoff.shape[1]
    Vc = sigma_g2 * Kloco + sigma_e2 * np.eye(N)
    Lv = np.linalg.cholesky(Vc)

    # LOCO PCs as fixed covariates
    w, V = np.linalg.eigh(Kloco)
    pcs = V[:, ::-1][:, :n_pc]
    one = np.ones((N, 1))
    C_cov = np.column_stack((one, pcs))

    # ---- the three operators ----
    P_ord = projector(one)
    P_cov = projector(C_cov)
    Vinv = np.linalg.inv(Vc)
    # P_c with fixed design = intercept (LMM handles structure via the random effect)
    VinvC = Vinv @ one
    P_c = Vinv - VinvC @ np.linalg.pinv(one.T @ VinvC) @ VinvC.T

    def score_pieces(P):
        PXc = P @ Xc
        q = np.einsum("ij,ij->j", Xc, PXc)
        return PXc, q

    PXc_ord, q_ord = score_pieces(P_ord)
    PXc_cov, q_cov = score_pieces(P_cov)
    PXc_c, q_c = score_pieces(P_c)

    # oracle moment ingredients (Eq 3/4) for the P_c model
    A = Xc.T @ PXc_c                      # x_j^T P_c x_k
    L_pldsc_full = (A ** 2).sum(1) / q_c  # sum_k q_k r_jk^2  (full chromosome)
    win = np.abs(pos_c[:, None] - pos_c[None, :]) <= window
    L_pldsc_win = ((A ** 2) * win).sum(1) / q_c

    # choose tau^2 so mean full-sum inflation ~ target_h2_infl
    tau2 = target_h2_infl / np.mean(L_pldsc_full)
    mu_eq3 = 1.0 + tau2 * L_pldsc_full    # analytic oracle E[z^2] (Eq 3)

    # ---- Monte Carlo ----
    # With signal:  y = Xc @ beta + eps,  beta~N(0,tau2), eps~N(0,Vc)  (Readout B)
    # Null:         y = eps only (tau2=0)                              (Readout A1)
    models = (("ord", PXc_ord, q_ord), ("cov", PXc_cov, q_cov), ("pldsc", PXc_c, q_c))
    z2, z2_null = {}, {}
    B = rng.standard_normal((Mc, n_rep)) * np.sqrt(tau2)
    E = Lv @ rng.standard_normal((N, n_rep))
    Y = Xc @ B + E
    for key, PX, q in models:
        z2[key] = np.mean((PX.T @ Y / np.sqrt(q)[:, None]) ** 2, axis=1)
        z2_null[key] = np.mean((PX.T @ E / np.sqrt(q)[:, None]) ** 2, axis=1)

    # ---- Readout B: MC vs Eq 3 for the oracle model ----
    var_z2 = np.var((PXc_c.T @ Y) / np.sqrt(q_c)[:, None], axis=1)
    mcse = np.sqrt(var_z2 / n_rep)
    resid_sigma = (z2["pldsc"] - mu_eq3) / mcse
    trunc = np.abs(L_pldsc_win - L_pldsc_full) / L_pldsc_full

    # ---- Readout A ----
    L_ord = windowed_r2_score(standardize(G[:, tgt]), pos_c, window)  # raw-geno LDSC
    Gres = P_cov @ G[:, tgt]
    L_cov = windowed_r2_score(Gres, pos_c, window)                    # cov-LDSC

    # A(1) confounding under a NULL target chromosome (tau2=0): mean z^2 IS the
    # LDSC-intercept-style confounding. Only the mixed-model (P_c) statistic is
    # calibrated to exactly 1; ordinary/cov OLS retain structure confounding.
    # This is the Section 1.1 use-case demonstrated directly.
    null_conf = {k: float(np.mean(v)) for k, v in z2_null.items()}

    # A(2) regressor bias: regress the SAME calibrated oracle z^2 on each LD
    # score. Only the correct (P_c) score should give intercept ~1 with high R^2;
    # ordinary/cov LD scores are the wrong regressor and bias the fit.
    def fit_full(x, y):
        a, b = fit_intercept_slope(x, y)
        yhat = a + b * x
        r2 = 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - np.mean(y)) ** 2)
        return a, b, r2

    reg_on_oracle = {
        "L_pldsc (correct)": fit_full(L_pldsc_full, z2["pldsc"]),
        "L_ordinary": fit_full(L_ord, z2["pldsc"]),
        "L_cov-LDSC": fit_full(L_cov, z2["pldsc"]),
    }

    # ---- long-range LD illustration: mean cross-chromosome r^2, raw vs P_c ----
    Xall = standardize(G)
    Xn = Xall / np.linalg.norm(Xall, axis=0)
    cross = np.outer(chrom == 0, chrom == 2)
    raw_r2 = (Xn.T @ Xn) ** 2
    # P_c-metric cross-correlation between target and chr-2 SNPs
    Xc2 = standardize(G[:, chrom == 2])
    PcXc2 = P_c @ Xc2
    q2 = np.einsum("ij,ij->j", Xc2, PcXc2)
    A02 = Xc.T @ PcXc2
    pc_r2 = (A02 ** 2) / np.outer(q_c, q2)
    lr_raw = raw_r2[cross].mean()
    lr_pc = pc_r2.mean()

    # ---- report ----
    print(f"\n{'='*74}\nSCENARIO: {name}")
    print(f"  N={N}  M={M}  target-chr SNPs={Mc}  fst={fst}  rho={rho}  "
          f"sib/pop={n_sib}  reps={n_rep}")
    print(f"  sigma_g2={sigma_g2} sigma_e2={sigma_e2}  tau2={tau2:.4g}  "
          f"mean full L={np.mean(L_pldsc_full):.3f}")
    print(f"\n  [Readout B] oracle moment Eq(3) vs Monte Carlo")
    print(f"    max |E[z^2]_MC - Eq3| / MCSE   = {np.max(np.abs(resid_sigma)):.2f} sigma "
          f"(should be O(1-3) if correct)")
    print(f"    RMS standardized residual      = {np.sqrt(np.mean(resid_sigma**2)):.2f} sigma")
    print(f"    mean |E[z^2]_MC - Eq3|         = {np.mean(np.abs(z2['pldsc']-mu_eq3)):.4f}")
    print(f"    window truncation |Lwin-Lfull|/Lfull: median={np.median(trunc):.3f} "
          f"max={np.max(trunc):.3f}")
    print(f"\n  [Readout A1] confounding under a NULL target chr (tau2=0): mean z^2")
    print(f"               (calibrated intercept = 1; excess = residual confounding)")
    print(f"    ordinary OLS    mean z^2 = {null_conf['ord']:7.3f}   "
          f"(confounding = {null_conf['ord']-1:+.3f})")
    print(f"    cov-adj OLS     mean z^2 = {null_conf['cov']:7.3f}   "
          f"(confounding = {null_conf['cov']-1:+.3f})")
    print(f"    mixed-model Pc  mean z^2 = {null_conf['pldsc']:7.3f}   "
          f"(confounding = {null_conf['pldsc']-1:+.3f})")
    print(f"\n  [Readout A2] regressor bias: regress the SAME calibrated oracle z^2")
    print(f"               on each LD score  (want intercept a=1, slope b=tau2={tau2:.4g})")
    for label, (a, b, r2) in reg_on_oracle.items():
        print(f"    {label:20s}  a={a:7.3f}   b={b:.4g}   R^2={r2:.3f}")
    print(f"\n  [long-range LD] mean chr0-chr2 r^2:  raw={lr_raw:.4f}   "
          f"P_c-metric={lr_pc:.4f}   (P_c should shrink structure-induced LD)")
    return {"null_confounding": null_conf, "reg_on_oracle": reg_on_oracle}


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    chrom_sizes = [80, 80, 80]
    print("PLDSC M-spike  (kill-or-continue)  seed=%d  bitgen=PCG64" % SEED)

    run_scenario("STRUCTURED + RELATED (two pops, fst=0.08, sib pairs)", rng,
                 chrom_sizes, n_found=150, n_sib=50, fst=0.08, rho=0.6,
                 sigma_g2=0.5, sigma_e2=0.5, target_h2_infl=0.5,
                 window=1.5, n_rep=6000)

    run_scenario("HOMOGENEOUS control (single pop, fst~0, no relatives)", rng,
                 chrom_sizes, n_found=200, n_sib=0, fst=1e-4, rho=0.6,
                 sigma_g2=0.5, sigma_e2=0.5, target_h2_infl=0.5,
                 window=1.5, n_rep=6000)

    print(f"\n{'='*74}")
    print("Kill-or-continue reading:")
    print("  CONTINUE if Readout B matches (small sigma) AND ordinary/cov-LDSC")
    print("  intercepts inflate above 1 under structure while PLDSC stays ~1,")
    print("  AND all three agree ~1 in the homogeneous control.")


if __name__ == "__main__":
    main()
