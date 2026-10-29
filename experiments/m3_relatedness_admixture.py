"""M3 scenario: does PLDSC beat cov-LDSC (not just ordinary LDSC)?

cov-LDSC corrects population structure with a few fixed principal components.
That works when the confounding covariance is low-rank (a continuous admixture
axis) but NOT when it is high-rank (cryptic relatedness: each family is its own
block). The LOCO mixed-model operator P_c uses the full GRM V_c, so it should
calibrate where fixed PCs cannot.

This experiment isolates that claim across three cohorts:
  A. cryptic relatedness (homogeneous pop, many sibships)   -> PCs should fail
  B. continuous admixture (no families)                     -> PCs should work
  C. admixture + relatedness (the realistic biobank case)   -> PCs should fail

Headline metric: mean null chi^2 = mean_j (x_j^T P V_c P x_j)/(x_j^T P x_j),
the LDSC-intercept-style confounding under a null target chromosome. It is
EXACTLY 1 for the oracle P_c (property P2: P_c V_c P_c = P_c) and inflated for
ordinary / cov-LDSC OLS operators. Computed analytically (no Monte Carlo); one
MC pass cross-checks the Eq (3) moment for PLDSC.

Run:  python experiments/m3_relatedness_admixture.py
"""
from __future__ import annotations

import numpy as np

SEED = 20260718


# --------------------------------------------------------------------------- #
# helpers (self-contained; mirror experiments/m_spike_kill_or_continue.py)
# --------------------------------------------------------------------------- #
def ndtri(p):
    """Inverse standard-normal CDF (Acklam), vectorized over an array."""
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
    if lo.any():
        q = np.sqrt(-2 * np.log(p[lo]))
        x[lo] = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if hi.any():
        q = np.sqrt(-2 * np.log(1 - p[hi]))
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
    return np.eye(C.shape[0]) - C @ np.linalg.pinv(C.T @ C) @ C.T


def fit_intercept_slope(x, y):
    A = np.column_stack((np.ones_like(x), x))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef[0], coef[1]


def ar1_haps(rng, freqs_mat, rho):
    """(n_hap x m) haplotypes: AR(1) latent LD, per-haplotype marginal freqs."""
    n_hap, m = freqs_mat.shape
    z = np.empty((n_hap, m))
    z[:, 0] = rng.standard_normal(n_hap)
    root = np.sqrt(1 - rho * rho)
    for s in range(1, m):
        z[:, s] = rho * z[:, s - 1] + root * rng.standard_normal(n_hap)
    return (z < ndtri(freqs_mat)).astype(np.float64)


def bn_freqs(rng, m, fst):
    anc = rng.uniform(0.1, 0.9, size=m)
    a, b = anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst
    return rng.beta(a, b), rng.beta(a, b)


# --------------------------------------------------------------------------- #
# cohort builders
# --------------------------------------------------------------------------- #
def _assemble(chrom_sizes):
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(chrom_sizes)])
    pos = np.concatenate([np.arange(m) * 0.15 for m in chrom_sizes])
    return chrom, pos


def build_relatedness(rng, chrom_sizes, n_fam, sibs, rho):
    """Homogeneous population, n_fam sibships of `sibs` full siblings."""
    freqs = [rng.uniform(0.1, 0.9, m) for m in chrom_sizes]
    rows = []
    for _ in range(n_fam):
        # 4 parental haplotypes per chromosome (2 parents x 2 haps)
        phap = [ar1_haps(rng, np.tile(freqs[c], (4, 1)), rho)
                for c in range(len(chrom_sizes))]
        for _s in range(sibs):
            child = [phap[c][rng.integers(2)] + phap[c][2 + rng.integers(2)]
                     for c in range(len(chrom_sizes))]
            rows.append(np.concatenate(child))
    return np.array(rows), *_assemble(chrom_sizes)


def build_admixture(rng, chrom_sizes, n_indiv, rho, fst, families=False, sibs=4):
    """Two ancestral pops, continuous per-individual admixture proportion.

    families=False: unrelated admixed individuals.
    families=True : sibships whose parents share a family admixture proportion
                    (endogamy) -> admixture axis PLUS within-family relatedness.
    """
    p1, p2 = zip(*[bn_freqs(rng, m, fst) for m in chrom_sizes])
    rows = []
    if not families:
        for _ in range(n_indiv):
            ai = rng.uniform(0.0, 1.0)
            ind = []
            for c in range(len(chrom_sizes)):
                freq = ai * p1[c] + (1 - ai) * p2[c]
                h = ar1_haps(rng, np.tile(freq, (2, 1)), rho)
                ind.append(h[0] + h[1])
            rows.append(np.concatenate(ind))
    else:
        n_fam = n_indiv // sibs
        for _ in range(n_fam):
            ai = rng.uniform(0.0, 1.0)                 # family admixture proportion
            phap = []
            for c in range(len(chrom_sizes)):
                freq = ai * p1[c] + (1 - ai) * p2[c]
                phap.append(ar1_haps(rng, np.tile(freq, (4, 1)), rho))
            for _s in range(sibs):
                child = [phap[c][rng.integers(2)] + phap[c][2 + rng.integers(2)]
                         for c in range(len(chrom_sizes))]
                rows.append(np.concatenate(child))
    return np.array(rows), *_assemble(chrom_sizes)


# --------------------------------------------------------------------------- #
# analysis
# --------------------------------------------------------------------------- #
def mean_relatedness_offdiag(K, thresh=0.35):
    """Count/mean of clearly-related pairs (GRM off-diag above the marker-noise
    floor). thresh=0.35 catches true sibs/parent-offspring (~0.5) while ignoring
    the ~1/sqrt(m) finite-marker noise in a small GRM."""
    n = K.shape[0]
    upper = K[np.triu_indices(n, 1)]
    hits = upper > thresh
    return float(np.mean(upper[hits]) if np.any(hits) else 0.0), int(np.sum(hits))


def _analyze_impl(name, G, chrom, pos, sigma_g2, sigma_e2, target_h2_infl,
                  window, pc_list, rng, mc_check=False, n_rep=4000):
    keep = G.var(0) > 0
    G, chrom, pos = G[:, keep], chrom[keep], pos[keep]
    N = G.shape[0]
    tgt = chrom == 0
    off = ~tgt
    Xc = standardize(G[:, tgt])
    Mc = Xc.shape[1]
    Zoff = standardize(G[:, off])
    Kloco = (Zoff @ Zoff.T) / Zoff.shape[1]
    Vc = sigma_g2 * Kloco + sigma_e2 * np.eye(N)

    one = np.ones((N, 1))
    W = np.linalg.inv(Vc)
    VinvC = W @ one
    P_c = W - VinvC @ np.linalg.pinv(one.T @ VinvC) @ VinvC.T

    ev, evec = np.linalg.eigh(Kloco)
    pcs_all = evec[:, ::-1]  # descending

    def null_chi2(P):
        PXc = P @ Xc
        q = np.einsum("ij,ij->j", Xc, PXc)
        num = np.einsum("ij,ij->j", PXc, Vc @ PXc)  # x_j^T P Vc P x_j
        return float(np.mean(num / q))

    rel_mean, n_rel = mean_relatedness_offdiag(Kloco)
    print(f"\n{'='*74}\nSCENARIO: {name}")
    print(f"  N={N}  target-chr SNPs={Mc}  sigma_g2={sigma_g2} sigma_e2={sigma_e2}")
    print(f"  clearly-related pairs (GRM off-diag>0.35): {n_rel}  mean={rel_mean:.2f}")
    print("\n  mean null chi^2 (LDSC-intercept confounding; want 1.000):")
    print(f"    ordinary  ( 0 PC)  : {null_chi2(projector(one)):.3f}")
    for k in pc_list:
        Ck = np.column_stack((one, pcs_all[:, :k]))
        print(f"    cov-LDSC  ({k:2d} PC)  : {null_chi2(projector(Ck)):.3f}")
    print(f"    PLDSC (full GRM)   : {null_chi2(P_c):.3f}   <-- mixed-model oracle")

    if mc_check:
        PXc = P_c @ Xc
        q = np.einsum("ij,ij->j", Xc, PXc)
        A = Xc.T @ PXc
        L_full = (A ** 2).sum(1) / q
        tau2 = target_h2_infl / np.mean(L_full)
        Lv = np.linalg.cholesky(Vc)
        Bmat = rng.standard_normal((Mc, n_rep)) * np.sqrt(tau2)
        E = Lv @ rng.standard_normal((N, n_rep))
        Y = Xc @ Bmat + E
        z2 = np.mean((PXc.T @ Y / np.sqrt(q)[:, None]) ** 2, axis=1)
        a, b = fit_intercept_slope(L_full, z2)
        print(f"  [Eq (3) MC check, PLDSC] regress z^2 on full L: "
              f"(intercept, slope)=({a:.3f}, {b:.4g})  vs tau2={tau2:.4g}")


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))
    cs = [80, 80, 80]
    pc_list = [2, 5, 10, 20, 50]
    print("PLDSC M3: relatedness / admixture -- PLDSC vs cov-LDSC")
    print(f"seed={SEED} bitgen=PCG64")

    G, chrom, pos = build_relatedness(rng, cs, n_fam=100, sibs=4, rho=0.6)
    _analyze_impl("A. cryptic relatedness (100 sibships x 4, homogeneous)",
                  G, chrom, pos, 0.7, 0.3, 0.5, 1.5, pc_list, rng,
                  mc_check=True, n_rep=15000)

    G, chrom, pos = build_admixture(rng, cs, n_indiv=400, rho=0.6, fst=0.15,
                                    families=False)
    _analyze_impl("B. continuous admixture (unrelated, fst=0.15)",
                  G, chrom, pos, 0.7, 0.3, 0.5, 1.5, pc_list, rng, mc_check=False)

    G, chrom, pos = build_admixture(rng, cs, n_indiv=400, rho=0.6, fst=0.15,
                                    families=True, sibs=4)
    _analyze_impl("C. admixture + relatedness (biobank-like, fst=0.15)",
                  G, chrom, pos, 0.7, 0.3, 0.5, 1.5, pc_list, rng, mc_check=False)

    print(f"\n{'='*74}")
    print("Reading: PLDSC null chi^2 should be ~1.000 in all three. cov-LDSC")
    print("should approach 1 with more PCs under admixture (B) but stay inflated")
    print("under relatedness (A, C) for any feasible number of PCs.")


if __name__ == "__main__":
    main()
