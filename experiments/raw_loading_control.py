"""The RAW (un-smeared) ancestry loading as a structure control for S-LDSC.

Follow-up to `identifiability_aliasing.py`, `local_adaptation_channel.py` and
`sldsc_pc_annotation.py`. Those established (i) that the rank-1 directional
component is exactly aliased, (ii) that the surviving heteroscedastic channel is
detectable only via a `corr^2 x ell^local` interaction whose naive test is badly
anti-conservative, and (iii) that a *smeared* PC-loading annotation removes
structure-driven false enrichment but cannibalises genuine enrichment.

The observation tested here: environmental confounding along an axis enters
`E[z_j^2]` through the SNP's OWN loading `a_j = corr(x_j,u_1)^2` and is NOT
propagated by LD, whereas an ancestry-differentiated genetic architecture enters
through the LD-SMEARED annotation score `ell_j(a) = sum_{k in W_j} r_jk^2 a_k`.
Those are separate regressors -- but only once the self term is removed, because
`r_jj^2 = 1` puts a whole copy of `a_j` inside `ell_j(a)`. S-LDSC's annotation LD
scores include that self term, which is precisely why a smeared PC annotation
cannot separate the two channels.

Part 1 -- collinearity and null calibration. Compares three statistics for the
heteroscedastic channel on identical phenotypes.
Part 2 -- the practical case. Does adding the raw loading fix structure-driven
false enrichment for a functional annotation without cannibalising a real one?

Run:  python experiments/raw_loading_control.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.genome import grm_eigh, unit_columns  # noqa: E402
from pldsc.simulate import simulate_admixture, standardize  # noqa: E402

SEED = 20260725
CHROMS = [400] * 5          # 5 chromosomes x 400 variants, 0.15 cM spacing
WINDOW = 1.0                # cM


# --------------------------------------------------------------- primitives

def annot_score(Xn, chrom, pos, window, annot=None, self_term=True):
    """`sum_{k in W_j} r_jk^2 * annot_k` (annot=None -> the plain LD score)."""
    out = np.zeros(Xn.shape[1])
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        Xc, pc = Xn[:, idx], pos[idx]
        ac = None if annot is None else annot[idx]
        for jj, j in enumerate(idx):
            lo = np.searchsorted(pc, pc[jj] - window, side="left")
            hi = np.searchsorted(pc, pc[jj] + window, side="right")
            r2 = (Xc[:, lo:hi].T @ Xc[:, jj]) ** 2
            if not self_term:
                r2[jj - lo] = 0.0
            out[j] = r2.sum() if ac is None else r2 @ ac[lo:hi]
    return out


def cohort(rng, n_indiv):
    """Admixed cohort -> standardized/unit genotypes, map, top axis, loadings."""
    G, chrom, pos = simulate_admixture(rng, CHROMS, n_indiv, fst=0.15, rho=0.9)
    Z = standardize(G)
    Xn, dead = unit_columns(Z)
    keep = ~dead
    Z, Xn, chrom, pos = Z[:, keep], Xn[:, keep], chrom[keep], pos[keep]
    _, U = grm_eigh(Z)
    u1 = U[:, -1]
    a = (Z.T @ u1 / np.sqrt(Z.shape[0])) ** 2
    return Z, Xn, chrom, pos, u1, a


def ols(X, y, k, xtx_inv=None, want_t=True):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    if not want_t:
        return beta[k]
    n, p = X.shape
    resid = y - X @ beta
    s2 = (resid @ resid) / (n - p)
    if xtx_inv is None:
        xtx_inv = np.linalg.pinv(X.T @ X)
    return beta[k] / np.sqrt(s2 * np.diag(xtx_inv)[k])


def _tune(f, target, hi=60.0, iters=40):
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------- part 1: collinearity + the null

def part1(rng, n_indiv, kappa, target_chi2, reps, h2=0.5):
    """ENV / DIR / HET arms; SPLIT vs NAIVE vs S-LDSC-convention statistics."""
    Z, Xn, chrom, pos, u1, a = cohort(rng, n_indiv)
    N, M = Z.shape
    l_loc = annot_score(Xn, chrom, pos, WINDOW)
    l_a = annot_score(Xn, chrom, pos, WINDOW, a, self_term=False)

    designs = {
        # the split: raw loading and LD-smeared loading as separate regressors
        "SPLIT": (np.column_stack([np.ones(M), l_loc, a, l_a]), 3),
        # the naive interaction test (identifiability_aliasing.py / local_adaptation)
        "NAIVE": (np.column_stack([np.ones(M), a, a * (l_loc - l_loc.mean())]), 2),
        # the S-LDSC convention: annotation LD score WITH its self term, no raw term
        "SLDSC": (np.column_stack([np.ones(M), l_loc, l_a + a]), 2),
    }
    inv = {k: np.linalg.pinv(X.T @ X) for k, (X, _) in designs.items()}

    res = {arm: {k: np.empty(reps) for k in designs} for arm in ("ENV", "DIR", "HET")}
    chi2 = {arm: np.empty(reps) for arm in res}
    for arm in res:
        for r in range(reps):
            e = rng.normal(size=N)
            e -= e.mean()
            e *= np.sqrt(1 - h2) / e.std()
            ze = Z.T @ e / np.sqrt(N)
            if arm == "HET":
                var = h2 / M * (1 + kappa * a / a.mean())
                g = Z @ rng.normal(0.0, np.sqrt(var))
                g -= g.mean()
                g /= g.std()
                zg = Z.T @ g / np.sqrt(N)
                cov = float(g @ e) / N
                s = _tune(lambda s: (((ze + s * zg) ** 2)
                                     / ((1 - h2) + s * s + 2 * s * cov)).mean(),
                          target_chi2)
                zt, vy = ze + s * zg, (1 - h2) + s * s + 2 * s * cov
            else:
                g = Z @ rng.normal(0.0, np.sqrt(h2 / M), M)
                g -= g.mean()
                g *= np.sqrt(h2) / g.std()
                zg = Z.T @ g / np.sqrt(N)
                axis = u1.copy() if arm == "ENV" else Z @ (Z.T @ u1)
                axis -= axis.mean()
                axis /= axis.std()
                za = Z.T @ axis / np.sqrt(N)
                cov = float((g + e) @ axis) / N
                vge = float(np.var(g + e))
                s = _tune(lambda s: (((zg + ze + s * za) ** 2)
                                     / (vge + s * s + 2 * s * cov)).mean(),
                          target_chi2)
                zt, vy = zg + ze + s * za, vge + s * s + 2 * s * cov
            z2 = (zt ** 2) / vy
            chi2[arm][r] = z2.mean()
            for k, (X, col) in designs.items():
                res[arm][k][r] = ols(X, z2, col, inv[k])
    return res, chi2, M, float(np.corrcoef(a, l_a)[0, 1]), float(np.corrcoef(a, l_a + a)[0, 1])


# ------------------------------------------------ part 2: S-LDSC enrichment

def part2(rng, n_indiv, reps, ancestry_correlated, h2=0.5, enrich=4.0, shift=0.8):
    Z, Xn, chrom, pos, u1, a = cohort(rng, n_indiv)
    N, M = Z.shape
    nfunc = M // 10
    if ancestry_correlated:
        c_idx = np.argsort(-a)[:nfunc]                 # the top-loading decile
    else:
        c_idx = rng.choice(M, size=nfunc, replace=False)
    c = np.zeros(M)
    c[c_idx] = 1.0

    l_loc = annot_score(Xn, chrom, pos, WINDOW)
    l_c = annot_score(Xn, chrom, pos, WINDOW, c)       # S-LDSC convention (self in)
    l_a = annot_score(Xn, chrom, pos, WINDOW, a)

    designs = {
        "A standard": np.column_stack([np.ones(M), l_loc, l_c]),
        "B +smeared": np.column_stack([np.ones(M), l_loc, l_c, l_a]),
        "C +raw":     np.column_stack([np.ones(M), l_loc, l_c, a]),
        "D +both":    np.column_stack([np.ones(M), l_loc, l_c, l_a, a]),
    }
    out = {s: {k: np.empty(reps) for k in designs} for s in ("CONF", "ENRICH")}
    for scen in out:
        for r in range(reps):
            if scen == "CONF":
                var = np.full(M, h2 / M)               # NO true enrichment
            else:
                w = 1.0 + (enrich - 1.0) * c           # true enrichment
                var = h2 / M * w / w.mean()
            g = Z @ rng.normal(0.0, np.sqrt(var))
            g -= g.mean()
            g *= np.sqrt(h2) / g.std()
            e = rng.normal(size=N)
            e -= e.mean()
            e *= np.sqrt(1 - h2) / e.std()
            y = g + e
            if scen == "CONF":                          # environmental axis shift
                ax = u1 - u1.mean()
                y = y + shift * ax / ax.std()
            z2 = ((Z.T @ y / np.sqrt(N)) / y.std()) ** 2
            for k, X in designs.items():
                out[scen][k][r] = ols(X, z2, 2, want_t=False)   # tau on ell(c)
    return out, M


# ------------------------------------------------------------------- driver

def main():
    rng = np.random.default_rng(SEED)

    print("PART 1 -- collinearity and null calibration "
          "(4 realizations x 60 phenotype reps; thresholds = ENV 95th pct)")
    for n_indiv, kappa in ((600, 1.0), (600, 2.0), (2400, 1.0)):
        acc = {arm: {} for arm in ("ENV", "DIR", "HET")}
        chi = {arm: [] for arm in acc}
        cols, cols_self = [], []
        for _ in range(4):
            res, chi2, M, col, col_self = part1(rng, n_indiv, kappa, 1.8, 60)
            cols.append(col)
            cols_self.append(col_self)
            for arm in res:
                chi[arm].append(chi2[arm])
                for k, v in res[arm].items():
                    acc[arm].setdefault(k, []).append(v)
        acc = {arm: {k: np.concatenate(v) for k, v in d.items()} for arm, d in acc.items()}
        chi = {arm: np.concatenate(v).mean() for arm, v in chi.items()}
        print(f"\n  N={n_indiv} kappa={kappa} M={M} | "
              f"corr(a, ell(a)) no-self {np.mean(cols):+.3f}, with-self {np.mean(cols_self):+.3f}"
              f" | mean chi2 ENV/DIR/HET {chi['ENV']:.2f}/{chi['DIR']:.2f}/{chi['HET']:.2f}")
        print(f"  {'stat':6s} {'med t ENV':>10s} {'med t DIR':>10s} {'med t HET':>10s}"
              f" {'thr95':>7s} {'pow HET':>8s} {'pow DIR':>8s}")
        for k in ("SPLIT", "NAIVE", "SLDSC"):
            thr = np.quantile(acc["ENV"][k], 0.95)
            print(f"  {k:6s} {np.median(acc['ENV'][k]):10.2f} {np.median(acc['DIR'][k]):10.2f}"
                  f" {np.median(acc['HET'][k]):10.2f} {thr:7.2f}"
                  f" {np.mean(acc['HET'][k] > thr):8.2f} {np.mean(acc['DIR'][k] > thr):8.2f}")

    print("\n\nPART 2 -- S-LDSC enrichment coefficient tau on ell(c) "
          "(3 realizations x 40 reps, N=600)")
    for flag in (True, False):
        acc = {s: {} for s in ("CONF", "ENRICH")}
        for _ in range(3):
            out, M = part2(rng, 600, 40, flag)
            for s in out:
                for k, v in out[s].items():
                    acc[s].setdefault(k, []).append(v)
        acc = {s: {k: np.concatenate(v).mean() for k, v in d.items()} for s, d in acc.items()}
        tag = "ancestry-CORRELATED" if flag else "ancestry-INDEPENDENT (control)"
        print(f"\n  functional annotation {tag}, M={M}")
        print(f"  {'design':11s} {'tau CONF':>11s} {'rel A':>7s} | "
              f"{'tau ENRICH':>11s} {'rel A':>7s}")
        for k in ("A standard", "B +smeared", "C +raw", "D +both"):
            print(f"  {k:11s} {acc['CONF'][k]:11.3e} {acc['CONF'][k]/acc['CONF']['A standard']:7.2f} | "
                  f"{acc['ENRICH'][k]:11.3e} {acc['ENRICH'][k]/acc['ENRICH']['A standard']:7.2f}")
        print(f"  (CONF wants tau = 0: no true enrichment, environmental axis shift only. "
              f"ENRICH wants tau > 0: {4.0:.0f}x true enrichment, no confounding.)")


if __name__ == "__main__":
    main()
