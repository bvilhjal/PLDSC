"""Monte Carlo check of the multi-ancestry meta-analysis LD score (MA2, MA6).

docs/MULTI_ANCESTRY.md derives the LD score a fixed-effect meta-analysis of
ancestry-specific GWAS follows, and the slope factor kappa that a
single-ancestry reference score induces. The unit tests check the in-sample
algebra exactly. This script checks the claims that matter in use: GWAS
samples and reference panels are *different* individuals, LD scores come
from finite reference panels (noise floor on within terms only, MA7), and
LDSC is an actual regression on noisy chi-square.

Two ancestries share ancestral allele frequencies (Balding-Nichols,
fst = 0.1) and have different LD lengths (AR(1) latent rho 0.8 "E" vs 0.55
"F"). GWAS: N_E = 3000, N_F = 1500, h2 = 0.5 in both. Effects are bivariate
normal per SNP on the within-ancestry standardized scale with correlation
rho_true. Meta-analysis: sample-size (Stouffer) or inverse-variance weights.
Reference panels: 400 independent individuals per ancestry. LDSC is OLS of
chi2 on the candidate score, and h2_hat = slope * M / N_total.

Pre-stated rule: PASS if (i) with the correct meta score, mean h2_hat is
within 5% of h2 in every arm, and (ii) with the E-only score, mean
h2_hat / h2 is within 0.05 of the panel-predicted kappa in every arm.

Outcome of the rule-bearing default run (N_REF = 400): **both FAIL.** Even
the single-ancestry arm (where the meta score IS the E score) returns
h2_hat / h2 = 0.80. The toy's LDSC is attenuated about 20% before any
multi-ancestry effect enters: errors-in-variables from a 400-person panel on
LD scores with a very narrow spread (sd 0.32 about a mean of 1.9). The
diagnostic in docs/MULTI_ANCESTRY.md §6 shows 0.85 / 0.97 / 1.02 at
N_REF = 400 / 1500 / 6000 on one realization. kappa is defined relative to
the slope on the correct score (MA6), so the ratio column
(E-only / meta) is the quantity MA6 predicts, and it matches in every arm.

Re-stated rules for the post-failure rerun `--n-ref 6000`, written before
that run: (i') with the correct meta score, mean h2_hat within 5% of h2 in
every arm; (iii) mean(E-only) / mean(meta) and mean(rho=1) / mean(meta)
within 0.03 of their predicted kappas in every arm. (iii) is also printed for
the default run, where it was observed after the fact, not pre-stated.

`--baseline-check` reproduces the diagnostic behind the failure: single-
ancestry LDSC on one realization with reference panels of 400 / 1500 / 6000
and with in-sample LD.

Run:  python experiments/multiancestry_meta_sim.py [--n-ref 6000]
      python experiments/multiancestry_meta_sim.py --baseline-check
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.multiancestry import (  # noqa: E402
    cross_ancestry_ld_scores,
    meta_ld_scores,
    meta_weights,
    single_reference_bias,
)
from pldsc.simulate import ar1_haps, bn_freqs  # noqa: E402

SEED = 20260925
N_CHROM, M_PER, SPACING = 4, 500, 0.1           # 2000 SNPs, 0.1 cM apart
LD_RHO = (0.8, 0.55)                            # E, F latent AR(1) correlation
N_GWAS = (3000, 1500)
N_REF = 400
H2 = 0.5
WINDOW = 1.0
REALIZATIONS, REPS = 3, 100


def genotypes(rng, n, freqs, rho):
    """Diploid genotypes, one independent AR(1) chain per chromosome."""
    cols = []
    for c in range(N_CHROM):
        f = freqs[c * M_PER:(c + 1) * M_PER]
        h = ar1_haps(rng, np.tile(f, (2 * n, 1)), rho)
        cols.append(h[:n] + h[n:])
    return np.hstack(cols)


def standardized(X):
    sd = X.std(0)
    out = np.zeros_like(X)
    live = sd > 0
    out[:, live] = (X[:, live] - X[:, live].mean(0)) / sd[live]
    return out, live


def realization(rng, arms, n_ref):
    M = N_CHROM * M_PER
    pos = (np.repeat(np.arange(N_CHROM) * 1000.0, M_PER)
           + np.tile(np.arange(M_PER) * SPACING, N_CHROM))  # chromosomes 1000 cM apart
    freqs = bn_freqs(rng, M, 0.1)
    raw = [genotypes(rng, n, f, r) for n, f, r in zip(N_GWAS, freqs, LD_RHO)]
    gwas = [standardized(G)[0] for G in raw]
    gwas_var = np.vstack([G.var(0) for G in raw])       # GWAS-side IVW information
    refs = [genotypes(rng, n_ref, f, r) for f, r in zip(freqs, LD_RHO)]
    s = cross_ancestry_ld_scores(refs, pos, window_cm=WINDOW)
    ell_E = s.scores[0, 0]

    out = {}
    for arm in arms:
        n_used, rho_true, scheme = arm
        n = np.asarray(n_used, dtype=float)
        # Panel-side prediction: weights from panel variances, availability from panels.
        avail = s.available & (n[:, None] > 0)
        om, _ = meta_weights(n, avail, variances=s.variances, scheme=scheme)
        R = [[1, rho_true], [rho_true, 1]]
        L_true = meta_ld_scores(s, om, rho=R)
        L_rho1 = meta_ld_scores(s, om)
        kappa = single_reference_bias(L_true, ell_E).kappa
        kappa_rho1 = single_reference_bias(L_true, L_rho1).kappa
        # GWAS-side meta z-weights, from their definitions.
        if scheme == "sample_size":
            c = np.sqrt(n[:, None] / n.sum()) * np.ones((1, M))
        else:
            w = n[:, None] * gwas_var
            c = np.sqrt(w / w.sum(0))
        est = {"meta": [], "E-only": [], "meta rho=1": []}
        chol = np.linalg.cholesky(np.array([[1.0, rho_true], [rho_true, 1.0]]) * H2 / M
                                  + 1e-15 * np.eye(2))
        for _ in range(REPS):
            beta = rng.standard_normal((M, 2)) @ chol.T
            chi_z = np.zeros(M)
            for a in range(2):
                if n[a] == 0:
                    continue
                X = gwas[a]
                y = X @ beta[:, a] + rng.standard_normal(X.shape[0]) * np.sqrt(1 - H2)
                y = (y - y.mean()) / y.std()
                chi_z += c[a] * (X.T @ y / np.sqrt(X.shape[0]))
            chi2 = chi_z ** 2
            for key, score in (("meta", L_true), ("E-only", ell_E), ("meta rho=1", L_rho1)):
                slope = np.polyfit(score, chi2, 1)[0]
                est[key].append(slope * M / n.sum() / H2)
        out[arm] = ({k: np.asarray(v) for k, v in est.items()}, kappa, kappa_rho1)
    return out


def baseline_check():
    """Single-ancestry LDSC attenuation vs reference-panel size (one realization)."""
    rng = np.random.default_rng(3)
    M = N_CHROM * M_PER
    pos = (np.repeat(np.arange(N_CHROM) * 1000.0, M_PER)
           + np.tile(np.arange(M_PER) * SPACING, N_CHROM))
    fE, _ = bn_freqs(rng, M, 0.1)
    G = genotypes(rng, 3000, fE, LD_RHO[0])
    X, _ = standardized(G)
    cands = {}
    for nref in (400, 1500, 6000):
        R = genotypes(rng, nref, fE, LD_RHO[0])
        cands[f"ref n={nref}"] = cross_ancestry_ld_scores([R], pos, window_cm=WINDOW).scores[0, 0]
    cands["in-sample, uncorrected"] = cross_ancestry_ld_scores(
        [G], pos, window_cm=WINDOW, finite_sample_correction=False).scores[0, 0]
    for k, v in cands.items():
        print(f"{k:24s} LD score mean {v.mean():6.3f} sd {v.std():5.3f}")
    est = {k: [] for k in cands}
    for _ in range(200):
        beta = rng.standard_normal(M) * np.sqrt(H2 / M)
        y = X @ beta + rng.standard_normal(3000) * np.sqrt(1 - H2)
        y = (y - y.mean()) / y.std()
        chi2 = (X.T @ y / np.sqrt(3000)) ** 2
        for k, v in cands.items():
            est[k].append(np.polyfit(v, chi2, 1)[0] * M / 3000 / H2)
    for k, v in est.items():
        v = np.asarray(v)
        print(f"{k:24s} h2_hat/h2 {v.mean():.3f} (se {v.std() / np.sqrt(v.size):.3f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-ref", type=int, default=N_REF)
    ap.add_argument("--baseline-check", action="store_true")
    args = ap.parse_args()
    if args.baseline_check:
        baseline_check()
        return
    rng = np.random.default_rng(SEED)
    both = N_GWAS
    arms = [((N_GWAS[0], 0), 1.0, "sample_size"), (both, 1.0, "sample_size"),
            (both, 0.6, "sample_size"), (both, 1.0, "ivw"), (both, 0.6, "ivw")]
    acc = {a: {"meta": [], "E-only": [], "meta rho=1": [], "k": [], "k1": []} for a in arms}
    for _ in range(REALIZATIONS):
        res = realization(rng, arms, args.n_ref)
        for a, (est, k, k1) in res.items():
            for key in ("meta", "E-only", "meta rho=1"):
                acc[a][key].append(est[key])
            acc[a]["k"].append(k)
            acc[a]["k1"].append(k1)
    print(f"M={N_CHROM * M_PER}, N_GWAS={N_GWAS}, N_REF={args.n_ref}, h2={H2}, "
          f"{REALIZATIONS} realizations x {REPS} reps; entries are mean h2_hat/h2 (MC se)")
    hdr = (f"{'N_E share':>9s} {'rho':>4s} {'weights':>11s} | {'meta score':>12s} | "
           f"{'E-only':>12s} {'ratio':>6s} {'pred kappa':>10s} | {'rho=1 score':>12s} {'ratio':>6s} {'pred':>6s}")
    print(hdr)
    print("-" * len(hdr))
    ok1 = ok2 = ok3 = True
    for a in arms:
        n_used, rho, scheme = a
        share = n_used[0] / sum(n_used)
        row = {k: np.concatenate(acc[a][k]) for k in ("meta", "E-only", "meta rho=1")}
        se = {k: v.std() / np.sqrt(v.size) for k, v in row.items()}
        kp, k1 = np.mean(acc[a]["k"]), np.mean(acc[a]["k1"])
        ok1 &= abs(row["meta"].mean() - 1) <= 0.05
        ok2 &= abs(row["E-only"].mean() - kp) <= 0.05
        rE = row["E-only"].mean() / row["meta"].mean()
        r1 = row["meta rho=1"].mean() / row["meta"].mean()
        ok3 &= abs(rE - kp) <= 0.03 and abs(r1 - k1) <= 0.03
        print(f"{share:9.2f} {rho:4.1f} {scheme:>11s} | {row['meta'].mean():6.3f} ({se['meta']:.3f}) | "
              f"{row['E-only'].mean():6.3f} ({se['E-only']:.3f}) {rE:6.3f} {kp:10.3f} | "
              f"{row['meta rho=1'].mean():6.3f} ({se['meta rho=1']:.3f}) {r1:6.3f} {k1:6.3f}")
    print(f"\n(i) correct meta score recovers h2 within 5%: {'PASS' if ok1 else 'FAIL'}")
    print(f"(ii) E-only bias matches panel-predicted kappa within 0.05: {'PASS' if ok2 else 'FAIL'}")
    print(f"(iii) slope ratios match predicted kappas within 0.03: {'PASS' if ok3 else 'FAIL'}")


if __name__ == "__main__":
    main()
