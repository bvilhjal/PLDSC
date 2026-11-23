#!/usr/bin/env python3
"""GBMI asthma: what Cochran's Q adds to a meta-analysis statistic (exploratory).

Part A rebuilds a two-ancestry meta-analysis from the EUR and EAS releases,
Z = cE zE + cF zF with inverse-variance weights on the per-allele scale, and
keeps only Z, Cochran's Q = (cF zE - cE zF)^2 and the weights. It asks how
much of the cross-ancestry genetic covariance g survives that compression,
against the fit that uses both strata and the fit that uses Z^2 alone.

Part B uses each ancestry's released heterogeneity p-value. With K datasets
of effective-size fractions f, shared genetic variance tau and dataset
heritability h2, the model is

    E[Z_a^2] = 1 + (N l / M) (tau + (h2 - tau) S2),      S2 = sum f^2,
    E[Q_a]   = (K - 1) + (N l / M) (h2 - tau) (1 - S2).

The ceiling tau / h2 is the genetic correlation between two datasets of one
ancestry when their effects are exchangeable. The corrected cross-ancestry
correlation divides g by sqrt(tau_EUR tau_EAS); the naive one divides by the
two metas' own heritabilities. Dataset sizes are solved from the per-variant
case and control counts and the direction strings.

Status. No decision rule was stated before this analysis. A first pass on
2026-11-20 used sample-size weights in Part A and read a score cache without
its provenance check; the analysis was then fixed in this driver and rerun on
rebuilt scores. The companion known-truth experiment with pre-stated rules is
q_heterogeneity_sim.py.

Inputs, filters, allele rules and effective N are those of
gbmi_asthma_casualty.py. Run from the repository root after
multiancestry_genome_kappa.py:
  python experiments/gbmi_asthma_heterogeneity.py --out report/asthma_heterogeneity.json

Outcome (2026-11-20, 937,544 SNPs, 92 s; report/methods.md Section 11).
  A. g is 0.1250 (SE 0.0106) from both strata. From Z and Q it is 0.0977
     (0.0161) under inverse-variance weights, 0.027 lower (paired SE 0.012),
     and 0.1236 (0.0230) under sample-size weights. Z^2 alone gives 0.049
     (0.057). Mean Q is 1.084 against 1.049 expected with identical effects.
  B. 14 EUR datasets: ceiling 0.910 (0.042). 5 EAS datasets, one dominant:
     0.84 (0.20). EUR-EAS correlation 0.738 (0.043) uncorrected, 0.744 with
     the EUR ceiling and 0.793 (0.070) with both.

Review, 2026-11-23. The Z and Q fit above has one free intercept per
equation. Under inverse-variance weights the weight share varies across SNPs,
so unequal EUR and EAS intercepts leak into g (q_model.py): on the chr22
example LD with these weights and intercepts 1.076 and 0.983 the leak is
-9.4% of g in expectation, against -22% observed here. The driver now also
fits the three stratum intercepts ("Z_and_Q_three_intercepts"); that refit has
not been run, because the GBMI files were not available. Until it is, the
inverse-variance shortfall is provisional and is not attributed to a cause.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time

import allelign
import numpy as np
from scipy.optimize import nnls
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gbmi_asthma_casualty as G  # noqa: E402
from q_model import zq_design  # noqa: E402

N_JK, GWS = 200, 29.716785


def extras(path, index, J):
    """Heterogeneity p, dataset presence and counts at the first row of each reference rsID."""
    hp, nd = np.full(J, np.nan), np.zeros(J, dtype=np.int64)
    nc, nk, pres = np.zeros(J), np.zeros(J), np.full(J, "", dtype=object)
    with gzip.open(path, "rt") as fh:
        f = fh.readline().rstrip("\r\n").split("\t")
        c = {k: f.index(k) for k in ("rsid", "inv_var_het_p", "direction", "N_case", "N_ctrl", "n_dataset")}
        for line in fh:
            p = line.rstrip("\r\n").split("\t")
            j = index.get(p[c["rsid"]])
            if j is None or pres[j] != "":
                continue
            try:
                hp[j] = float(p[c["inv_var_het_p"]])
            except ValueError:
                pass
            pres[j] = "".join("0" if ch == "?" else "1" for ch in p[c["direction"]])
            nc[j], nk[j], nd[j] = float(p[c["N_case"]]), float(p[c["N_ctrl"]]), int(p[c["n_dataset"]])
    return hp, pres, nc, nk, nd


def dataset_sizes(pres, nc, nk, ok):
    """Per-dataset cases, controls and effective N, and each SNP's sum of squared fractions."""
    pats, first, inv, cnt = np.unique(pres[ok].astype(str), return_index=True, return_inverse=True, return_counts=True)
    A = np.array([[ch == "1" for ch in p] for p in pats], dtype=float)
    w = np.sqrt(cnt)[:, None]
    case, _ = nnls(A * w, nc[ok][first] * w[:, 0])
    ctrl, _ = nnls(A * w, nk[ok][first] * w[:, 0])
    neff = np.where((case >= 1) & (ctrl >= 1), 4.0 / (1.0 / np.maximum(case, 1e-9) + 1.0 / np.maximum(ctrl, 1e-9)), 0.0)
    P = A[inv]
    S2 = (P @ neff ** 2) / (P @ neff) ** 2
    resid = max(np.abs(P @ case - nc[ok]).max(), np.abs(P @ ctrl - nk[ok]).max())
    return case, ctrl, neff, S2, float(resid)


def fit(groups, starts):
    """Stacked WLS over row groups that share SNP order; full and delete-block fits."""
    p = groups[0][0].shape[1]
    A, b = np.zeros((starts.size, p, p)), np.zeros((starts.size, p))
    for X, y, w in groups:
        Xw = X * w[:, None]
        A += np.add.reduceat(np.einsum("ni,nj->nij", Xw, X), starts, axis=0)
        b += np.add.reduceat(Xw * y[:, None], starts, axis=0)
    At, bt = A.sum(0), b.sum(0)
    return np.linalg.solve(At, bt), np.linalg.solve(At[None] - A, (bt[None] - b)[..., None])[..., 0]


def se(dele):
    k = dele.shape[0]
    return float(np.sqrt((k - 1) / k * np.sum((dele - dele.mean()) ** 2)))


def est(full, dele):
    return {"estimate": float(full), "se": se(np.asarray(dele))}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--score-dir", help="directory with genome-wide chromosome score caches")
    ap.add_argument("--out", help="write the estimates to JSON")
    args = ap.parse_args()
    if args.score_dir:
        G.GK.OUT = args.score_dir
    t0 = time.time()
    names = ("EUR", "EAS")
    paths = {a: os.path.join(G.GBMI_DIR, G.PATTERN.format(G.STRATA[a])) for a in names}
    for a in names:
        if G.sha256(paths[a]) != G.SHA256[a]:
            raise ValueError(f"{a} differs from the recorded asthma input")
    D = G.load_scores()
    M, S, J = D["M"], D["scores"], D["rsid"].size
    index = {rs: j for j, rs in enumerate(D["rsid"])}
    pa = [G.ANC.index(a) for a in names]

    Z, N, W = np.full((2, J), np.nan), np.zeros((2, J)), np.zeros((2, J))
    Q, K, S2 = np.full((2, J), np.nan), np.zeros((2, J)), np.full((2, J), np.nan)
    datasets = {}
    for s, a in enumerate(names):
        rows, ss, _ = G.read_gbmi(paths[a], index)
        if s == 0:
            ref0, alt0 = np.full(J, "", dtype=object), np.full(J, "", dtype=object)
            ref0[rows], alt0[rows] = ss.oa, ss.ea
            reference = ss.as_variants()
            reference.pos = np.zeros(len(ss), dtype=np.int64)
            reference.freeze()
            rows0 = rows
            Z[0, rows], N[0, rows], W[0, rows] = ss.z, ss.n_eff, 1.0 / ss.se ** 2
        else:
            ok = ((ref0[rows] == ss.oa) & (alt0[rows] == ss.ea)) | ((ref0[rows] == ss.ea) & (alt0[rows] == ss.oa))
            h = allelign.harmonize(ss.subset(ok), reference, drop_ambiguous=False)
            r = rows0[h.var_index]
            Z[1, r], N[1, r], W[1, r] = h.beta / h.se, h.n_eff, 1.0 / h.se ** 2
        hp, pres, nc, nk, nd = extras(paths[a], index, J)
        have = np.isfinite(Z[s]) & (pres != "")
        case, ctrl, neff, s2, resid = dataset_sizes(pres, nc, nk, have)
        S2[s, have], K[s] = s2, nd
        with np.errstate(invalid="ignore"):
            Q[s] = np.where((nd >= 2) & np.isfinite(hp) & (hp > 0), chi2.isf(np.clip(hp, 1e-300, 1.0), np.maximum(nd - 1, 1)), np.nan)
        datasets[a] = {"cases": np.round(case).astype(int).tolist(), "controls": np.round(ctrl).astype(int).tolist(),
                       "n_eff": np.round(neff).astype(int).tolist(), "max_count_residual": resid}
        print(f"[{time.time() - t0:4.0f}s] {a}: {int(have.sum())} variants, {len(case)} datasets, "
              f"largest count residual {resid:.2f}; dataset effective N {datasets[a]['n_eff']}", flush=True)

    sel = G.GK.regression_set(D) & D["available"][pa[0]] & D["available"][pa[1]]
    for s in (0, 1):
        sel &= np.isfinite(Z[s]) & (N[s] > 0) & np.isfinite(Q[s]) & np.isfinite(S2[s])
    for s in (0, 1):
        sel &= G.chi2_ok(Z, N, s, sel)
    idx = np.flatnonzero(sel)
    n = idx.size
    starts = G.block_starts(n, N_JK)
    zE, zF, NE, NF = Z[0, idx], Z[1, idx], N[0, idx], N[1, idx]
    lEE, lFF, lEF = S[pa[0], pa[0]][idx], S[pa[1], pa[1]][idx], S[pa[0], pa[1]][idx]
    xE, xF, xX = NE * lEE / M, NF * lFF / M, np.sqrt(NE * NF) * lEF / M
    dW = np.maximum(lEE, 1.0)
    one, zero = np.ones(n), np.zeros(n)
    out = {"n_snps": int(n), "M": int(M), "median_n_eff": [float(np.median(NE)), float(np.median(NF))],
           "mean_scores": {"EUR": float(lEE.mean()), "EAS": float(lFF.mean()), "cross": float(lEF.mean())},
           "datasets": datasets}
    print(f"{n} SNPs, M = {M}; median effective N {np.median(NE):.0f} and {np.median(NF):.0f}")

    # ---- both strata
    XE, XF, XX = np.column_stack((one, xE)), np.column_stack((one, xF)), np.column_stack((one, xX))
    muE = np.maximum(XE @ fit([(XE, zE ** 2, 1 / dW)], starts)[0], 1.0)
    muF = np.maximum(XF @ fit([(XF, zF ** 2, 1 / dW)], starts)[0], 1.0)
    q0 = XX @ fit([(XX, zE * zF, 1 / dW)], starts)[0]
    hE, hEd = fit([(XE, zE ** 2, 1 / (dW * muE ** 2))], starts)
    hF, hFd = fit([(XF, zF ** 2, 1 / (dW * muF ** 2))], starts)
    g, gd = fit([(XX, zE * zF, 1 / (dW * (muE * muF + q0 ** 2)))], starts)
    naive, naived = g[1] / np.sqrt(hE[1] * hF[1]), gd[:, 1] / np.sqrt(hEd[:, 1] * hFd[:, 1])
    out["strata"] = {"h2_EUR": est(hE[1], hEd[:, 1]), "h2_EAS": est(hF[1], hFd[:, 1]), "g": est(g[1], gd[:, 1]),
                     "correlation": est(naive, naived)}
    print(f"A. both strata: h2 {hE[1]:.4f} ({se(hEd[:, 1]):.4f}) and {hF[1]:.4f} ({se(hFd[:, 1]):.4f}); "
          f"g {g[1]:.4f} ({se(gd[:, 1]):.4f}); correlation {naive:.3f} ({se(naived):.3f})")

    # ---- Part A: the rebuilt two-ancestry meta, under two weightings
    for label, wE, wF in (("inverse_variance", W[0, idx], W[1, idx]), ("sample_size", NE, NF)):
        cE, cF = np.sqrt(wE / (wE + wF)), np.sqrt(wF / (wE + wF))
        Zm, Dm = cE * zE + cF * zF, cF * zE - cE * zF
        Qm = Dm ** 2
        XZ, XQ, _ = zq_design(cE, cF, xE, xF, xX, "shared")         # the recorded 2026-11-20 form
        f0, _ = fit([(XZ, Zm ** 2, 1 / dW), (XQ, Qm, 1 / dW)], starts)
        f, d = fit([(XZ, Zm ** 2, 1 / (dW * np.maximum(XZ @ f0, 1) ** 2)), (XQ, Qm, 1 / (dW * np.maximum(XQ @ f0, 1) ** 2))], starts)
        XZ3, XQ3, form3 = zq_design(cE, cF, xE, xF, xX, "three")
        f30, _ = fit([(XZ3, Zm ** 2, 1 / dW), (XQ3, Qm, 1 / dW)], starts)
        f3, d3 = fit([(XZ3, Zm ** 2, 1 / (dW * np.maximum(XZ3 @ f30, 1) ** 2)), (XQ3, Qm, 1 / (dW * np.maximum(XQ3 @ f30, 1) ** 2))], starts)
        Xz = XZ[:, [0, 2, 3, 4]]
        f0, _ = fit([(Xz, Zm ** 2, 1 / dW)], starts)
        fz, dz = fit([(Xz, Zm ** 2, 1 / (dW * np.maximum(Xz @ f0, 1) ** 2))], starts)

        def expected_q(a, b, c, given_z=False):
            vZ = 1 + cE ** 2 * xE * a + cF ** 2 * xF * b + 2 * cE * cF * xX * c
            vD = 1 + cF ** 2 * xE * a + cE ** 2 * xF * b - 2 * cE * cF * xX * c
            if not given_z:
                return vD
            cv = cE * cF * (xE * a - xF * b) + (cF ** 2 - cE ** 2) * xX * c
            return vD - cv ** 2 / vZ + (cv / vZ) ** 2 * Zm ** 2

        top = Zm ** 2 > GWS
        res = {"EUR_weight_share_median": float(np.median(cE ** 2)),
               "Z_and_Q": {"intercepts": [float(f[0]), float(f[1])], "h2_EUR": est(f[2], d[:, 2]), "h2_EAS": est(f[3], d[:, 3]), "g": est(f[4], d[:, 4]),
                           "correlation": est(f[4] / np.sqrt(f[2] * f[3]), d[:, 4] / np.sqrt(d[:, 2] * d[:, 3])),
                           "g_minus_strata": est(f[4] - g[1], d[:, 4] - gd[:, 1])},
               "Z_and_Q_three_intercepts": {"intercept_form": form3, "intercepts": [float(v) for v in f3[:-3]],
                                            "h2_EUR": est(f3[-3], d3[:, -3]), "h2_EAS": est(f3[-2], d3[:, -2]), "g": est(f3[-1], d3[:, -1]),
                                            "g_minus_strata": est(f3[-1] - g[1], d3[:, -1] - gd[:, 1]),
                                            "g_minus_shared_form": est(f3[-1] - f[4], d3[:, -1] - d[:, 4])},
               "Z_only": {"h2_EUR": est(fz[1], dz[:, 1]), "h2_EAS": est(fz[2], dz[:, 2]), "g": est(fz[3], dz[:, 3])},
               "mean_Z2": float(np.mean(Zm ** 2)), "mean_Q": float(Qm.mean()), "share_Q_above_3.84": float(np.mean(Qm > 3.841459)),
               "mean_Q_expected_identical_effects": float(expected_q(hE[1], hE[1], hE[1]).mean()),
               "mean_Q_expected_from_strata": float(expected_q(hE[1], hF[1], g[1]).mean()),
               "hits": {"n": int(top.sum()), "mean_Q": float(Qm[top].mean()), "share_Q_above_3.84": float(np.mean(Qm[top] > 3.841459)),
                        "mean_Q_expected_identical_effects": float(expected_q(hE[1], hE[1], hE[1], True)[top].mean()),
                        "mean_Q_expected_from_strata": float(expected_q(hE[1], hF[1], g[1], True)[top].mean())}}
        out[f"meta_{label}"] = res
        zq, zo, z3 = res["Z_and_Q"], res["Z_only"], res["Z_and_Q_three_intercepts"]
        print(f"A. {label} meta, three stratum intercepts ({z3['intercept_form']} form): g {z3['g']['estimate']:.4f} ({z3['g']['se']:.4f}), "
              f"difference from both strata {z3['g_minus_strata']['estimate']:.4f} ({z3['g_minus_strata']['se']:.4f})")
        print(f"A. {label} meta (median EUR weight {res['EUR_weight_share_median']:.3f}): Z and Q give g {zq['g']['estimate']:.4f} ({zq['g']['se']:.4f}), "
              f"difference from both strata {zq['g_minus_strata']['estimate']:.4f} ({zq['g_minus_strata']['se']:.4f}), "
              f"correlation {zq['correlation']['estimate']:.3f} ({zq['correlation']['se']:.3f}); Z only gives g {zo['g']['estimate']:.4f} ({zo['g']['se']:.4f}), "
              f"h2_EAS {zo['h2_EAS']['estimate']:.3f} ({zo['h2_EAS']['se']:.3f})")
        print(f"   mean Q {res['mean_Q']:.4f}; expected {res['mean_Q_expected_identical_effects']:.4f} with identical effects and "
              f"{res['mean_Q_expected_from_strata']:.4f} from the strata estimates; at {res['hits']['n']} pooled hits mean Q {res['hits']['mean_Q']:.3f} "
              f"against {res['hits']['mean_Q_expected_identical_effects']:.3f} and {res['hits']['mean_Q_expected_from_strata']:.3f}")

    # ---- Part B: within-ancestry Q as a ceiling
    tau, taud = {}, {}
    out["within"] = {}
    for s, (a, x, ell) in enumerate(zip(names, (xE, xF), (lEE, lFF))):
        k1, s2, z2, q = K[s, idx] - 1.0, S2[s, idx], Z[s, idx] ** 2, Q[s, idx]
        dA = np.maximum(ell, 1.0)
        XZ = np.column_stack((one, zero, x, x * s2))                # aZ, aQ per df, tau, h2 - tau
        XQ = np.column_stack((zero, k1, zero, x * (1 - s2)))
        f0, _ = fit([(XZ, z2, 1 / dA), (XQ, q, 1 / (dA * k1))], starts)
        f, d = fit([(XZ, z2, 1 / (dA * np.maximum(XZ @ f0, 1) ** 2)), (XQ, q, k1 / (dA * np.maximum(XQ @ f0, 1) ** 2))], starts)
        tau[a], taud[a] = f[2], d[:, 2]
        kk = k1.astype(int) + 1
        out["within"][a] = {"datasets_mode": int(np.bincount(kk).argmax()), "share_at_mode": float(np.mean(kk == np.bincount(kk).argmax())),
                            "S2_median": float(np.median(s2)), "mean_Q_per_df": float(np.mean(q / k1)), "Q_intercept_per_df": float(f[1]),
                            "tau": est(f[2], d[:, 2]), "h2_minus_tau": est(f[3], d[:, 3]),
                            "ceiling": est(f[2] / (f[2] + f[3]), d[:, 2] / (d[:, 2] + d[:, 3]))}
        w = out["within"][a]
        print(f"B. {a}: {w['datasets_mode']} datasets at {w['share_at_mode']:.3f} of SNPs, median S2 {w['S2_median']:.3f}; mean Q per df {w['mean_Q_per_df']:.4f}; "
              f"tau {w['tau']['estimate']:.4f} ({w['tau']['se']:.4f}); h2 - tau {w['h2_minus_tau']['estimate']:.4f} ({w['h2_minus_tau']['se']:.4f}); "
              f"ceiling {w['ceiling']['estimate']:.3f} ({w['ceiling']['se']:.3f})")
    both, bothd = g[1] / np.sqrt(tau["EUR"] * tau["EAS"]), gd[:, 1] / np.sqrt(taud["EUR"] * taud["EAS"])
    eur, eurd = g[1] / np.sqrt(tau["EUR"] * hF[1]), gd[:, 1] / np.sqrt(taud["EUR"] * hFd[:, 1])
    out["correlation"] = {"naive": est(naive, naived), "corrected_EUR_ceiling": est(eur, eurd), "corrected_both_ceilings": est(both, bothd),
                          "corrected_minus_naive": est(both - naive, bothd - naived)}
    c = out["correlation"]
    print(f"B. EUR-EAS correlation: naive {naive:.3f} ({c['naive']['se']:.3f}); EUR ceiling only {eur:.3f} ({c['corrected_EUR_ceiling']['se']:.3f}); "
          f"both ceilings {both:.3f} ({c['corrected_both_ceilings']['se']:.3f}); corrected minus naive {both - naive:.3f} ({c['corrected_minus_naive']['se']:.3f})")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
