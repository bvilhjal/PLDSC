#!/usr/bin/env python3
"""Cochran's Q as a second GWAS summary: known-truth simulation on 1000 Genomes LD.

Question. A two-ancestry meta-analysis reports Z = cE zE + cF zF. The
heterogeneity contrast D = cF zE - cE zF is the other half of the same
rotation and Cochran's Q = D^2. Part 1 asks whether a release holding only Z,
Q and the weights recovers the cross-ancestry genetic covariance g, and
whether it can tell different LD from different effects. Part 2 asks whether
each ancestry's own within-ancestry Q can serve as a ceiling: cohorts of one
ancestry are replicate measurements of a phenotype, so their heterogeneity
estimates how far a cross-ancestry correlation falls below one for reasons
unrelated to ancestry.

Design. Genotypes are the bigsnpr 1000 Genomes set used by
multiancestry_genome_kappa.py. SNPs with MAF >= 0.05 in both groups are cut
into consecutive 250-SNP blocks; LD is the in-panel correlation within a
block and zero between blocks, so every score is known exactly. Standardized
effects have variance h2/M (the PLDSC baseline). Association statistics are
exact Gaussian draws, z = sqrt(N) R b + e with cov(e) = R. The meta-analysis
is inverse-variance weighted on the per-allele scale. All fits are two-step
weighted least squares with free intercepts and a 200-group jackknife over
LD blocks. The full run uses chromosomes 14-22. Total N is set so that an
all-EUR GWAS of that size would have mean chi-square near 2.

Rules, stated before the full run (2026-11-20).
  R1 bias. With exact scores, in the infinitesimal and sparse arms, the mean
     Z+Q estimate of g lies within max(3 Monte Carlo SE, 5% of g) of truth.
  R2 calibration. At rho = 1 with equal h2, exact scores and infinitesimal
     effects, the one-sided 5% jackknife test of rho < 1 from Z+Q rejects in
     at most 9% of 200 replicates, and mean jackknife SE of g divided by its
     empirical SD lies in [0.8, 1.25].
  R3 ceiling. With exchangeable cohorts the median corrected correlation is
     within 0.05 of the true ancestry correlation, and the median naive
     correlation is within 0.03 of its analytic attenuated value.
Rule added 2026-11-23, before the unequal-intercept arm was first run.
  R4 intercepts. In the arm with stratum intercepts 1.08 and 1.0, the mean
     three-intercept (ZQ3) estimate of g lies within max(3 Monte Carlo SE,
     5% of g) of truth, as in R1.
Everything else is descriptive: power, the Z^2-only fit, shared per-allele
effects, unequal heritabilities, the EUR-AFR pair, cohorts that share a
definition, cohorts whose effects are proportional but on different scales,
and reference scores from an independent half panel. The
half-panel arm is a stress test: both the generating LD and the reference
come from about 250 people, so the mismatch variance is roughly four times
that of a 500-person panel against a large GWAS.

Run from the repository root:
  python experiments/q_heterogeneity_sim.py --out report/q_heterogeneity_simulation.json
  python experiments/q_heterogeneity_sim.py --quick        # chr22, 10 replicates

Outcome (2026-11-20; 234,000 SNPs in 936 blocks, total N 46,475, about 15
minutes; report/methods.md Section 11). The rules were in this docstring
before the full run but were not committed separately first.
  R1 PASS. The Z+Q mean of g is within 2% of truth in all ten cells.
  R2 PASS. At rho = 1 the test rejects in 4.5% (85% EUR) and 4.0% (50%);
     the SE ratio is 0.96 and 1.01.
  R3 PASS. Corrected medians are within 0.042 of truth and naive medians
     within 0.021 of their expectation.
At 85% EUR the SD of g from Z+Q is 1.7 times that from both strata, and Z^2
alone 6.7 times; at 50% EUR Z+Q equals the strata. A mean-Q test rejects in
98-100% of replicates with identical effects. Shared per-allele effects bias
g by -9% (Z+Q) and +7% (strata). The two scale cells were added after the
first full run and carry their own seed keys, so the earlier draws are fixed.

Review, 2026-11-23. Every cell above has both stratum intercepts at one, so
it cannot show that the Z+Q fit with one free intercept per equation lets
unequal stratum intercepts leak into g when the inverse-variance weights vary
across SNPs (q_model.py). Each cell also reports ZQ3, the fit with three
stratum intercepts, and one arm has intercepts 1.08 and 1.0 (rule R4). Both
were first run on 2026-11-23; the outcome is in experiments/RESULTS.md.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import multiancestry_genome_kappa as GK  # noqa: E402
from q_model import zq_design  # noqa: E402

BLOCK, N_JK, H2, SEED, PI_SPARSE = 250, 200, 0.3, 20261120, 0.02
CRIT = 1.6448536269514722


# ------------------------------------------------------------------ genotypes
def load_groups(groups, chroms):
    """Standardized genotypes (blocks, n, BLOCK) per group, and genotype SDs."""
    chrom, _, _ = GK.read_bim(GK.BFILE)
    n, m_total = len(next(iter(groups.values()))[1]), chrom.size
    X = {g: [] for g in groups}
    sd = {g: [] for g in groups}
    for c in chroms:
        j = np.flatnonzero(chrom == c)
        G = GK.read_bed(GK.BFILE, n, m_total, int(j[0]), int(j[-1]) + 1)
        sub = {g: G[rows] for g, (rows, _) in groups.items()}
        p = {g: s.mean(axis=0, dtype=np.float64) / 2.0 for g, s in sub.items()}
        keep = np.ones(j.size, dtype=bool)
        for g in groups:
            keep &= np.minimum(p[g], 1.0 - p[g]) >= 0.05
        idx = np.flatnonzero(keep)
        nb = idx.size // BLOCK
        idx = idx[:nb * BLOCK]
        for g in groups:
            A = sub[g][:, idx].astype(np.float32)
            s = A.std(axis=0)
            A = (A - A.mean(axis=0)) / s
            X[g].append(np.ascontiguousarray(A.reshape(A.shape[0], nb, BLOCK).transpose(1, 0, 2)))
            sd[g].append(np.sqrt(2.0 * p[g][idx] * (1.0 - p[g][idx])))
        del G, sub
    return {g: np.concatenate(X[g]) for g in groups}, {g: np.concatenate(sd[g]) for g in groups}


def ld_scores(Xa, Xb, *, correct=False, chunk=64):
    """Within-a, within-b and signed cross scores, exact for the block LD.

    With correct=True the within-panel squared correlations get PLDSC's
    (1 - r^2)/(n - 2) adjustment, as for a finite reference panel."""
    na, nb = Xa.shape[1], Xb.shape[1]
    out = [np.empty(Xa.shape[0] * BLOCK) for _ in range(3)]
    for s in range(0, Xa.shape[0], chunk):
        Ra = np.matmul(Xa[s:s + chunk].transpose(0, 2, 1), Xa[s:s + chunk]).astype(np.float64) / na
        Rb = np.matmul(Xb[s:s + chunk].transpose(0, 2, 1), Xb[s:s + chunk]).astype(np.float64) / nb
        A, B = Ra * Ra, Rb * Rb
        if correct:
            A, B = A - (1.0 - A) / (na - 2), B - (1.0 - B) / (nb - 2)
        sl = slice(s * BLOCK, (s + Ra.shape[0]) * BLOCK)
        out[0][sl], out[1][sl], out[2][sl] = A.sum(-1).ravel(), B.sum(-1).ravel(), (Ra * Rb).sum(-1).ravel()
    return out


def draw_z(X, beta, N, rng):
    """Exact z for standardized effects beta (blocks, BLOCK) at sample size N."""
    n = X.shape[1]
    s = np.matmul(X, beta[..., None])[..., 0]
    e = rng.standard_normal(s.shape, dtype=np.float32)
    v = (np.sqrt(N) / n) * s + e / np.sqrt(n)
    return np.matmul(X.transpose(0, 2, 1), v[..., None])[..., 0].ravel().astype(np.float64)


# ----------------------------------------------------------------- regression
def fit(groups, starts):
    """Stacked WLS over row groups that share SNP order; full and delete-group fits."""
    p = groups[0][0].shape[1]
    A, b = np.zeros((starts.size, p, p)), np.zeros((starts.size, p))
    for X, y, w in groups:
        Xw = X * w[:, None]
        A += np.add.reduceat(np.einsum("ni,nj->nij", Xw, X), starts, axis=0)
        b += np.add.reduceat(Xw * y[:, None], starts, axis=0)
    At, bt = A.sum(0), b.sum(0)
    return np.linalg.solve(At, bt), np.linalg.solve(At[None] - A, (bt[None] - b)[..., None])[..., 0]


def jk_se(dele):
    k = dele.shape[0]
    return np.sqrt((k - 1) / k * np.sum((dele - dele.mean(0)) ** 2, axis=0))


def slope(x, y, dW, starts, cross=None):
    """Two-step weighted slope of y on x with a free intercept."""
    X = np.column_stack((np.ones(x.size), x))
    full, _ = fit([(X, y, 1.0 / dW)], starts)
    mu = X @ full
    w = 1.0 / (dW * (np.maximum(mu, 1.0) ** 2 if cross is None else np.maximum(cross + mu * mu, 1e-6)))
    full, dele = fit([(X, y, w)], starts)
    return full, dele, X


def derived(hE, hF, g):
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where(hE * hF > 0, g / np.sqrt(np.abs(hE * hF)), np.nan)
    return rho, hE + hF - 2.0 * g


def summarise(est, truth_g):
    """est: replicate rows (g, se_g, rho, se_rho, v, se_v)."""
    e = np.asarray(est)
    ok = np.isfinite(e[:, 2]) & np.isfinite(e[:, 3])
    q = np.nanpercentile(e[ok, 2], [25, 50, 75]) if ok.any() else [np.nan] * 3
    return {
        "g_mean": float(e[:, 0].mean()), "g_sd": float(e[:, 0].std(ddof=1)), "g_mcse": float(e[:, 0].std(ddof=1) / np.sqrt(len(e))),
        "g_jk_se_mean": float(e[:, 1].mean()), "g_truth": float(truth_g),
        "rho_median": float(q[1]), "rho_robust_sd": float((q[2] - q[0]) / 1.349), "rho_invalid": int((~ok).sum()),
        "reject_rho_lt_1": float(np.mean(((1.0 - e[ok, 2]) / e[ok, 3]) > CRIT)) if ok.any() else float("nan"),
        "reject_v_gt_0": float(np.mean(e[:, 4] / e[:, 5] > CRIT)),
    }


# ------------------------------------------------------------- part 1: Z and Q
def part1_cell(P, shareE, rho, arch, reps, rng, *, hE=H2, hF=H2, ref=None, intE=1.0, intF=1.0):
    """One Part 1 cell. ``intE``/``intF`` >= 1 add independent noise so the stratum
    intercepts differ (intE = 1.08 is the unequal-intercept arm). ``ZQ`` is the
    shared-intercept fit recorded on 2026-11-20; ``ZQ3`` fits the three
    stratum intercepts (q_model.py)."""
    if intE < 1.0 or intF < 1.0:
        raise ValueError("stratum intercepts below one are not simulated")
    XE, XF, sdE, sdF = P["XE"], P["XF"], P["sdE"], P["sdF"]
    B = XE.shape[0]
    M = B * BLOCK
    NE, NF = shareE * P["N"], (1.0 - shareE) * P["N"]
    rEE, rFF, rEF = P["scores"] if ref is None else ref
    wE, wF = NE * sdE ** 2, NF * sdF ** 2
    cE, cF = np.sqrt(wE / (wE + wF)), np.sqrt(wF / (wE + wF))
    xE, xF, xX = NE * rEE / M, NF * rFF / M, np.sqrt(NE * NF) * rEF / M
    dW = np.maximum(rEE, 1.0)
    starts = (np.arange(N_JK) * B // N_JK) * BLOCK
    one = np.ones(M)
    XZ, XQ, _ = zq_design(cE, cF, xE, xF, xX, "shared")
    XZ3, XQ3, form3 = zq_design(cE, cF, xE, xF, xX, "three")
    sE3, sF3 = sdE.reshape(B, BLOCK).astype(np.float32), sdF.reshape(B, BLOCK).astype(np.float32)
    rows = {"strata": [], "ZQ": [], "ZQ3": [], "Z": []}
    meanQ, meanZ2, naive = [], [], []
    for _ in range(reps):
        u, v = rng.standard_normal((2, B, BLOCK), dtype=np.float32)
        aE, aF = u, rho * u + np.sqrt(max(1.0 - rho * rho, 0.0)) * v
        if arch == "sparse":
            mask = (rng.random((B, BLOCK)) < PI_SPARSE) / np.sqrt(PI_SPARSE)
            bE, bF = np.sqrt(hE / M) * aE * mask, np.sqrt(hF / M) * aF * mask
        elif arch == "allele":          # correlated per-allele effects, constant per-allele variance
            bE, bF = aE * sE3 * np.sqrt(hE / np.sum(sdE ** 2)), aF * sF3 * np.sqrt(hF / np.sum(sdF ** 2))
        else:
            bE, bF = np.sqrt(hE / M) * aE, np.sqrt(hF / M) * aF
        zE, zF = draw_z(XE, bE.astype(np.float32), NE, rng), draw_z(XF, bF.astype(np.float32), NF, rng)
        if intE > 1.0 or intF > 1.0:        # stratum-specific inflation, independent across SNPs and strata
            zE = zE + np.sqrt(intE - 1.0) * rng.standard_normal(M)
            zF = zF + np.sqrt(intF - 1.0) * rng.standard_normal(M)
        Z, D = cE * zE + cF * zF, cF * zE - cE * zF
        Q = D * D
        meanQ.append(Q.mean()); meanZ2.append(np.mean(Z * Z))
        bq = np.add.reduceat(Q - 1.0, starts) / np.add.reduceat(one, starts)
        naive.append(bq.mean() / (bq.std(ddof=1) / np.sqrt(N_JK)) > CRIT)
        # both strata
        fE, dE, X1 = slope(xE, zE * zE, dW, starts)
        fF, dF, X2 = slope(xF, zF * zF, dW, starts)
        fX, dX, _ = slope(xX, zE * zF, dW, starts, cross=np.maximum(X1 @ fE, 1.0) * np.maximum(X2 @ fF, 1.0))
        r, vv = derived(fE[1], fF[1], fX[1]); rd, vd = derived(dE[:, 1], dF[:, 1], dX[:, 1])
        rows["strata"].append((fX[1], jk_se(dX[:, 1]), r, jk_se(rd), vv, jk_se(vd)))
        # Z^2 and Q
        f0, _ = fit([(XZ, Z * Z, 1.0 / dW), (XQ, Q, 1.0 / dW)], starts)
        f, d = fit([(XZ, Z * Z, 1.0 / (dW * np.maximum(XZ @ f0, 1.0) ** 2)), (XQ, Q, 1.0 / (dW * np.maximum(XQ @ f0, 1.0) ** 2))], starts)
        r, vv = derived(f[2], f[3], f[4]); rd, vd = derived(d[:, 2], d[:, 3], d[:, 4])
        rows["ZQ"].append((f[4], jk_se(d[:, 4]), r, jk_se(rd), vv, jk_se(vd)))
        # Z^2 and Q with three stratum intercepts
        f0, _ = fit([(XZ3, Z * Z, 1.0 / dW), (XQ3, Q, 1.0 / dW)], starts)
        f, d = fit([(XZ3, Z * Z, 1.0 / (dW * np.maximum(XZ3 @ f0, 1.0) ** 2)), (XQ3, Q, 1.0 / (dW * np.maximum(XQ3 @ f0, 1.0) ** 2))], starts)
        r, vv = derived(f[-3], f[-2], f[-1]); rd, vd = derived(d[:, -3], d[:, -2], d[:, -1])
        rows["ZQ3"].append((f[-1], jk_se(d[:, -1]), r, jk_se(rd), vv, jk_se(vd)))
        # Z^2 only
        Xz = XZ[:, [0, 2, 3, 4]]
        f0, _ = fit([(Xz, Z * Z, 1.0 / dW)], starts)
        f, d = fit([(Xz, Z * Z, 1.0 / (dW * np.maximum(Xz @ f0, 1.0) ** 2))], starts)
        r, vv = derived(f[1], f[2], f[3]); rd, vd = derived(d[:, 1], d[:, 2], d[:, 3])
        rows["Z"].append((f[3], jk_se(d[:, 3]), r, jk_se(rd), vv, jk_se(vd)))
    g = rho * np.sqrt(hE * hF) if arch != "allele" else rho * np.sqrt(hE * hF) * np.sum(sdE * sdF) / np.sqrt(np.sum(sdE ** 2) * np.sum(sdF ** 2))
    out = {k: summarise(v, g) for k, v in rows.items()}
    out.update(mean_Q=float(np.mean(meanQ)), mean_Z2=float(np.mean(meanZ2)), naive_Q_reject=float(np.mean(naive)),
               N_E=float(NE), N_F=float(NF), reps=reps,
               ZQ3_intercept_form=form3, stratum_intercepts=[float(intE), float(intF)])
    return out


# --------------------------------------------------- part 2: within-ancestry Q
def part2_cell(P, shareE, fE, fF, rw, ranc, reps, rng, *, shared_pair=False, scales=None):
    XE, XF = P["XE"], P["XF"]
    B = XE.shape[0]
    M = B * BLOCK
    NE, NF = shareE * P["N"], (1.0 - shareE) * P["N"]
    lEE, lFF, lEF = P["scores"]
    xE, xF, xX = NE * lEE / M, NF * lFF / M, np.sqrt(NE * NF) * lEF / M
    dW = np.maximum(lEE, 1.0)
    starts = (np.arange(N_JK) * B // N_JK) * BLOCK
    S2E, S2F = float(np.sum(np.square(fE))), float(np.sum(np.square(fF)))
    sc = np.float32(np.sqrt(H2 / M))

    sE = np.ones(len(fE)) if scales is None else np.sqrt(np.asarray(scales, dtype=float))

    def cohorts(X, base, f, Na, tie, s=None):
        Zs, SS, shared = 0.0, 0.0, None
        for c, fc in enumerate(f):
            if tie and c < 2:
                shared = rng.standard_normal((B, BLOCK), dtype=np.float32) if shared is None else shared
                dev = shared
            else:
                dev = rng.standard_normal((B, BLOCK), dtype=np.float32)
            sc_c = sc if s is None else sc * np.float32(s[c])
            z = draw_z(X, sc_c * (np.float32(np.sqrt(rw)) * base + np.float32(np.sqrt(1.0 - rw)) * dev), fc * Na, rng)
            Zs, SS = Zs + np.sqrt(fc) * z, SS + z * z
        return Zs, SS - Zs * Zs

    rows = []
    for _ in range(reps):
        u, v = rng.standard_normal((2, B, BLOCK), dtype=np.float32)
        ZE, QE = cohorts(XE, u, fE, NE, shared_pair, sE)
        ZF, QF = cohorts(XF, np.float32(ranc) * u + np.float32(np.sqrt(max(1.0 - ranc * ranc, 0.0))) * v, fF, NF, False)
        fmE, dmE, X1 = slope(xE, ZE * ZE, dW, starts)
        fmF, dmF, X2 = slope(xF, ZF * ZF, dW, starts)
        fHE, dHE, _ = slope(xE, QE, dW, starts)
        fHF, dHF, _ = slope(xF, QF, dW, starts)
        fG, dG, _ = slope(xX, ZE * ZF, dW, starts, cross=np.maximum(X1 @ fmE, 1.0) * np.maximum(X2 @ fmF, 1.0))

        def combine(mE, mF, HE, HF, G):
            tE, tF = mE - HE * S2E / (1.0 - S2E), mF - HF * S2F / (1.0 - S2F)
            with np.errstate(invalid="ignore", divide="ignore"):
                naive = G / np.sqrt(mE * mF)
                corr = np.where(tE * tF > 0, G / np.sqrt(np.abs(tE * tF)), np.nan)
                ceilE, ceilF = tE / (tE + HE / (1.0 - S2E)), tF / (tF + HF / (1.0 - S2F))
            return naive, corr, ceilE, ceilF

        full = combine(fmE[1], fmF[1], fHE[1], fHF[1], fG[1])
        dele = combine(dmE[:, 1], dmF[:, 1], dHE[:, 1], dHF[:, 1], dG[:, 1])
        rows.append((*full, jk_se(dele[1]), jk_se(dele[0])))
    e = np.asarray(rows)
    tau = rw * H2
    fs = float(np.dot(fE, sE))
    hmE = tau * fs ** 2 + (H2 - tau) * float(np.dot(np.square(fE), sE ** 2))
    naive_truth = fs * tau * ranc / np.sqrt(hmE * (tau + (H2 - tau) * S2F))
    med = lambda x: float(np.nanmedian(x))
    rsd = lambda x: float(np.subtract(*np.nanpercentile(x, [75, 25])) / 1.349)
    return {"naive_median": med(e[:, 0]), "naive_truth": float(naive_truth), "naive_robust_sd": rsd(e[:, 0]),
            "corrected_median": med(e[:, 1]), "corrected_truth": float(ranc), "corrected_robust_sd": rsd(e[:, 1]),
            "corrected_jk_se_median": med(e[:, 4]), "corrected_invalid": int(np.isnan(e[:, 1]).sum()),
            "ceiling_E_median": med(e[:, 2]), "ceiling_F_median": med(e[:, 3]), "ceiling_truth": float(rw),
            "reject_corrected_lt_1": float(np.nanmean((1.0 - e[:, 1]) / e[:, 4] > CRIT)),
            "reject_naive_lt_1": float(np.nanmean((1.0 - e[:, 0]) / e[:, 5] > CRIT)), "reps": reps}


# ----------------------------------------------------------------------- main
def panel(sp, a, b, chroms, rng=None):
    groups = {}
    for name, lab in (("E", a), ("F", b)):
        rows = np.flatnonzero(sp == lab)
        if rng is None:
            groups[name] = (rows, sp)
        else:
            perm = rng.permutation(rows)
            groups[name + "gen"], groups[name + "ref"] = (np.sort(perm[:rows.size // 2]), sp), (np.sort(perm[rows.size // 2:]), sp)
    return load_groups(groups, chroms)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true", help="chr22 and 10 replicates per cell")
    ap.add_argument("--out", help="write all cell summaries to JSON")
    ap.add_argument("--cells", help="run only cells whose label matches this regular expression")
    args = ap.parse_args()
    chroms = [22] if args.quick else list(range(14, 23))
    reps_main, reps = (10, 10) if args.quick else (200, 100)
    sp = GK.read_superpops(GK.BFILE)
    t0 = time.time()
    res, cell = {"config": {"chroms": chroms, "block": BLOCK, "h2": H2, "seed": SEED, "jackknife_groups": N_JK}}, 0

    def run(label, fn, *a, key=None, **k):
        nonlocal cell
        if key is None:                 # later additions carry their own key, so earlier draws stay fixed
            cell += 1
            key = cell
        if args.cells and not re.search(args.cells, label):
            return
        res[label] = fn(*a, np.random.default_rng([SEED, key]), **k)
        print(f"[{time.time() - t0:5.0f}s] {label}: {json.dumps(res[label])}", flush=True)

    def setup(a, b, split_rng=None):
        X, sd = panel(sp, a, b, chroms, split_rng)
        if split_rng is None:
            P = {"XE": X["E"], "XF": X["F"], "sdE": sd["E"], "sdF": sd["F"]}
            P["scores"] = ld_scores(P["XE"], P["XF"])
            ref = None
        else:
            P = {"XE": X["Egen"], "XF": X["Fgen"], "sdE": sd["Egen"], "sdF": sd["Fgen"]}
            P["scores"] = ld_scores(P["XE"], P["XF"])
            ref = ld_scores(X["Eref"], X["Fref"], correct=True)
        M = P["XE"].shape[0] * BLOCK
        P["N"] = 1.0 * M / (H2 * P["scores"][0].mean())
        l = P["scores"]
        print(f"[{time.time() - t0:5.0f}s] {a}-{b}{' half panels' if ref else ''}: {M} SNPs in {M // BLOCK} blocks, n = {P['XE'].shape[1]} and {P['XF'].shape[1]}, "
              f"total N {P['N']:.0f}; mean scores {l[0].mean():.2f} {l[1].mean():.2f} {l[2].mean():.2f}; "
              f"corr(EE,EF) {np.corrcoef(l[0], l[2])[0, 1]:.3f}, corr(FF,EF) {np.corrcoef(l[1], l[2])[0, 1]:.3f}", flush=True)
        res[f"panel {a}-{b}{' half' if ref else ''}"] = {"M": M, "N": P["N"], "mean_scores": [float(x.mean()) for x in l]}
        return P, ref

    # Part 1, EUR-EAS, exact scores
    P, _ = setup("EUR", "EAS")
    for share in (0.85, 0.5):
        for rho in (1.0, 0.9, 0.74):
            run(f"P1 EUR-EAS share {share} rho {rho} inf", part1_cell, P, share, rho, "inf", reps_main)
    for rho in (1.0, 0.74):
        run(f"P1 EUR-EAS share 0.85 rho {rho} sparse", part1_cell, P, 0.85, rho, "sparse", reps)
    run("P1 EUR-EAS share 0.85 rho 1.0 allele", part1_cell, P, 0.85, 1.0, "allele", reps)
    run("P1 EUR-EAS share 0.85 rho 1.0 inf unequal h2", part1_cell, P, 0.85, 1.0, "inf", reps, hE=0.25, hF=0.35)
    # Unequal stratum intercepts (added 2026-11-23 after review; own seed key; rule R4).
    run("P1 EUR-EAS share 0.85 rho 0.74 inf intercepts 1.08 1.0", part1_cell, P, 0.85, 0.74, "inf", reps_main, key=1003, intE=1.08, intF=1.0)
    # Part 2, same panels
    fE, fF = (0.55, 0.25, 0.12, 0.08), (0.7, 0.3)
    for share in (0.85, 0.5):
        for rw, ranc in ((1.0, 1.0), (0.8, 1.0), (0.8, 0.85), (0.6, 1.0)):
            if share == 0.5 and (rw, ranc) not in ((0.8, 1.0), (0.8, 0.85)):
                continue
            run(f"P2 share {share} within {rw} ancestry {ranc}", part2_cell, P, share, fE, fF, rw, ranc, reps)
    run("P2 share 0.5 within 0.8 ancestry 1.0 two EUR cohorts share a definition", part2_cell, P, 0.5, fE, fF, 0.8, 1.0, reps, shared_pair=True)
    # Cohorts with proportional effects on different scales (added after the first full run).
    run("P2 share 0.5 within 1.0 ancestry 1.0 EUR cohort h2 scales 1.4 0.8 0.6 1.0", part2_cell, P, 0.5, fE, fF, 1.0, 1.0, reps, key=1001, scales=(1.4, 0.8, 0.6, 1.0))
    run("P2 share 0.5 within 1.0 ancestry 0.85 EUR cohort h2 scales 1.4 0.8 0.6 1.0", part2_cell, P, 0.5, fE, fF, 1.0, 0.85, reps, key=1002, scales=(1.4, 0.8, 0.6, 1.0))
    del P
    # Part 1, reference scores from an independent half panel
    P, ref = setup("EUR", "EAS", np.random.default_rng([SEED, 999]))
    for rho in (1.0, 0.74):
        run(f"P1 EUR-EAS share 0.85 rho {rho} inf half-panel reference", part1_cell, P, 0.85, rho, "inf", reps, ref=ref)
    del P, ref
    # Part 1, EUR-AFR
    P, _ = setup("EUR", "AFR")
    for rho in (1.0, 0.74):
        run(f"P1 EUR-AFR share 0.85 rho {rho} inf", part1_cell, P, 0.85, rho, "inf", reps)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
