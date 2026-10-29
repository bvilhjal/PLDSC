"""When is LDSC miscalibrated for mixed-model GWAS summary statistics?

Paper angle #1. A mixed-model (LMM) GWAS produces per-SNP statistics whose null
covariance is the LMM null, so the exact moment is
    E[z_j^2] = 1 + tau^2 * L_j^(P_c)          (docs/METHOD.md Eq M6, oracle V_c),
with L^(P_c) the P_c-metric LD score. A practitioner instead regresses chi^2 on
a CONVENTIONAL LD score (ordinary or covariate-adjusted, in-sample or from an
external reference panel). This maps the systematic bias in the LDSC intercept
and slope (-> h^2) as the confounding regime varies, and asks whether the
P_c-matched score removes it.

Design choices that keep the map honest:
- The bias is SYSTEMATIC (large-M), so we regress the *exact* analytic E[z_j^2]
  (no Monte Carlo) on each candidate LD score.
- N is held FIXED (=400) across every regime, and true h^2 (tau^2) is fixed
  within each sweep, so drift reflects the LD-score mismatch, not N or the model.
- The GWAS sample and the external reference share the SAME loci (same allele
  frequencies), differing only in individuals -- so the external panel is a fair
  reference that simply cannot contain the GWAS sample's relatedness.
- Two unambiguous readouts: the intercept (absolute, target 1) and the slope
  drift = slope(regime)/slope(anchor) within each sweep (target 1; <1 means h^2
  underestimated).

Run:  python experiments/map_ldsc_miscalibration.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from pldsc._lmm import ld_scores, residual_precision, score_pieces  # noqa: E402
from pldsc.regression import ld_score_regression  # noqa: E402
from pldsc.simulate import ar1_haps, bn_freqs, grm, standardize  # noqa: E402
from m4_fitted_vc import reml_fit  # noqa: E402

SEED = 20260718
CHROMS = [150, 250, 250]      # target chr 0; others build the LOCO GRM
NCHR = len(CHROMS)
N = 400                       # held fixed across ALL regimes
NREF = 400
RHO = 0.6
SG2, SE2 = 0.6, 0.4          # background variance components (h2_bg = 0.6)
WINDOW_CM = 2.0


# --------------------------------------------------------------------------- #
# matched-loci genotype generators (shared allele frequencies -> same SNPs)
# --------------------------------------------------------------------------- #
def _haps(rng, freqs_by_chr):
    """Full-genome haplotypes; freqs_by_chr[c] is an (n_hap, m_chr) matrix.
    AR(1) LD resets at each chromosome boundary. Batched over all haplotypes."""
    return np.hstack([ar1_haps(rng, freqs_by_chr[c], RHO) for c in range(NCHR)])


def _unrelated_homog(rng, freqs, n):
    if n <= 0:
        return np.empty((0, sum(CHROMS)))
    h = _haps(rng, [np.tile(freqs[c], (2 * n, 1)) for c in range(NCHR)])
    return h[0::2] + h[1::2]


def _unrelated_admixed(rng, p1, p2, n):
    a2 = np.repeat(rng.uniform(0.0, 1.0, n), 2)            # per-haplotype ancestry
    fbc = [a2[:, None] * p1[c][None, :] + (1 - a2)[:, None] * p2[c][None, :]
           for c in range(NCHR)]
    h = _haps(rng, fbc)
    return h[0::2] + h[1::2]


def _sibpairs_homog(rng, freqs, n_pairs):
    if n_pairs <= 0:
        return np.empty((0, sum(CHROMS)))
    # 4 parental haplotypes per pair, per chromosome (batched)
    phap = [ar1_haps(rng, np.tile(freqs[c], (4 * n_pairs, 1)), RHO) for c in range(NCHR)]
    rows = []
    for p in range(n_pairs):
        for _s in range(2):  # two full siblings; per-chromosome independent assortment
            rows.append(np.concatenate(
                [phap[c][4 * p + rng.integers(2)] + phap[c][4 * p + 2 + rng.integers(2)]
                 for c in range(NCHR)]))
    return np.array(rows)


def _positions():
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(CHROMS)])
    pos = np.concatenate([np.arange(m) * 0.15 for m in CHROMS])
    return chrom, pos


def build(regime, seed_off=0):
    rng = np.random.default_rng(np.random.PCG64(SEED + regime["salt"] + 1000 * seed_off))
    if regime["kind"] == "related":
        freqs = [rng.uniform(0.1, 0.9, m) for m in CHROMS]
        n_pairs = regime["n_pairs"]
        blocks = []
        if n_pairs > 0:
            blocks.append(_sibpairs_homog(rng, freqs, n_pairs))
        if N - 2 * n_pairs > 0:
            blocks.append(_unrelated_homog(rng, freqs, N - 2 * n_pairs))
        G = np.vstack(blocks)
        Gref = _unrelated_homog(rng, freqs, NREF)   # same loci, unrelated
    else:
        p1, p2 = zip(*[bn_freqs(rng, m, regime["fst"]) for m in CHROMS])
        G = _unrelated_admixed(rng, p1, p2, N)
        Gref = _unrelated_admixed(rng, p1, p2, NREF)  # same loci + structure
    chrom, pos = _positions()
    keep = (G.var(0) > 0) & (Gref.var(0) > 0)
    return G[:, keep], Gref[:, keep], chrom[keep], pos[keep]


# --------------------------------------------------------------------------- #
def corr_ldscore(Xstd, mask):
    Xn = Xstd / np.linalg.norm(Xstd, axis=0)
    R = np.clip(Xn.T @ Xn, -1, 1) ** 2
    return (R * mask).sum(1)


def analyze(regime, seed_off=0, tau2=None):
    G, Gref, chrom, pos = build(regime, seed_off)
    tgt = chrom == 0
    n = G.shape[0]
    Xc = standardize(G[:, tgt])
    pos_c = pos[tgt]
    mask = (np.abs(pos_c[:, None] - pos_c[None, :]) <= WINDOW_CM).astype(np.float64)
    K = grm(G[:, ~tgt])
    Vc = SG2 * K + SE2 * np.eye(n)
    one = np.ones((n, 1))

    P_c = residual_precision(Vc, one)
    PXc, q = score_pieces(P_c, Xc)
    A = Xc.T @ PXc
    L_full = (A ** 2).sum(1) / q
    if tau2 is None:
        tau2 = 0.30 / float(np.mean(L_full))
    Ez2 = 1.0 + tau2 * L_full            # exact LMM-GWAS moment (oracle)

    scores = {"ordinary in-sample": corr_ldscore(Xc, mask),
              "ordinary external": corr_ldscore(standardize(Gref[:, tgt]), mask)}
    w, V = np.linalg.eigh(K)
    for k in (2, 10):
        C = np.column_stack((one, V[:, ::-1][:, :k]))
        Xres = (np.eye(n) - C @ np.linalg.pinv(C.T @ C) @ C.T) @ G[:, tgt]
        scores[f"cov-LDSC {k}PC"] = corr_ldscore(Xres, mask)
    scores["PLDSC oracle"] = ld_scores(P_c, Xc, window_mask=mask, PX=PXc, q=q)
    lam, U = np.linalg.eigh(K)
    ynull = np.linalg.cholesky(Vc) @ np.random.default_rng(
        np.random.PCG64(SEED + regime["salt"] + 999 + 1000 * seed_off)).standard_normal(n)
    sg, se = reml_fit(U.T @ ynull, U.T @ one, lam)
    Ph = residual_precision(sg * K + se * np.eye(n), one)
    PXh, qh = score_pieces(Ph, Xc)
    scores["PLDSC fitted"] = ld_scores(Ph, Xc, window_mask=mask, PX=PXh, q=qh)

    res = {name: ld_score_regression(Ez2, L, n_blocks=Xc.shape[1])
           for name, L in scores.items()}
    return res, tau2


METHODS = ["ordinary in-sample", "ordinary external", "cov-LDSC 2PC",
           "cov-LDSC 10PC", "PLDSC oracle", "PLDSC fitted"]


def sweep_point(regime, tau2, n_geno):
    """Mean intercept and mean slope per method over n_geno genotype realizations."""
    ic = {m: [] for m in METHODS}
    sl = {m: [] for m in METHODS}
    for s in range(n_geno):
        res, _ = analyze(regime, seed_off=s, tau2=tau2)
        for m in METHODS:
            ic[m].append(res[m].intercept)
            sl[m].append(res[m].slope)
    return ({m: float(np.mean(ic[m])) for m in METHODS},
            {m: float(np.mean(sl[m])) for m in METHODS})


def run_sweep(title, anchor, regimes, n_geno=6):
    _, tau2 = analyze(anchor, 0)               # fix true model on the anchor
    a_ic, a_sl = sweep_point(anchor, tau2, n_geno)
    print(f"\n{'='*100}\n{title}")
    print(f"  anchor (clean): tau2={tau2:.4g}; mean over {n_geno} genotype realizations")
    print(f"\n  {'regime':16s} " + " ".join(
        f"{m.replace('ordinary','ord').replace('-LDSC','')[:12]:>13s}" for m in METHODS))
    print(f"  {'clean anchor':16s} " +
          " ".join(f"{a_ic[m]:>13.3f}" for m in METHODS) + "   [intercept]")
    for label, reg in regimes:
        r_ic, r_sl = sweep_point(reg, tau2, n_geno)
        ic = " ".join(f"{r_ic[m]:>13.3f}" for m in METHODS)
        dr = " ".join(f"{r_sl[m]/a_sl[m]:>13.2f}" for m in METHODS)
        print(f"  {label:16s} " + ic + "   [intercept]")
        print(f"  {'':16s} " + dr + "   [h2 drift]")


def main():
    print("MAP: LDSC miscalibration on mixed-model GWAS summary statistics")
    print(f"seed={SEED}  N={N} (fixed)  window={WINDOW_CM}cM  background h2={SG2/(SG2+SE2):.1f}")
    print("methods:", ", ".join(METHODS))

    run_sweep(
        "SWEEP A -- cryptic relatedness (homogeneous; vary # sib-pairs at fixed N)",
        {"kind": "related", "n_pairs": 0, "salt": 0},
        [("50 sib-pairs", {"kind": "related", "n_pairs": 50, "salt": 1}),
         ("100 sib-pairs", {"kind": "related", "n_pairs": 100, "salt": 2}),
         ("180 sib-pairs", {"kind": "related", "n_pairs": 180, "salt": 3})],
    )
    run_sweep(
        "SWEEP B -- population structure / admixture (unrelated; vary fst)",
        {"kind": "structured", "fst": 1e-4, "salt": 10},
        [("fst=0.05", {"kind": "structured", "fst": 0.05, "salt": 11}),
         ("fst=0.10", {"kind": "structured", "fst": 0.10, "salt": 12}),
         ("fst=0.15", {"kind": "structured", "fst": 0.15, "salt": 13})],
    )

    print(f"\n{'='*100}")
    print("Read: intercept target 1; h2 drift target 1.00 (<1 = h2 underestimated,")
    print(">1 = overestimated). PLDSC should hold both; conventional scores reveal")
    print("WHEN and HOW MUCH they break -- the map.")


if __name__ == "__main__":
    main()
