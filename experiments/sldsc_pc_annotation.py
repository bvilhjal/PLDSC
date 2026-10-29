"""Does a PC-loading annotation de-bias functional enrichment in S-LDSC?

Idea: population structure inflates chi^2 by ~ b^2 * gamma_{j,1}^2 (squared PC
loading; Bhatia et al.). If a functional annotation correlates with the ancestry
axis, standard S-LDSC misreads that inflation as functional enrichment. Adding a
squared-PC-loading annotation should soak the structure into its own coefficient
(tau_pc) and de-bias the functional coefficient (tau_func).

Design: two-population cohort. A functional annotation is enriched for
heritability AND correlated with ancestry (flagged preferentially at
high-Fst SNPs). Phenotype = polygenic (with functional enrichment) + a tunable
structure component. Ordinary (OLS) GWAS -> chi^2 (structure confounding is
present, as it is for real sumstats S-LDSC consumes). PC loadings come from a
LOCO GRM (other chromosomes) to avoid proximal contamination. Compare S-LDSC
   A: intercept + base + func
   B: intercept + base + func + pc-loading
by how well each recovers the TRUE functional coefficient, with structure off
and on, averaged over replicates.

Run:  python experiments/sldsc_pc_annotation.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc.simulate import ar1_haps, bn_freqs, grm, standardize  # noqa: E402

SEED = 20260718
CHROMS = [400, 500]          # target chr 0 (annotated); chr 1 builds the LOCO GRM/PCs
N = 500
FST = 0.12
RHO = 0.6
WINDOW_CM = 2.0
FRAC_FUNC = 0.20
TAU_BASE = 1.0               # per-SNP effect variance (arbitrary units)
TAU_FUNC = 2.0              # functional SNPs carry TAU_BASE+TAU_FUNC -> 3x -> enriched


def simulate(rng):
    # two equal populations
    p1, p2 = zip(*[bn_freqs(rng, m, FST) for m in CHROMS])
    half = N // 2
    def pop_geno(pfreq):
        h = np.hstack([ar1_haps(rng, np.tile(pfreq[c], (2 * half, 1)), RHO)
                       for c in range(len(CHROMS))])
        return h[0::2] + h[1::2]
    G = np.vstack([pop_geno(p1), pop_geno(p2)])
    anc = np.concatenate([np.full(half, -1.0), np.full(half, 1.0)])  # ancestry axis
    chrom = np.concatenate([np.full(m, c) for c, m in enumerate(CHROMS)])
    keep = G.var(0) > 0
    G, chrom = G[:, keep], chrom[keep]
    # allele-frequency differentiation per SNP (for correlating func with ancestry)
    fdiff = np.abs(G[anc < 0].mean(0) - G[anc > 0].mean(0)) / 2.0
    return G, chrom, anc, fdiff


def windowed_strat_ldscore(Xn, pos, annot, mask):
    """Sum_k annot[k] r_jk^2 within the window (Xn = unit-norm standardized)."""
    R2 = np.clip(Xn.T @ Xn, -1, 1) ** 2 * mask
    return R2 @ annot


def sldsc_fit(chi2, ld_cols, ell_base):
    """S-LDSC-style WLS of chi2 on [1, *ld_cols].

    Two-step weights ~ 1 / (ell_base * max(pred,1)^2): the 1/ell_base factor is
    the LD-overlap weight, the 1/pred^2 factor the chi^2 heteroscedasticity
    weight. Stabilizes the fit that plain OLS renders degenerate under structure.
    """
    A = np.column_stack([np.ones_like(chi2)] + ld_cols)
    coef, *_ = np.linalg.lstsq(A, chi2, rcond=None)
    pred = np.maximum(A @ coef, 1.0)
    w = 1.0 / (np.maximum(ell_base, 1.0) * pred ** 2)
    Aw = A * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(Aw, chi2 * np.sqrt(w), rcond=None)
    return coef


def run(rng, structure, func_correlated=True):
    G, chrom, anc, fdiff = simulate(rng)
    tgt = chrom == 0
    Xc = standardize(G[:, tgt])
    n, M = Xc.shape
    pos_c = np.arange(M) * 0.15
    mask = (np.abs(pos_c[:, None] - pos_c[None, :]) <= WINDOW_CM).astype(float)

    # functional annotation: enriched. Correlated with ancestry (prob rises with
    # Fst) OR ancestry-independent (constant prob) -- the control.
    if func_correlated:
        fd = fdiff[tgt]
        pfunc = np.clip(FRAC_FUNC * (fd / fd.mean()), 0, 1)
    else:
        pfunc = np.full(M, FRAC_FUNC)
    func = (rng.uniform(size=M) < pfunc).astype(float)

    # true per-SNP effect variance and phenotype
    var_k = TAU_BASE + TAU_FUNC * func
    beta = rng.standard_normal(M) * np.sqrt(var_k)
    y = Xc @ beta
    y = y / y.std()
    if structure > 0:
        y = y + structure * anc          # ancestry-driven confounding component
    y = y + rng.standard_normal(n)       # environmental noise

    # ordinary (OLS) GWAS: marginal z per SNP -> chi^2
    Xn = Xc / np.linalg.norm(Xc, axis=0)
    yc = (y - y.mean())
    z = (Xn.T @ yc) / yc.std()           # ~ standardized marginal score
    chi2 = z ** 2

    # LOCO PC1 loadings for the target chromosome, from the OTHER chromosome's GRM
    K = grm(G[:, ~tgt])
    w, V = np.linalg.eigh(K)
    pc1 = V[:, -1]
    loading2 = ((Xn.T @ (pc1 - pc1.mean())) / np.linalg.norm(pc1 - pc1.mean())) ** 2

    base = np.ones(M)
    ell_base = windowed_strat_ldscore(Xn, pos_c, base, mask)
    ell_func = windowed_strat_ldscore(Xn, pos_c, func, mask)
    ell_pc = windowed_strat_ldscore(Xn, pos_c, loading2, mask)

    coef_A = sldsc_fit(chi2, [ell_base, ell_func], ell_base)          # intercept, base, func
    coef_B = sldsc_fit(chi2, [ell_base, ell_func, ell_pc], ell_base)  # + pc-loading
    return {"func_A": coef_A[2], "func_B": coef_B[2], "pc_B": coef_B[3],
            "int_A": coef_A[0], "int_B": coef_B[0]}


def scenario(name, func_correlated, reps=40):
    print(f"\n{'='*80}\n{name}")
    # paired over seeds: each seed contributes run(b=0) and run(b) so ratios are
    # within-realization (removes cohort noise). Ratio = func coef(b)/func coef(b=0).
    print(f"  {'b':>5s}  {'A func /A(b0)':>14s}  {'B func /B(b0)':>14s}  "
          f"{'pc coef':>10s}  {'A intercept':>12s}")
    refs = [run(np.random.default_rng(np.random.PCG64(SEED + 200 + r)), 0.0, func_correlated)
            for r in range(reps)]
    fA0 = np.mean([x["func_A"] for x in refs])
    fB0 = np.mean([x["func_B"] for x in refs])
    for b in [0.0, 0.5, 1.0]:
        rows = [run(np.random.default_rng(np.random.PCG64(SEED + 200 + r)), b, func_correlated)
                for r in range(reps)]
        fA = np.mean([x["func_A"] for x in rows])
        fB = np.mean([x["func_B"] for x in rows])
        pcB = np.mean([x["pc_B"] for x in rows])
        iA = np.mean([x["int_A"] for x in rows])
        print(f"  {b:>5.1f}  {fA/fA0:>14.2f}  {fB/fB0:>14.2f}  {pcB:>10.4g}  {iA:>12.3f}")
    print(f"  (model-A func coef at b=0 = {fA0:.3g}; model-B func coef at b=0 = {fB0:.3g}"
          f"  -> B/A at b=0 = {fB0/fA0:.2f} shows the collinearity cost)")


def main():
    print("S-LDSC with a PC-loading annotation: does it de-bias functional enrichment?")
    print(f"seed={SEED} N={N} fst={FST} frac_func={FRAC_FUNC} "
          f"true tau_func/tau_base={TAU_FUNC}/{TAU_BASE}  (WLS, paired seeds)")
    scenario("SCENARIO 1: functional annotation CORRELATED with ancestry", True)
    scenario("SCENARIO 2: functional annotation INDEPENDENT of ancestry (control)", False)
    print("\n" + "=" * 80)
    print("Read: 'A func /A(b0)' >1 means structure inflates standard-S-LDSC")
    print("functional enrichment; 'B func /B(b0)' ~1 means the pc-loading annotation")
    print("stabilizes it. The B/A-at-b0 line is the price: if <1, the pc annotation")
    print("also steals genuine enrichment (collinearity) -- expected only in Scenario 1.")


if __name__ == "__main__":
    main()
