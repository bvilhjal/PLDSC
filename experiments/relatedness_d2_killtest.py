"""D2 kill-test: is the panel-relatedness noise-floor correction a paper or a paragraph?

The one decision-gating measurement from docs/PROJECT_DIRECTIONS.md §4 (D2):
compute the effective pair degrees of freedom ``tr(Atilde)^2 / tr(Atilde^2)``
(METHOD.md §10.1, Eq. M22) against the naive residual df on real LD panels.
If the ratio is ~1 on the standard panels, the under-subtraction is negligible
and D2 collapses to a one-paragraph caveat; if a standard panel shows a
materially inflated unlinked-pair r^2 noise floor, the "LDSC intercept is not a
clean confounding barometer on related in-sample panels" claim is alive.

Panels: 1000G EUR, AFR, AMR, FULL (real, this script), plus simulated full-sib
panels as positive controls (calibrated against
experiments/relatedness_effective_df.py, which reported N_eff 399 -> 84 for a
16-sib design). UKB-unrelated is omitted: it is unrelated *by construction*
(KING kinship < 0.0442), so its ratio is ~1 by definition of the panel, not by
measurement.

Two estimation notes, learned from a first pass of this script:

- ``Atilde = P_C A P_C`` is evaluated with C = intercept + top-10 GRM PCs (the
  covariates a GWAS actually adjusts for). Intercept-only numbers are printed
  for contrast, but they conflate ancestry sharing with relatedness: in a
  multi-ancestry panel same-continent pairs carry F_st-scale GRM off-diagonals,
  which is structure (cov-LDSC territory), not the cryptic relatedness D2 is
  about.
- ``tr(Atilde^2)`` is estimated by the split-half cross-product
  ``tr(Atilde_a Atilde_b)`` from disjoint SNP halves. A plain ``tr(K^2)`` of a
  finite-SNP GRM is inflated by ~N^2/M of estimation noise, which biases N_eff
  downward — materially so at biobank N.

Reads the 1000G phase-3 VCFs through the benchmark script's cached parser
(experiments/realdata_1000g_benchmarks.py), thinned 1-in-25 with MAF >= 0.05.

Run:  python experiments/relatedness_d2_killtest.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import realdata_1000g_benchmarks as bench  # noqa: E402
from pldsc.simulate import grm, simulate_relatedness  # noqa: E402

THIN = 25
MAF = 0.05
K_PC = 10
SEED = 20260811
S_WINDOW = 2000  # window-accumulated SNP count used in the METHOD §10.1 illustration


def _project(K: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """P_C K P_C with P_C = I - Q Q^T, Q an orthonormal basis of the design."""
    KQ = K @ Q
    return K - Q @ (Q.T @ K) - KQ @ Q.T + Q @ (Q.T @ KQ) @ Q.T


def panel_report(name: str, G: np.ndarray) -> dict:
    """Split-half N_eff estimates (intercept-only and PC-adjusted) per cohort."""
    p = G.mean(0) / 2.0
    keep = (np.minimum(p, 1 - p) >= MAF) & (G.var(0) > 0)
    Z = bench.std(G[:, keep])
    n, m = Z.shape
    # Disjoint SNP halves (interleaved so LD is shared, noise independent).
    Za = np.ascontiguousarray(Z[:, 0::2], dtype=np.float32)
    Zb = np.ascontiguousarray(Z[:, 1::2], dtype=np.float32)
    Ka = (Za @ Za.T) / Za.shape[1]
    Kb = (Zb @ Zb.T) / Zb.shape[1]
    K = 0.5 * (Ka + Kb)
    pcs = np.linalg.eigh(K.astype(np.float64))[1][:, -K_PC:]

    out = {}
    for label, C in (("intercept", None), (f"{K_PC} PCs", pcs)):
        Q = np.linalg.qr(C if C is not None else np.ones((n, 1)))[0]
        Aa, Ab = _project(Ka, Q), _project(Kb, Q)
        tr_a = 0.5 * (np.trace(Aa) + np.trace(Ab))
        tr_a2 = float(np.sum(Aa * Ab))  # unbiased for tr(Atilde^2)
        n_eff = float(tr_a**2 / tr_a2)
        naive = n - Q.shape[1] - 1  # the package default: residual_df - 1
        ratio = n_eff / naive
        # Excess unlinked-pair noise floor accumulated over a window of
        # S_WINDOW SNPs, in LD-score units (typical LDSC intercepts ~1.0-1.1).
        excess = S_WINDOW * (1.0 / n_eff - 1.0 / naive)
        Am = 0.5 * (Aa + Ab)
        off = Am[np.triu_indices(n, k=1)]
        census = ", ".join(f">{t}: {int((off > t).sum())}" for t in (0.05, 0.1))
        print(f"  {name:12s} {label:10s} N={n:5d} M={m:6d}  N_eff={n_eff:9.1f}  "
              f"naive={naive:6d}  ratio={ratio:.4f}  window excess={excess:8.3f}  "
              f"pairs {census}")
        out[label] = (ratio, excess)
    return out


def main():
    t0 = time.time()
    print("D2 kill-test: effective pair df vs naive residual df on real panels")
    bench.THIN.update({22: THIN, 21: THIN})
    samples, G22, _ = bench.load_chrom(22)
    s21, G21, _ = bench.load_chrom(21)
    assert s21 == samples, "sample order mismatch across chromosomes"
    G = np.hstack([G22, G21])
    print(f"  chr21+chr22 parsed: {G.shape[1]} SNP lines x {G.shape[0]} samples "
          f"({time.time() - t0:.0f}s)")

    pop, superpop = bench.load_panel()
    masks = bench.cohort_masks(samples, pop, superpop)
    verdicts = {}
    for cname in ("EUR", "AFR", "AMR", "FULL"):
        verdicts[cname] = panel_report(cname, G[masks[cname]])

    print("  positive control (simulated full-sib panels):")
    rng = np.random.default_rng(SEED)
    for n_fam, sibs in ((100, 4), (25, 16)):
        Gs, _, _ = simulate_relatedness(rng, [3000, 3000], n_fam, sibs)
        # Pedigree-structured but homogeneous: no PCs needed, intercept only.
        Gs = Gs.astype(np.float64)
        K = grm(Gs)
        n = K.shape[0]
        Q = np.linalg.qr(np.ones((n, 1)))[0]
        A = _project(K, Q)
        n_eff = float(np.trace(A) ** 2 / np.sum(A * A))
        ratio = n_eff / (n - 2)
        excess = S_WINDOW * (1.0 / n_eff - 1.0 / (n - 2))
        print(f"  sibs {n_fam}x{sibs:<3d} intercept  N={n:5d}        "
              f"N_eff={n_eff:9.1f}  naive={n - 2:6d}  ratio={ratio:.4f}  "
              f"window excess={excess:8.3f}")

    real = {c: v[f"{K_PC} PCs"][0] for c, v in verdicts.items()}
    print()
    if all(r > 0.99 for r in real.values()):
        print("VERDICT: PC-adjusted ratio ~1 on every real panel -> the "
              "under-subtraction is negligible; D2 is a one-paragraph caveat, "
              "not a paper.")
    else:
        worst = min(real, key=real.get)
        excess = verdicts[worst][f"{K_PC} PCs"][1]
        print(f"VERDICT: {worst} departs from 1 (PC-adjusted ratio "
              f"{real[worst]:.4f}, window excess {excess:.3f} LD-score units); "
              "judge against typical intercepts (~1.0-1.1) and their SEs "
              "before declaring D2 alive.")
    print(f"\nTotal runtime {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
