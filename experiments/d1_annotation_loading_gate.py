"""D1(a) gate: do baseline-LD annotations correlate with ancestry PC loadings?

The first kill-test of direction D1 (docs/PROJECT_DIRECTIONS.md §4): the S-LDSC
self-term collinearity only matters if real annotations actually carry ancestry
signal. Compute, for the real baseline-LD v2.2 annotation set (Gazal et al.
2017, 1000G phase 3), the correlation of every annotation column with the
top-10 ancestry PC loadings a_j = corr(z_j, u_p) of the multi-ancestry 1000G
FULL panel. The record's decision rule: "if essentially none, D1 collapses to a
footnote."

Discipline:
- LOCO loadings: the GRM eigenspace is built from chr21 only; chr22 SNP
  loadings are evaluated against it, so no SNP's loading sits inside the GRM
  that defined the axes.
- Annotations are the real baselineLD.22.annot.gz (downloaded from
  broad-alkesgroup-public/LDSCORE), matched to the cached chr22 SNP set by
  position. No annotation is recomputed or approximated.
- Per annotation: max |corr| over the 10 loadings, and R^2 of the annotation
  regressed on all 10 jointly (they are near-orthogonal, but the joint number
  is the one S-LDSC's design matrix sees).

Run:  python experiments/d1_annotation_loading_gate.py
"""
from __future__ import annotations

import gzip
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import realdata_1000g_benchmarks as bench  # noqa: E402

ANNOT = os.path.join(HERE, "..", ".scratch_1000g", "baselineLD.22.annot.gz")
K_PC = 10
SEED = 20260811


def load_annotations(pos_wanted: np.ndarray):
    """Return (names, A) with A rows matched (by bp) to pos_wanted order."""
    wanted = {int(p): i for i, p in enumerate(pos_wanted)}
    with gzip.open(ANNOT, "rt") as fh:
        header = fh.readline().split()
        rows: dict[int, list[float]] = {}
        for line in fh:
            parts = line.split()
            bp = int(parts[1])
            if bp in wanted:
                rows[bp] = [float(x) for x in parts[4:]]
    names = header[4:]
    idx = [wanted[bp] for bp in rows]
    A = np.full((len(pos_wanted), len(names)), np.nan)
    for r, (bp, vals) in enumerate(rows.items()):
        A[wanted[bp]] = vals
    matched = ~np.isnan(A[:, 0])
    print(f"  annotations: {matched.sum()} of {len(pos_wanted)} chr22 SNPs "
          f"matched ({len(names)} columns)")
    return names, A, matched


def main():
    t0 = time.time()
    print("D1(a) gate: baseline-LD annotations vs ancestry PC loadings (LOCO)")
    samples, G22, pos22 = bench.load_chrom(22)   # t10 cache
    s21, G21, _ = bench.load_chrom(21)           # t25 cache
    assert s21 == samples
    p21 = G21.mean(0) / 2.0
    keep21 = (np.minimum(p21, 1 - p21) >= 0.05) & (G21.var(0) > 0)
    Z21 = bench.std(G21[:, keep21]).astype(np.float32)
    K = ((Z21 @ Z21.T) / Z21.shape[1]).astype(np.float64)  # LOCO GRM (chr21)
    U = np.linalg.eigh(K)[1][:, -K_PC:]
    n = Z21.shape[0]

    p22 = G22.mean(0) / 2.0
    keep22 = (np.minimum(p22, 1 - p22) >= 0.05) & (G22.var(0) > 0)
    Z22 = bench.std(G22[:, keep22])
    pos22 = pos22[keep22]
    # loadings a_jp = corr(z_j, u_p); std() columns have norm sqrt(N)
    loadings = (Z22.T @ U) / np.sqrt(n)
    print(f"  GRM chr21 ({keep21.sum()} SNPs), loadings chr22 ({keep22.sum()} SNPs) "
          f"({time.time() - t0:.0f}s)")

    names, A, matched = load_annotations(pos22)
    L, A = loadings[matched], A[matched]
    Lc = L - L.mean(0)

    rows = []
    for j, nm in enumerate(names):
        c = A[:, j]
        if np.nanstd(c) == 0 or np.isnan(c).any():
            continue
        cc = c - c.mean()
        denom = np.sqrt((cc @ cc) * (Lc * Lc).sum(0))
        with np.errstate(invalid="ignore", divide="ignore"):
            corrs = np.where(denom > 0, (cc @ Lc) / denom, 0.0)
        X = np.column_stack([np.ones(L.shape[0]), L])
        beta, res, *_ = np.linalg.lstsq(X, c, rcond=None)
        tss = float(((c - c.mean()) ** 2).sum())
        rss = float(res[0]) if res.size else float(((c - X @ beta) ** 2).sum())
        r2 = 1.0 - rss / tss if tss > 0 else 0.0
        rows.append((nm, float(np.max(np.abs(corrs))), float(max(r2, 0.0))))

    rows.sort(key=lambda r: -r[1])
    print(f"\n  {'annotation':<40}{'max|corr|':>10}{'R2 (10 PCs)':>12}")
    for nm, mc, r2 in rows[:15]:
        print(f"  {nm[:40]:<40}{mc:>10.3f}{r2:>12.3f}")
    mc_all = np.array([r[1] for r in rows])
    r2_all = np.array([r[2] for r in rows])
    print(f"\n  {len(rows)} non-constant annotations: "
          f"|corr|>0.05: {(mc_all > 0.05).sum()}, "
          f"|corr|>0.1: {(mc_all > 0.1).sum()}, "
          f"|corr|>0.2: {(mc_all > 0.2).sum()}; "
          f"R2>0.1: {(r2_all > 0.1).sum()}, R2>0.3: {(r2_all > 0.3).sum()}")
    print(f"  median max|corr| {np.median(mc_all):.4f}, max R2 {r2_all.max():.3f}")
    print()
    if (mc_all > 0.1).sum() == 0 and (r2_all > 0.1).sum() == 0:
        print("VERDICT: essentially no annotation-loading correlation -> "
              "D1 collapses to a footnote.")
    else:
        print("VERDICT: real baseline-LD annotations carry ancestry-loading "
              "signal -> D1 survives its gate; proceed to the at-scale test "
              "(weighted S-LDSC, jackknife, M >= 1e5) and the scoped prior-art "
              "sweep.")
    print(f"\nTotal runtime {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
