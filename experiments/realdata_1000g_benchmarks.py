"""Thorough 1000 Genomes benchmarks for the LD-score decomposition.

Follow-up to the chr22 pilot (`realdata_1000g_structure.py`), whose caveats this
suite is designed to close. Three chromosomes (22, 21, 10), a caching parser
(each VCF is parsed once into `.scratch_1000g/cache/`), and six benchmarks:

  B1  DENSITY x PHYSICAL WINDOW (chr22) -- was the pilot's 19-62% windowed
      structure share just thinning density? Fixed physical windows
      (+-100/250/1000 kb) at dense (~6 kb) vs pilot-matched sparse (~13 kb)
      spacing, EUR and FULL.
  B2  LOCO vs IN-CHROMOSOME GRM (chr22 target) -- the pilot used an
      in-chromosome GRM (no other chromosome existed). GRM from chr21+chr10
      gives a proper LOCO operator and quantifies the bias.
  B3  PC-PROJECTION vs P_c (FULL, chr22) -- how many hard-projected PCs
      (cov-LDSC style, K = 2..50) remove the same structure share as the GLS
      operator; plus the exact spectral capture (how much of the genome-wide
      LD score lives in the top-k eigenvectors, via the METHOD sec.14
      identity), i.e. how many PCs cov-LDSC needs on real data.
  B4  METRIC SENSITIVITY (chr22) -- P_c at h2 in {0.1, 0.5, 0.9}; the pilot
      claimed insensitivity without checking.
  B5  MAF GRADIENT (chr22) -- per-SNP structure share in MAF bins
      [0.01-0.05, 0.05-0.2, 0.2-0.5], EUR and FULL.
  B6  ANCESTRY GRADIENT (chr22) -- EUR -> admixed American (CLM/MXL/PUR) ->
      AFR -> FULL at a matched SNP set: the survey's question 1 done properly.

Caveats carried from the pilot: phase-3 calls are essentially complete, so any
non-0/1 genotype is treated as homozygous reference; SNP sets are thinned
(1-in-10 for chr22, 1-in-25 for chr21/10 -- GRM fodder only); P_c uses an
intercept-only design (no covariates) throughout.

Run:  conda run -n pldsc python experiments/realdata_1000g_benchmarks.py
"""
from __future__ import annotations

import gzip
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.join(HERE, "..", ".scratch_1000g")
CACHE = os.path.join(SCRATCH, "cache")
PANEL = os.path.join(SCRATCH, "samples.panel")

MAF_CACHE = 0.01          # MAF floor stored in the cache (full panel)
MAF_COMMON = 0.05         # common-SNP set for cross-cohort comparisons
THIN = {22: 10, 21: 25, 10: 25}
WINDOWS_KB = (100, 250, 1000)


# ---------------------------------------------------------------- data layer

def vcf_path(chrom):
    return os.path.join(SCRATCH, f"ALL.chr{chrom}.vcf.gz")


def cache_paths(chrom):
    t = THIN[chrom]
    return (os.path.join(CACHE, f"chr{chrom}_t{t}_geno.npy"),
            os.path.join(CACHE, f"chr{chrom}_t{t}_pos.npy"),
            os.path.join(CACHE, f"chr{chrom}_t{t}_samples.txt"))


def parse_and_cache(chrom):
    gp, pp, sp = cache_paths(chrom)
    if all(os.path.exists(p) for p in (gp, pp, sp)):
        return
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    thin_every = THIN[chrom]
    kept_lines, pos = [], []
    with gzip.open(vcf_path(chrom), "rt") as f:
        samples = None
        n_seen = 0
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            tab = line.find("\t")
            p1 = line.find("\t", tab + 1)
            p2 = line.find("\t", p1 + 1)
            p3 = line.find("\t", p2 + 1)
            p4 = line.find("\t", p3 + 1)
            ref = line[p2 + 1:p3]
            alt = line[p3 + 1:p4]
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue
            n_seen += 1
            if n_seen % thin_every:
                continue
            kept_lines.append(line)
            pos.append(int(line[tab + 1:p1]))
    n, m = len(samples), len(kept_lines)
    G = np.empty((n, m), dtype=np.float32)
    for j, line in enumerate(kept_lines):
        fields = line.rstrip("\n").split("\t")
        col = np.empty(n, dtype=np.float32)
        for i, gt in enumerate(fields[9:]):
            c = gt[0:3]
            col[i] = 2.0 if c == "1|1" else (1.0 if c in ("0|1", "1|0") else 0.0)
        G[:, j] = col
    p = G.mean(0) / 2.0
    maf = np.minimum(p, 1 - p)
    keep = maf >= MAF_CACHE
    G, pos = G[:, keep], np.array(pos, dtype=np.int64)[keep]
    np.save(gp, G)
    np.save(pp, pos)
    with open(sp, "w") as fh:
        fh.write("\n".join(samples))
    print(f"  [cache] chr{chrom}: {n_seen} SNP lines, kept {G.shape[1]} "
          f"(1-in-{thin_every}, MAF>={MAF_CACHE}), {time.time()-t0:.0f}s")


def load_chrom(chrom):
    parse_and_cache(chrom)
    gp, pp, sp = cache_paths(chrom)
    with open(sp) as fh:
        samples = fh.read().split("\n")
    return samples, np.load(gp), np.load(pp)


def load_panel():
    pop, superpop = {}, {}
    with open(PANEL) as f:
        cols = f.readline().strip().split("\t")
        i_s, i_p, i_sp = cols.index("sample"), cols.index("pop"), cols.index("super_pop")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            pop[parts[i_s]] = parts[i_p]
            superpop[parts[i_s]] = parts[i_sp]
    return pop, superpop


# ---------------------------------------------------------------- operators

def std(G):
    p = G.mean(0) / 2.0
    s = np.sqrt(2.0 * p * (1.0 - p))
    s = np.where(s == 0, 1.0, s)
    return (G.astype(np.float64) - 2.0 * p) / s


def eigh_grm(X):
    K = (X @ X.T) / X.shape[1]
    return np.linalg.eigh(K)


def pcov(w, U, h2):
    """P_c (intercept-only) from the GRM eigendecomposition; V = h2 K + (1-h2) I."""
    d = 1.0 / (h2 * w + (1.0 - h2))
    W1 = (U * d) @ (U.T @ np.ones(U.shape[0]))
    W = (U * d) @ U.T
    return W - np.outer(W1, W1) / W1.sum()


def project(X, U, k):
    Uk = U[:, ::-1][:, :k]
    return X - Uk @ (Uk.T @ X)


def ld_scores(X, pos, windows_kb=(), block=1000):
    """Per-SNP sums of squared CORRELATIONS: genome-wide ('full') and windowed.

    Columns are normalized to unit norm first, so both ordinary standardized
    genotypes and whitened/projected genotypes yield correlations (the
    cov-LDSC convention), not covariances. Variants with numerically zero
    residual variance are excluded (reported in the returned dict).
    """
    N, M = X.shape
    norms = np.linalg.norm(X, axis=0)
    tol = np.finfo(np.float64).eps * max(X.shape) * max(norms.max(), 1.0)
    dead = norms <= tol
    Xn = X / np.where(dead, 1.0, norms)
    out = {"full": np.zeros(M), "n_excluded": int(dead.sum())}
    for w in windows_kb:
        out[w] = np.zeros(M)
    for lo in range(0, M, block):
        hi = min(M, lo + block)
        R2 = (Xn[:, lo:hi].T @ Xn) ** 2
        out["full"][lo:hi] = R2.sum(1)
        for jj in range(lo, hi):
            for w in windows_kb:
                w0 = np.searchsorted(pos, pos[jj] - w * 1000, "left")
                w1 = np.searchsorted(pos, pos[jj] + w * 1000, "right")
                out[w][jj] += R2[jj - lo, w0:w1].sum()
    return out


def _cmask(Gc, min_maf=0.01):
    """Boolean mask of cohort-polymorphic variants (guards zero-variance cols)."""
    p = Gc.mean(0) / 2.0
    maf = np.minimum(p, 1 - p)
    return maf >= min_maf


def pc_score(X, w, U, k=10):
    Uk = U[:, ::-1][:, :k]
    lam = w[::-1][:k]
    c2 = ((X.T @ Uk) ** 2) / X.shape[0]
    return (lam[None, :] / lam.sum() * c2).sum(1)


def shares(scores_o, scores_w, key):
    so, sw = scores_o[key], scores_w[key]
    fo, fw = scores_o["full"], scores_w["full"]
    return (fo.mean(), fw.mean(), (fo - fw).mean() / fo.mean() * 100,
            (so.mean() - sw.mean()) / so.mean() * 100)


def cohort_masks(samples, pop, superpop):
    pops = np.array([pop[s] for s in samples])
    sps = np.array([superpop[s] for s in samples])
    return {
        "EUR": sps == "EUR",
        "AMR": np.isin(pops, ["CLM", "MXL", "PUR"]),
        "AFR": sps == "AFR",
        "FULL": np.ones(len(samples), bool),
    }


def sparse_subset(pos, min_gap_bp=13000):
    keep, last = [], -10**18
    for j, p in enumerate(pos):
        if p - last >= min_gap_bp:
            keep.append(j)
            last = p
    return np.array(keep)


# ---------------------------------------------------------------- benchmarks

def bench_B1(G, pos, masks):
    print("\n== B1  DENSITY x PHYSICAL WINDOW (chr22; SNP set MAF>=0.05 FULL) ==")
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= MAF_COMMON)
    sparse = keep[sparse_subset(pos[keep])]
    print(f"  dense: {len(keep)} SNPs (~{int(np.median(np.diff(pos[keep]))/1000)} kb); "
          f"sparse: {len(sparse)} SNPs (~{int(np.median(np.diff(pos[sparse]))/1000)} kb)")
    print(f"  {'cohort':6s} {'density':7s} {'full_ord':>9s} {'full_wht':>9s} "
          f"{'struct%full':>11s} {'struct%250kb':>12s}")
    for cname in ("EUR", "FULL"):
        for dname, idx in (("dense", keep), ("sparse", sparse)):
            Gsub = G[masks[cname]][:, idx]
            cm = _cmask(Gsub)
            Gc, pc = Gsub[:, cm], pos[idx][cm]
            X = std(Gc)
            w, U = eigh_grm(X)
            P = pcov(w, U, 0.5)
            so = ld_scores(X, pc, WINDOWS_KB)
            sw = ld_scores(P @ X, pc, WINDOWS_KB)
            fo, fw, sfull, s250 = shares(so, sw, 250)
            print(f"  {cname:6s} {dname:7s} {fo:9.2f} {fw:9.2f} {sfull:10.0f}% {s250:11.0f}%")


def bench_B2(G, pos, masks, loco_grm):
    print("\n== B2  LOCO (chr21+10 GRM) vs IN-CHROMOSOME GRM (chr22 target) ==")
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= MAF_COMMON)
    print(f"  {'cohort':6s} {'GRM':12s} {'full_wht':>9s} {'struct%full':>11s} "
          f"{'struct%250kb':>12s}")
    for cname in ("EUR", "FULL"):
        Gsub = G[masks[cname]][:, keep]
        cm = _cmask(Gsub)
        Gc, pc = Gsub[:, cm], pos[keep][cm]
        X = std(Gc)
        so = ld_scores(X, pc, (250,))
        for gname, (w, U) in (("in-chrom", eigh_grm(X)), ("LOCO", loco_grm[cname])):
            P = pcov(w, U, 0.5)
            sw = ld_scores(P @ X, pc, (250,))
            fo, fw, sfull, s250 = shares(so, sw, 250)
            print(f"  {cname:6s} {gname:12s} {fw:9.2f} {sfull:10.0f}% {s250:11.0f}%")


def bench_B3(G, pos, masks):
    print("\n== B3  PC-PROJECTION vs P_c (FULL, chr22) ==")
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= MAF_COMMON)
    Gc, pc = G[masks["FULL"]][:, keep], pos[keep]
    X = std(Gc)
    w, U = eigh_grm(X)
    so = ld_scores(X, pc, (250,))
    P = pcov(w, U, 0.5)
    sw = ld_scores(P @ X, pc, (250,))
    fo, fw = so["full"].mean(), sw["full"].mean()
    wo, ww = so[250].mean(), sw[250].mean()
    # exact spectral capture of the genome-wide score by the top-k eigenvectors
    c2all = ((X.T @ U) ** 2) / X.shape[0]       # (M, N)
    per_pc = w[None, :] * c2all                  # lambda_p * sum_j corr^2
    tot = per_pc.sum()
    capture = np.cumsum(per_pc.sum(0)[::-1]) / tot   # top-k cumulative fraction
    ks = (1, 2, 5, 10, 20, 50, 100)
    print("  operator            full_mean  struct%full  struct%250kb  (vs P_c anchor)")
    for k in (2, 5, 10, 20, 50):
        Xp = project(X, U, k)
        sp = ld_scores(Xp, pc, (250,))
        rem_full = (sp["full"].mean() - fw) / (fo - fw) * 100
        rem_win = (sp[250].mean() - ww) / (wo - ww) * 100
        print(f"  project top-{k:<3d} PCs {sp['full'].mean():9.2f} {rem_full:10.0f}%"
              f" {rem_win:11.0f}%  structure REMAINING")
    print(f"  P_c (GLS, h2=0.5)  {fw:9.2f}   anchor (0% remaining)")
    cap = {k: float(capture[k - 1]) for k in ks}
    k50 = next(k for k in ks if cap[k] >= 0.5)
    k90 = next((k for k in ks if cap[k] >= 0.9), None)
    print("  spectral capture of the genome-wide score (exact, METHOD sec.14):")
    print("   " + "  ".join(f"top-{k}: {cap[k]*100:.0f}%" for k in ks))
    print(f"   -> 50% of the genome-wide score lives in top-{k50}, "
          f"90% in top-{k90} eigenvectors")


def bench_B4(G, pos, masks):
    print("\n== B4  METRIC SENSITIVITY: P_c at h2 in {0.1, 0.5, 0.9} (chr22) ==")
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= MAF_COMMON)
    for cname in ("EUR", "FULL"):
        Gsub = G[masks[cname]][:, keep]
        cm = _cmask(Gsub)
        Gc, pc = Gsub[:, cm], pos[keep][cm]
        X = std(Gc)
        w, U = eigh_grm(X)
        so = ld_scores(X, pc, (250,))
        line = f"  {cname:6s}"
        for h2 in (0.1, 0.5, 0.9):
            sw = ld_scores(pcov(w, U, h2) @ X, pc, (250,))
            _, _, sfull, s250 = shares(so, sw, 250)
            line += f"   h2={h2}: full {sfull:.0f}% / win {s250:.0f}%"
        print(line)


def bench_B5(G, pos, masks):
    print("\n== B5  MAF GRADIENT of the structure share (chr22) ==")
    bins = ((0.01, 0.05), (0.05, 0.2), (0.2, 0.5))
    for cname in ("EUR", "FULL"):
        Gc = G[masks[cname]]
        p = Gc.mean(0) / 2.0
        maf = np.minimum(p, 1 - p)
        keep = maf >= MAF_CACHE
        Gk, mafk = Gc[:, keep], maf[keep]
        X = std(Gk)
        w, U = eigh_grm(X)
        so = ld_scores(X, pos[keep], ())
        sw = ld_scores(pcov(w, U, 0.5) @ X, pos[keep], ())
        share_j = (so["full"] - sw["full"]) / so["full"]
        line = f"  {cname:6s}"
        for lo, hi in bins:
            b = (mafk >= lo) & (mafk < hi)
            line += f"   MAF[{lo},{hi}): {np.mean(share_j[b])*100:.0f}% (n={b.sum()})"
        print(line)


def bench_B6(G, pos, masks):
    print("\n== B6  ANCESTRY GRADIENT at a matched SNP set (chr22, MAF>=0.05 FULL) ==")
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= MAF_COMMON)
    print(f"  {'cohort':6s} {'N':>5s} {'lambda1':>8s} {'struct%full':>11s} "
          f"{'struct%250kb':>12s} {'corr(struct,PC10)':>17s}")
    for cname in ("EUR", "AMR", "AFR", "FULL"):
        Gsub = G[masks[cname]][:, keep]
        cm = _cmask(Gsub)
        Gc, pc = Gsub[:, cm], pos[keep][cm]
        X = std(Gc)
        w, U = eigh_grm(X)
        so = ld_scores(X, pc, (250,))
        sw = ld_scores(pcov(w, U, 0.5) @ X, pc, (250,))
        _, _, sfull, s250 = shares(so, sw, 250)
        struct_j = so["full"] - sw["full"]
        c = np.corrcoef(struct_j, pc_score(X, w, U))[0, 1]
        print(f"  {cname:6s} {X.shape[0]:5d} {w[-1]:8.1f} {sfull:10.0f}% "
              f"{s250:11.0f}% {c:17.3f}")


def main():
    t0 = time.time()
    print("1000G decomposition benchmarks (follow-up to the chr22 pilot)")
    pop, superpop = load_panel()
    samples, G22, pos22 = load_chrom(22)
    assert np.all(np.diff(pos22) > 0)
    masks = cohort_masks(samples, pop, superpop)
    print(f"  chr22 loaded: {G22.shape[1]} SNPs x {G22.shape[0]} samples "
          f"({time.time()-t0:.0f}s)")

    # LOCO GRMs from chr21+chr10 (sample order reconciled by name)
    s21, G21, _ = load_chrom(21)
    s10, G10, _ = load_chrom(10)
    assert s21 == s10 == samples, "sample order mismatch across chromosomes"
    loco_grm = {}
    for cname in ("EUR", "FULL"):
        m = masks[cname]
        Xl = std(np.hstack([G21[m], G10[m]]))
        loco_grm[cname] = eigh_grm(Xl)
        del Xl
    print(f"  LOCO GRMs built ({time.time()-t0:.0f}s)")

    bench_B1(G22, pos22, masks)
    bench_B2(G22, pos22, masks, loco_grm)
    bench_B3(G22, pos22, masks)
    bench_B4(G22, pos22, masks)
    bench_B5(G22, pos22, masks)
    bench_B6(G22, pos22, masks)
    print(f"\nTotal runtime {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
