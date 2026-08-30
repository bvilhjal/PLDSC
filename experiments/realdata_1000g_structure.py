"""1000 Genomes pilot: the real structure fraction of the LD score.

REAL_DATA_SURVEY.md question 1, smallest version: how much of a real cohort's
LD score is structure (non-local) vs linkage (local), across a
homogeneous -> multi-ancestry gradient? The simulations gave ~50% genome-wide /
~5% within a local window at fst=0.1 (ld_decomposition.py); this measures the
real value on chr22 of 1000 Genomes phase 3 (20130502), for EUR-only (503),
AFR-only (661), and the full panel (2504).

Mirrors ld_decomposition.py: ordinary LD score ell_ord (sum of r^2), ell_local
with structure-whitened correlations (P_c-metric with V_c = sigma_g2 K +
sigma_e2 I from the LOCO-free GRM; chr22 is the only chromosome here, so the
GRM is in-chromosome -- stated as a limit), ell_struct = ell_ord - ell_local
(the split is additive BY CONSTRUCTION; the empirical content is the
correlation with the PC score), and the top-10 PC score. Windowed arm: +-50
SNP window. Genome-wide arm: blocked full-chromosome sums.

Inputs (not committed): .scratch_1000g/ALL.chr22.vcf.gz + samples.panel.
Parsing is dependency-free: biallelic SNPs only, MAF >= 0.05 in the full
panel, deterministic thinning to ~THIN_EVERY-th variant.

Run:  conda run -n pldsc python experiments/realdata_1000g_structure.py
"""
from __future__ import annotations

import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pldsc._lmm import residual_precision  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "..", ".scratch_1000g", "ALL.chr22.vcf.gz")
PANEL = os.path.join(HERE, "..", ".scratch_1000g", "samples.panel")

MAF_MIN = 0.05
THIN_EVERY = 25           # keep every 25th common biallelic SNP (~20k SNPs)
WINDOW = 50               # +-50 SNP "local" window (no genetic map: SNP-index)
K_PC = 10
SG2, SE2 = 0.5, 0.5       # P_c metric shape; results are insensitive (checked)


def load_panel(path):
    pop, superpop = {}, {}
    with open(path) as f:
        header = f.readline()
        cols = header.strip().split("\t")
        i_s, i_p, i_sp = cols.index("sample"), cols.index("pop"), cols.index("super_pop")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            pop[parts[i_s]] = parts[i_p]
            superpop[parts[i_s]] = parts[i_sp]
    return pop, superpop


def parse_vcf(path, maf_min, thin_every):
    """Stream the VCF; return (samples, G float32 (n,m), positions)."""
    kept_lines, pos = [], []
    with gzip.open(path, "rt") as f:
        samples = None
        n_seen = 0
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            # cheap pre-filter before any real parsing
            tab = line.find("\t")
            p1 = line.find("\t", tab + 1)
            p2 = line.find("\t", p1 + 1)
            p3 = line.find("\t", p2 + 1)
            p4 = line.find("\t", p3 + 1)
            ref = line[p2 + 1:p3]
            alt = line[p3 + 1:p4]
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue          # SNPs only, biallelic only
            # AF from the INFO field would need parsing; filter by AC later.
            n_seen += 1
            if n_seen % thin_every:
                continue
            kept_lines.append(line)
            pos.append(int(line[tab + 1:p1]))
    n = len(samples)
    m = len(kept_lines)
    print(f"  parsed chr22: {n_seen} biallelic SNP lines seen, kept {m} after "
          f"1-in-{thin_every} thinning; {n} samples")
    G = np.empty((n, m), dtype=np.float32)
    for j, line in enumerate(kept_lines):
        fields = line.rstrip("\n").split("\t")
        col = np.empty(n, dtype=np.float32)
        for i, gt in enumerate(fields[9:]):
            c = gt[0:3]
            col[i] = 2.0 if c == "1|1" else (1.0 if c in ("0|1", "1|0") else 0.0)
        G[:, j] = col
    # MAF filter in the full panel
    p = G.mean(0) / 2.0
    maf = np.minimum(p, 1 - p)
    keep = maf >= maf_min
    return samples, G[:, keep], np.array(pos)[keep]


def ldscores_blocked(X, window):
    """ell_j = sum_k r_jk^2, windowed (SNP-index) and full, in blocks of rows."""
    N, M = X.shape
    ell_win = np.zeros(M)
    ell_full = np.zeros(M)
    B = 2000
    for lo in range(0, M, B):
        hi = min(M, lo + B)
        R = (X[:, lo:hi].T @ X) / N            # (rows, M)
        R2 = R * R
        ell_full[lo:hi] = R2.sum(1)
        for jj in range(lo, hi):
            w0, w1 = max(0, jj - window), min(M, jj + window + 1)
            ell_win[jj] = R2[jj - lo, w0:w1].sum()
    return ell_win, ell_full


def decomposition(G, label):
    N, M = G.shape
    print(f"\n[{label}] N={N} M={M}")
    # standardize (expected-variance convention, as pldsc.simulate.standardize)
    p = G.mean(0) / 2.0
    s = np.sqrt(2.0 * p * (1.0 - p))
    s = np.where(s == 0, 1.0, s)
    X = ((G - 2.0 * p) / s).astype(np.float64)

    K = (X @ X.T) / M
    w, U = np.linalg.eigh(K)
    Vc = SG2 * K + SE2 * np.eye(N)
    P = residual_precision(Vc, np.ones((N, 1)))
    PX = P @ X
    q = np.einsum("ij,ij->j", X, PX)
    A = (X.T @ PX) / np.sqrt(np.outer(q, q))     # P_c-metric correlations

    # whitened LD scores (blocked)
    ell_win_w = np.zeros(M)
    ell_full_w = np.zeros(M)
    B = 2000
    for lo in range(0, M, B):
        hi = min(M, lo + B)
        Rb = A[lo:hi, :]
        R2 = Rb * Rb
        ell_full_w[lo:hi] = R2.sum(1)
        for jj in range(lo, hi):
            w0, w1 = max(0, jj - WINDOW), min(M, jj + WINDOW + 1)
            ell_win_w[jj] = R2[jj - lo, w0:w1].sum()

    ell_win_o, ell_full_o = ldscores_blocked(X, WINDOW)

    struct_win = ell_win_o - ell_win_w
    struct_full = ell_full_o - ell_full_w

    # top-K PC score (LOCO-free: chr22 is the whole dataset here)
    lam = w[::-1][:K_PC]
    Utop = U[:, ::-1][:, :K_PC]
    corr2 = ((X.T @ Utop) ** 2) / N             # (M, K_PC)
    pc_score = (lam[None, :] / lam.sum() * corr2).sum(1)

    c_win = np.corrcoef(struct_win, pc_score)[0, 1]
    c_full = np.corrcoef(struct_full, pc_score)[0, 1]
    print(f"  mean scores: ordinary full {ell_full_o.mean():.3f}  "
          f"whitened full {ell_full_w.mean():.3f}  struct full {struct_full.mean():.3f} "
          f"({struct_full.mean()/ell_full_o.mean()*100:.0f}% of full)")
    print(f"               ordinary win  {ell_win_o.mean():.3f}  "
          f"whitened win  {ell_win_w.mean():.3f}  struct win  {struct_win.mean():.3f} "
          f"({struct_win.mean()/ell_win_o.mean()*100:.0f}% of windowed)")
    print(f"  corr(ell_struct, top-{K_PC} PC score): full {c_full:.3f}  windowed {c_win:.3f}")
    print(f"  top GRM eigenvalue lambda1 = {w[-1]:.1f}  (N={N})")


def main():
    print("1000G chr22 pilot: real structure fraction of the LD score")
    pop, superpop = load_panel(PANEL)
    samples, G, pos = parse_vcf(VCF, MAF_MIN, THIN_EVERY)
    sp = np.array([superpop[s] for s in samples])
    G = G.astype(np.float64)

    for label, mask in [
        ("EUR", sp == "EUR"),
        ("AFR", sp == "AFR"),
        ("FULL", np.ones(len(samples), bool)),
    ]:
        # re-thin lightly within cohort to keep eigh/GRM cost sane is NOT done:
        # same SNP set, cohort-specific MAF filter
        Gc = G[mask]
        p = Gc.mean(0) / 2.0
        maf = np.minimum(p, 1 - p)
        keep = maf >= 0.01
        decomposition(Gc[:, keep], f"{label} (cohort-MAF>=0.01: {keep.sum()} SNPs)")


if __name__ == "__main__":
    main()
