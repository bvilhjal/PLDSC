"""What does a EUR LD score cost a multi-ancestry meta-analysis? (1000G chr22)

The first real-data step of the multi-ancestry direction (docs/STRATEGY.md,
docs/MULTI_ANCESTRY.md). From the 1000 Genomes phase-3 superpopulation panels
it computes the within- and cross-ancestry LD scores (MA1, MA7), the LD score
a fixed-effect meta-analysis actually follows (MA2), and the slope factor
kappa = h2_hat / h2 that LDSC incurs when the meta is regressed on the 1000G
EUR score instead (MA6). No GWAS is needed: kappa is a property of the panels,
the ancestry mix and the assumed trans-ancestry genetic correlation rho.

Data: chr22 from the repo's t10 cache (every 10th biallelic SNP, full-panel
MAF >= 0.01; ~17k variants, near HapMap3 density), the S3 mirror of the
phase-3 release (v5a file), and the SHAPEIT4 GRCh37 genetic map. Window
+/-1 cM (LDSC's default). Panels are superpopulations with no within-panel
covariates (standard reference practice): EUR 503, AFR 661, EAS 504,
SAS 489, AMR 347. Regression SNPs: full-panel MAF >= 0.05. A variant enters
ancestry a's GWAS if its MAF in that panel is >= 0.01.

Pre-stated kill rule (from the direction proposal, made precise before
running): at the "realistic" setting -- EUR share 0.75 of total N, the rest
from ONE other superpopulation X, rho = 0.8 for every pair, sample-size
weights -- the direction is KEPT if |kappa - 1| >= 0.10 for at least one
X in {AFR, EAS, SAS, AMR}, and otherwise shrinks to a practice note.
Everything else printed is descriptive.

Data (not committed; .scratch_1000g/ is gitignored). The EBI FTP named in
realdata_genome_loco.download_plan() was unreachable from the 2026-09-25
environment, so the files came from mirrors:

  S3=https://1000genomes.s3.amazonaws.com/release/20130502
  curl -o .scratch_1000g/ALL.chr22.v5a.vcf.gz \
       $S3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  ln -s ALL.chr22.v5a.vcf.gz .scratch_1000g/ALL.chr22.vcf.gz
  curl -o .scratch_1000g/samples.panel $S3/integrated_call_samples_v3.20130502.ALL.panel
  curl -o .scratch_1000g/maps/genetic_maps.b37.tar.gz \
       https://raw.githubusercontent.com/odelaneau/shapeit4/master/maps/genetic_maps.b37.tar.gz
  tar -xzf .scratch_1000g/maps/genetic_maps.b37.tar.gz -C .scratch_1000g/maps chr22.b37.gmap.gz

Run:  python experiments/multiancestry_1000g_meta.py [--thin 3]
      (`--thin 3` is a density-robustness rerun; it parses a denser cache.)

Outcome (2026-09-25; docs/MULTI_ANCESTRY.md Tables 2-5). The kill rule reads
KEEP as pre-stated (Table B). A robustness check under LDSC-style weights,
run after the rule, found two defects in the Table B computation: it
ignored per-SNP N_j, and it kept EUR-rare regression SNPs. Table D corrects
both and still reads KEEP (rule-cell kappa 0.750 AFR / 0.841 EAS /
0.868 SAS / 0.907 AMR under OLS).
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from pldsc.genome import interpolate_cm  # noqa: E402
from pldsc.multiancestry import (  # noqa: E402
    cross_ancestry_ld_scores,
    meta_ld_scores,
    meta_weights,
    single_reference_bias,
)
import realdata_1000g_benchmarks as B  # noqa: E402

CHROM = 22
MAP = os.path.join(HERE, "..", ".scratch_1000g", "maps", "chr22.b37.gmap.gz")
ANCESTRIES = ("EUR", "AFR", "EAS", "SAS", "AMR")
WINDOW_CM = 1.0
MAF_REGRESSION = 0.05
MAF_GWAS = 0.01
EUR_SHARES = (0.9, 0.75, 0.5)
RHOS = (1.0, 0.8, 0.6)


def read_gmap(path):
    """SHAPEIT4 map: header `pos chr cM`; returns strictly increasing bp, cM."""
    bp, cm = [], []
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            p, _, c = line.split()
            bp.append(float(p))
            cm.append(float(c))
    bp, cm = np.asarray(bp), np.asarray(cm)
    keep = np.concatenate(([True], np.diff(bp) > 0))
    return bp[keep], cm[keep]


def load(thin):
    B.THIN[CHROM] = thin
    samples, G, pos = B.load_chrom(CHROM)
    _, superpop = B.load_panel()
    sp = np.array([superpop[s] for s in samples])
    map_bp, map_cm = read_gmap(MAP)
    cm = interpolate_cm(pos.astype(float), map_bp, map_cm)
    panels = [G[sp == a].astype(np.float64) for a in ANCESTRIES]
    return panels, G, cm


def maf(X):
    p = X.mean(0) / 2.0
    return np.minimum(p, 1.0 - p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thin", type=int, default=10)
    args = ap.parse_args()
    t0 = time.time()
    panels, G, cm = load(args.thin)
    reg = np.flatnonzero(maf(G) >= MAF_REGRESSION)
    s = cross_ancestry_ld_scores(panels, cm, window_cm=WINDOW_CM, regression=reg)
    gwas_ok = np.stack([maf(P)[reg] >= MAF_GWAS for P in panels]) & s.available
    A = len(ANCESTRIES)
    E = ANCESTRIES.index("EUR")
    ell_E = s.scores[E, E]
    print(f"chr{CHROM}, thin 1-in-{args.thin}: {G.shape[1]} reference variants, "
          f"{reg.size} regression variants, {cm[-1] - cm[0]:.1f} cM, window +/-{WINDOW_CM} cM "
          f"({time.time() - t0:.0f}s)")

    print("\nTable A. Mean LD scores (diagonal: within-ancestry, noise-floor corrected; "
          "off-diagonal: cross-ancestry)")
    print("      " + " ".join(f"{a:>7s}" for a in ANCESTRIES))
    for a in range(A):
        print(f"{ANCESTRIES[a]:>5s} " + " ".join(f"{s.scores[a, b].mean():7.2f}" for b in range(A)))
    print("  corr(ell^EE, ell^ab) over SNPs: "
          + ", ".join(f"{ANCESTRIES[a]}-{ANCESTRIES[b]} {np.corrcoef(ell_E, s.scores[a, b])[0, 1]:.2f}"
                      for a in range(A) for b in range(a, A) if (a, b) != (E, E)))

    def kappa_for(n, rho, scheme):
        om, _ = meta_weights(n, gwas_ok, variances=s.variances, scheme=scheme)
        R = np.full((A, A), rho)
        np.fill_diagonal(R, 1.0)
        L = meta_ld_scores(s, om, rho=R)
        b = single_reference_bias(L, ell_E)
        within = meta_ld_scores(s, om, rho=np.eye(A))
        return b, L, within

    print("\nTable B. kappa = h2_hat/h2 when the meta is regressed on the EUR score; "
          "offset/mean(L) = intercept shift as a share of the LD-driven excess chi2")
    print(f"{'mix':>16s} {'weights':>11s} " + " ".join(f"{'rho=' + str(r):>13s}" for r in RHOS))
    verdict = {}
    for X in ANCESTRIES[1:]:
        for fE in EUR_SHARES:
            n = np.zeros(A)
            n[E], n[ANCESTRIES.index(X)] = fE, 1.0 - fE
            for scheme in ("sample_size", "ivw"):
                cells = []
                for rho in RHOS:
                    b, L, _ = kappa_for(n, rho, scheme)
                    cells.append(f"{b.kappa:6.3f} ({b.offset / L.mean():+.2f})")
                    if fE == 0.75 and rho == 0.8 and scheme == "sample_size":
                        verdict[X] = b.kappa
                print(f"{f'EUR {fE:.2f}+{X}':>16s} {scheme:>11s} " + " ".join(f"{c:>13s}" for c in cells))

    n4 = np.array([0.6, 0.1, 0.2, 0.0, 0.1])     # illustrative 4-way mix: EUR/AFR/EAS/-/AMR
    print("\n  illustrative 4-way mix EUR 0.6 / EAS 0.2 / AFR 0.1 / AMR 0.1, sample-size weights: "
          + ", ".join(f"rho={r}: {kappa_for(n4, r, 'sample_size')[0].kappa:.3f}" for r in RHOS))

    # Table D: the corrected computation. Two defects of Table B surfaced in a
    # weighted-regression robustness check run after the rule: (1) LDSC's
    # regressor is N_j * ell_j, and a meta reports a smaller N_j where fewer
    # ancestries contribute; (2) Table B's regression set keeps variants that
    # are rare or absent in EUR, which no EUR LD-score file contains (their
    # EUR score is ~0, and LDSC-style weights make them dominate). Table D
    # scales both regressors by N_j and restricts to EUR MAF >= 0.05.
    eur_common = maf(panels[E])[reg] >= MAF_REGRESSION
    print(f"\nTable D (corrected). N_j-scaled regressors, regression set restricted to "
          f"EUR MAF >= {MAF_REGRESSION} ({eur_common.sum()} of {reg.size} variants); "
          "kappa under OLS / LDSC-style weights (offset/mean share, OLS)")
    print(f"{'mix':>16s} " + " ".join(f"{'rho=' + str(r):>22s}" for r in RHOS))
    corrected = {}
    for X in ANCESTRIES[1:]:
        for fE in EUR_SHARES:
            n = np.zeros(A)
            n[E], n[ANCESTRIES.index(X)] = fE, 1.0 - fE
            cells = []
            for rho in RHOS:
                om, n_tot = meta_weights(n, gwas_ok, variances=s.variances, scheme="sample_size")
                R = np.full((A, A), rho)
                np.fill_diagonal(R, 1.0)
                L = meta_ld_scores(s, om, rho=R)[eur_common]
                nj, lE = n_tot[eur_common], ell_E[eur_common]
                ols = single_reference_bias(L, lE, sample_sizes=nj)
                x = lE * nj / nj.mean()
                w = 1.0 / (np.maximum(lE, 1.0) * (1.0 + 0.5 / x.mean() * x) ** 2)
                wtd = single_reference_bias(L, lE, sample_sizes=nj, weights=w)
                share = ols.offset / float(np.mean(L * nj / nj.mean()))
                cells.append(f"{ols.kappa:.3f} / {wtd.kappa:.3f} ({share:+.2f})")
                if fE == 0.75 and rho == 0.8:
                    corrected[X] = (ols.kappa, wtd.kappa)
            print(f"{f'EUR {fE:.2f}+{X}':>16s} " + " ".join(f"{c:>22s}" for c in cells))

    print("\nTable C. Is rho identifiable from the meta chi2 alone? corr of the N_j-scaled "
          "within part and cross part (MA9), on Table D's regression set")
    for X in ANCESTRIES[1:]:
        n = np.zeros(A)
        n[E], n[ANCESTRIES.index(X)] = 0.75, 0.25
        om, n_tot = meta_weights(n, gwas_ok, variances=s.variances, scheme="sample_size")
        full = meta_ld_scores(s, om)[eur_common] * n_tot[eur_common]
        within = meta_ld_scores(s, om, rho=np.eye(A))[eur_common] * n_tot[eur_common]
        print(f"  EUR 0.75 + {X}: corr = {np.corrcoef(within, full - within)[0, 1]:.3f}")

    worst = max(verdict, key=lambda k: abs(verdict[k] - 1))
    keep = abs(verdict[worst] - 1) >= 0.10
    print("\nKill rule as pre-stated (Table B cell: EUR 0.75 + X, rho = 0.8, sample-size weights): "
          + ", ".join(f"{k} kappa {v:.3f}" for k, v in verdict.items())
          + f" -> {'KEEP' if keep else 'KILL (practice note)'} (largest |kappa-1| = "
          f"{abs(verdict[worst] - 1):.3f}, {worst})")
    ok_ols = max(abs(v[0] - 1) for v in corrected.values()) >= 0.10
    ok_w = max(abs(v[1] - 1) for v in corrected.values()) >= 0.10
    print("Same cell, corrected (Table D), OLS / LDSC-weighted: "
          + ", ".join(f"{k} {v[0]:.3f} / {v[1]:.3f}" for k, v in corrected.items())
          + f" -> rule would read {'KEEP' if ok_ols else 'KILL'} / {'KEEP' if ok_w else 'KILL'}")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
