"""Genome-wide LOCO structure-share run, built on the tested `pldsc.genome` engine.

This is the driver that scales the B2/B6 benchmarks from three thinned
chromosomes to a genuine leave-one-chromosome-out sweep with genetic-distance
(cM) windows and the O(M N^2) spectral genome-wide score (no M x M matrix). It is
the shared infrastructure for `docs/WHEN_STRUCTURE_MATTERS.md` (#6) and the
per-axis study (#4).

For each target chromosome c and each cohort, it:
  1. builds a LOCO GRM from a thinned SNP set pooled over every OTHER cached
     chromosome (never chromosome c), and forms the intercept-only P_c;
  2. computes ordinary and P_c-whitened per-SNP LD scores on chromosome c's
     common SNPs -- windowed in cM (via the Beagle GRCh37 map) and genome-wide
     via the spectral identity;
  3. reports the per-chromosome structure share, then the SNP-weighted
     genome-wide aggregate across the ancestry gradient.

Run on the cached three chromosomes (rotating LOCO: target 22 <- 21+10 fodder,
target 21 <- 22+10, target 10 <- 22+21) it is a complete end-to-end
demonstration; pointed at all 22 autosomes (see `--- SCALING ---` and
`docs/REAL_DATA_SURVEY.md`) it is the genome-wide run itself. Every numerical
primitive it calls is unit-tested in `tests/test_genome.py` and reconciled with
the brute-force benchmark in the B2/B6 suite.

Run:  conda run -n pldsc python experiments/realdata_genome_loco.py

  --- SCALING to all 22 autosomes ---
Memory is the only constraint, and this demo driver does NOT yet stream: `run()`
loads every requested chromosome's genotypes into `geno` up front. That is fine
for the cached, thinned three chromosomes (~0.5 GB) but would exhaust memory on
22 full autosomes (~10 GB per chromosome at common-SNP density). The genome-wide
run therefore needs one refactor, noted here rather than silently shipped:

  * hold only a THINNED common-SNP fodder subset resident for all 22 chromosomes
    (~50-100k SNPs genome-wide is ample for the GRM; N x N eigendecomposition is
    trivial at N=2504), and
  * stream each target chromosome's full common SNPs in for its LD-score pass,
    then release it (one chromosome resident at a time).

The download is ~15 GB of phase-3 VCFs (chr1 1.1 GB ... chr22 0.2 GB); the
parser/cacher in `realdata_1000g_benchmarks.py` already streams and thins them.
The genetic maps are the 22.8 MB Beagle GRCh37 set already in
`.scratch_1000g/maps/`. See `download_plan()` for the concrete manifest.
(Reproducibility note, 2026-08-30: the scratch tree was rebuilt on 2026-08-11
with chr21/chr22 + baselineLD only -- the maps, the chr10 VCF and the chr10
cache are no longer present, so regenerating the 3-chromosome figures below
requires re-fetching them per the download plan. The D1/D2 gate runs are
unaffected.)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from pldsc.genome import (  # noqa: E402
    grm_eigh,
    interpolate_cm,
    pc_operator,
    read_plink_genetic_map,
    spectral_full_scores,
    standardize,
    structure_share,
    unit_columns,
    windowed_scores,
)
import realdata_1000g_benchmarks as B  # noqa: E402

MAPS = os.path.join(HERE, "..", ".scratch_1000g", "maps")
COHORTS = ("EUR", "AMR", "AFR", "FULL")
WINDOW_CM = 1.0            # LD-score window in genetic distance
GRM_H2 = 0.5              # assumed h2 in the P_c operator
GRM_THIN = 4              # extra thinning of fodder SNPs for the LOCO GRM


def map_path(chrom):
    return os.path.join(MAPS, f"plink.chr{chrom}.GRCh37.map")


def cohort_common(G, mask, pos):
    """Common (MAF>=0.05 in FULL), cohort-polymorphic SNP columns for a cohort."""
    keep = np.flatnonzero(np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2) >= B.MAF_COMMON)
    Gsub = G[mask][:, keep]
    cm = B._cmask(Gsub)
    return Gsub[:, cm], pos[keep][cm]


def loco_grm_eig(geno, mask, targets, this_chrom):
    """Eigendecomposition of a LOCO GRM: thinned fodder from all but this_chrom."""
    fodder = []
    for c in targets:
        if c == this_chrom:
            continue
        G = geno[c][mask]
        fodder.append(G[:, ::GRM_THIN])
    X = standardize(np.hstack(fodder))
    return grm_eigh(X)


def run(chroms):
    t0 = time.time()
    print("Genome-wide LOCO structure-share run (pldsc.genome engine)")
    print(f"  chromosomes: {chroms}   window: {WINDOW_CM} cM   P_c h2: {GRM_H2}")
    pop, superpop = B.load_panel()

    geno, posbp, cmpos, samples = {}, {}, {}, None
    for c in chroms:
        s, G, p = B.load_chrom(c)
        if samples is None:
            samples = s
        else:
            assert s == samples, f"sample order mismatch on chr{c}"
        geno[c], posbp[c] = G, p
        mbp, mcm = read_plink_genetic_map(map_path(c))
        cmpos[c] = interpolate_cm(p.astype(float), mbp, mcm)
    masks = B.cohort_masks(samples, pop, superpop)
    print(f"  loaded {len(chroms)} chromosomes, {len(samples)} samples "
          f"({time.time() - t0:.0f}s)\n")

    # accumulate SNP-weighted structure shares per cohort
    agg = {cn: {"win_o": 0.0, "win_w": 0.0, "full_o": 0.0, "full_w": 0.0, "m": 0}
           for cn in COHORTS}

    header = (f"  {'chr':>4s} {'cohort':>6s} {'M':>7s} {'ord_win':>8s} "
              f"{'struct%win':>11s} {'struct%full':>12s}")
    for target in chroms:
        print(header)
        for cn in COHORTS:
            mask = masks[cn]
            Gc, bp = cohort_common(geno[target], mask, posbp[target])
            cm = interpolate_cm(bp.astype(float),
                                *read_plink_genetic_map(map_path(target)))
            X = standardize(Gc)
            Xn, _ = unit_columns(X)

            w, U = loco_grm_eig(geno, mask, chroms, target)
            P = pc_operator(w, U, GRM_H2)
            Xw_n, _ = unit_columns(P @ X)

            o_win = windowed_scores(Xn, cm, WINDOW_CM)
            w_win = windowed_scores(Xw_n, cm, WINDOW_CM)
            o_full = spectral_full_scores(Xn)
            w_full = spectral_full_scores(Xw_n)

            s_win = structure_share(o_win, w_win) * 100
            s_full = structure_share(o_full, w_full) * 100
            m = X.shape[1]
            a = agg[cn]
            a["win_o"] += o_win.sum(); a["win_w"] += w_win.sum()
            a["full_o"] += o_full.sum(); a["full_w"] += w_full.sum()
            a["m"] += m
            print(f"  {target:>4d} {cn:>6s} {m:>7d} {o_win.mean():>8.3f} "
                  f"{s_win:>10.1f}% {s_full:>11.1f}%")
        print()

    print("  === SNP-weighted genome-wide aggregate ===")
    print(f"  {'cohort':>6s} {'M':>8s} {'struct%win':>11s} {'struct%full':>12s}")
    for cn in COHORTS:
        a = agg[cn]
        sw = (a["win_o"] - a["win_w"]) / a["win_o"] * 100
        sf = (a["full_o"] - a["full_w"]) / a["full_o"] * 100
        print(f"  {cn:>6s} {a['m']:>8d} {sw:>10.1f}% {sf:>11.1f}%")
    print(f"\n  runtime {time.time() - t0:.0f}s")
    print("  NOTE: three chromosomes is a demonstration of the machinery; the")
    print("  genome-wide headline needs all 22 autosomes (see download_plan()).")


def download_plan():
    base = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
            "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b."
            "20130502.genotypes.vcf.gz")
    print("Manifest for the full 22-autosome run:")
    print("  Everything below is fetched by scripts/fetch_1000g.py (S3 `v5a`")
    print("  primary, sha256-verified, maps converted to the plink format this")
    print("  driver reads). One command:")
    print("    python scripts/fetch_1000g.py --chrom 1 ... 22 --maps all \\")
    print("        --panel --build-cache 1 ... 22")
    print("  The per-file EBI `v5b` URLs the plan used to print (a mirror that")
    print("  was unreachable from this environment and serves a different")
    print("  re-release than the recorded results used):")
    for c in range(1, 23):
        print(f"    chr{c}: {base.format(c=c)}")
    print("  genetic maps: .scratch_1000g/maps/, both formats, from the same")
    print("  script (the historical Beagle plink.GRCh37 host is gone; the")
    print("  script converts the SHAPEIT4 b37 map marker-for-marker)")
    print("  panel: .scratch_1000g/samples.panel, from the same script")
    print("  GRM fodder: thin each chromosome to a common-SNP subset; ~50-100k")
    print("  SNPs genome-wide is ample. Per-target LD scores hold one chromosome.")


if __name__ == "__main__":
    if "--plan" in sys.argv:
        download_plan()
    else:
        cached = [c for c in (22, 21, 10)
                  if os.path.exists(B.cache_paths(c)[0])
                  or os.path.exists(B.vcf_path(c))]
        run(cached)
