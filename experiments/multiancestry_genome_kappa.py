"""Genome-wide kappa: does the chr22 headline survive all 22 autosomes? (STRATEGY step 1)

STRATEGY Table 3 step 1. The kappa table of docs/MULTI_ANCESTRY.md (Tables 3-4)
was measured on chr22 alone. This script computes the same quantities on every
autosome and pools the regression SNPs genome-wide.

Engine. Per chromosome it calls the tested `cross_ancestry_ld_scores` on
overlapping segments. Each segment carries every reference variant within
+/-1 cM of its regression variants, so its scores equal one whole-chromosome
call (`--check-segments` measures the difference on chr22) while peak memory
stays below about 1 GB on chr2, where one whole-chromosome call would need
several float64 copies of a 2,490 x 139k matrix. Per-chromosome scores are
cached under .scratch_1000g/kappa_gw/ with their configuration, so a rerun
resumes and a changed configuration recomputes.

Data. bigsnpr's 1000G_phase3_common_norel (bigsnpr::download_1000G): the
1000 Genomes phase-3 genotypes of 2,490 of the 2,504 samples (bigsnpr drops
remaining relatives) at 1,664,852 autosomal SNPs, GRCh37, with super-population
labels in the .fam2 file. bigsnpr documents the variant set as common variants
shared with HapMap3 or the UK Biobank. sha256 of the copy used: bed 4a33207e...,
bim 3184ad9c..., fam d16a66fe..., fam2 833478cb... (full values in CHECKSUMS
below). Expected at .scratch_1000g/1000G_phase3_common_norel.{bed,bim,fam,fam2}
(gitignored; a symlink is fine). Genetic map: SHAPEIT4 GRCh37, as for chr22.
Superpopulations: EUR 503, AFR 652, EAS 504, SAS 484, AMR 347.

Settings carried over from multiancestry_1000g_meta.py: window +/-1 cM;
reference variants with full-panel MAF >= 0.01; regression variants with
full-panel MAF >= 0.05; a variant enters ancestry a's GWAS if its MAF in panel
a is >= 0.01; no within-panel covariates.

Pre-stated rules (2026-10-08, written before the first genome-wide run):

  G1 (the STRATEGY step-1 gate). PASS if the rule-cell AFR kappa is <= 0.90
     genome-wide under both OLS and LDSC-style weights. The rule cell is
     MULTI_ANCESTRY Table 4's corrected computation: EUR 0.75 of N plus AFR
     0.25, rho = 0.8 for every pair, sample-size weights, N_j-scaled
     regressors, regression SNPs with full-panel and EUR MAF >= 0.05. Genome-
     wide it also drops regression SNPs in the MHC (chr6:25-34 Mb), as ldsc's
     regression-SNP files do. FAIL means chr22 was unrepresentative and the
     headline shrinks to the genome-wide number.

  G2 (bridge; the data deviate from step 1's protocol, a 1-in-10 thinning of
     the full phase-3 VCFs, of which only chr10/21/22 are local). On chr22
     (VCF cache 1-in-10), chr21 and chr10 (1-in-25), the four rule-cell kappas
     (OLS) are computed from this variant set and from the VCF caches through
     the same code. If any |delta kappa| > 0.03 (twice the density effect in
     MULTI_ANCESTRY Table 3), the genome-wide table is labelled
     variant-set-specific and the thinned-VCF run is required before the
     headline changes. Otherwise the genome-wide numbers may replace it.

Everything else printed is descriptive: the full grid, intercept shares, the
per-chromosome profile, corr(D, X), and the sensitivity rows (MHC kept,
HapMap3-only regression SNPs, IVW weights, eta_X in {0.8, 1.25}). The
200-block jackknife SE is over contiguous blocks of regression SNPs in genome
order. It measures how much kappa varies along the genome, not panel sampling
error: kappa is a property of these 2,490 genomes.

Run:  python experiments/multiancestry_genome_kappa.py [--chroms 1 ... 22]
          [--check-segments] [--figure kappa_profile.pdf]

Outcome (2026-10-08, 309 s; docs/MULTI_ANCESTRY.md Table 6). 1,190,723 rule
regression SNPs on 22 autosomes. Segmented and whole-chromosome scores agree
to 1.2e-12 on chr22, and the VCF-cache route reproduces chr22's recorded
0.750 / 0.841 / 0.868 / 0.907 exactly.

  G1 PASS. Rule-cell AFR kappa 0.737 (OLS; jackknife SE 0.005) / 0.726
     (weighted). EAS 0.846 / 0.837, SAS 0.877 / 0.866, AMR 0.896 / 0.891.
  G2 FAILED its tolerance. The largest bridge difference is 0.039 (chr10 AFR:
     0.726 here, 0.765 on the 1-in-25 cache), so this table is
     variant-set-specific and the headline stays chr22's until the
     thinned-VCF run. The miss cannot reverse G1: AFR would need +0.16.

Per chromosome, OLS kappa spans 0.71-0.82 (AFR). chr12 and chr14 are outliers
because a few extreme-LD regions dominate the OLS fit (on chr14, 0.1% of SNPs
carry 10% of the EUR-score variance). Keeping the MHC lowers OLS kappa by
0.01-0.03 but leaves the weighted fit unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import pldsc  # noqa: E402
from pldsc.genome import interpolate_cm  # noqa: E402
from pldsc.multiancestry import (  # noqa: E402
    cross_ancestry_ld_scores,
    meta_ld_scores,
    meta_weights,
    single_reference_bias,
)
import realdata_1000g_benchmarks as B  # noqa: E402
from multiancestry_1000g_meta import read_gmap  # noqa: E402

SCRATCH = os.path.join(ROOT, ".scratch_1000g")
BFILE = os.path.join(SCRATCH, "1000G_phase3_common_norel")
MAPS = os.path.join(SCRATCH, "maps")
OUT = os.path.join(SCRATCH, "kappa_gw")
HM3 = os.path.join(ROOT, "..", "reference-data", "1000g-phase3", "panels",
                   "kgp_phase3_20130502_superpops_hm3_grch37.npz")
CHECKSUMS = {
    "bed": "4a33207e17ad1e7861db2dc7d849acba7b267e6f649afbcf0dd6603553331647",
    "bim": "3184ad9c8afe03487f6c1ab0abd5ccf2c6665dce3a8cf67977331b1018a8ea2f",
    "fam": "d16a66fee04a6dcfe8844c3100d6aa01df658b56e8849de3aa48451871af4750",
    "fam2": "833478cbc4a10db395a433ba035064fcff8563cbe15b71956e26518b880c4f7c",
}

ANCESTRIES = ("EUR", "AFR", "EAS", "SAS", "AMR")
E = 0
WINDOW_CM = 1.0
MAF_REF = 0.01
MAF_REGRESSION = 0.05
MAF_GWAS = 0.01
MHC = (6, 25_000_000, 34_000_000)
SEGMENT = 4096                      # regression variants per engine call
EUR_SHARES = (0.9, 0.75, 0.5)
RHOS = (1.0, 0.8, 0.6)
GATE = 0.90
BRIDGE_TOL = 0.03
BRIDGE = {22: 10, 21: 25, 10: 25}   # chromosome -> VCF-cache thinning
N_BLOCKS = 200


# ------------------------------------------------------------------- data

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_bim(prefix):
    chrom, rsid, pos = [], [], []
    with open(prefix + ".bim") as fh:
        for line in fh:
            f = line.split()
            chrom.append(int(f[0]))
            rsid.append(f[1])
            pos.append(int(f[3]))
    return np.array(chrom, dtype=np.int64), np.array(rsid), np.array(pos, dtype=np.int64)


def read_superpops(prefix):
    """Super-population per .fam row, from bigsnpr's tab-separated .fam2."""
    with open(prefix + ".fam") as fh:
        iids = [line.split()[1] for line in fh]
    sp = {}
    with open(prefix + ".fam2") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id, i_sp = header.index("sample.ID"), header.index("Super Population")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sp[parts[i_id]] = parts[i_sp]
    return np.array([sp[i] for i in iids])


# 2-bit PLINK codes, low bits first: 00 hom A1 (2 copies of A1), 01 missing,
# 10 het, 11 hom A2. Genotypes count the bim A1 allele; every panel shares it.
_LUT = np.array([[(2, -1, 1, 0)[(b >> s) & 3] for s in (0, 2, 4, 6)] for b in range(256)],
                dtype=np.int8)


def read_bed(prefix, n, m_total, start, stop):
    """(n, stop - start) int8 A1 counts for variants [start, stop)."""
    with open(prefix + ".bed", "rb") as fh:
        if fh.read(3) != b"\x6c\x1b\x01":
            raise ValueError("not a SNP-major PLINK .bed")
    bpv = (n + 3) // 4
    mm = np.memmap(prefix + ".bed", dtype=np.uint8, mode="r", offset=3, shape=(m_total, bpv))
    G = _LUT[np.asarray(mm[start:stop])].reshape(stop - start, 4 * bpv)[:, :n]
    if np.any(G < 0):
        raise ValueError("missing genotypes are not supported")
    return np.ascontiguousarray(G.T)


def chrom_cm(chrom, pos):
    map_bp, map_cm = read_gmap(os.path.join(MAPS, f"chr{chrom}.b37.gmap.gz"))
    return interpolate_cm(pos.astype(float), map_bp, map_cm)


# ----------------------------------------------------------------- scoring

def score_chromosome(panels, cm, maf_full, *, segment=SEGMENT):
    """cross_ancestry_ld_scores over exact-window segments of one chromosome."""
    reg = np.flatnonzero(maf_full >= MAF_REGRESSION)
    parts = []
    for s in range(0, reg.size, segment):
        chunk = reg[s:s + segment]
        lo = int(np.searchsorted(cm, cm[chunk[0]] - WINDOW_CM, side="left"))
        hi = int(np.searchsorted(cm, cm[chunk[-1]] + WINDOW_CM, side="right"))
        parts.append(cross_ancestry_ld_scores([P[:, lo:hi] for P in panels], cm[lo:hi],
                                              window_cm=WINDOW_CM, regression=chunk - lo))
    maf = np.stack([np.minimum(P[:, reg].mean(0) / 2.0, 1.0 - P[:, reg].mean(0) / 2.0)
                    for P in panels])
    return dict(scores=np.concatenate([p.scores for p in parts], axis=2),
                available=np.concatenate([p.available for p in parts], axis=1),
                variances=np.concatenate([p.variances for p in parts], axis=1),
                maf=maf, reg=reg)


def config(data_sha):
    return json.dumps(dict(data=data_sha, window_cm=WINDOW_CM, maf_ref=MAF_REF,
                           maf_regression=MAF_REGRESSION, pldsc=pldsc.__version__,
                           engine="cross_ancestry_ld_scores/segments"), sort_keys=True)


def save(path, cfg, **arrays):
    tmp = path[:-4] + ".tmp.npz"
    np.savez(tmp, config=np.array(cfg), **arrays)
    os.replace(tmp, path)


def cached(path, cfg):
    if not os.path.exists(path):
        return None
    d = np.load(path, allow_pickle=False)
    return d if str(d["config"]) == cfg else None


def bigsnpr_chromosome(c, bim, sp, data_sha, *, check=False):
    path = os.path.join(OUT, f"chr{c}.npz")
    cfg = config(data_sha)
    d = cached(path, cfg)
    if d is not None and not check:
        return d
    chrom, rsid, pos = bim
    rows = np.flatnonzero(chrom == c)
    start, stop = int(rows[0]), int(rows[-1]) + 1
    if rows.size != stop - start or np.any(np.diff(pos[start:stop]) < 0):
        raise ValueError(f"chr{c} variants must be contiguous and position-sorted")
    t0 = time.time()
    G = read_bed(BFILE, sp.size, chrom.size, start, stop)
    p = G.mean(0) / 2.0
    keep = np.minimum(p, 1.0 - p) >= MAF_REF
    G = G[:, keep]
    maf_full = np.minimum(p, 1.0 - p)[keep]
    cm = chrom_cm(c, pos[start:stop][keep])
    panels = [G[sp == a] for a in ANCESTRIES]
    del G
    r = score_chromosome(panels, cm, maf_full)
    if check:
        whole = cross_ancestry_ld_scores(panels, cm, window_cm=WINDOW_CM, regression=r["reg"])
        diff = np.abs(whole.scores - r["scores"]).max()
        print(f"  [check] chr{c}: max |segmented - whole-chromosome| score = {diff:.2e}")
        if d is not None:
            return d
    reg = r["reg"]
    p_kept = pos[start:stop][keep]
    save(path, cfg, scores=r["scores"], available=r["available"], variances=r["variances"],
         maf=r["maf"], maf_full=maf_full[reg], pos=p_kept[reg], cm=cm[reg],
         rsid=rsid[start:stop][keep][reg], n_reference=np.array(keep.sum()))
    print(f"  [scores] chr{c}: {keep.sum()} reference, {reg.size} regression variants "
          f"({time.time() - t0:.0f}s)")
    return np.load(path, allow_pickle=False)


def vcf_chromosome(c, thin):
    """The same scoring on the repo's VCF cache (2,504 samples, 1-in-thin SNPs)."""
    path = os.path.join(OUT, f"vcf_chr{c}_t{thin}.npz")
    cfg = config(f"vcf-cache-t{thin}")
    d = cached(path, cfg)
    if d is not None:
        return d
    B.THIN[c] = thin
    samples, G, pos = B.load_chrom(c)
    _, superpop = B.load_panel()
    sp = np.array([superpop[s] for s in samples])
    p = G.mean(0) / 2.0
    maf_full = np.minimum(p, 1.0 - p)
    cm = chrom_cm(c, pos)
    r = score_chromosome([G[sp == a] for a in ANCESTRIES], cm, maf_full)
    save(path, cfg, scores=r["scores"], available=r["available"], variances=r["variances"],
         maf=r["maf"], maf_full=maf_full[r["reg"]], pos=pos[r["reg"]], cm=cm[r["reg"]])
    return np.load(path, allow_pickle=False)


# ------------------------------------------------------------------- kappa

def pooled(ds, chroms):
    keys = ("scores", "available", "variances", "maf", "maf_full", "pos")
    D = {k: np.concatenate([d[k] for d in ds], axis=-1) for k in keys}
    D["chrom"] = np.concatenate([np.full(d["pos"].size, c) for d, c in zip(ds, chroms)])
    if all("rsid" in d.files for d in ds):
        D["rsid"] = np.concatenate([d["rsid"] for d in ds])
    D["gwas_ok"] = (D["maf"] >= MAF_GWAS) & D["available"]
    return D


def regression_set(D, *, mhc=False, hm3=None):
    sel = D["maf"][E] >= MAF_REGRESSION
    if not mhc:
        c, lo, hi = MHC
        sel &= ~((D["chrom"] == c) & (D["pos"] >= lo) & (D["pos"] < hi))
    if hm3 is not None:
        sel &= np.isin(D["rsid"], hm3)
    return sel


def meta_score(D, X, fE, rho, *, scheme="sample_size", eta=None):
    A = len(ANCESTRIES)
    n = np.zeros(A)
    n[E], n[X] = fE, 1.0 - fE
    om, n_tot = meta_weights(n, D["gwas_ok"], variances=D["variances"], scheme=scheme)
    R = np.full((A, A), rho)
    np.fill_diagonal(R, 1.0)
    return meta_ld_scores(D["scores"], om, rho=R, h2_ratio=eta), n_tot, om


def fit(L, lE, nj):
    """Table D's two fits: OLS and LDSC-style weights; plus the OLS offset share."""
    ols = single_reference_bias(L, lE, sample_sizes=nj)
    x = lE * nj / nj.mean()
    w = 1.0 / (np.maximum(lE, 1.0) * (1.0 + 0.5 / x.mean() * x) ** 2)
    wtd = single_reference_bias(L, lE, sample_sizes=nj, weights=w)
    return ols.kappa, wtd.kappa, ols.offset / float(np.mean(L * nj / nj.mean()))


def kappa_cell(D, sel, X, fE=0.75, rho=0.8, **kw):
    L, n_tot, _ = meta_score(D, X, fE, rho, **kw)
    return fit(L[sel], D["scores"][E, E][sel], n_tot[sel])


def jackknife(D, sel, X, fE=0.75, rho=0.8, n_blocks=N_BLOCKS):
    L, n_tot, _ = meta_score(D, X, fE, rho)
    L, lE, nj = L[sel], D["scores"][E, E][sel], n_tot[sel]
    edges = np.linspace(0, L.size, n_blocks + 1).astype(int)
    th = np.empty((n_blocks, 2))
    for b in range(n_blocks):
        keep = np.ones(L.size, dtype=bool)
        keep[edges[b]:edges[b + 1]] = False
        th[b] = fit(L[keep], lE[keep], nj[keep])[:2]
    return np.sqrt((n_blocks - 1) / n_blocks * ((th - th.mean(0)) ** 2).sum(0))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chroms", type=int, nargs="+", default=list(range(1, 23)))
    ap.add_argument("--check-segments", action="store_true",
                    help="compare segmented vs one whole-chromosome call on chr22")
    ap.add_argument("--figure", default=None, help="write the per-chromosome kappa profile")
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    got = {k: sha256(f"{BFILE}.{k}") for k in CHECKSUMS}
    for k, v in got.items():
        if v != CHECKSUMS[k]:
            print(f"  [warn] {k} sha256 {v[:12]}... differs from the recorded copy")
    data_sha = hashlib.sha256("".join(got[k] for k in sorted(got)).encode()).hexdigest()
    bim = read_bim(BFILE)
    sp = read_superpops(BFILE)
    print(f"data: {sp.size} samples ("
          + ", ".join(f"{a} {(sp == a).sum()}" for a in ANCESTRIES)
          + f"), {bim[0].size} variants")

    if args.check_segments:
        bigsnpr_chromosome(22, bim, sp, data_sha, check=True)
    ds = [bigsnpr_chromosome(c, bim, sp, data_sha) for c in args.chroms]
    D = pooled(ds, args.chroms)
    sel = regression_set(D)
    A = len(ANCESTRIES)
    print(f"\n{len(args.chroms)} chromosomes: {sum(int(d['n_reference']) for d in ds)} reference, "
          f"{D['pos'].size} regression variants; {sel.sum()} in the rule regression set "
          f"(EUR MAF >= {MAF_REGRESSION}, MHC excluded) ({time.time() - t0:.0f}s)")

    print("\nTable GW1. Mean LD scores over the rule regression set (diagonal MA7-corrected)")
    print("      " + " ".join(f"{a:>7s}" for a in ANCESTRIES))
    for a in range(A):
        print(f"{ANCESTRIES[a]:>5s} " + " ".join(f"{D['scores'][a, b][sel].mean():7.2f}" for b in range(A)))

    print("\nTable GW2. kappa, genome-wide: OLS / LDSC-style weights (OLS offset share)")
    print(f"{'mix':>16s} " + " ".join(f"{'rho=' + str(r):>22s}" for r in RHOS))
    rule = {}
    for X in range(1, A):
        for fE in EUR_SHARES:
            cells = []
            for rho in RHOS:
                ko, kw, share = kappa_cell(D, sel, X, fE, rho)
                cells.append(f"{ko:.3f} / {kw:.3f} ({share:+.2f})")
                if fE == 0.75 and rho == 0.8:
                    rule[ANCESTRIES[X]] = (ko, kw)
            print(f"{f'EUR {fE:.2f}+{ANCESTRIES[X]}':>16s} " + " ".join(f"{c:>22s}" for c in cells))

    print(f"\nTable GW3. Rule cell (EUR 0.75 + X, rho = 0.8) with {N_BLOCKS}-block jackknife SE "
          "(variation along the genome)")
    for X in range(1, A):
        se = jackknife(D, sel, X)
        ko, kw = rule[ANCESTRIES[X]]
        print(f"  {ANCESTRIES[X]}: OLS {ko:.3f} (SE {se[0]:.3f}), weighted {kw:.3f} (SE {se[1]:.3f})")

    print("\nTable GW4. Sensitivity at the rule cell, OLS / weighted")
    hm3 = None
    if os.path.exists(HM3) and "rsid" in D:
        hm3 = np.load(HM3, allow_pickle=False)["ids"]
    rows = [("as ruled", sel, {}), ("MHC kept", regression_set(D, mhc=True), {}),
            ("IVW weights", sel, {"scheme": "ivw"})]
    if hm3 is not None:
        hs = regression_set(D, hm3=hm3)
        rows.append((f"HapMap3 regression SNPs only ({hs.sum()})", hs, {}))
    for label, s, kw_ in rows:
        print(f"  {label:>34s}: " + ", ".join(
            "{} {:.3f} / {:.3f}".format(ANCESTRIES[X], *kappa_cell(D, s, X, **kw_)[:2]) for X in range(1, A)))
    for eta_x in (0.8, 1.25):
        cells = []
        for X in range(1, A):
            eta = np.ones(A)
            eta[X] = eta_x
            cells.append("{} {:.3f} / {:.3f}".format(ANCESTRIES[X], *kappa_cell(D, sel, X, eta=eta)[:2]))
        print(f"  {f'eta_X = {eta_x}':>34s}: " + ", ".join(cells))

    print("\nTable GW5. corr(D, X) of the N_j-scaled within and cross parts (MA9), EUR 0.75")
    for X in range(1, A):
        n = np.zeros(A)
        n[E], n[X] = 0.75, 0.25
        om, n_tot = meta_weights(n, D["gwas_ok"], variances=D["variances"], scheme="sample_size")
        full = meta_ld_scores(D["scores"], om)[sel] * n_tot[sel]
        within = meta_ld_scores(D["scores"], om, rho=np.eye(A))[sel] * n_tot[sel]
        print(f"  EUR 0.75 + {ANCESTRIES[X]}: corr = {np.corrcoef(within, full - within)[0, 1]:.3f}")

    print("\nTable GW6. Rule-cell kappa per chromosome, OLS")
    print("  chr " + " ".join(f"{a:>6s}" for a in ANCESTRIES[1:]))
    prof = np.empty((len(args.chroms), A - 1))
    for i, c in enumerate(args.chroms):
        sc = sel & (D["chrom"] == c)
        prof[i] = [kappa_cell(D, sc, X)[0] for X in range(1, A)]
        print(f"  {c:>3d} " + " ".join(f"{k:6.3f}" for k in prof[i]))
    print("  range: " + ", ".join(f"{ANCESTRIES[X + 1]} {prof[:, X].min():.3f}-{prof[:, X].max():.3f}"
                                  for X in range(A - 1)))
    if args.figure:
        from pldsc.plotting import kappa_profile
        kappa_profile(prof, args.chroms, ancestries=ANCESTRIES[1:], path=args.figure,
                      title="1000G rule-cell kappa per chromosome")

    print("\nTable GW7. Bridge (G2): rule-cell kappa (OLS), this variant set vs the VCF cache")
    worst = 0.0
    for c, thin in BRIDGE.items():
        if c not in args.chroms:
            continue
        V = pooled([vcf_chromosome(c, thin)], [c])
        vs = regression_set(V)
        cs = sel & (D["chrom"] == c)
        cells = []
        for X in range(1, A):
            kb, kv = kappa_cell(D, cs, X)[0], kappa_cell(V, vs, X)[0]
            worst = max(worst, abs(kb - kv))
            cells.append(f"{ANCESTRIES[X]} {kb:.3f} vs {kv:.3f}")
        print(f"  chr{c} (VCF 1-in-{thin}, {vs.sum()} vs {cs.sum()} regression SNPs): " + ", ".join(cells))

    ko, kw = rule["AFR"]
    g1 = ko <= GATE and kw <= GATE
    print(f"\nG1: genome-wide rule-cell AFR kappa {ko:.3f} (OLS) / {kw:.3f} (weighted) -> "
          f"{'PASS' if g1 else 'FAIL'} (gate <= {GATE})")
    print(f"G2: largest bridge |delta kappa| = {worst:.3f} -> "
          + ("variant set does not matter at this tolerance" if worst <= BRIDGE_TOL
             else "VARIANT-SET-SPECIFIC; the thinned-VCF run is required")
          + f" (tolerance {BRIDGE_TOL})")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
