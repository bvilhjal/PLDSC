"""The GBMI asthma casualty: rho from the ancestry metas, then h2 on two scores (STRATEGY step 3)

GBMI estimated SNP-heritability from its all-biobank meta-analyses with LDSC
and UK Biobank (European) LD scores. Asthma: h2_SNP = 0.085 on the liability
scale, with population prevalence set to the pooled sample prevalence
(Tsuo et al. 2022, doi 10.1016/j.xgen.2022.100212, citing the GBMI flagship
paper's Table S9, doi 10.1016/j.xgen.2022.100192). MULTI_ANCESTRY MA6 says a
EUR score understates h2 for a multi-ancestry meta. This script measures how
much on the real asthma meta. It estimates rho from GBMI's single-ancestry
metas, as STRATEGY step 3 asks, instead of bracketing it.

Data. GBMI release 052021, asthma, six files: the all-biobank meta (ALL) and
the EUR, AFR, EAS, AMR and SAS metas, from gbmi-sumstats.s3.amazonaws.com, in
.scratch_gbmi/asthma/ (gitignored; sha256 values in SHA256 below). They are
GRCh38 and joined to the 1000G scores by rsID only. LD scores: the genome-wide
1000G cache that multiancestry_genome_kappa.py writes (run it first). It holds
1.36M regression variants and the 5 x 5 within- and cross-ancestry scores,
over a reference of M = 1,664,534 variants.

Method, fixed before the GBMI files were read (2026-10-08):

  Variants. Biallelic SNPs whose rsID occurs once in the file, dropping GBMI's
  is_strand_flip and is_diff_AF_gnomAD flags. A stratum's (REF, ALT) must equal
  the ALL file's or be its swap, in which case z is negated; anything else is
  dropped. Regression variants: EUR MAF >= 0.05 in 1000G, outside the MHC
  (chr6:25-34 Mb), polymorphic in every panel the regression uses, and
  chi2 <= max(80, 0.001 N) in every file it uses (ldsc's rule).

  Sample size. Per variant and file, N is the effective N implied by the
  standard error, N = 4 / (se^2 * 2 p (1 - p)), with p the file's own meta
  allele frequency. Totals of N_case and N_ctrl would overstate the N of a
  meta whose biobanks differ in case fraction. h2 is therefore on the
  observed scale at P = 0.5. It is converted to the liability scale only for
  the comparison with GBMI, using K = 153,763 / 1,800,785, GBMI's choice.

  Regressions. chi2_j ~ a + b (N_j l_j) gives h2 = b M. z_aj z_bj ~ a + b
  (sqrt(N_aj N_bj) l^ab_j) gives the genetic covariance g_ab = b M. Then
  rho_ab = g_ab / sqrt(h2_a h2_b), on the standardized scale MA2 uses. Weights
  are 1 / max(l, 1), with l the regression's own within score, or the
  geometric mean of the two for a cross regression. The intercept is free.
  SEs come from a 200-block jackknife over contiguous variants in genome
  order. All quantities that share a selection are deleted together, so
  rho and the casualty difference have paired SEs.

  rho used in MA2. A pair's estimate is used, clipped to [-1, 1], when both
  strata have an h2 z-score > 4. Otherwise that pair is bracketed at
  rho = 0.6 / 0.8 / 1.0 and the casualty is reported for each.

  The casualty. The ALL file's chi2 is regressed on N_j l^EUR_j (what GBMI
  did, with a 1000G EUR score in place of the UKB one) and on N_j L_j. Here
  L_j is MA2 with inverse-variance weights from the per-variant stratum N,
  the 1000G variances, the rho above and eta = 1. N_j is the ALL file's N.
  Both fits use the same variants, weights and blocks. Regression variants
  must also appear in the EUR stratum file.

Pre-stated rules:

  C1 (STRATEGY step 3's gate). CASUALTY if |h2_meta - h2_EUR| > 2 SE(h2_EUR),
     with the jackknife SE of the EUR-score estimate, at the rho actually used.
     Otherwise this becomes a practice note (STRATEGY section 7).
  C2 (does MA6 predict it?). Descriptive check: the observed ratio
     h2_EUR / h2_meta should lie within 2 paired-jackknife SEs of the kappa
     that single_reference_bias predicts from the same L, l^EUR, N_j and
     weights.

Everything else is descriptive: per-stratum h2, the full rho matrix, the
rho brackets, the ALL-file coverage check (strata N summed over the ALL N),
the N_j = sum-of-strata row, and the liability-scale comparison with GBMI's
0.085. That comparison cannot be exact. GBMI's LD scores were UKB-based and
denser, and its N convention is not stated.

Run:  python -I experiments/gbmi_asthma_casualty.py

Outcome (2026-10-08, 110 s; docs/MULTI_ANCESTRY.md Table 7). 953,596
regression SNPs. Strata N sum to the ALL N (median ratio 1.011).

  Per stratum, h2 z is 10.7 for EUR and 9.4 for EAS. AFR (2.1), SAS (3.4)
  and AMR (0.4) are too small, so only EUR-EAS is estimated: rho = 0.755
  (SE 0.046). The other pairs are bracketed. Their point estimates
  (EUR-AFR 0.73, EUR-SAS 0.83) fall inside the bracket.
  C1 NOT MET ROBUSTLY. h2_EUR is 0.133 (SE 0.012). h2_meta is 0.159 / 0.156 /
     0.153 at brackets 0.6 / 0.8 / 1.0. The shift is 2.2 / 1.9 / 1.6 SE(h2_EUR),
     so at the central bracket asthma is a practice note. The shift itself is
     precise: +0.023 with a paired SE of 0.002. The gate's SE is mostly
     genome-wide noise that the two estimates share.
  C2 consistent (ratio 0.852 vs kappa 0.850). Table C4 shows the check is
     nearly mechanical: the two regressors correlate 0.9965.

Post hoc (Table C4): the composition-weighted mean of the per-stratum h2,
each on its own ancestry's score, is 0.154 (SE 0.012). That matches h2_meta
(0.156) and is 14% above GBMI's method on the same scores (0.133). GBMI's
0.085 is not comparable to 0.133: the LD scores differ (UKB, denser) and
GBMI's N convention is unstated. If the ratio carries over to GBMI's scores,
their value would be about 0.10 on the meta score. That is an extrapolation,
not a measurement.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import sys
import time
from statistics import NormalDist

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from pldsc.multiancestry import meta_ld_scores, meta_weights, single_reference_bias  # noqa: E402
import multiancestry_genome_kappa as GK  # noqa: E402

GBMI_DIR = os.path.join(ROOT, ".scratch_gbmi", "asthma")
PATTERN = "Asthma_Bothsex_{}inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
STRATA = {"ALL": "", "EUR": "eur_", "AFR": "afr_", "EAS": "eas_", "SAS": "sas_", "AMR": "amr_"}
SHA256 = {                          # the 2026-10-08 download, recorded before the first run
    "ALL": "eba2007927274725e05d8f77049fab14e9f29b6d308df52049e8c7b873499591",
    "EUR": "83fec73908cebc20ec896597b36c9aee47527bd90aceaabf82154d18df682a88",
    "AFR": "80d59a041dc78eefba9576e145e86394abddbc8b23eb0177ebf7a73e8c9fec7b",
    "EAS": "3c24d2df1fb0149502c472b047f937e8a4f88d3b2f3116551989c5436f516f86",
    "SAS": "f29d8d8c51cdf0c3b714363fb980aaeca082001733f0ed346ef2d9d04179028e",
    "AMR": "aa3ea725321f3d22a354b2e42a5caeaec1f6ea7b803f181e81a4e1e1918dbe98",
}
ANC = GK.ANCESTRIES                 # EUR AFR EAS SAS AMR (score-cache order)
E = 0
N_CASES, N_TOTAL = 153_763, 153_763 + 1_647_022
H2_GBMI = 0.085
Z_MIN = 4.0
RHO_BRACKET = (0.6, 0.8, 1.0)
N_BLOCKS = 200


# -------------------------------------------------------------------- data

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_scores():
    chroms = list(range(1, 23))
    ds = [np.load(os.path.join(GK.OUT, f"chr{c}.npz"), allow_pickle=False) for c in chroms]
    D = GK.pooled(ds, chroms)
    D["M"] = sum(int(d["n_reference"]) for d in ds)
    return D


def read_gbmi(path, index):
    """Rows of one GBMI file at the cached variants: index, alleles, z, SE-implied N."""
    rows, rsids, ref, alt, z, n = [], [], [], [], [], []
    count = {}
    with gzip.open(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        c = {k: h.index(k) for k in ("REF", "ALT", "rsid", "all_meta_AF", "inv_var_meta_beta",
                                     "inv_var_meta_sebeta", "is_strand_flip", "is_diff_AF_gnomAD")}
        for line in fh:
            rs = line.split("\t", c["rsid"] + 1)[c["rsid"]]
            j = index.get(rs)
            if j is None:
                continue
            count[rs] = count.get(rs, 0) + 1
            f = line.rstrip("\n").split("\t")
            if (len(f[c["REF"]]) != 1 or len(f[c["ALT"]]) != 1 or f[c["is_strand_flip"]] == "yes"
                    or f[c["is_diff_AF_gnomAD"]] == "yes"):
                continue
            se, p = float(f[c["inv_var_meta_sebeta"]]), float(f[c["all_meta_AF"]])
            if not (se > 0 and 0 < p < 1):
                continue
            rows.append(j)
            rsids.append(rs)
            ref.append(f[c["REF"]])
            alt.append(f[c["ALT"]])
            z.append(float(f[c["inv_var_meta_beta"]]) / se)
            n.append(4.0 / (se * se * 2.0 * p * (1.0 - p)))
    dup = {rs for rs, k in count.items() if k > 1}
    keep = np.array([rs not in dup for rs in rsids], dtype=bool)
    return (np.array(rows, dtype=np.int64)[keep], np.array(ref)[keep], np.array(alt)[keep],
            np.array(z)[keep], np.array(n)[keep], len(dup))


def load_gbmi(D):
    index = {rs: j for j, rs in enumerate(D["rsid"])}
    J = D["rsid"].size
    Z = np.full((len(STRATA), J), np.nan)
    N = np.zeros((len(STRATA), J))
    ref0 = alt0 = None
    for s, (name, tok) in enumerate(STRATA.items()):
        path = os.path.join(GBMI_DIR, PATTERN.format(tok))
        t0 = time.time()
        rows, ref, alt, z, n, ndup = read_gbmi(path, index)
        if name == "ALL":
            ref0, alt0 = np.full(J, ""), np.full(J, "")
            ref0[rows], alt0[rows] = ref, alt
            sign = np.ones(rows.size)
        else:
            same = (ref0[rows] == ref) & (alt0[rows] == alt)
            swap = (ref0[rows] == alt) & (alt0[rows] == ref)
            ok = same | swap
            rows, z, n = rows[ok], z[ok], n[ok]
            sign = np.where(swap[ok], -1.0, 1.0)
        Z[s, rows], N[s, rows] = sign * z, n
        print(f"  [gbmi] {name}: {rows.size} variants at the 1000G regression set "
              f"({ndup} duplicated rsIDs dropped), median N {np.median(n):.0f} ({time.time() - t0:.0f}s)")
    return Z, N


# ------------------------------------------------------------- regressions

def block_stats(x, y, w, blocks):
    """Per-block WLS sufficient statistics for y ~ a + b x: (B, 5)."""
    cols = np.stack([w, w * x, w * y, w * x * x, w * x * y], axis=1)
    return np.add.reduceat(cols, blocks, axis=0)


def slope_full_and_deleted(S):
    """Slope and intercept from the full sums and from each delete-one-block sum."""
    def solve(T):
        sw, sx, sy, sxx, sxy = T[..., 0], T[..., 1], T[..., 2], T[..., 3], T[..., 4]
        b = (sw * sxy - sx * sy) / (sw * sxx - sx * sx)
        return b, (sy - b * sx) / sw
    tot = S.sum(axis=0)
    return solve(tot), solve(tot[None, :] - S)


def jk_se(full, deleted):
    B = deleted.shape[0]
    return float(np.sqrt((B - 1) / B * np.sum((deleted - deleted.mean()) ** 2)))


def block_starts(n, n_blocks=N_BLOCKS):
    return np.linspace(0, n, min(n_blocks, n) + 1).astype(int)[:-1]


def chi2_ok(Z, N, s, sel):
    nbar = np.median(N[s][sel]) if sel.any() else 0.0
    return np.nan_to_num(Z[s] ** 2, nan=np.inf) <= max(80.0, 0.001 * nbar)


def liability(h2_obs50, K):
    t = NormalDist().inv_cdf(1.0 - K)
    phi = NormalDist().pdf(t)
    return h2_obs50 * K ** 2 * (1 - K) ** 2 / (0.25 * phi ** 2)


def main():
    t0 = time.time()
    for name, tok in STRATA.items():
        path = os.path.join(GBMI_DIR, PATTERN.format(tok))
        got = sha256(path)
        if SHA256 and got != SHA256.get(name):
            print(f"  [warn] {name} sha256 {got[:12]}... differs from the recorded download")
    D = load_scores()
    M = D["M"]
    base = GK.regression_set(D)
    print(f"1000G scores: {D['rsid'].size} regression variants, M = {M}, "
          f"{base.sum()} with EUR MAF >= {GK.MAF_REGRESSION} outside the MHC")
    Z, N = load_gbmi(D)
    S = D["scores"]
    names = list(STRATA)                     # ALL EUR AFR EAS SAS AMR
    panel = {name: ANC.index(name) for name in names[1:]}

    def present(*ss):
        m = base.copy()
        for s in ss:
            si = names.index(s)
            m &= np.isfinite(Z[si]) & (N[si] > 0)
            if s != "ALL":
                m &= D["available"][panel[s]]
        for s in ss:
            m &= chi2_ok(Z, N, names.index(s), m)
        return m

    def h2_and_cov(pairs, sel):
        """Joint jackknife of h2 per stratum and g per pair on one selection."""
        idx = np.flatnonzero(sel)
        bl = block_starts(idx.size)
        out = {}
        for key in pairs:
            a, b = key
            ia, ib = names.index(a), names.index(b)
            pa, pb = panel.get(a, E), panel.get(b, E)
            ell = S[pa, pb][idx]
            if a == b:
                x, y, w = N[ia][idx] * ell, Z[ia][idx] ** 2, 1.0 / np.maximum(ell, 1.0)
            else:
                x = np.sqrt(N[ia][idx] * N[ib][idx]) * ell
                y = Z[ia][idx] * Z[ib][idx]
                w = 1.0 / np.maximum(np.sqrt(np.maximum(S[pa, pa][idx], 1.0) * np.maximum(S[pb, pb][idx], 1.0)), 1.0)
            (b_f, a_f), (b_d, _) = slope_full_and_deleted(block_stats(x, y, w, bl))
            out[key] = (b_f * M, b_d * M, a_f)
        return out

    # Table C1. Per-stratum h2.
    print("\nTable C1. Per-stratum h2 (observed scale at P = 0.5; liability at GBMI's K for reference)")
    K = N_CASES / N_TOTAL
    h2z, h2_of = {}, {}
    for s in names:
        sel = present(s)
        r = h2_and_cov([(s, s)], sel)[(s, s)]
        se = jk_se(r[0], r[1])
        h2z[s] = r[0] / se
        h2_of[s] = (r[0], se, np.median(N[names.index(s)][sel]))
        print(f"  {s:>4s}: {sel.sum():7d} variants, median N {np.median(N[names.index(s)][sel]):8.0f}, "
              f"mean chi2 {np.mean(Z[names.index(s)][sel] ** 2):.3f}, h2 {r[0]:.4f} (SE {se:.4f}, z {h2z[s]:.1f}), "
              f"intercept {r[2]:.3f}, liability {liability(r[0], K):.3f}")

    # Table C2. rho for every pair of ancestry strata.
    print(f"\nTable C2. rho_ab (SE); used in MA2 when both strata have h2 z > {Z_MIN:g}")
    A = len(ANC)
    rho_hat = np.eye(A)
    usable = np.eye(A, dtype=bool)
    for i in range(A):
        for k in range(i + 1, A):
            a, b = ANC[i], ANC[k]
            sel = present(a, b)
            if sel.sum() < 10 * N_BLOCKS:
                print(f"  {a}-{b}: {sel.sum()} shared variants, not estimated")
                continue
            r = h2_and_cov([(a, a), (b, b), (a, b)], sel)
            with np.errstate(invalid="ignore"):
                full = r[(a, b)][0] / np.sqrt(r[(a, a)][0] * r[(b, b)][0])
                dele = r[(a, b)][1] / np.sqrt(r[(a, a)][1] * r[(b, b)][1])
            se = jk_se(full, dele) if np.all(np.isfinite(dele)) else np.nan
            ok = h2z[a] > Z_MIN and h2z[b] > Z_MIN and np.isfinite(full)
            rho_hat[i, k] = rho_hat[k, i] = np.clip(full, -1.0, 1.0) if ok else np.nan
            usable[i, k] = usable[k, i] = ok
            print(f"  {a}-{b}: rho {full:.3f} (SE {se:.3f}) on {sel.sum()} variants"
                  f"{'' if ok else '  [bracketed]'}")

    # Table C3. The casualty.
    sel = present("ALL", "EUR")
    idx = np.flatnonzero(sel)
    iA = names.index("ALL")
    strata_n = N[1:, idx]                                       # stratum order = names[1:]
    order = [names[1:].index(a) for a in ANC]                   # -> score-cache order
    n_anc = strata_n[order]
    avail = (n_anc > 0) & D["available"][:, idx]
    cover = n_anc.sum(axis=0) / N[iA][idx]
    print(f"\nTable C3. The casualty on the ALL meta: {idx.size} variants; strata N / ALL N: "
          f"median {np.median(cover):.3f} (IQR {np.quantile(cover, 0.25):.3f}-{np.quantile(cover, 0.75):.3f})")
    sub = {k: S[..., idx] for k in ("scores",)}
    var = D["variances"][:, idx]
    om, _ = meta_weights(n_anc, avail, variances=var, scheme="ivw")
    lE = S[E, E][idx]
    w = 1.0 / np.maximum(lE, 1.0)
    bl = block_starts(idx.size)
    y = Z[iA][idx] ** 2

    def casualty(rho, n_j, label):
        L = meta_ld_scores(sub["scores"], om, rho=rho)
        (bE, aE), (bEd, _) = slope_full_and_deleted(block_stats(n_j * lE, y, w, bl))
        (bL, aL), (bLd, _) = slope_full_and_deleted(block_stats(n_j * L, y, w, bl))
        hE, hL = bE * M, bL * M
        seE, seL = jk_se(hE, bEd * M), jk_se(hL, bLd * M)
        se_diff = jk_se(hL - hE, (bLd - bEd) * M)
        ratio, ratio_d = hE / hL, bEd / bLd
        kap = single_reference_bias(L, lE, sample_sizes=n_j, weights=w).kappa
        print(f"  {label}: h2_EUR {hE:.4f} (SE {seE:.4f}, int {aE:.3f}) | h2_meta {hL:.4f} "
              f"(SE {seL:.4f}, int {aL:.3f}) | diff {hL - hE:+.4f} = {(hL - hE) / seE:+.1f} SE(h2_EUR) "
              f"(paired SE {se_diff:.4f}) | ratio {ratio:.3f} (SE {jk_se(ratio, ratio_d):.3f}) vs kappa {kap:.3f}")
        return hE, hL, seE, ratio, jk_se(ratio, ratio_d), kap

    n_all = N[iA][idx]
    if usable.all():
        primary = casualty(rho_hat, n_all, "rho estimated")
        brackets = []
    else:
        primary = None
        brackets = []
        for r in RHO_BRACKET:
            R = np.where(usable, rho_hat, r)
            brackets.append((r, casualty(R, n_all, f"rho estimated where usable, else {r}")))
    for r in RHO_BRACKET:
        R = np.full((A, A), r)
        np.fill_diagonal(R, 1.0)
        casualty(R, n_all, f"all pairs rho = {r} (bracket)")
    R_used = rho_hat if usable.all() else np.where(usable, rho_hat, 0.8)
    casualty(R_used, n_anc.sum(axis=0), "N_j = sum of strata N (sensitivity)")

    # Table C4 was added after the first run. It is post hoc and descriptive,
    # not a rule.
    L_used = meta_ld_scores(sub["scores"], om, rho=R_used)
    f = np.array([h2_of[a][2] for a in ANC])
    f = f / f.sum()
    mix = float(sum(fa * h2_of[a][0] for fa, a in zip(f, ANC)))
    mix_se = float(np.sqrt(sum((fa * h2_of[a][1]) ** 2 for fa, a in zip(f, ANC))))
    print("\nTable C4 (post hoc). corr(N l_EUR, N L) on the casualty set = "
          f"{np.corrcoef(n_all * lE, n_all * L_used)[0, 1]:.4f}, which makes C2 nearly "
          "mechanical. The composition-weighted mean of the per-stratum h2, each on its own "
          "ancestry's score (f = " + ", ".join(f"{a} {fa:.3f}" for a, fa in zip(ANC, f))
          + f"), is {mix:.4f} (SE {mix_se:.4f}, strata treated as independent)")

    print("\nVerdicts")
    cells = [("rho estimated", primary)] if primary else [(f"bracket {r}", c) for r, c in brackets]
    for label, (hE, hL, seE, ratio, se_ratio, kap) in cells:
        c1 = abs(hL - hE) > 2 * seE
        c2 = abs(ratio - kap) <= 2 * se_ratio
        print(f"  {label}: C1 {'CASUALTY' if c1 else 'practice note'} (|diff| = {abs(hL - hE) / seE:.1f} SE); "
              f"C2 observed ratio {ratio:.3f} vs kappa {kap:.3f} -> {'consistent' if c2 else 'NOT consistent'}")
        print(f"     liability scale (K = {K:.4f}): h2_EUR {liability(hE, K):.3f}, h2_meta {liability(hL, K):.3f}; "
              f"GBMI reported {H2_GBMI}")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
