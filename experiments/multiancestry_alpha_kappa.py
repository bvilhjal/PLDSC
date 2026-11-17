"""Does the kappa headline depend on effect model A1? (MULTI_ANCESTRY section 7)

MULTI_ANCESTRY section 7 lists the effect model as an unexplored caveat: A1
gives every variant the same standardized effect variance within each
ancestry, and alpha != -1 or eta != 1 would change kappa. Skoglund & Mathieson
(Nat Rev Genet 27:935, 2026) summarise the case that stabilizing selection
shapes complex-trait architecture; its standard parameterisation (Schoech et
al., Nat Commun 10:790, 2019) gives per-allele effects variance proportional to
[p(1 - p)]^alpha, with a best-fit mean alpha of -0.38 (SE 0.02) over 25 UK
Biobank traits.

Effect models. A1 is `effect_scales=None`. Model S(alpha): per-allele
effects are shared across ancestries (correlation rho on the per-allele scale)
with variance proportional to [pbar(1 - pbar)]^alpha, where pbar is the pooled
1000G frequency. The standardized effect in ancestry a then has standard
deviation v_ak = sqrt(2 p_ak q_ak) [pbar_k qbar_k]^(alpha/2), and the scores
are sum_k v_ak v_bk r^a_jk r^b_jk with eta = 1 (ancestry heritability
differences live in v). v is divided by one constant so that the EUR mean of
v^2 over the reference variants is 1, putting every arm on the EUR
heritability scale. S(-1) isolates "shared per-allele effects" from "standardized
within each ancestry"; S(-0.38) adds the frequency dependence.

The analyst regresses on the ordinary (A1) EUR score in every arm. kappa_meta
is the slope ratio for the meta-analysis; kappa_EUR is the same quantity for an
EUR-only GWAS under the same arm, which carries the familiar single-ancestry
alpha bias and has nothing to do with multiple ancestries. The gated quantity
is the multi-ancestry penalty pi = kappa_meta / kappa_EUR, which does not
depend on the normalisation of v. Under A1, kappa_EUR = 1 and pi is the
existing kappa.

Data, map, filters, windows, regression set and fits are those of
multiancestry_genome_kappa.py (bigsnpr 1000G, SHAPEIT4 GRCh37, +/-1 cM,
reference MAF >= 0.01, regression MAF >= 0.05, GWAS MAF >= 0.01, MHC dropped
genome-wide, OLS and LDSC-style weights).

Pre-stated rules (2026-11-17, written before the first run):

  R0 (implementation). On chr22, the A1 arm's scores equal the cached
     kappa_gw/chr22.npz scores to within 1e-10. Otherwise stop; nothing below
     is interpreted.

  R1 (the decision). Rule cell: EUR 0.75 + AFR, rho = 0.8, sample-size
     weights, N_j-scaled regressors, regression SNPs with full-panel and EUR
     MAF >= 0.05. If |pi(S(alpha)) - pi(A1)| <= 0.03 for both alpha in
     {-1, -0.38} under both OLS and LDSC-style weights, the chr22 headline is
     A1-robust: section 7's caveat is closed for chr22 with these numbers and no
     genome-wide run is needed. If any difference exceeds 0.03, the headline is
     reported as conditional on A1, a sensitivity row is added, and the
     genome-wide S(-0.38) run (`--chroms 1 ... 22`) is required before any
     headline number changes. 0.03 is the G2 bridge tolerance (twice the
     density effect in MULTI_ANCESTRY Table 3); it is not revised after the run.

Everything else printed is descriptive: the implied eta per ancestry, the
other ancestries, the rho grid, kappa_EUR and kappa_meta separately.

Run:  python experiments/multiancestry_alpha_kappa.py [--chroms 22 ...]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from pldsc.multiancestry import cross_ancestry_ld_scores, meta_ld_scores, meta_weights  # noqa: E402
import multiancestry_genome_kappa as GK  # noqa: E402

ANC = GK.ANCESTRIES
E = GK.E
ALPHAS = (-1.0, -0.38)
RHOS = (1.0, 0.8, 0.6)
TOL = 0.03
R0_TOL = 1e-10


def effect_scales(panels, p_full, alpha):
    """v_ak = sqrt(2 p_ak q_ak) (pbar_k qbar_k)^(alpha/2), unnormalised."""
    shared = (p_full * (1.0 - p_full)) ** (alpha / 2.0)
    out = []
    for P in panels:
        p = P.mean(0) / 2.0
        out.append(np.sqrt(2.0 * p * (1.0 - p)) * shared)
    return out


def score(panels, cm, maf_full, scales=None, segment=GK.SEGMENT):
    """GK.score_chromosome with optional effect scales sliced per segment."""
    reg = np.flatnonzero(maf_full >= GK.MAF_REGRESSION)
    parts = []
    for s in range(0, reg.size, segment):
        chunk = reg[s:s + segment]
        lo = int(np.searchsorted(cm, cm[chunk[0]] - GK.WINDOW_CM, side="left"))
        hi = int(np.searchsorted(cm, cm[chunk[-1]] + GK.WINDOW_CM, side="right"))
        sc = None if scales is None else [v[lo:hi] for v in scales]
        parts.append(cross_ancestry_ld_scores([P[:, lo:hi] for P in panels], cm[lo:hi],
                                              window_cm=GK.WINDOW_CM, regression=chunk - lo,
                                              effect_scales=sc))
    maf = np.stack([np.minimum(P[:, reg].mean(0) / 2.0, 1.0 - P[:, reg].mean(0) / 2.0)
                    for P in panels])
    return dict(scores=np.concatenate([p.scores for p in parts], axis=2),
                available=np.concatenate([p.available for p in parts], axis=1),
                variances=np.concatenate([p.variances for p in parts], axis=1),
                maf=maf, reg=reg)


def chromosome(c, bim, sp):
    chrom, _, pos = bim
    rows = np.flatnonzero(chrom == c)
    start, stop = int(rows[0]), int(rows[-1]) + 1
    G = GK.read_bed(GK.BFILE, sp.size, chrom.size, start, stop)
    p = G.mean(0) / 2.0
    keep = np.minimum(p, 1.0 - p) >= GK.MAF_REF
    G, p = G[:, keep], p[keep]
    cm = GK.chrom_cm(c, pos[start:stop][keep])
    panels = [G[sp == a] for a in ANC]
    del G
    maf_full = np.minimum(p, 1.0 - p)
    arms = {"A1": score(panels, cm, maf_full)}
    sumsq = {}
    for alpha in ALPHAS:
        v = effect_scales(panels, p, alpha)
        sumsq[alpha] = np.array([float(np.sum(x ** 2)) for x in v])
        arms[alpha] = score(panels, cm, maf_full, scales=v)
    reg = arms["A1"]["reg"]
    return dict(arms=arms, sumsq=sumsq, n_reference=int(keep.sum()),
                pos=pos[start:stop][keep][reg], chrom=np.full(reg.size, c))


def pooled(per_chrom):
    """Concatenate chromosomes; normalise each S arm by the genome EUR mean of v^2."""
    m = sum(d["n_reference"] for d in per_chrom)
    out = {}
    for arm in ("A1", *ALPHAS):
        D = {k: np.concatenate([d["arms"][arm][k] for d in per_chrom], axis=-1)
             for k in ("scores", "available", "variances", "maf")}
        D["pos"] = np.concatenate([d["pos"] for d in per_chrom])
        D["chrom"] = np.concatenate([d["chrom"] for d in per_chrom])
        D["gwas_ok"] = (D["maf"] >= GK.MAF_GWAS) & D["available"]
        if arm != "A1":
            ss = sum(d["sumsq"][arm] for d in per_chrom)      # (A,) sum of v_a^2
            D["scores"] = D["scores"] / (ss[E] / m)             # EUR mean v^2 -> 1
            D["eta"] = ss / ss[E]
        else:
            D["eta"] = np.ones(len(ANC))
        out[arm] = D
    return out


def meta_L(D, X, fE, rho):
    A = len(ANC)
    n = np.zeros(A)
    n[E], n[X] = fE, 1.0 - fE
    om, n_tot = meta_weights(n, D["gwas_ok"], variances=D["variances"], scheme="sample_size")
    R = np.full((A, A), rho)
    np.fill_diagonal(R, 1.0)
    return meta_ld_scores(D["scores"], om, rho=R), n_tot


def penalty(arms, arm, sel, X, fE=0.75, rho=0.8):
    """(pi, kappa_meta, kappa_EUR), each as (OLS, weighted); regressor = A1 EUR score."""
    D, lR = arms[arm], arms["A1"]["scores"][E, E][sel]
    L, n_tot = meta_L(D, X, fE, rho)
    km = GK.fit(L[sel], lR, n_tot[sel])[:2]
    ke = GK.fit(D["scores"][E, E][sel], lR, np.ones(lR.size))[:2]
    return tuple(a / b for a, b in zip(km, ke)), km, ke


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chroms", type=int, nargs="+", default=[22])
    args = ap.parse_args()
    t0 = time.time()

    got = {k: GK.sha256(f"{GK.BFILE}.{k}") for k in GK.CHECKSUMS}
    for k, v in got.items():
        if v != GK.CHECKSUMS[k]:
            print(f"  [warn] {k} sha256 {v[:12]}... differs from the recorded copy")
    bim = GK.read_bim(GK.BFILE)
    sp = GK.read_superpops(GK.BFILE)
    per = []
    for c in args.chroms:
        t = time.time()
        per.append(chromosome(c, bim, sp))
        print(f"  chr{c}: {per[-1]['n_reference']} reference, {per[-1]['pos'].size} regression "
              f"variants, three arms scored ({time.time() - t:.0f}s)")

    if 22 in args.chroms:
        cache = os.path.join(GK.OUT, "chr22.npz")
        if os.path.exists(cache):
            ref = np.load(cache, allow_pickle=False)["scores"]
            mine = per[args.chroms.index(22)]["arms"]["A1"]["scores"]
            diff = float(np.abs(ref - mine).max()) if ref.shape == mine.shape else np.inf
            ok = diff <= R0_TOL
            print(f"\nR0: chr22 A1 scores vs cached kappa_gw scores, max |diff| = {diff:.2e} "
                  f"-> {'PASS' if ok else 'FAIL'}")
            if not ok:
                print("R0 FAILED: stopping; nothing below would be interpretable.")
                return
        else:
            print(f"\nR0: no cache at {cache}; the control cannot be checked.")

    arms = pooled(per)
    sel = GK.regression_set(arms["A1"])
    print(f"\n{sel.sum()} rule regression SNPs on chr {', '.join(map(str, args.chroms))}")

    print("\nTable AK1. Implied eta_a = h_a^2 / h_EUR^2 (mean v^2 relative to EUR); A1 is 1")
    print("        " + " ".join(f"{a:>6s}" for a in ANC))
    for alpha in ALPHAS:
        print(f"  S({alpha:+.2f}) " + " ".join(f"{e:6.3f}" for e in arms[alpha]["eta"]))

    X = ANC.index("AFR")
    base, _, _ = penalty(arms, "A1", sel, X)
    print("\nTable AK2 (R1). Rule cell EUR 0.75 + AFR, rho = 0.8: pi = kappa_meta / kappa_EUR, "
          "OLS / weighted")
    print(f"  {'A1':>8s}: {base[0]:.3f} / {base[1]:.3f}")
    worst = 0.0
    for alpha in ALPHAS:
        pi, _, _ = penalty(arms, alpha, sel, X)
        d = (pi[0] - base[0], pi[1] - base[1])
        worst = max(worst, abs(d[0]), abs(d[1]))
        print(f"  S({alpha:+.2f}): {pi[0]:.3f} / {pi[1]:.3f}   (difference from A1 "
              f"{d[0]:+.3f} / {d[1]:+.3f})")
    verdict = "A1-ROBUST" if worst <= TOL else "A1-CONDITIONAL"
    print(f"  R1: largest |difference| {worst:.3f} vs tolerance {TOL} -> {verdict}")

    print("\nTable AK3. pi for every ancestry and rho, EUR 0.75 + X, OLS / weighted (descriptive)")
    print(f"  {'arm':>8s} {'X':>4s} " + " ".join(f"{'rho=' + str(r):>15s}" for r in RHOS))
    for arm in ("A1", *ALPHAS):
        label = arm if arm == "A1" else f"S({arm:+.2f})"
        for x in range(1, len(ANC)):
            cells = [penalty(arms, arm, sel, x, rho=r)[0] for r in RHOS]
            print(f"  {label:>8s} {ANC[x]:>4s} " + " ".join(f"{c[0]:.3f} / {c[1]:.3f}".rjust(15)
                                                          for c in cells))

    print("\nTable AK4. Rule cell decomposed: kappa_meta and kappa_EUR, OLS / weighted (descriptive)")
    for arm in ("A1", *ALPHAS):
        label = arm if arm == "A1" else f"S({arm:+.2f})"
        _, km, ke = penalty(arms, arm, sel, X)
        print(f"  {label:>8s}: kappa_meta {km[0]:.3f} / {km[1]:.3f}, kappa_EUR {ke[0]:.3f} / {ke[1]:.3f}")
    print(f"\ntotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
