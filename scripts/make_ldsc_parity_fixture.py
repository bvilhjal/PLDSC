#!/usr/bin/env python3
"""Write tests/data/ldsc_parity.npz: inputs and reference ldsc estimates.

The expected values come from the estimator code of bulik/ldsc (commit
2fdeeb3: ldscore/regressions.py, irwls.py, jackknife.py). That code is
Python 2; the only changes needed to import it under Python 3 are
``xrange`` -> ``range``, ``print x`` -> ``print(x)`` in irwls.py, and
relative imports (``from . import jackknife as jk``, ``from .irwls import
IRWLS``), and, for NumPy 2, ``float(np.ravel(v)[0])`` in place of
``float(v)`` for the three (1, 1) arrays at the end of ``RG.__init__``.
No arithmetic is changed. Usage:

    git clone https://github.com/bulik/ldsc /tmp/ldsc   # then apply the port above
    python scripts/make_ldsc_parity_fixture.py --ldsc /tmp/ldsc --out tests/data/ldsc_parity.npz

The inputs are synthetic: 3,000 SNPs with per-SNP N, over-counting scores
w_ld that differ from the regression scores, a few large chi-squares so the
two-step estimator drops rows, and correlated z scores for rg.

The cases map to ldsc's command line as follows. ``ldsc --h2`` runs the
``twostep`` h2 case. ``ldsc --rg`` runs the ``onestep`` rg case, because
``estimate_rg`` turns ``intercept_h2`` into a list before testing it against
None, so its default cutoff never applies; ``ldsc --rg --two-step 30`` runs
the ``twostep`` case. Checked through the ldsc CLI on these inputs on
2026-11-23: r_g 0.3924 (0.0900) and 0.3284 (0.0895).
"""
import argparse
import sys

import numpy as np


def inputs(seed=20261123, n=3000):
    rng = np.random.default_rng(seed)
    ld = rng.gamma(2.0, 20.0, n) + 0.5
    w_ld = np.maximum(ld * rng.uniform(0.6, 1.0, n), 0.8)
    N1 = rng.uniform(40_000, 60_000, n)
    N2 = rng.uniform(20_000, 30_000, n)
    M = 2.0e5
    h1, h2, rho, ns_over = 0.3, 0.2, 0.5, 0.05
    v1 = 1.02 + N1 * h1 * ld / M
    v2 = 1.00 + N2 * h2 * ld / M
    cv = ns_over + np.sqrt(N1 * N2) * rho * np.sqrt(h1 * h2) * ld / M
    r = cv / np.sqrt(v1 * v2)
    a, b = rng.standard_normal((2, n))
    z1 = np.sqrt(v1) * a
    z2 = np.sqrt(v2) * (r * a + np.sqrt(1 - r * r) * b)
    big = rng.choice(n, 12, replace=False)
    z1[big] += rng.choice([-1, 1], 12) * rng.uniform(6, 9, 12)
    return dict(ld=ld, w_ld=w_ld, N1=N1, N2=N2, z1=z1, z2=z2, M=np.array(M))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ldsc", required=True, help="directory holding the Python 3 port of ldscore/")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    sys.path.insert(0, args.ldsc)
    from ldscore import regressions as reg

    d = inputs()
    n = d["ld"].size
    col = lambda v: np.asarray(v, dtype=float).reshape(n, 1)           # noqa: E731
    Mm = np.array([[float(d["M"])]])
    out = dict(d)
    for name, kw in (("twostep", dict(twostep=30)), ("onestep", dict()), ("fixed", dict(intercept=1.0))):
        h = reg.Hsq(col(d["z1"] ** 2), col(d["ld"]), col(d["w_ld"]), col(d["N1"]), Mm, n_blocks=200, **kw)
        isd = 0.0 if h.constrain_intercept else float(h.intercept_se)
        out[f"h2_{name}"] = np.array([float(h.tot), float(h.tot_se), float(h.intercept), isd])
    for name, kw in (("twostep", dict(twostep=30)), ("onestep", dict()), ("fixed", dict(intercept_gencov=0.0))):
        g = reg.RG(col(d["z1"]), col(d["z2"]), col(d["ld"]), col(d["w_ld"]), col(d["N1"]), col(d["N2"]), Mm,
                   n_blocks=200, **kw)
        gsd = 0.0 if g.gencov.constrain_intercept else float(g.gencov.intercept_se)
        out[f"rg_{name}"] = np.array([float(g.rg_ratio), float(g.rg_se), float(g.gencov.tot), float(g.gencov.intercept), gsd,
                                      float(g.hsq1.tot), float(g.hsq2.tot), float(g.gencov.tot_se),
                                      float(g.hsq1.tot_se), float(g.hsq2.tot_se)])
        print(name, out[f"h2_{name}"], out[f"rg_{name}"])
    np.savez_compressed(args.out, **out)


if __name__ == "__main__":
    main()
