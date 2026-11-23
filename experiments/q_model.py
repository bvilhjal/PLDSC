"""Moment model of a pooled statistic and its Cochran's Q, shared by the Q drivers.

A two-stratum fixed-effect meta-analysis reports Z = cE zE + cF zF and
Q = D^2 with D = cF zE - cE zF, where cE^2 + cF^2 = 1 per SNP. With stratum
moments

    E[zE^2] = aE + xE hE,  E[zF^2] = aF + xF hF,  E[zE zF] = aEF + xX g,

the rotation gives

    E[Z^2] = cE^2 aE + cF^2 aF + 2 cE cF aEF + cE^2 xE hE + cF^2 xF hF + 2 cE cF xX g,
    E[Q]   = cF^2 aE + cE^2 aF - 2 cE cF aEF + cF^2 xE hE + cE^2 xF hF - 2 cE cF xX g.

``"three"`` fits the three intercepts aE, aF and aEF through the columns
(cE^2, cF^2, 2 cE cF) and (cF^2, cE^2, -2 cE cF). ``"shared"`` is the form
recorded on 2026-11-20 (one free constant per equation). The two agree when
the weights are proportional across SNPs (sample-size weights), because the
intercept parts of E[Z^2] and E[Q] are then constants. Under inverse-variance
weights c varies with allele frequency, and unequal stratum intercepts leak
into the slopes of the shared form. On the chr22 example bundle, with
EUR/EAS inverse-variance weights, aE = 1.076 and aF = 0.983, the shared form
biases g by -9.4% in expectation and the three-intercept form is exact
(tests/test_q_heterogeneity.py).

With constant c the three intercept columns span only two dimensions, so
:func:`zq_design` returns the shared form, which then spans the same space.
"""
from __future__ import annotations

import numpy as np

INTERCEPT_FORMS = ("three", "shared")


def weights_vary(cE, rtol=1e-10):
    """True when the stratum weight share cE^2 varies across SNPs."""
    c2 = np.asarray(cE, dtype=float) ** 2
    return bool(np.ptp(c2) > rtol * max(1.0, float(np.max(np.abs(c2)))))


def zq_design(cE, cF, xE, xF, xX, intercepts="three"):
    """Design matrices for the stacked Z^2 and Q regressions.

    Returns ``(XZ, XQ, form)``. The last three columns are always the
    slopes (hE, hF, g); the leading columns are intercepts. ``form`` is the
    intercept form actually used: ``"three"`` falls back to ``"shared"``
    when the weights do not vary, since only two intercept combinations are
    then identified.
    """
    if intercepts not in INTERCEPT_FORMS:
        raise ValueError(f"intercepts must be one of {INTERCEPT_FORMS}")
    cE, cF = np.broadcast_arrays(np.asarray(cE, dtype=float), np.asarray(cF, dtype=float))
    n = np.broadcast_shapes(cE.shape, np.shape(xE), np.shape(xF), np.shape(xX))[0]
    cE, cF = np.broadcast_to(cE, (n,)), np.broadcast_to(cF, (n,))
    sZ = (cE ** 2 * xE, cF ** 2 * xF, 2 * cE * cF * xX)
    sQ = (cF ** 2 * xE, cE ** 2 * xF, -2 * cE * cF * xX)
    if intercepts == "three" and weights_vary(cE):
        iZ = (cE ** 2, cF ** 2, 2 * cE * cF)
        iQ = (cF ** 2, cE ** 2, -2 * cE * cF)
        form = "three"
    else:
        one, zero = np.ones(n), np.zeros(n)
        iZ, iQ = (one, zero), (zero, one)
        form = "shared"
    return np.column_stack(iZ + sZ), np.column_stack(iQ + sQ), form


def zq_moments(cE, cF, xE, xF, xX, hE, hF, g, aE=1.0, aF=1.0, aEF=0.0):
    """Exact E[Z^2] and E[Q] under the stratum moment model above."""
    vE, vF, vX = aE + xE * hE, aF + xF * hF, aEF + xX * g
    return (cE ** 2 * vE + cF ** 2 * vF + 2 * cE * cF * vX,
            cF ** 2 * vE + cE ** 2 * vF - 2 * cE * cF * vX)
