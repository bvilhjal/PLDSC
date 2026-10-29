"""Allocation-free delete-a-block IRLS for cross-trait LDSC.

Each replicate refits both marginal regressions and the cross regression with
its own weights. The NumPy implementation in :mod:`pldsc.bivariate` remains the
oracle and handles poorly conditioned systems. No fastmath: changing a small
heritability's sign changes whether the jackknife standard error is defined.
"""

import numpy as np

from ._numba import _jit, _jit_parallel, prange


@_jit
def _solve_sums(sw, swx, swxx, swy, swxy):
    """Solve two normal equations, or request the reference solve with NaNs.

    Relative determinant and slope-numerator checks avoid cancellation near
    a constant predictor or zero slope. This is a routing threshold, not a
    singularity definition: the NumPy reference solve decides those cases.
    """
    denom = sw * swxx
    det = denom - swx * swx
    numerator = sw * swxy - swx * swy
    if (not np.isfinite(det) or not np.isfinite(numerator)
            or det <= 1e-8 * abs(denom)
            or abs(numerator) <= 1e-10 * (abs(sw * swxy) + abs(swx * swy))):
        return np.nan, np.nan
    slope = numerator / det
    intercept = (swxx * swy - swx * swxy) / det
    return slope, intercept


@_jit
def _marginal_fit(x, y, ell, lo, hi, n_iter):
    slope = intercept = 0.0
    for it in range(n_iter + 1):
        sw = swx = swxx = swy = swxy = 0.0
        for part in range(2):
            a, e = (0, lo) if part == 0 else (hi, len(x))
            for i in range(a, e):
                pred = 1.0 if it == 0 else max(intercept + slope * x[i], 1.0)
                w = (1.0 / (2.0 * pred * pred)) * (1.0 / ell[i])
                sw += w
                swx += w * x[i]
                swxx += w * x[i] * x[i]
                swy += w * y[i]
                swxy += w * x[i] * y[i]
        slope, intercept = _solve_sums(sw, swx, swxx, swy, swxy)
        if not np.isfinite(slope) or not np.isfinite(intercept):
            break
    return slope, intercept


def _jackknife_fits_impl(x1, x2, xc, y1, y2, yc, ell, starts, stops,
                         n_iter, constrained, c):
    """Return (h1, h2, gcov, intercept) per deletion; NaNs request fallback.

    Deletions are independent: each writes only its own row of ``out`` and
    every reduction stays serial within a replicate, so the parallel twin
    below is bit-identical to the serial one at any thread count.
    """
    out = np.full((len(starts), 4), np.nan)
    for b in prange(len(starts)):
        lo, hi = starts[b], stops[b]
        h1, i1 = _marginal_fit(x1, y1, ell, lo, hi, n_iter)
        h2, i2 = _marginal_fit(x2, y2, ell, lo, hi, n_iter)
        if not (np.isfinite(h1) and np.isfinite(i1)
                and np.isfinite(h2) and np.isfinite(i2)):
            continue
        slope = 0.0
        intercept = c if constrained else 0.0
        for _ in range(n_iter + 1):
            sw = swx = swxx = swy = swxy = 0.0
            for part in range(2):
                a, e = (0, lo) if part == 0 else (hi, len(xc))
                for i in range(a, e):
                    p1 = max(i1 + h1 * x1[i], 1.0)
                    p2 = max(i2 + h2 * x2[i], 1.0)
                    pc = intercept + slope * xc[i]
                    w = 1.0 / (ell[i] * max(p1 * p2 + pc * pc, 1e-6))
                    if constrained:
                        swxx += w * xc[i] * xc[i]
                        swxy += w * xc[i] * (yc[i] - c)
                    else:
                        sw += w
                        swx += w * xc[i]
                        swxx += w * xc[i] * xc[i]
                        swy += w * yc[i]
                        swxy += w * xc[i] * yc[i]
            if constrained:
                if swxx <= 0.0 or not np.isfinite(swxx):
                    slope = np.nan
                else:
                    slope = swxy / swxx
                intercept = c
            else:
                slope, intercept = _solve_sums(sw, swx, swxx, swy, swxy)
            if not np.isfinite(slope) or not np.isfinite(intercept):
                break
        out[b] = h1, h2, slope, intercept
    return out


# Serial and parallel twins of one source, as in bipred._bivar_kernels: the
# serial compile treats ``prange`` as ``range`` and never enters a parallel
# region, so an ``ncores=1`` caller pays no threading cost.
_jackknife_fits = _jit(_jackknife_fits_impl)
_jackknife_fits_par = _jit_parallel(_jackknife_fits_impl)
