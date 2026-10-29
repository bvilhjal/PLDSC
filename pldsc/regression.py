"""LD Score regression with a provisional block-jackknife helper (milestone M5).

Weighted least-squares regression of a per-variant statistic (chi^2 or z^2) on
an LD score, with a delete-a-block jackknife for provisional standard errors.
Kept separate from score construction: it consumes an LD score and a statistic
and knows nothing about how either was built.

The intercept SE reported here is the *regression* SE at a fixed operator. When
the LD score / statistic come from a fitted V_c, that operator is itself
uncertain; ``experiments/m5_inference.py`` illustrates that omitted uncertainty
but does not provide a validated propagation correction.
The caller must supply chromosome/position-sorted observations and choose blocks
that cover the dependence range. Two-step weights are estimated once from the
full data and held fixed across replicates, so this helper is not a complete
genomic block-jackknife implementation; see ``docs/METHOD.md`` Section 9.1.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Float = NDArray[np.float64]


@dataclass
class RegressionResult:
    intercept: float
    slope: float
    intercept_se: float
    slope_se: float
    n_blocks: int
    weighting: str


def _wls_beta(x: Float, y: Float, w: Float) -> Float:
    """WLS of ``y ~ 1 + x`` with weights ``w``; returns ``[intercept, slope]``."""
    # Center after shifting to avoid cancellation for nearly constant scores.
    x0, y0 = x[np.argmin(np.abs(x))], y[np.argmin(np.abs(y))]
    dx, dy = x - x0, y - y0
    sw = np.sum(w)
    mx, my = np.dot(w, dx) / sw, np.dot(w, dy) / sw
    dx, dy = dx - mx, dy - my
    wx = w * dx
    var = np.dot(wx, dx)
    if var <= 0.0 or not np.isfinite(var):
        raise ValueError("singular weighted design (constant ldscore, or a jackknife "
                         "replicate whose kept points carry no weight)")
    slope = np.dot(wx, dy) / var
    return np.array([y0 + my - slope * (x0 + mx), slope])


def ld_score_regression(
    chi2: Float,
    ldscore: Float,
    *,
    weights: Float | None = None,
    two_step: bool = False,
    n_blocks: int = 200,
) -> RegressionResult:
    """Regress ``chi2`` on ``ldscore`` and return estimates with jackknife SEs.

    Parameters
    ----------
    chi2, ldscore
        1-D arrays of equal length (per-variant statistic and LD score).
    weights
        Optional per-variant regression weights. ``None`` is ordinary LS.
    two_step
        If true, apply a *simplified* LDSC-style heteroscedasticity weighting: a
        first pass predicts ``E[chi2]`` and reweights by ``1 / (2 E[chi2]^2)``
        (clipped at ``E[chi2] >= 1``). Note this is the heteroscedasticity term
        only; full LDSC additionally down-weights redundant correlated outcomes
        using an LD-overlap approximation (Bulik-Sullivan 2015). That factor is
        not automatically valid for a changed metric. Composes with ``weights``
        if both are given. The preliminary weights are not recomputed inside
        delete-block replicates.
    n_blocks
        Number of contiguous array-order jackknife blocks (capped at the number
        of points). Input genomic order is assumed and is not checked.

    Returns
    -------
    RegressionResult
        Point estimates and delete-a-block jackknife standard errors.
    """
    chi2 = np.asarray(chi2, dtype=np.float64)
    ld = np.asarray(ldscore, dtype=np.float64)
    if chi2.ndim != 1 or chi2.shape != ld.shape:
        raise ValueError("chi2 and ldscore must be 1-D arrays of equal length")
    if not np.all(np.isfinite(chi2)) or not np.all(np.isfinite(ld)):
        raise ValueError("chi2 and ldscore must contain only finite values")
    n = chi2.size
    if n < 4:
        raise ValueError("need at least 4 points for a block jackknife")

    if weights is None:
        w = np.ones(n)
        tag = "ols"
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape != chi2.shape:
            raise ValueError("weights must match chi2 shape")
        if not np.all(np.isfinite(w)):
            raise ValueError("weights must contain only finite values")
        if np.any(w <= 0.0):
            raise ValueError("weights must be strictly positive")
        tag = "user"
    if two_step:
        b0 = _wls_beta(ld, chi2, w)
        pred = np.maximum(b0[0] + b0[1] * ld, 1.0)
        w = w / (2.0 * pred ** 2)
        tag = "two-step" if tag == "ols" else tag + "+two-step"

    beta = _wls_beta(ld, chi2, w)

    nb = int(min(n_blocks, n))
    blocks = np.array_split(np.arange(n), nb)
    ests = np.empty((nb, 2))
    for i, blk in enumerate(blocks):
        keep = np.ones(n, dtype=bool)
        keep[blk] = False
        ests[i] = _wls_beta(ld[keep], chi2[keep], w[keep])
    mean = ests.mean(axis=0)
    se = np.sqrt((nb - 1) / nb * np.sum((ests - mean) ** 2, axis=0))

    return RegressionResult(
        intercept=float(beta[0]),
        slope=float(beta[1]),
        intercept_se=float(se[0]),
        slope_se=float(se[1]),
        n_blocks=nb,
        weighting=tag,
    )
