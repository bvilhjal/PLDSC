"""LD Score regression with block-jackknife inference (milestone M5).

Weighted least-squares regression of a per-variant statistic (chi^2 or z^2) on
an LD score, with a delete-a-block jackknife for standard errors -- the standard
LDSC inference scheme. Kept separate from score construction: it consumes an LD
score and a statistic and knows nothing about how either was built.

The intercept SE reported here is the *regression* SE at a fixed operator. When
the LD score / statistic come from a fitted V_c, that operator is itself
uncertain; see ``experiments/m5_inference.py`` for propagating variance-
component uncertainty into the intercept SE (the requirement M4 surfaced).
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
    X = np.column_stack((np.ones_like(x), x))
    XtW = X.T * w
    return np.linalg.solve(XtW @ X, XtW @ y)


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
        only; full LDSC additionally down-weights by the sampling variance of the
        LD score itself (Bulik-Sullivan 2015). Composes with ``weights`` if both
        are given.
    n_blocks
        Number of contiguous jackknife blocks (capped at the number of points).

    Returns
    -------
    RegressionResult
        Point estimates and delete-a-block jackknife standard errors.
    """
    chi2 = np.asarray(chi2, dtype=np.float64)
    ld = np.asarray(ldscore, dtype=np.float64)
    if chi2.ndim != 1 or chi2.shape != ld.shape:
        raise ValueError("chi2 and ldscore must be 1-D arrays of equal length")
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
