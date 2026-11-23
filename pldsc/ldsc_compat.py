"""LD Score regression with the estimator of the reference ``ldsc`` software.

This is an opt-in alternative to the default one-step fits of
:func:`pldsc.ldsc_h2` and :func:`pldsc.ldsc_rg`. It follows
``bulik/ldsc`` (``ldscore/regressions.py``, ``irwls.py``, ``jackknife.py``,
commit 2fdeeb3) for a single LD-score annotation:

* initial weights from the aggregate estimator, then two weight updates;
* a separate over-counting score ``w_ld`` (``1 / max(w_ld, 1)``);
* heteroscedasticity weights with the heritability clipped to [0, 1] and the
  fitted intercept used without a floor;
* the block jackknife with the *final weights held fixed*
  (``LstsqJackknifeFast``), 200 blocks at ``floor(linspace(0, n, B + 1))``;
* the two-step estimator for free intercepts (intercept from rows with
  chi-square below 30, slope on all rows with that intercept fixed), with
  ldsc's adjustment of the step-two delete values;
* no high-chi-square cap; h2 and gcov reported from the full-data fit;
* rg with ldsc's ratio jackknife.

Two quirks are reproduced on purpose: step one of the two-step estimator
computes its weight updates from ``N * ell / mean(N)`` rather than ``ell``
(identical when N is constant), and the genetic covariance is clipped to
[-1, 1] inside the cross-trait weights. ``M`` must be given explicitly; the
published ldsc results use the MAF 5-50% count (``M_5_50``).

Parity is tested against the reference code's output in
``tests/test_ldsc_compat.py``.
"""
from __future__ import annotations

import numpy as np

from ._common import _as_n_vector, _finite_control, _integer_at_least
from .ldsc import LDSCResult

TWOSTEP_CUTOFF = 30.0


def _separators(n, n_blocks):
    return np.floor(np.linspace(0, n, n_blocks + 1)).astype(np.int64)


def _jackknife(x, y, sep):
    """Fixed-weight delete-a-block least squares (ldsc LstsqJackknifeFast).

    ``x`` and ``y`` are already multiplied by the square-root weights.
    Returns the full estimate, the delete values and the jackknife covariance.
    """
    nb = sep.size - 1
    p = x.shape[1]
    xtx = np.empty((nb, p, p))
    xty = np.empty((nb, p))
    for b in range(nb):
        xs, ys = x[sep[b]:sep[b + 1]], y[sep[b]:sep[b + 1]]
        xtx[b] = xs.T @ xs
        xty[b] = xs.T @ ys
    XtX, Xty = xtx.sum(0), xty.sum(0)
    est = np.linalg.solve(XtX, Xty)
    dele = np.linalg.solve(XtX[None] - xtx, (Xty[None] - xty)[..., None])[..., 0]
    return est, dele


def _pseudo_cov(est, dele):
    nb = dele.shape[0]
    pseudo = nb * est[None, :] - (nb - 1) * dele
    return np.atleast_2d(np.cov(pseudo.T, ddof=1)) / nb


def _irwls(x, y, update, w0, sep):
    """ldsc IRWLS: two weight updates, then the fixed-weight jackknife."""
    w = np.sqrt(w0)
    for _ in range(2):
        coef = np.linalg.lstsq(x * w[:, None], y * w, rcond=None)[0]
        w = np.sqrt(update(coef))
    return _jackknife(x * w[:, None], y * w, sep)


def _hsq_weights(ld, w_ld, N, M, hsq, intercept):
    hsq = min(max(hsq, 0.0), 1.0)
    ld = np.fmax(ld, 1.0)
    c = hsq * N / M
    return 1.0 / (2.0 * np.square(intercept + c * ld)) / np.fmax(w_ld, 1.0)


def _gencov_weights(ld, w_ld, N1, N2, M, h1, h2, rho_g, int_gc, int1, int2):
    h1, h2 = min(max(h1, 0.0), 1.0), min(max(h2, 0.0), 1.0)
    rho_g = max(min(rho_g, 1.0), -1.0)
    ld = np.fmax(ld, 1.0)
    a = N1 * h1 * ld / M + int1
    b = N2 * h2 * ld / M + int2
    c = np.sqrt(N1 * N2) * rho_g * ld / M + int_gc
    return 1.0 / (a * b + c * c) / np.fmax(w_ld, 1.0)


def _regression(y, ld, w_ld, N, M, n_blocks, intercept, step1, null_intercept, weights):
    """ldsc's LD_Score_Regression for one annotation.

    ``weights(ld_used, w_ld, rows, coef_slope_scaled, intercept)`` returns
    regression weights for the selected ``rows`` (a boolean mask or None).
    Returns (slope, intercept, slope delete values, intercept delete values,
    jackknife covariance of (slope, intercept), n step-two filtered), with
    the slope already converted to h2 or gcov units.
    """
    n = y.size
    agg = M * (np.mean(y) - (null_intercept if intercept is None else intercept)) / np.mean(ld * N)
    w0 = weights(ld, None, agg, intercept if intercept is not None else None)
    Nbar = float(np.mean(N))
    xs = N * ld / Nbar
    sep = _separators(n, n_blocks)
    if intercept is not None:
        yp = y - intercept
        X = xs[:, None]
        update = lambda a: weights(ld, None, M * a[0] / Nbar, intercept)   # noqa: E731
        est, dele = _irwls(X, yp, update, w0, sep)
        cov = _pseudo_cov(est, dele)
        return (M * est[0] / Nbar, float(intercept), M * dele[:, 0] / Nbar, None,
                float(np.sqrt(cov[0, 0])) * M / Nbar, 0.0, 0)
    X = np.column_stack((xs, np.ones(n)))
    if step1 is None:
        update = lambda a: weights(ld, None, M * a[0] / Nbar, a[1])          # noqa: E731
        est, dele = _irwls(X, y, update, w0, sep)
        cov = _pseudo_cov(est, dele)
        return (M * est[0] / Nbar, float(est[1]), M * dele[:, 0] / Nbar, dele[:, 1],
                float(np.sqrt(cov[0, 0])) * M / Nbar, float(np.sqrt(cov[1, 1])), 0)
    # Two-step estimator. Step one: free intercept on rows below the cutoff,
    # with ldsc's weight updates computed from N*ell/Nbar on those rows.
    n1 = int(step1.sum())
    X1 = X[step1]
    sep1 = _separators(n1, n_blocks)
    update1 = lambda a: weights(xs[step1], step1, M * a[0] / Nbar, a[1])     # noqa: E731
    est1, dele1 = _irwls(X1, y[step1], update1, w0[step1], sep1)
    int1 = float(est1[1])
    # Step two: slope on all rows with that intercept fixed, on blocks whose
    # boundaries are the step-one boundaries mapped back to all rows.
    rows1 = np.flatnonzero(step1)
    sep2 = np.concatenate(([0], rows1[sep1[1:-1]], [n]))
    update2 = lambda a: weights(ld, None, M * a[0] / Nbar, int1)             # noqa: E731
    est2, dele2 = _irwls(xs[:, None], y - int1, update2, w0, sep2)
    c = np.sum(w0 * xs) / np.sum(w0 * xs * xs)
    dele = np.column_stack((dele2[:, 0] - c * (dele1[:, 1] - int1), dele1[:, 1]))
    est = np.array([est2[0], int1])
    cov = _pseudo_cov(est, dele)
    return (M * est[0] / Nbar, int1, M * dele[:, 0] / Nbar, dele[:, 1],
            float(np.sqrt(cov[0, 0])) * M / Nbar, float(np.sqrt(cov[1, 1])), n - n1)


def _check(chisq_or_z, ld_scores, w_ld, name):
    y = np.asarray(chisq_or_z, dtype=float)
    ell = np.asarray(ld_scores, dtype=float)
    if y.ndim != 1 or ell.shape != y.shape or y.size < 2:
        raise ValueError(f"{name} and ld_scores must be one-dimensional with equal length >= 2")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(ell)) or np.any(ell <= 0.0):
        raise ValueError(f"{name} and ld_scores must be finite; ld_scores positive")
    wl = ell if w_ld is None else np.asarray(w_ld, dtype=float)
    if wl.shape != ell.shape or not np.all(np.isfinite(wl)):
        raise ValueError("w_ld must be finite with one value per row")
    return y, ell, wl


def _options(m_snps, n_blocks, twostep, constrain_intercept, m):
    if m_snps is None:
        raise ValueError("ldsc-compatible mode requires an explicit m_snps (ldsc uses M_5_50)")
    M = _finite_control("m_snps", m_snps, lower=1.0)
    nb = min(_integer_at_least("n_blocks", n_blocks, 2), m)
    if constrain_intercept is not None:
        constrain_intercept = _finite_control("constrain_intercept", constrain_intercept)
        twostep = None                     # ldsc: twostep is incompatible with a fixed intercept
    elif twostep is not None:
        twostep = _finite_control("twostep", twostep, lower=0.0)
    return M, nb, twostep, constrain_intercept


def ldsc_h2_compat(chisq, ld_scores, n_eff, *, m_snps, w_ld=None, n_blocks=200,
                   twostep=TWOSTEP_CUTOFF, constrain_intercept=None):
    """SNP heritability with the reference ldsc estimator (see module docstring).

    ``w_ld`` are the over-counting LD scores (ldsc ``--w-ld``); they default
    to ``ld_scores``, which is what ldsc does when the same files are passed
    to ``--ref-ld`` and ``--w-ld``. ``twostep=None`` disables the two-step
    estimator. ``m_snps`` is required.
    """
    y, ell, wl = _check(chisq, ld_scores, w_ld, "chisq")
    if np.any(y < 0.0):
        raise ValueError("chisq must be non-negative")
    N = _as_n_vector(n_eff, y.size)
    M, nb, twostep, fixed = _options(m_snps, n_blocks, twostep, constrain_intercept, y.size)

    def weights(ld, rows, hsq, intercept):
        Nr, wr = (N, wl) if rows is None else (N[rows], wl[rows])
        return _hsq_weights(ld, wr, Nr, M, hsq, 1.0 if intercept is None else intercept)

    step1 = None if twostep is None else y < twostep
    h2, intercept, _, _, h2_se, int_se, _ = _regression(y, ell, wl, N, M, nb, fixed, step1, 1.0, weights)
    mean_chisq = float(y.mean())
    ratio = (intercept - 1.0) / (mean_chisq - 1.0) if mean_chisq > 1.0 else float("nan")
    return LDSCResult(h2=float(h2), h2_se=h2_se, intercept=intercept, intercept_se=int_se,
                      mean_chisq=mean_chisq, ratio=ratio)


def _h2_parts(y, ell, wl, N, M, nb, fixed, twostep):
    def weights(ld, rows, hsq, intercept):
        Nr, wr = (N, wl) if rows is None else (N[rows], wl[rows])
        return _hsq_weights(ld, wr, Nr, M, hsq, 1.0 if intercept is None else intercept)
    step1 = None if twostep is None else y < twostep
    return _regression(y, ell, wl, N, M, nb, fixed, step1, 1.0, weights)


def ldsc_rg_compat(z1, z2, ld_scores, n_eff1, n_eff2, *, m_snps, w_ld=None, n_blocks=200,
                   twostep=TWOSTEP_CUTOFF, constrain_intercept=None):
    """Genetic correlation with the reference ldsc estimator (see module docstring).

    ``z1`` and ``z2`` are signed z scores aligned to the same allele and in
    genomic order. Both marginal intercepts are free; ``constrain_intercept``
    fixes the cross-trait intercept. ldsc raises an error when a fixed
    cross-trait intercept meets its default two-step estimator; here the
    cross-trait regression then runs without the two-step estimator and the
    marginal fits keep it. Pass ``twostep=None`` to match ``ldsc
    --intercept-gencov ... --two-step`` left unset for every fit.
    """
    from .bivariate import LDSCRgResult

    z1, ell, wl = _check(z1, ld_scores, w_ld, "z1")
    z2 = np.asarray(z2, dtype=float)
    if z2.shape != z1.shape or not np.all(np.isfinite(z2)):
        raise ValueError("z2 must be finite with one value per row")
    N1, N2 = _as_n_vector(n_eff1, z1.size), _as_n_vector(n_eff2, z1.size)
    M, nb, twostep_gc, fixed = _options(m_snps, n_blocks, twostep, constrain_intercept, z1.size)
    _, _, twostep_h2, _ = _options(m_snps, n_blocks, twostep, None, z1.size)
    h1, i1, d1, _, h1_se, _, _ = _h2_parts(z1 * z1, ell, wl, N1, M, nb, None, twostep_h2)
    h2, i2, d2, _, h2_se, _, _ = _h2_parts(z2 * z2, ell, wl, N2, M, nb, None, twostep_h2)

    def weights(ld, rows, rho_g, intercept):
        n1, n2, wr = (N1, N2, wl) if rows is None else (N1[rows], N2[rows], wl[rows])
        return _gencov_weights(ld, wr, n1, n2, M, h1, h2, rho_g,
                               0.0 if intercept is None else intercept, i1, i2)

    step1 = None if twostep_gc is None else (z1 * z1 < twostep_gc) & (z2 * z2 < twostep_gc)
    gc, ic, dg, di, gc_se, ic_se, _ = _regression(z1 * z2, ell, wl, np.sqrt(N1 * N2), M, nb, fixed,
                                              step1, 0.0, weights)
    if h1 > 0.0 and h2 > 0.0:
        rg = gc / np.sqrt(h1 * h2)
        with np.errstate(invalid="ignore"):
            ratio_dele = dg / np.sqrt(d1 * d2)
        pseudo = nb * rg - (nb - 1) * ratio_dele
        rg_se = float(np.sqrt(np.var(pseudo, ddof=1) / nb)) if np.all(np.isfinite(pseudo)) else float("nan")
    else:
        rg, rg_se = float("nan"), float("nan")
    return LDSCRgResult(rg=float(rg), rg_se=rg_se, gcov=float(gc), gcov_intercept=float(ic),
                        h2=(float(h1), float(h2)), constrain_intercept=fixed,
                        gcov_intercept_se=ic_se, n_jackknife_blocks=nb, gcov_se=gc_se,
                        h2_se=(h1_se, h2_se))
