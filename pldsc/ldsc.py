"""Efficient univariate LDSC, adapted from LDpred3 (commit 2f5d108).

The estimand and regression weights match the original implementation. Scores
must be built over the same M variants used in the regression scaling.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ._numba import HAVE_NUMBA, _jit
from ._common import _as_n_vector, _integer_at_least, _finite_control

@dataclass
class LDSCResult:
    """Output of :func:`ldsc_h2`."""

    h2: float                   # SNP heritability (regression slope on x = N*ell/M)
    h2_se: float                # block-jackknife standard error
    intercept: float            # ~1 under a compatible, unconfounded model
    intercept_se: float
    mean_chisq: float
    ratio: float                # (intercept - 1) / (mean_chisq - 1), model-dependent

    @property
    def h2_ci(self):
        """Approximate 95% CI from the jackknife SE."""
        return (self.h2 - 1.96 * self.h2_se, self.h2 + 1.96 * self.h2_se)

    def __repr__(self):
        return (f"LDSCResult(h2={self.h2:.3f} ± {self.h2_se:.3f}, "
                f"intercept={self.intercept:.3f}, mean_chi2={self.mean_chisq:.2f})")


def _wls(x, y, w, constrain_intercept):
    """Return WLS slope/intercept, using minimum-norm fit for singular subsets."""
    if constrain_intercept is None:
        # Shift before centering, using values nearest zero so a downweighted
        # outlier cannot erase small differences among the other rows.
        x0, y0 = x[np.argmin(np.abs(x))], y[np.argmin(np.abs(y))]
        dx = np.asarray(x, dtype=float) - x0
        dy = np.asarray(y, dtype=float) - y0
        sw = np.sum(w)
        mx, my = np.dot(w, dx) / sw, np.dot(w, dy) / sw
        dx -= mx
        dy -= my
        wx = w * dx
        variance = np.dot(wx, dx)
        if variance == 0.0:
            intercept = (y0 + my) / (1.0 + x0 * x0)
            return float(x0 * intercept), float(intercept)
        slope = np.dot(wx, dy) / variance
        return float(slope), float(y0 + my - slope * (x0 + mx))
    c = float(constrain_intercept)
    slope = np.sum(w * x * (y - c)) / np.sum(w * x * x)
    return float(slope), c


def _weights(pred_mean, ell_w):
    """Heteroscedasticity (1 / 2 mu^2) x overcounting (1 / ell) regression weights."""
    het = 1.0 / (2.0 * np.maximum(pred_mean, 1e-6) ** 2)
    over = 1.0 / np.maximum(ell_w, 1.0)
    return het * over


def _jackknife_fits(x, y, ell_w, starts, stops, n_iter, constrained, c,
                    h2_out, int_out, singular):
    """Every delete-a-block IRLS fit of :func:`ldsc_h2`, without masked copies.

    Free-intercept fits cache weights, then center in a second pass to avoid
    cancellation. Fixed-intercept fits need only one pass. Exactly singular
    subsets are flagged for the caller's minimum-norm route.
    """
    m = x.shape[0]
    weights = np.empty(0 if constrained else m)
    x0, y0 = np.min(x), np.min(y)       # both are non-negative in ldsc_h2
    for b in range(starts.shape[0]):
        lo = starts[b]
        hi = stops[b]
        first_x = x[0 if lo > 0 else hi]
        varies = False
        slope = 0.0
        intercept = 0.0
        for it in range(n_iter + 1):
            sw = 0.0
            swx = 0.0
            swxx = 0.0
            swy = 0.0
            swxy = 0.0
            # The kept rows are the two stretches either side of [lo, hi).
            for part in range(2):
                a = 0 if part == 0 else hi
                e = lo if part == 0 else m
                for i in range(a, e):
                    xi = x[i]
                    if it == 0:
                        pred = 1.0
                    else:
                        pred = intercept + slope * xi
                        if pred < 1.0:
                            pred = 1.0
                    li = ell_w[i]
                    if li < 1.0:
                        li = 1.0
                    w = (1.0 / (2.0 * pred * pred)) * (1.0 / li)
                    yi = y[i]
                    if constrained:
                        swxx += w * xi * xi
                        swxy += w * xi * (yi - c)
                    else:
                        varies = varies or xi != first_x
                        sw += w
                        swx += w * (xi - x0)
                        swy += w * (yi - y0)
                        weights[i] = w
            if constrained:
                slope = swxy / swxx
                intercept = c
            else:
                mx, my = swx / sw, swy / sw
                for part in range(2):
                    a = 0 if part == 0 else hi
                    e = lo if part == 0 else m
                    for i in range(a, e):
                        dx = (x[i] - x0) - mx
                        wx = weights[i] * dx
                        swxx += wx * dx
                        swxy += wx * ((y[i] - y0) - my)
                if not varies or swxx == 0.0:
                    singular[b] = True
                    break
                slope = swxy / swxx
                intercept = y0 + my - slope * (x0 + mx)
        h2_out[b] = slope
        int_out[b] = intercept


_jackknife_fits_jit = _jit(_jackknife_fits)


#: High-chi-square exclusion for LD-score regression: drop a variant whose
#: chi-square exceeds ``max(LDSC_CHI2_CAP, LDSC_CHI2_N_SCALE * N)``. These are
#: the reference implementation's own numbers (``max(80, N/1000)``), and the
#: ones ``smartpred._runner_trait_qc._QC_LDSC_PARAMETERS`` uses, so a CLI regression and
#: a SMARTpred regression exclude the same rows. The cap matters because a
#: handful of very large chi-squares at high LD scores dominate the slope, and
#: the slope is the heritability.
LDSC_CHI2_CAP = 80.0
LDSC_CHI2_N_SCALE = 0.001


def ldsc_chi2_keep(chisq, n_eff):
    """Rows below the standard high-chi-square cap, and the cap itself.

    Returns ``(keep, cap)`` with ``cap`` per variant when ``n_eff`` is a
    vector. A caller that regresses every row is free to ignore this, but
    should say so: it is the largest single lever on the point estimate.
    """
    chisq = np.asarray(chisq, dtype=float)
    cap = np.maximum(LDSC_CHI2_CAP,
                     LDSC_CHI2_N_SCALE * np.asarray(n_eff, dtype=float))
    return np.isfinite(chisq) & (chisq <= cap), cap


def ldsc_h2(chisq, ld_scores, n_eff, *, m_snps=None, n_blocks=100,
            constrain_intercept=None, n_iter=2):
    """Estimate SNP heritability by LD Score regression.

    Parameters
    ----------
    chisq : array_like (m,)
        ``z_j^2`` from original effects and SEs. For standardized effects,
        ``z^2 = N*beta_std^2/(1-beta_std^2)``; ``N*beta_std^2`` is approximate.
    ld_scores : array_like (m,)
        LD scores from :func:`ld_scores`.
    n_eff : float or array_like
        GWAS sample size (scalar or per-variant).
    m_snps : int, optional
        Number of SNPs the heritability is defined over (default: ``len(chisq)``,
        the regression rows). It only rescales ``x``, hence ``h2``; pass the
        reference-panel count when the scores come from that panel, and the
        pre-exclusion count when rows were dropped from the regression only.
    n_blocks : int, default 100
        Contiguous delete-a-block groups. This affects SE precision and cost,
        not the full-data point estimate.
    constrain_intercept : float, optional
        Fix the intercept (e.g. ``1.0`` to assume no confounding) instead of
        estimating it.
    n_iter : int, default 2
        Number of weight-update iterations (weights depend on the fitted mean).

    Returns
    -------
    LDSCResult

    Notes
    -----
    Regress ``chisq_j`` on ``x_j = N_j ell_j / M`` so the slope is ``h2``
    and the intercept can absorb confounding, genomic-control deflation,
    or LD/model mismatch.
    Weights are ``1 / (2 mu_j² max(ell_j, 1))``: the heteroscedasticity of
    a chi-square with mean ``mu_j`` times LD over-counting. The first fit
    takes ``mu_j = 1``; each of ``n_iter`` refits sets ``mu_j =
    max(intercept + h2 x_j, 1)``. Standard errors are the delete-a-block
    jackknife ``sqrt((B-1)/B sum_b (theta_b - mean theta)²)`` over ``B``
    contiguous blocks in input order, each refit with the same iterations.
    ``ratio = (intercept - 1)/(mean chisq - 1)``. The high-chi-square row cap
    (:func:`ldsc_chi2_keep`) is the caller's to apply.

    Free-intercept fits use weighted centering to avoid cancellation for
    nearly constant predictors. A constant full-data predictor with at least
    two rows is rejected. Single-row fits and constant-predictor jackknife
    subsets use minimum-norm coefficients; these do not identify slope and
    intercept separately. With fewer than two jackknife groups the slope SE
    is unavailable. A constrained intercept has SE zero by definition.
    """
    # Validate before a malformed N or control can resemble a valid null fit.
    chisq = np.asarray(chisq, dtype=float)
    ell = np.asarray(ld_scores, dtype=float)
    if chisq.ndim != 1 or ell.ndim != 1:
        raise ValueError("chisq and ld_scores must be one-dimensional")
    m = chisq.shape[0]
    if m == 0:
        raise ValueError("chisq must contain at least one variant")
    if ell.shape[0] != m:
        raise ValueError("chisq and ld_scores must have the same length")
    if not np.all(np.isfinite(chisq)) or np.any(chisq < 0.0):
        raise ValueError("chisq must contain only finite, non-negative values")
    if not np.all(np.isfinite(ell)) or np.any(ell <= 0.0):
        raise ValueError("ld_scores must contain only finite, positive values")
    N = _as_n_vector(n_eff, m)
    n_iter = _integer_at_least("n_iter", n_iter, 0)
    n_blocks = _integer_at_least("n_blocks", n_blocks, 1)
    if m_snps is not None:
        m_snps = _finite_control("m_snps", m_snps, lower=1.0)
    if constrain_intercept is not None:
        constrain_intercept = _finite_control("constrain_intercept",
                                              constrain_intercept)
    M = float(m_snps if m_snps is not None else m)

    x = N * ell / M               # slope on x is h2 directly
    if constrain_intercept is None and m >= 2 and np.ptp(x) == 0.0:
        raise ValueError(
            "LD scores (times N) have no variation across variants, so the "
            "LDSC slope and intercept are not jointly identified; provide "
            "variants spanning a range of LD scores or fix the intercept with "
            "constrain_intercept")
    ell_w = np.maximum(ell, 1.0)

    def fit(mask=None):
        xi, yi, li = (x, chisq, ell_w) if mask is None else (x[mask], chisq[mask], ell_w[mask])
        pred = np.ones_like(yi)
        slope = intercept = 0.0
        for _ in range(n_iter + 1):
            w = _weights(pred, li)
            slope, intercept = _wls(xi, yi, w, constrain_intercept)
            pred = np.maximum(intercept + slope * xi, 1.0)
        return slope, intercept

    h2, intercept = fit()

    nb = int(min(n_blocks, m))
    if nb < 2:
        # One block has no delete-a-block uncertainty estimate.
        h2_se = float("nan")
        int_se = 0.0 if constrain_intercept is not None else float("nan")
    else:
        # np.array_split's contiguous ranges: the first m % nb get one extra.
        sizes = np.full(nb, m // nb, dtype=np.int64)
        sizes[:m % nb] += 1
        stops = np.cumsum(sizes)
        starts = stops - sizes
        h2_jk = np.empty(nb)
        int_jk = np.empty(nb)
        singular = np.zeros(nb, dtype=np.bool_)
        if HAVE_NUMBA:
            # Visit retained ranges directly, reusing one weights buffer.
            constrained = constrain_intercept is not None
            _jackknife_fits_jit(
                x, chisq, ell_w, starts, stops, n_iter, constrained,
                float(constrain_intercept) if constrained else 0.0,
                h2_jk, int_jk, singular)
        else:
            singular[:] = True
        for b in np.flatnonzero(singular):
            keep = np.ones(m, dtype=bool)
            keep[starts[b]:stops[b]] = False
            h2_jk[b], int_jk[b] = fit(keep)
        fac = (nb - 1) / nb
        h2_se = float(np.sqrt(fac * np.sum((h2_jk - h2_jk.mean()) ** 2)))
        int_se = (0.0 if constrain_intercept is not None
                  else float(np.sqrt(fac * np.sum((int_jk - int_jk.mean()) ** 2))))

    mean_chisq = float(chisq.mean())
    ratio = (intercept - 1.0) / (mean_chisq - 1.0) if mean_chisq > 1.0 else float("nan")
    return LDSCResult(h2=float(h2), h2_se=h2_se, intercept=float(intercept),
                      intercept_se=int_se, mean_chisq=mean_chisq, ratio=ratio)
