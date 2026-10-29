"""Bivariate LDSC, adapted from bipred (commit fa211a1).

The fast delete-block fits and their conservative singularity fallback match
bipred. The two traits must be harmonized and sorted in genomic order.
"""
from __future__ import annotations
from contextlib import contextmanager
from dataclasses import dataclass
import numpy as np
from ._numba import HAVE_NUMBA, _set_threads
from ._common import _as_n_vector, _finite_control, _integer_at_least
from .ldsc import _weights, _wls
from ._ldsc_kernels import _jackknife_fits, _jackknife_fits_par

@contextmanager
def _pinned_threads(ncores):
    """Pin Numba's thread mask for one parallel launch, then restore it."""
    if not HAVE_NUMBA or ncores <= 1:
        yield
        return
    from numba import get_num_threads, set_num_threads

    previous = get_num_threads()
    _set_threads(ncores)
    try:
        yield
    finally:
        set_num_threads(previous)


def _require_slope_information(x, constrain):
    """Raise when a slope is not identified by the selected design."""
    if x.size == 0:
        raise np.linalg.LinAlgError("LDSC regression has no observations")
    if constrain is None:
        if x.size < 2 or not np.any(x != x[0]):
            raise np.linalg.LinAlgError(
                "LDSC slope and intercept require variation in LD score times N")
    elif not np.any(x != 0.0):
        raise np.linalg.LinAlgError(
            "constrained-intercept LDSC requires a nonzero slope predictor")


def _fit_slope(y, x, ell_w, n_iter, constrain):
    """Iterated WLS of y on x with LDSC heteroscedasticity/overcounting weights."""
    _require_slope_information(x, constrain)
    pred = np.ones_like(y)
    slope = intercept = 0.0
    for _ in range(n_iter + 1):
        slope, intercept = _wls(x, y, _weights(pred, ell_w), constrain)
        pred = np.maximum(intercept + slope * x, 1.0)
    return slope, intercept


def _as_finite_vector(value, name):
    """Return a nonempty, one-dimensional finite float array."""
    try:
        value = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a one-dimensional numeric array") from None
    if value.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if value.size == 0:
        raise ValueError(f"{name} must be nonempty")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain only finite values")
    return value


def _as_sample_size(value, name, m):
    """Strict ``n_eff``: reject strings/booleans outright, then the seam's check.

    ldpred3's ``_as_n_vector`` coerces numeric strings and booleans; bipred
    rejects them at the boundary, so the strict pre-check stays here and the
    scalar-or-length-m mechanics are delegated to the seam.
    """
    raw = np.asarray(value)
    invalid = raw.dtype.kind in "bSU" or (raw.dtype.kind == "O" and
        any(isinstance(x, (bool, np.bool_, str, bytes)) for x in raw.flat))
    if invalid:
        raise ValueError(f"{name} must be a positive finite scalar or length-m vector")
    return _as_n_vector(value, m)


def _z_from_standardized(beta_std, n_eff, name):
    """Recover exact signed z scores from LDpred3-standardized effects."""
    if np.any(np.abs(beta_std) >= 1.0):
        raise ValueError(
            f"{name} must contain standardized effects with absolute value < 1")
    beta2 = beta_std * beta_std
    np.subtract(1.0, beta2, out=beta2)
    np.sqrt(beta2, out=beta2)
    z = np.sqrt(n_eff)
    z *= beta_std
    z /= beta2
    return z


def _as_finite_scalar(value, name, *, positive=False):
    """Return a finite float, rejecting booleans and arrays.

    ``positive=True`` also requires a value above zero.
    """
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite scalar")
    try:
        value = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite scalar") from None
    if value.ndim != 0 or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite scalar")
    value = float(value)
    if positive and value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


@dataclass
class LDSCRgResult:
    """Output of :func:`ldsc_rg`."""

    rg: float                   # genetic correlation (can fall outside [-1,1] when noisy)
    rg_se: float                # block-jackknife standard error
    gcov: float                 # genetic covariance (cross-trait slope)
    gcov_intercept: float       # sample overlap and/or correlated confounding
    h2: tuple                   # (h2_1, h2_2) marginal heritabilities
    constrain_intercept: float | None = None  # value passed to ldsc_rg, or None if free
    gcov_intercept_se: float = float("nan")  # delete-a-block SE, independent of rg validity
    n_jackknife_blocks: int = 0

    def __repr__(self):
        return (f"LDSCRgResult(rg={self.rg:+.3f} ± {self.rg_se:.3f}, "
                f"gcov={self.gcov:+.4f}, h2=({self.h2[0]:.3f}, {self.h2[1]:.3f}))")


def ldsc_chi2_mask(beta_hat, n_eff, *, cap=80.0, n_scale=0.001):
    """LDSC chi-square exclusion mask for standardized effects.

    Retain rows whose exact z² is at most max(cap, n_scale * N). Apply this
    only to the regression rows; keep m_snps at the full score-map count.
    """
    b = _as_finite_vector(beta_hat, "beta_hat")
    n = _as_sample_size(n_eff, "n_eff", b.size)
    cap = _as_finite_scalar(cap, "cap", positive=True)
    n_scale = _as_finite_scalar(n_scale, "n_scale")
    if n_scale < 0.0:
        raise ValueError("n_scale must be finite and >= 0")
    z = _z_from_standardized(b, n, "beta_hat")
    return (z * z) <= np.maximum(cap, n_scale * n)


def ldsc_rg(beta_hat1, beta_hat2, ld_scores, n_eff1, n_eff2, *, m_snps=None,
            n_blocks=200, n_iter=2, constrain_intercept=None, ncores=1,
            _input_is_z=False):
    """Genetic correlation by cross-trait LD Score regression.

    ``beta_hat1`` and ``beta_hat2`` are LDpred3-standardized marginal effects,
    not original per-allele effects. Their exact signed z scores are
    ``sqrt(N) * beta_hat / sqrt(1 - beta_hat**2)``. To provide signed z scores
    directly, use :func:`ldsc_rg_z`.

    Each LD score must be positive and refer to the same full SNP map as
    ``m_snps``. If regression rows are a subset, retain the full score-map
    count in ``m_snps``; the default is the number of supplied rows. Sample
    sizes can be scalar or per variant. All arrays must be allele-aligned and
    in genomic order: contiguous rows define the jackknife blocks.

    The fit regresses z1*z2 on sqrt(N1*N2)*ell/M and fits both marginal
    heritabilities on the common rows. The slope is genetic covariance and
    rg = gcov/sqrt(h2_1*h2_2). The cross intercept is free unless fixed by
    ``constrain_intercept``. Each of ``n_blocks`` deletions refits both
    marginals and cross-trait weights. With Numba, ``ncores`` parallelizes
    independent deletions; the reference NumPy path remains available.

    This one-step fit does not filter rows. Apply an appropriate long-range LD
    exclusion and the high-chi-square cap before fitting where warranted.
    The rg or its SE is NaN if the relevant marginal heritabilities are
    non-positive, or a required delete-block fit is invalid. The cross
    intercept SE is reported separately and is zero for a valid fixed
    intercept. Neither the intercept nor the fit alone identifies the cause of
    confounding, overlap, or LD/reference mismatch.
    """
    b1 = _as_finite_vector(beta_hat1, "beta_hat1")
    b2 = _as_finite_vector(beta_hat2, "beta_hat2")
    ell = _as_finite_vector(ld_scores, "ld_scores")
    m = b1.shape[0]
    if b2.shape != (m,) or ell.shape != (m,):
        raise ValueError("beta_hat1, beta_hat2, and ld_scores must have equal length")
    if np.any(ell <= 0.0):
        raise ValueError("ld_scores must contain only positive values")
    N1 = _as_sample_size(n_eff1, "n_eff1", m)
    N2 = _as_sample_size(n_eff2, "n_eff2", m)
    M = float(m) if m_snps is None else _as_finite_scalar(
        m_snps, "m_snps", positive=True)
    n_blocks = _integer_at_least("n_blocks", n_blocks, 1)
    n_iter = _integer_at_least("n_iter", n_iter, 0)
    ncores = _integer_at_least("ncores", ncores, 1)
    if constrain_intercept is not None:
        constrain_intercept = _finite_control(
            "constrain_intercept", constrain_intercept)

    z1 = b1.copy() if _input_is_z else _z_from_standardized(b1, N1, "beta_hat1")
    chi1 = z1 * z1
    z2 = b2.copy() if _input_is_z else _z_from_standardized(b2, N2, "beta_hat2")
    chi2 = z2 * z2
    cross = z1 * z2
    if not (np.all(np.isfinite(chi1)) and np.all(np.isfinite(chi2))
            and np.all(np.isfinite(cross))):
        raise ValueError("z scores are too large to form finite regression statistics")
    del z1, z2
    sqrt_n1n2 = np.sqrt(N1) * np.sqrt(N2)
    x1 = N1 * ell / M
    x2 = N2 * ell / M
    xc = sqrt_n1n2 * ell / M
    ell_w = np.maximum(ell, 1.0)

    def fit(sel):
        # Gather each selected column exactly once. The jackknife calls this nb
        # times, so repeating ``x1[sel]`` / ``ell_w[sel]`` inline would allocate
        # and copy several extra length-m arrays per delete-a-block replicate.
        ell_s = ell_w[sel]
        x1s, x2s, xcs, cross_s = x1[sel], x2[sel], xc[sel], cross[sel]
        h1, i1 = _fit_slope(chi1[sel], x1s, ell_s, n_iter, None)
        h2, i2 = _fit_slope(chi2[sel], x2s, ell_s, n_iter, None)
        _require_slope_information(xcs, constrain_intercept)
        pred1 = np.maximum(i1 + h1 * x1s, 1.0)
        pred2 = np.maximum(i2 + h2 * x2s, 1.0)
        # For approximately bivariate-normal z scores,
        # Var(z1*z2) = E[z1^2] E[z2^2] + E[z1*z2]^2. This is the
        # Gencov.weights formula in the reference LDSC implementation; ell_w
        # supplies its LD-overcounting factor.
        pred_cross = np.full_like(
            cross_s, 0.0 if constrain_intercept is None else constrain_intercept)
        for _ in range(n_iter + 1):
            variance = pred1 * pred2 + pred_cross * pred_cross
            w = 1.0 / (ell_s * np.maximum(variance, 1e-6))
            gcov, ic = _wls(xcs, cross_s, w, constrain_intercept)
            pred_cross = ic + gcov * xcs
        return h1, h2, gcov, ic

    try:
        # A view keeps the reference arithmetic without copying every full
        # column merely to select all rows. Deleted subsets still gather below.
        h1, h2, gcov, ic = fit(slice(None))
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "LDSC regression is singular; LD scores must vary when intercepts "
            "are estimated") from exc

    rg = gcov / np.sqrt(h1 * h2) if h1 > 0.0 and h2 > 0.0 else float("nan")

    nb = int(min(n_blocks, m))
    rg_jk = []
    intercept_jk = []
    jackknife_valid = nb >= 2
    if nb >= 2:
        full = np.arange(m)
        # ``array_split`` yields contiguous ranges, so delete-a-block is two
        # contiguous slices of the index vector rather than a length-m boolean
        # mask rebuilt (and re-scanned by every gather) for each of nb blocks.
        splits = np.array_split(full, nb)
        starts = np.asarray([s[0] for s in splits], dtype=np.int64)
        stops = np.asarray([s[-1] + 1 for s in splits], dtype=np.int64)
        # Keep the full-data point estimate on the reference route. Only the
        # repeated jackknife fits enter the compiled, allocation-free loop.
        native = None
        if HAVE_NUMBA:
            kernel = _jackknife_fits_par if ncores > 1 else _jackknife_fits
            with _pinned_threads(ncores):
                native = kernel(
                    x1, x2, xc, chi1, chi2, cross, ell_w, starts, stops,
                    n_iter, constrain_intercept is not None,
                    0.0 if constrain_intercept is None else constrain_intercept)
        for b, (start, stop) in enumerate(zip(starts, stops)):
            if native is not None and np.all(np.isfinite(native[b])):
                hb1, hb2, gb, ib = native[b]
            else:
                # The reference solve owns near-singular and near-zero-slope
                # decisions. It is also the no-Numba implementation.
                keep = np.concatenate((full[:start], full[stop:]))
                try:
                    hb1, hb2, gb, ib = fit(keep)
                except np.linalg.LinAlgError:
                    jackknife_valid = False
                    break
            intercept_jk.append(ib)
            value = gb / np.sqrt(hb1 * hb2) if hb1 > 0.0 and hb2 > 0.0 else float("nan")
            rg_jk.append(value)

    def jackknife_se(values):
        # Never condition an SE on the subset of deletions that happened to fit.
        if not jackknife_valid or len(values) != nb or not np.all(np.isfinite(values)):
            return float("nan")
        values = np.asarray(values)
        return float(np.sqrt((nb - 1) / nb * np.sum((values - values.mean()) ** 2)))

    rg_se = jackknife_se(rg_jk) if np.isfinite(rg) else float("nan")

    return LDSCRgResult(rg=float(rg), rg_se=rg_se, gcov=float(gcov),
                        gcov_intercept=float(ic), h2=(float(h1), float(h2)),
                        constrain_intercept=constrain_intercept,
                        gcov_intercept_se=jackknife_se(intercept_jk),
                        n_jackknife_blocks=nb)



def ldsc_rg_z(z1, z2, ld_scores, n_eff1, n_eff2, **kwargs):
    """Cross-trait LDSC on signed z scores, with the same options as ldsc_rg.

    Unlike ldsc_rg, these are original signed z scores, not LDpred3-standardized
    effects. Both arrays must refer to the same effect allele and genomic order.
    """
    return ldsc_rg(z1, z2, ld_scores, n_eff1, n_eff2, _input_is_z=True, **kwargs)
