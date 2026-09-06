"""Covariate-adjusted linkage-disequilibrium scores."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatMatrix = NDArray[np.float64]
FloatVector = NDArray[np.float64]


def _as_genotype_matrix(genotypes: ArrayLike) -> FloatMatrix:
    values = np.asarray(genotypes, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("genotypes must have shape (samples, variants)")
    if values.shape[0] < 3 or values.shape[1] == 0:
        raise ValueError("genotypes must contain at least 3 samples and 1 variant")
    if not np.all(np.isfinite(values)):
        raise ValueError("genotypes must contain only finite values")
    return values


def _design_matrix(n_samples: int, covariates: ArrayLike | None) -> FloatMatrix:
    intercept = np.ones((n_samples, 1), dtype=np.float64)
    if covariates is None:
        return intercept

    values = np.asarray(covariates, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != n_samples:
        raise ValueError("covariates must have shape (samples, covariates)")
    if not np.all(np.isfinite(values)):
        raise ValueError("covariates must contain only finite values")
    return np.column_stack((intercept, values))


def _residualize(
    genotypes: ArrayLike, covariates: ArrayLike | None
) -> tuple[FloatMatrix, int]:
    values = _as_genotype_matrix(genotypes)
    design = _design_matrix(values.shape[0], covariates)
    coefficients, _, rank, _ = np.linalg.lstsq(design, values, rcond=None)
    residuals = values - design @ coefficients
    residual_df = values.shape[0] - rank
    if residual_df <= 1:
        raise ValueError("covariates leave fewer than 2 residual degrees of freedom")
    return residuals, residual_df


def residualize_genotypes(
    genotypes: ArrayLike, covariates: ArrayLike | None = None
) -> FloatMatrix:
    """Project an intercept and fixed covariates out of each variant.

    Parameters
    ----------
    genotypes
        Numeric matrix with samples in rows and variants in columns.
    covariates
        Optional sample-aligned vector or matrix. An intercept is always added;
        redundant columns are handled by the least-squares rank calculation.

    Returns
    -------
    numpy.ndarray
        Residual genotypes with the same shape as ``genotypes``.
    """

    residuals, _ = _residualize(genotypes, covariates)
    return residuals


def effective_pair_df(
    relationship: ArrayLike, covariates: ArrayLike | None = None
) -> float:
    """Noise-floor denominator for unlinked-variant squared correlations.

    For two unlinked variants in a panel with relationship matrix ``A`` and LD
    covariates ``C``, the finite-sample noise floor of the squared correlation
    is ``E[r^2] ~= tr(Atilde^2) / tr(Atilde)^2`` with ``Atilde = P_C A P_C``.
    This returns its reciprocal ``tr(Atilde)^2 / tr(Atilde^2)``, the effective
    denominator to pass as ``pair_df`` to :func:`local_ld_scores`.

    For unrelated individuals (``A = I``, intercept only) this is ``N - 1``,
    reproducing the usual convention up to the customary O(1) offset. With
    relatives it is strictly smaller, so the default ``residual_df - 1``
    under-subtracts. Verified in ``experiments/relatedness_effective_df.py``.

    Parameters
    ----------
    relationship
        Symmetric sample-by-sample relationship matrix (a GRM, or a pedigree
        additive relationship matrix).
    covariates
        Optional sample-aligned fixed covariates. An intercept is always added,
        matching :func:`local_ld_scores`.

    Returns
    -------
    float
        The effective pair degrees of freedom.
    """

    A = np.asarray(relationship, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("relationship must be a square matrix")
    if not np.all(np.isfinite(A)):
        raise ValueError("relationship must contain only finite values")
    if not np.allclose(A, A.T):
        raise ValueError("relationship must be symmetric")

    design = _design_matrix(A.shape[0], covariates)
    q, _ = np.linalg.qr(design)
    # Atilde = P A P with P = I - q q^T, formed without an explicit N x N P.
    qa = q.T @ A
    At = A - q @ qa - (A @ q) @ q.T + q @ (qa @ q) @ q.T

    trace = float(np.trace(At))
    trace_sq = float(np.sum(At * At))
    if trace_sq <= 0.0:
        raise ValueError("projected relationship matrix has no variance")
    return trace**2 / trace_sq


def local_ld_scores(
    genotypes: ArrayLike,
    positions_cm: ArrayLike,
    *,
    covariates: ArrayLike | None = None,
    window_cm: float = 1.0,
    finite_sample_correction: bool = True,
    pair_df: float | None = None,
    include_self: bool = True,
) -> FloatVector:
    """Calculate covariate-adjusted LD scores for one chromosome.

    Correlations are computed from genotypes after projecting out an intercept
    and the supplied fixed covariates. The finite-sample adjustment uses the
    residual rank as specified by Equation (7) in ``docs/IDEA_REVIEW.md``.
    Negative adjusted terms are retained.

    This fixed-covariate calculation does not implement LOCO or mixed-model
    adjustment. If PCs were estimated from the same genotype matrix, the
    fixed-design finite-sample argument is only approximate.

    Parameters
    ----------
    genotypes
        Numeric matrix with samples in rows and variants in columns.
    positions_cm
        Finite, nondecreasing genetic positions for a single chromosome.
    covariates
        Optional sample-aligned fixed covariates used by the GWAS.
    window_cm
        Inclusive distance on either side of each variant, in centimorgans.
    finite_sample_correction
        Apply the residual-degrees-of-freedom adjustment to squared
        correlations.
    pair_df
        Denominator used in that adjustment, ``r2 - (1 - r2) / pair_df``.
        Defaults to ``residual_df - 1``, which assumes independent samples. A
        panel containing relatives has a higher unlinked-pair noise floor and
        needs the smaller :func:`effective_pair_df`; see ``docs/METHOD.md``
        Section 10.1.
    include_self
        Include the standard LD Score self-term, fixed exactly to one.

    Returns
    -------
    numpy.ndarray
        One LD score per input variant.
    """

    values = _as_genotype_matrix(genotypes)
    positions = np.asarray(positions_cm, dtype=np.float64)
    if positions.ndim != 1 or positions.shape[0] != values.shape[1]:
        raise ValueError("positions_cm must contain one value per variant")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_cm must contain only finite values")
    if np.any(np.diff(positions) < 0):
        raise ValueError("positions_cm must be nondecreasing")
    if not np.isfinite(window_cm) or window_cm < 0:
        raise ValueError("window_cm must be finite and nonnegative")
    if pair_df is not None and (not np.isfinite(pair_df) or pair_df <= 0.0):
        raise ValueError("pair_df must be finite and positive")

    residuals, residual_df = _residualize(values, covariates)
    denominator = float(residual_df - 1) if pair_df is None else float(pair_df)
    norms = np.linalg.norm(residuals, axis=0)
    input_norms = np.linalg.norm(values, axis=0)
    tolerance = (
        np.finfo(np.float64).eps
        * max(values.shape)
        * np.maximum(input_norms, 1.0)
    )
    invalid = np.flatnonzero(norms <= tolerance)
    if invalid.size:
        indices = ", ".join(map(str, invalid.tolist()))
        raise ValueError(f"variants have zero residual variance at indices: {indices}")

    normalized = residuals / norms
    scores = np.empty(values.shape[1], dtype=np.float64)

    for variant in range(values.shape[1]):
        left = np.searchsorted(
            positions, positions[variant] - window_cm, side="left"
        )
        right = np.searchsorted(
            positions, positions[variant] + window_cm, side="right"
        )
        correlations = normalized[:, left:right].T @ normalized[:, variant]
        squared = np.clip(correlations, -1.0, 1.0) ** 2
        squared[variant - left] = 1.0
        if finite_sample_correction:
            squared = squared - (1.0 - squared) / denominator
        if not include_self:
            squared[variant - left] = 0.0
        scores[variant] = np.sum(squared)

    return scores
