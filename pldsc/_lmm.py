"""Private dense mixed-model LD-score oracle (milestone M2).

Reference-quality float64 dense implementation of the generalized residual
precision ``P_c`` and the oracle conditional moment derived in
``docs/METHOD.md``. This is deliberately private (underscore module): it makes
no public estimator claim and is used by the test suite and the ``experiments/``
scripts to validate the theory before any scalable or fitted-component path
exists. It never constructs an ``M x M`` LD matrix and permits dense ``N x N``
matrices only because it is the reference path.

All operators take the covariance and design as arrays, so the same code serves
the oracle-``V_c`` case (Sections M2-M3) and, by passing a different ``V_true``
to :func:`null_inflation`, the misspecified-``V_c`` diagnostic (Section M4).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Float = NDArray[np.float64]


def residual_precision(V: Float, C: Float) -> Float:
    r"""Generalized residual precision (Eq. M4).

    ``P = V^{-1} - V^{-1} C (C^T V^{-1} C)^+ C^T V^{-1}``.

    Parameters
    ----------
    V
        ``(N, N)`` symmetric positive-definite null covariance.
    C
        ``(N, K)`` fixed-effect design; include the intercept column. A 1-D
        array is treated as a single column. Rank-deficient ``C`` is handled by
        the pseudoinverse.

    Returns
    -------
    numpy.ndarray
        ``(N, N)`` symmetric PSD operator satisfying ``P C = 0`` and
        ``P V P = P``.
    """
    V = np.asarray(V, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    if C.ndim == 1:
        C = C[:, None]
    W = np.linalg.inv(V)
    WC = W @ C
    P = W - WC @ np.linalg.pinv(C.T @ WC) @ WC.T
    return 0.5 * (P + P.T)  # symmetrize against round-off


def score_pieces(P: Float, X: Float) -> tuple[Float, Float]:
    """Return ``PX = P @ X`` and the self-weights ``q_j = x_j^T P x_j``."""
    PX = P @ X
    q = np.einsum("ij,ij->j", X, PX)
    return PX, q


def ld_scores(
    P: Float,
    X: Float,
    window_mask: Float | None = None,
    *,
    PX: Float | None = None,
    q: Float | None = None,
) -> Float:
    r"""PLDSC LD score ``L_j = sum_k q_k {r_jk^{(P)}}^2`` (Eq. M7).

    Equals ``(1 / q_j) sum_k (x_j^T P x_k)^2``. With ``window_mask`` (an
    ``(M, M)`` boolean array) the sum is restricted to ``True`` entries;
    otherwise it is the full sum over all columns of ``X``.
    """
    if PX is None or q is None:
        PX, q = score_pieces(P, X)
    A = X.T @ PX  # a_jk = x_j^T P x_k
    A2 = A * A
    if window_mask is not None:
        A2 = A2 * window_mask
    return A2.sum(axis=1) / q


def null_inflation(
    P: Float,
    X: Float,
    V_true: Float,
    *,
    PX: Float | None = None,
    q: Float | None = None,
) -> Float:
    r"""Null (``tau^2 = 0``) moment ``E[z_j^2] = x_j^T P V_true P x_j / q_j``.

    Equals 1 exactly when ``V_true`` is the covariance used to build ``P``
    (the oracle case, by the identity ``P V P = P``). Passing a different
    ``V_true`` quantifies the confounding that survives a misspecified or fitted
    ``V_c`` (the intercept-as-diagnostic use-case, Section M4/M9).
    """
    if PX is None or q is None:
        PX, q = score_pieces(P, X)
    return np.einsum("ij,ij->j", PX, V_true @ PX) / q


def oracle_moment(
    P: Float,
    X: Float,
    tau2: float,
    V_true: Float | None = None,
    window_mask: Float | None = None,
) -> Float:
    r"""Homoscedastic conditional moment ``E[z_j^2 | X]`` (Eq. M6/M7).

    ``null_inflation + tau2 * L_j``. With ``V_true=None`` the intercept is the
    oracle value 1; pass ``V_true`` to use a (possibly misspecified) covariance.
    """
    PX, q = score_pieces(P, X)
    null = (
        np.ones_like(q)
        if V_true is None
        else null_inflation(P, X, V_true, PX=PX, q=q)
    )
    return null + tau2 * ld_scores(P, X, window_mask, PX=PX, q=q)
