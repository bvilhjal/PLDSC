"""LD scores for multi-ancestry GWAS: meta-analysed and pooled designs.

Derivations, assumptions and numbered equations are in
``docs/MULTI_ANCESTRY.md``. In brief, for a fixed-effect meta-analysis of
ancestry-specific GWAS with z-weights ``c_aj`` (``sum_a c_aj^2 = 1``), MA2 is

    E[chi2_j] = 1 + (N_j h^2 / M) L_j,
    L_j = sum_{a,b} omega_{ab,j} rho_ab sqrt(eta_a eta_b) ell^{ab}_j,

with ``omega_{ab,j} = c_aj c_bj sqrt(N_a N_b) / N_j`` (MA3), and ``ell^{ab}_j``
the cross-ancestry LD score ``sum_k r^a_jk r^b_jk`` (MA1). ``rho_ab`` is the
standardized-scale trans-ancestry genetic correlation, and
``eta_a = h_a^2 / h^2`` is ancestry a's heritability relative to the one
reported. A single-ancestry LD score used in place of ``L`` biases the LDSC
slope by the factor returned by :func:`single_reference_bias` (MA6).

Estimation from reference panels (MA7): the within-ancestry terms
``ell^{aa}`` carry the usual squared-correlation noise floor, so the default
finite-sample correction is applied to them. The cross terms ``ell^{ab}``
(``a != b``) are products of correlations estimated in *independent* panels,
so they carry no noise floor and are left uncorrected.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .ld import _as_genotype_matrix, _residualize

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]
IntArray = NDArray[np.intp]

__all__ = [
    "CrossAncestryScores",
    "ReferenceBias",
    "cross_ancestry_ld_scores",
    "meta_weights",
    "meta_ld_scores",
    "pooled_ld_scores",
    "single_reference_bias",
]


@dataclass(frozen=True)
class CrossAncestryScores:
    """Per-pair LD scores at the regression variants.

    Attributes
    ----------
    scores
        ``(A, A, J)`` symmetric array. ``scores[a, a]`` is ancestry a's LD
        score (noise-floor corrected unless disabled); ``scores[a, b]`` is the
        uncorrected cross-ancestry score ``sum_k r^a_jk r^b_jk``.
    available
        ``(A, J)``: the regression variant is polymorphic in panel a. An
        unavailable variant contributes nothing to any score.
    variances
        ``(A, J)``: residual genotype variance of each regression variant in
        panel a, the per-allele information an inverse-variance meta-analysis
        weights by (``ddof = 0``).
    index
        ``(J,)`` indices of the regression variants in the reference set.
    """

    scores: FloatArray
    available: BoolArray
    variances: FloatArray
    index: IntArray


@dataclass(frozen=True)
class ReferenceBias:
    """How a wrong LD score distorts LD Score regression (MA6).

    ``kappa`` is the ratio of the population slope on the used score to the
    slope on the correct one, i.e. ``h2_hat / h2``. ``offset`` is the intercept
    shift in LD-score units; multiply by ``N h^2 / M`` for chi-square units.
    """

    kappa: float
    offset: float


def _window_bounds(positions: FloatArray, idx: IntArray, window_cm: float):
    lo = np.searchsorted(positions, positions[idx] - window_cm, side="left")
    hi = np.searchsorted(positions, positions[idx] + window_cm, side="right")
    return lo, hi


def _check_positions(positions_cm: ArrayLike, m: int, window_cm: float) -> FloatArray:
    positions = np.asarray(positions_cm, dtype=np.float64)
    if positions.ndim != 1 or positions.shape[0] != m:
        raise ValueError("positions_cm must contain one value per variant")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_cm must contain only finite values")
    if np.any(np.diff(positions) < 0):
        raise ValueError("positions_cm must be nondecreasing")
    if not np.isfinite(window_cm) or window_cm < 0:
        raise ValueError("window_cm must be finite and nonnegative")
    return positions


def _check_regression(regression: ArrayLike | None, m: int) -> IntArray:
    if regression is None:
        return np.arange(m)
    idx = np.asarray(regression)
    if idx.ndim != 1 or idx.size == 0 or not np.issubdtype(idx.dtype, np.integer):
        raise ValueError("regression must be a non-empty 1-D integer index array")
    if np.any(np.diff(idx) <= 0) or idx[0] < 0 or idx[-1] >= m:
        raise ValueError("regression indices must be strictly increasing and in range")
    return idx.astype(np.intp)


def _unit_columns(residuals: FloatArray, raw: FloatArray) -> tuple[FloatArray, BoolArray]:
    """Unit-norm columns; zero-variance columns become zero and unavailable."""
    norms = np.linalg.norm(residuals, axis=0)
    tolerance = (np.finfo(np.float64).eps * max(raw.shape)
                 * np.maximum(np.linalg.norm(raw, axis=0), 1.0))
    ok = norms > tolerance
    out = np.zeros_like(residuals)
    out[:, ok] = residuals[:, ok] / norms[ok]
    return out, ok


def cross_ancestry_ld_scores(
    panels: Sequence[ArrayLike],
    positions_cm: ArrayLike,
    *,
    window_cm: float = 1.0,
    regression: ArrayLike | None = None,
    covariates: Sequence[ArrayLike | None] | None = None,
    finite_sample_correction: bool = True,
    pair_df: Sequence[float] | None = None,
    block_size: int = 256,
) -> CrossAncestryScores:
    """Within- and cross-ancestry LD scores from per-ancestry reference panels.

    Parameters
    ----------
    panels
        One genotype matrix per ancestry (samples in rows), all with the same
        reference variants in the same order and the same counted allele. The
        sign of every cross term depends on consistent allele coding.
    positions_cm
        Shared nondecreasing genetic positions of the reference variants.
    window_cm
        Inclusive window half-width, in centimorgans.
    regression
        Strictly increasing indices of the variants to score (LDSC's
        regression SNPs). Defaults to every reference variant. Windows always
        run over the full reference set.
    covariates
        Optional per-panel within-ancestry covariates (e.g. panel PCs); an
        intercept is always included.
    finite_sample_correction
        Subtract the squared-correlation noise floor from the within-ancestry
        terms, ``r2 - (1 - r2) / pair_df`` (MA7). Cross terms are never
        corrected.
    pair_df
        Per-panel noise-floor denominators. Defaults to ``residual_df - 1``
        (independent samples); pass :func:`pldsc.effective_pair_df` values for
        panels containing relatives.
    block_size
        Regression variants processed per matrix product.

    Returns
    -------
    CrossAncestryScores
    """

    if len(panels) == 0:
        raise ValueError("need at least one panel")
    mats = [_as_genotype_matrix(p) for p in panels]
    m = mats[0].shape[1]
    if any(x.shape[1] != m for x in mats):
        raise ValueError("all panels must have the same variants (columns)")
    A = len(mats)
    positions = _check_positions(positions_cm, m, window_cm)
    idx = _check_regression(regression, m)
    if covariates is None:
        covariates = [None] * A
    if len(covariates) != A:
        raise ValueError("covariates must have one entry per panel")
    if pair_df is not None:
        if len(pair_df) != A:
            raise ValueError("pair_df must have one entry per panel")
        if not all(np.isfinite(d) and d > 0 for d in pair_df):
            raise ValueError("pair_df entries must be finite and positive")
    if block_size < 1:
        raise ValueError("block_size must be positive")

    unit, avail, dfs, var = [], [], [], []
    for a, (X, C) in enumerate(zip(mats, covariates)):
        resid, rdf = _residualize(X, C)
        Xn, ok = _unit_columns(resid, X)
        unit.append(Xn)
        avail.append(ok)
        dfs.append(float(rdf - 1) if pair_df is None else float(pair_df[a]))
        var.append(np.sum(resid[:, idx] ** 2, axis=0) / X.shape[0] * ok[idx])

    J = idx.size
    scores = np.zeros((A, A, J))
    lo_all, hi_all = _window_bounds(positions, idx, window_cm)
    for start in range(0, J, block_size):
        stop = min(start + block_size, J)
        cols = idx[start:stop]
        lo, hi = int(lo_all[start]), int(hi_all[stop - 1])
        rows = np.arange(lo, hi)
        inwin = ((rows[:, None] >= lo_all[None, start:stop])
                 & (rows[:, None] < hi_all[None, start:stop]))
        self_row = cols - lo
        bcols = np.arange(stop - start)
        R = []
        for a in range(A):
            r = np.clip(unit[a][:, lo:hi].T @ unit[a][:, cols], -1.0, 1.0)
            r[self_row, bcols] = avail[a][cols].astype(np.float64)
            R.append(r * inwin)
        for a in range(A):
            sq = R[a] ** 2
            if finite_sample_correction:
                live = inwin & avail[a][lo:hi, None] & avail[a][None, cols]
                sq = sq - live * (1.0 - sq) / dfs[a]
            scores[a, a, start:stop] = sq.sum(axis=0)
            for b in range(a + 1, A):
                cross = np.einsum("kj,kj->j", R[a], R[b])
                scores[a, b, start:stop] = cross
                scores[b, a, start:stop] = cross

    available = np.stack([ok[idx] for ok in avail])
    return CrossAncestryScores(scores=scores, available=available,
                               variances=np.stack(var), index=idx)


def meta_weights(
    n_samples: ArrayLike,
    available: ArrayLike,
    *,
    variances: ArrayLike | None = None,
    scheme: str = "sample_size",
) -> tuple[FloatArray, FloatArray]:
    """Meta-analysis pair weights ``omega_{ab,j}`` (MA3) and reported ``N_j``.

    Parameters
    ----------
    n_samples
        GWAS sample size per ancestry, shape ``(A,)`` or per variant ``(A, J)``.
    available
        ``(A, J)`` booleans: ancestry a contributes variant j to the meta.
    variances
        ``(A, J)`` genotype variances, required for ``scheme="ivw"``.
    scheme
        ``"sample_size"`` (Stouffer / METAL SAMPLESIZE: ``c_aj^2 = N_a/N_j``)
        or ``"ivw"`` (fixed-effect inverse variance on the per-allele scale:
        ``c_aj^2 = N_a s_aj^2 / sum_b N_b s_bj^2``).

    Returns
    -------
    omega, n_total
        ``(A, A, J)`` weights and the ``(J,)`` sample size the meta reports,
        ``N_j = sum`` of contributing ``N_a``.
    """

    avail = np.asarray(available, dtype=bool)
    if avail.ndim != 2:
        raise ValueError("available must have shape (ancestries, variants)")
    A, J = avail.shape
    n = np.asarray(n_samples, dtype=np.float64)
    if n.ndim == 1:
        n = np.broadcast_to(n[:, None], (A, J))
    if n.shape != (A, J):
        raise ValueError("n_samples must have shape (A,) or (A, J)")
    if not np.all(np.isfinite(n)) or np.any(n < 0):
        raise ValueError("n_samples must be finite and nonnegative")
    n_a = n * avail
    n_total = n_a.sum(axis=0)
    safe = np.where(n_total > 0, n_total, 1.0)

    if scheme == "sample_size":
        c2 = n_a / safe
    elif scheme == "ivw":
        if variances is None:
            raise ValueError("scheme='ivw' needs variances")
        v = np.asarray(variances, dtype=np.float64)
        if v.shape != (A, J) or not np.all(np.isfinite(v)) or np.any(v < 0):
            raise ValueError("variances must be finite, nonnegative, shape (A, J)")
        w = n_a * v
        wsum = w.sum(axis=0)
        c2 = w / np.where(wsum > 0, wsum, 1.0)
    else:
        raise ValueError("scheme must be 'sample_size' or 'ivw'")

    c_sqrt_n = np.sqrt(c2) * np.sqrt(n_a)                  # c_aj sqrt(N_a)
    omega = c_sqrt_n[:, None, :] * c_sqrt_n[None, :, :] / safe
    omega[:, :, n_total == 0] = 0.0
    return omega, n_total


def _pair_factors(A: int, rho: ArrayLike | None, h2_ratio: ArrayLike | None) -> FloatArray:
    R = np.ones((A, A)) if rho is None else np.asarray(rho, dtype=np.float64)
    if R.shape != (A, A) or not np.all(np.isfinite(R)):
        raise ValueError("rho must be a finite (A, A) matrix")
    if not np.allclose(R, R.T) or not np.allclose(np.diag(R), 1.0) or np.any(np.abs(R) > 1.0):
        raise ValueError("rho must be symmetric with unit diagonal and |rho| <= 1")
    # A correlation matrix must be positive semidefinite. Allow roundoff on
    # singular matrices such as the all-ones shared-effect model.
    eigenvalues = np.linalg.eigvalsh((R + R.T) * 0.5)
    if eigenvalues[0] < -1e-10 * max(1.0, float(eigenvalues[-1])):
        raise ValueError("rho must be positive semidefinite")
    eta = np.ones(A) if h2_ratio is None else np.asarray(h2_ratio, dtype=np.float64)
    if eta.shape != (A,) or not np.all(np.isfinite(eta)) or np.any(eta <= 0):
        raise ValueError("h2_ratio must be A positive finite values")
    root = np.sqrt(eta)
    return R * root[:, None] * root[None, :]


def meta_ld_scores(
    scores: CrossAncestryScores | ArrayLike,
    omega: ArrayLike,
    *,
    rho: ArrayLike | None = None,
    h2_ratio: ArrayLike | None = None,
) -> FloatArray:
    """The LD score a meta-analysed GWAS actually follows (MA2).

    ``L_j = sum_{a,b} omega_{ab,j} rho_ab sqrt(eta_a eta_b) ell^{ab}_j``.
    ``rho`` defaults to all ones (identical standardized effects);
    ``rho = I`` returns the within-ancestry part alone, and the difference
    from ``rho = ones`` is the cross-ancestry part (MA9).
    """

    S = scores.scores if isinstance(scores, CrossAncestryScores) else np.asarray(scores, dtype=np.float64)
    W = np.asarray(omega, dtype=np.float64)
    if S.ndim != 3 or S.shape[0] != S.shape[1] or W.shape != S.shape:
        raise ValueError("scores and omega must both have shape (A, A, J)")
    F = _pair_factors(S.shape[0], rho, h2_ratio)
    return np.einsum("abj,abj,ab->j", W, S, F)


def single_reference_bias(
    true_scores: ArrayLike,
    reference_scores: ArrayLike,
    *,
    weights: ArrayLike | None = None,
    sample_sizes: ArrayLike | None = None,
) -> ReferenceBias:
    """Slope ratio and intercept shift from regressing on the wrong LD score (MA6).

    If ``E[chi2_j] = 1 + s * (N_j / Nbar) * true_scores_j`` and LDSC regresses
    on ``(N_j / Nbar) * reference_scores_j`` (LDSC scales each LD score by the
    variant's sample size), the (weighted) population slope is ``kappa * s``
    and the intercept is ``1 + s * offset``. ``sample_sizes`` is the per-variant
    ``N_j`` a meta-analysis reports (``Nbar`` is its mean); omit it for a
    constant ``N``. The slope interpretation assumes ``true_scores`` really
    represents the GWAS expectation. If it is estimated from finite reference
    panels, the returned value is a plug-in sensitivity, not a calibrated
    estimate of the absolute heritability bias.
    """

    L = np.asarray(true_scores, dtype=np.float64)
    R = np.asarray(reference_scores, dtype=np.float64)
    if L.ndim != 1 or L.shape != R.shape or L.size < 2:
        raise ValueError("scores must be 1-D arrays of equal length >= 2")
    if not np.all(np.isfinite(L)) or not np.all(np.isfinite(R)):
        raise ValueError("scores must contain only finite values")
    if sample_sizes is not None:
        n = np.asarray(sample_sizes, dtype=np.float64)
        if n.shape != L.shape or not np.all(np.isfinite(n)) or np.any(n < 0) or n.mean() <= 0:
            raise ValueError("sample_sizes must be finite, nonnegative, one per variant")
        L, R = L * (n / n.mean()), R * (n / n.mean())
    w = np.ones_like(L) if weights is None else np.asarray(weights, dtype=np.float64)
    if w.shape != L.shape or not np.all(np.isfinite(w)) or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("weights must be finite, nonnegative and not all zero")
    w = w / w.sum()
    mL, mR = float(w @ L), float(w @ R)
    var = float(w @ (R - mR) ** 2)
    if var <= 0.0:
        raise ValueError("reference_scores has no weighted variance")
    kappa = float(w @ ((L - mL) * (R - mR))) / var
    return ReferenceBias(kappa=kappa, offset=mL - kappa * mR)


def pooled_ld_scores(
    genotypes: ArrayLike,
    groups: ArrayLike,
    positions_cm: ArrayLike,
    *,
    covariates: ArrayLike | None = None,
    window_cm: float = 1.0,
    regression: ArrayLike | None = None,
    rho: ArrayLike | None = None,
    h2_ratio: ArrayLike | None = None,
    block_size: int = 256,
) -> FloatArray:
    """The LD score of a pooled multi-ancestry GWAS, given its genotypes (MA11).

    ``L_j = sum_{a,b} f_a f_b rho_ab sqrt(eta_a eta_b) sum_{k in W_j} m^a_jk m^b_jk``
    with ``m^a_jk = (1 / (N_a s_ak)) sum_{i in a} xtilde_ij g_ik``, where
    ``xtilde_j`` is variant j residualized on the GWAS covariates and scaled to
    unit mean square over the pooled sample, and ``s_ak`` is variant k's
    standard deviation in group a. This is the exact conditional expectation
    given genotypes, with effects standardized within each group. It includes
    the structure channel: group mean differences left unadjusted by the
    covariates enter through ``sum_{i in a} xtilde_ij != 0``.

    No finite-sample correction is applied. Use it in-sample, or with a large
    composition-matched pooled panel. With group indicators among the
    covariates spanning exactly the group indicators it reproduces the
    inverse-variance meta score exactly (MA12). Additional covariates generally
    change the residualized genotypes and require matched meta adjustments.
    """

    G = _as_genotype_matrix(genotypes)
    N, m = G.shape
    g = np.asarray(groups)
    if g.shape != (N,):
        raise ValueError("groups must have one label per sample")
    labels = np.unique(g)
    A = labels.size
    positions = _check_positions(positions_cm, m, window_cm)
    idx = _check_regression(regression, m)
    F = _pair_factors(A, rho, h2_ratio)

    resid, _ = _residualize(G, covariates)
    unit, _ = _unit_columns(resid, G)
    xt = unit * np.sqrt(N)               # unit mean square over the pooled sample

    members = [np.flatnonzero(g == lab) for lab in labels]
    f = np.array([mem.size / N for mem in members])
    scaled = []                          # g_ik / (N_a s_ak), zero if s_ak = 0
    for mem in members:
        Ga = G[mem]
        sd = Ga.std(axis=0)
        inv = np.zeros(m)
        live = sd > 0
        inv[live] = 1.0 / (mem.size * sd[live])
        scaled.append(Ga * inv)

    J = idx.size
    out = np.zeros(J)
    lo_all, hi_all = _window_bounds(positions, idx, window_cm)
    for start in range(0, J, block_size):
        stop = min(start + block_size, J)
        cols = idx[start:stop]
        lo, hi = int(lo_all[start]), int(hi_all[stop - 1])
        rows = np.arange(lo, hi)
        inwin = ((rows[:, None] >= lo_all[None, start:stop])
                 & (rows[:, None] < hi_all[None, start:stop]))
        Mom = [(scaled[a][:, lo:hi].T @ xt[members[a]][:, cols]) * inwin
               for a in range(A)]
        for a in range(A):
            for b in range(A):
                out[start:stop] += (f[a] * f[b] * F[a, b]
                                    * np.einsum("kj,kj->j", Mom[a], Mom[b]))
    return out
