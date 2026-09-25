"""Genome-wide LOCO LD-score engine (numerical primitives).

Shared infrastructure for the real-data survey (`docs/REAL_DATA_SURVEY.md`) and
the two candidate notes (`docs/WHEN_STRUCTURE_MATTERS.md`, `docs/NOTE_OUTLINE.md`).
The `experiments/realdata_1000g_benchmarks.py` suite established the operators on
three thinned chromosomes with an in-chromosome / two-chromosome-LOCO bracket and
physical-kb windows. Scaling to all 22 autosomes needs three things that brute
force cannot supply, all provided here and unit-tested against the brute-force
definitions:

1. A **genome-wide LD score without the M x M correlation matrix.** The total
   per-SNP squared-correlation sum is exactly the GRM quadratic form
   (METHOD.md Section 14, Eq. D2):

       sum_k r_jk^2 = (M / N) sum_p lambda_p corr(xtilde_j, u_p)^2,

   where xtilde_k = x_k / ||x_k|| are unit-norm columns and (lambda_p, u_p) are
   the eigenpairs of K = Xtilde Xtilde^T / M. This is O(M N^2), not O(M^2), so it
   is feasible at millions of SNPs. `spectral_full_scores` computes it and is
   validated to reproduce `brute_full_scores` to floating-point tolerance.

2. **Genetic-distance (cM) windows** instead of physical kb (`windowed_scores`),
   fed by `interpolate_cm` from a plink-format genetic map.

3. **LOCO whitening operators** built from the other chromosomes' GRM
   (`grm_eigh`, `pc_operator`, `hard_projection`), matching the definitions used
   throughout the benchmarks.

Only numpy. VCF parsing and cohort bookkeeping stay in the experiment driver;
this module is the tested kernel.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatMatrix = NDArray[np.float64]
FloatVector = NDArray[np.float64]


# ------------------------------------------------------------------ scaling

def standardize(genotypes: ArrayLike) -> FloatMatrix:
    """HWE-standardize a 0/1/2 genotype matrix (samples x variants).

    Column ``j`` becomes ``(g_j - 2 p_j) / sqrt(2 p_j (1 - p_j))`` with
    ``p_j = mean(g_j) / 2``. Monomorphic columns (zero variance) are left
    centered-only to avoid division by zero; callers should mask them out.
    """
    G = np.asarray(genotypes, dtype=np.float64)
    if G.ndim != 2:
        raise ValueError("genotypes must be a 2-D samples-by-variants matrix")
    p = G.mean(axis=0) / 2.0
    scale = np.sqrt(2.0 * p * (1.0 - p))
    scale = np.where(scale == 0.0, 1.0, scale)
    return (G - 2.0 * p) / scale


def unit_columns(X: ArrayLike) -> tuple[FloatMatrix, NDArray[np.bool_]]:
    """Scale each column to unit L2 norm; return (scaled, dead_mask).

    ``dead_mask`` flags numerically zero-norm columns, which are returned
    unchanged (norm forced to 1) so downstream code can exclude them. Squared
    correlations are then plain column dot products.
    """
    X = np.asarray(X, dtype=np.float64)
    norms = np.linalg.norm(X, axis=0)
    tol = np.finfo(np.float64).eps * max(X.shape) * max(float(norms.max()), 1.0)
    dead = norms <= tol
    return X / np.where(dead, 1.0, norms), dead


# --------------------------------------------------------------- operators

def grm_eigh(X: ArrayLike) -> tuple[FloatVector, FloatMatrix]:
    """Eigendecomposition of the GRM ``K = X X^T / M`` (ascending eigenvalues).

    ``X`` is a standardized samples-by-variants matrix. Returns ``(w, U)`` from
    ``numpy.linalg.eigh`` (eigenvalues ascending, so the leading structure axes
    are the last columns of ``U``).
    """
    X = np.asarray(X, dtype=np.float64)
    K = (X @ X.T) / X.shape[1]
    return np.linalg.eigh(K)


def pc_operator(w: ArrayLike, U: ArrayLike, h2: float) -> FloatMatrix:
    """Intercept-only mixed-model residual-precision operator ``P_c``.

    Built from a GRM eigendecomposition ``(w, U)`` under ``V = h2 K + (1-h2) I``:
    ``P_c = V^{-1} - V^{-1} 1 (1^T V^{-1} 1)^{-1} 1^T V^{-1}`` (METHOD.md Eq. M4,
    intercept-only design). Matches ``pcov`` in the benchmark suite exactly.
    """
    w = np.asarray(w, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    if not 0.0 < h2 < 1.0:
        raise ValueError("h2 must lie strictly between 0 and 1")
    d = 1.0 / (h2 * w + (1.0 - h2))
    Vinv = (U * d) @ U.T
    v1 = Vinv @ np.ones(U.shape[0])
    return Vinv - np.outer(v1, v1) / v1.sum()


def hard_projection(X: ArrayLike, U: ArrayLike, k: int) -> FloatMatrix:
    """cov-LDSC-style residualization off the top-``k`` GRM eigenvectors."""
    X = np.asarray(X, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    if k <= 0:
        return X.copy()
    Uk = U[:, ::-1][:, :k]
    return X - Uk @ (Uk.T @ X)


# ------------------------------------------------------------- LD scores

def brute_full_scores(Xn: ArrayLike, block: int = 1000) -> FloatVector:
    """Genome-wide per-SNP sum of squared correlations, by brute force.

    ``Xn`` must have unit-norm columns (see :func:`unit_columns`). O(M^2);
    used only to *validate* :func:`spectral_full_scores` on small inputs, and as
    the definition the spectral form must reproduce.
    """
    Xn = np.asarray(Xn, dtype=np.float64)
    M = Xn.shape[1]
    out = np.empty(M)
    for lo in range(0, M, block):
        hi = min(M, lo + block)
        out[lo:hi] = ((Xn[:, lo:hi].T @ Xn) ** 2).sum(axis=1)
    return out


def spectral_full_scores(Xn: ArrayLike) -> FloatVector:
    """Genome-wide per-SNP sum of squared correlations via the GRM spectrum.

    Exact evaluation of ``sum_k r_jk^2`` using the identity (METHOD.md Eq. D2)

        sum_k r_jk^2 = xtilde_j^T (Xtilde Xtilde^T) xtilde_j
                     = || Xtilde^T xtilde_j ||^2,

    for unit-norm columns ``xtilde``. Costs O(M N^2) with no M x M matrix, so it
    scales to millions of SNPs where :func:`brute_full_scores` cannot. The two
    agree to floating-point tolerance (see ``tests/test_genome.py``).
    """
    Xn = np.asarray(Xn, dtype=np.float64)
    # score_j = ||Xtilde^T xtilde_j||^2 = xtilde_j^T (Xn Xn^T) xtilde_j.
    # Form the N x N sample Gram once, then a per-column quadratic form: this
    # stays O(M N^2) time and O(N^2) memory, never materializing M x M.
    S = Xn @ Xn.T                       # N x N
    return np.einsum("ij,ij->j", Xn, S @ Xn)


def windowed_scores(
    Xn: ArrayLike,
    coord: ArrayLike,
    window: float,
    *,
    include_self: bool = True,
) -> FloatVector:
    """Per-SNP sum of squared correlations within +/- ``window`` of ``coord``.

    ``Xn`` has unit-norm columns; ``coord`` is a per-variant, nondecreasing
    coordinate (physical bp or genetic cM) and ``window`` is in the same units.
    Banded: O(M * variants-per-window * N), feasible genome-wide. Variants are
    assumed sorted by ``coord`` (as they come off a single chromosome).

    Unlike :func:`pldsc.ld.local_ld_scores`, correlations are not clipped to
    [-1, 1] and the self term is not pinned to exactly 1: with unit-norm columns
    both effects are roundoff-level, and the spectral identity validation
    depends on this staying the raw quadratic form.
    """
    Xn = np.asarray(Xn, dtype=np.float64)
    coord = np.asarray(coord, dtype=np.float64)
    M = Xn.shape[1]
    if coord.shape != (M,):
        raise ValueError("coord must have one value per variant")
    if np.any(np.diff(coord) < 0):
        raise ValueError("coord must be nondecreasing")
    out = np.zeros(M)
    for j in range(M):
        lo = np.searchsorted(coord, coord[j] - window, side="left")
        hi = np.searchsorted(coord, coord[j] + window, side="right")
        r2 = (Xn[:, lo:hi].T @ Xn[:, j]) ** 2
        if not include_self:
            r2[j - lo] = 0.0
        out[j] = r2.sum()
    return out


# ---------------------------------------------------------- structure share

def structure_share(ordinary: ArrayLike, whitened: ArrayLike) -> float:
    """Fraction of the mean LD score removed by a whitening operator.

    ``(mean(ordinary) - mean(whitened)) / mean(ordinary)``. Returns the share as
    a fraction (not a percent). Both inputs are per-SNP LD scores on the same
    variant set (ordinary standardized vs operator-whitened).
    """
    ordinary = np.asarray(ordinary, dtype=np.float64)
    whitened = np.asarray(whitened, dtype=np.float64)
    if ordinary.shape != whitened.shape:
        raise ValueError("ordinary and whitened must have the same length")
    fo = ordinary.mean()
    if fo == 0.0:
        raise ValueError("mean ordinary LD score is zero")
    return float((fo - whitened.mean()) / fo)


# --------------------------------------------------------------- genetic map

def read_plink_genetic_map(path: str) -> tuple[FloatVector, FloatVector]:
    """Read a plink-format genetic map: columns ``chrom . cM bp``.

    Returns ``(map_bp, map_cm)`` sorted by physical position, with exact
    duplicate bp collapsed (keeping the first cM) so :func:`interpolate_cm`'s
    strict-increase precondition holds. This is the Beagle ``plink.GRCh37.map``
    layout (one file per chromosome).
    """
    bp_list: list[float] = []
    cm_list: list[float] = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4 or parts[0].startswith("#"):
                continue
            cm_list.append(float(parts[2]))
            bp_list.append(float(parts[3]))
    bp = np.asarray(bp_list, dtype=np.float64)
    cm = np.asarray(cm_list, dtype=np.float64)
    order = np.argsort(bp, kind="stable")
    bp, cm = bp[order], cm[order]
    keep = np.concatenate(([True], np.diff(bp) > 0))
    return bp[keep], cm[keep]


def interpolate_cm(
    bp: ArrayLike, map_bp: ArrayLike, map_cm: ArrayLike
) -> FloatVector:
    """Interpolate genetic position (cM) at physical positions ``bp``.

    ``map_bp`` / ``map_cm`` are a plink-format genetic map for one chromosome
    (strictly increasing ``map_bp``). Linear interpolation inside the map;
    outside the map range the rate of the nearest flanking interval is
    extrapolated (constant-rate shoulders), never clamped flat, so distinct SNPs
    beyond the last map marker keep distinct cM and windows stay well defined.
    """
    bp = np.asarray(bp, dtype=np.float64)
    map_bp = np.asarray(map_bp, dtype=np.float64)
    map_cm = np.asarray(map_cm, dtype=np.float64)
    if map_bp.ndim != 1 or map_bp.shape != map_cm.shape:
        raise ValueError("map_bp and map_cm must be 1-D and the same length")
    if map_bp.size < 2:
        raise ValueError("genetic map needs at least two markers")
    if np.any(np.diff(map_bp) <= 0):
        raise ValueError("map_bp must be strictly increasing")

    cm = np.interp(bp, map_bp, map_cm)
    # np.interp clamps flat outside the range; replace with constant-rate slopes.
    left_rate = (map_cm[1] - map_cm[0]) / (map_bp[1] - map_bp[0])
    right_rate = (map_cm[-1] - map_cm[-2]) / (map_bp[-1] - map_bp[-2])
    lo = bp < map_bp[0]
    hi = bp > map_bp[-1]
    cm[lo] = map_cm[0] + (bp[lo] - map_bp[0]) * left_rate
    cm[hi] = map_cm[-1] + (bp[hi] - map_bp[-1]) * right_rate
    return cm
