"""Bounded-memory LD scores, adapted from LDpred3 (commit 2f5d108)."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ._numba import HAVE_NUMBA, _jit

_LD_SCORE_WORK_BYTES = 32 << 20
_Q8 = 127.0

@_jit
def _int8_sumsq(a, out):
    for i in range(a.shape[0]):
        total = 0
        for j in range(a.shape[1]):
            x = int(a[i, j])
            total += x * x
        out[i] = total


def _row_sumsq(a):
    a = np.asarray(a)
    if HAVE_NUMBA and a.dtype == np.int8 and a.flags.c_contiguous:
        out = np.empty(a.shape[0], dtype=float)
        _int8_sumsq(a, out)
        return out
    return np.einsum("ij,ij->i", a, a, dtype=np.float64)


@dataclass
class LowRankLD:
    """Factor for R = scale² U U.T + diagonal, whose diagonal is one.

    Float U uses scale 1; int8 U uses a positive quantization scale. This
    representation preserves full within-block LD without materializing R.
    Rows must have scaled squared norm at most one.
    """
    U: np.ndarray
    scale: float = 1.0

    def __post_init__(self):
        self.U = np.asarray(self.U)

    @property
    def m(self):
        return self.U.shape[0]

    @property
    def is_int8(self):
        return np.asarray(self.U).dtype == np.int8

def _lowrank_ld_scores(R):
    """Float64 contractions, tiled when widening the whole factor is costly."""
    m, rank = R.U.shape
    scale = float(R.scale) if R.is_int8 else 1.0
    if 8 * (2 * m * rank + rank * rank) <= _LD_SCORE_WORK_BYTES:
        U = np.asarray(R.U, dtype=np.float64)
        if R.is_int8:
            U *= scale
        gram = U.T @ U
        score = np.einsum("ij,ij->i", U @ gram, U)
        diag = np.einsum("ij,ij->i", U, U)
        score += 1.0 - diag ** 2
        return score

    # A score contracts U_j with every Gram tile. Off-diagonal rank tiles
    # contribute twice; the diagonal correction restores the exact self term.
    # 512-row/rank windows bound scratch independently of both block axes.
    score = np.zeros(m)
    for a in range(0, rank, 512):
        for b in range(a, rank, 512):
            gram = np.zeros((min(512, rank-a), min(512, rank-b)))
            for start in range(0, m, 512):
                left = np.asarray(R.U[start:start+512, a:a+512], dtype=float)
                right = (left if a == b else np.asarray(
                    R.U[start:start+512, b:b+512], dtype=float))
                gram += left.T @ right
            gram *= scale**2
            for start in range(0, m, 512):
                left = np.asarray(R.U[start:start+512, a:a+512], dtype=float)
                right = (left if a == b else np.asarray(
                    R.U[start:start+512, b:b+512], dtype=float))
                term = np.einsum("ij,ij->i", left @ gram, right)
                score[start:start+512] += term * (scale**2 if a == b else 2 * scale**2)
    diag = _row_sumsq(R.U) * scale**2
    score += 1.0 - diag**2
    return score



def ld_scores(blocks, *, n_ref=None):
    """LD scores from a dense LD matrix or indexed, disjoint LD blocks.

    Blocks are ``(matrix_or_LowRankLD, indices)`` pairs covering 0..M-1.
    Dense int8 blocks represent correlations scaled by 127; floating matrices
    contain correlations directly. A LowRankLD or LDpred3 LowRankLD factor is
    contracted without forming its dense matrix. Scores omit between-block LD.
    ``n_ref`` optionally applies the conventional finite-reference correction.
    Calculate scores on the full variant map before subsetting regression rows.
    """
    if n_ref is not None:
        if isinstance(n_ref, (bool, np.bool_)) or not isinstance(n_ref, (int, np.integer)) or n_ref <= 2:
            raise ValueError("n_ref must be an integer > 2")
    if isinstance(blocks, np.ndarray):
        if blocks.ndim != 2:
            raise ValueError("LD matrix must be two-dimensional")
        blocks = [(blocks, np.arange(blocks.shape[0]))]
    else:
        blocks = list(blocks)
    if any(not isinstance(b, (list, tuple)) or len(b) != 2 for b in blocks):
        raise ValueError("blocks must contain (LD, index) pairs")
    m = sum(np.asarray(idx).size for _, idx in blocks)
    if not m:
        raise ValueError("at least one variant is required")
    result = np.empty(m)
    seen = np.zeros(m, dtype=bool)
    for R, index in blocks:
        idx = np.asarray(index)
        if idx.ndim != 1 or idx.dtype.kind not in "iu" or not idx.size:
            raise ValueError("block indices must be nonempty integer vectors")
        if np.any(idx < 0) or np.any(idx >= m) or np.unique(idx).size != idx.size or np.any(seen[idx]):
            raise ValueError("block indices must tile 0..M-1 without duplicates")
        seen[idx] = True
        k = len(idx)
        if hasattr(R, "U") and hasattr(R, "scale") and hasattr(R, "is_int8"):
            U = np.asarray(R.U)
            scale = float(R.scale)
            if U.ndim != 2 or U.shape[0] != k or not 1 <= U.shape[1] <= k:
                raise ValueError("low-rank factor has invalid shape")
            if U.dtype != np.int8 and not np.issubdtype(U.dtype, np.floating):
                raise ValueError("low-rank factor must be floating point or int8")
            if not np.isfinite(scale) or scale <= 0 or (U.dtype != np.int8 and scale != 1.0):
                raise ValueError("invalid low-rank factor scale")
            if U.dtype == np.int8:
                if U.min() == -128:
                    raise ValueError("invalid low-rank factor values")
            else:
                rows = max(1, min(512, (8 << 20) // max(1, U.shape[1] * U.itemsize)))
                for start in range(0, k, rows):
                    if not np.all(np.isfinite(U[start:start + rows])):
                        raise ValueError("invalid low-rank factor values")
            if np.max(_row_sumsq(U) * scale * scale) > 1.0 + 2e-5:
                raise ValueError("low-rank factor has a row norm above one")
            score = _lowrank_ld_scores(R)
        else:
            arr = np.asarray(R)
            if arr.shape != (k, k) or (arr.dtype != np.int8 and not np.issubdtype(arr.dtype, np.floating)):
                raise ValueError("dense LD block must be a square float or int8 matrix")
            scale = 1.0 / _Q8 if arr.dtype == np.int8 else 1.0
            row_count = max(1, min(512, (8 << 20) // (8 * k)))
            for start in range(0, k, row_count):
                rows = arr[start:start + row_count]
                if (not np.all(np.isfinite(rows)) or np.max(np.abs(rows.astype(float) * scale)) > 1.0 + 2e-5
                        or not np.allclose(rows, arr[:, start:start + row_count].T, atol=2e-5 / scale, rtol=2e-5)):
                    raise ValueError("dense LD block must be finite, symmetric and bounded")
            if not np.allclose(np.diag(arr) * scale, 1.0, atol=2e-5, rtol=2e-5):
                raise ValueError("dense LD block must have unit diagonal")
            score = _row_sumsq(arr) * scale**2
        if n_ref is not None:
            score = score - (k - score) / (n_ref - 2.0)
        result[idx] = score
    if not np.all(seen):
        raise ValueError("block indices must tile 0..M-1")
    return result
