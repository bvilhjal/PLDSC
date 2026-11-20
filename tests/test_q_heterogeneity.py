"""Checks of the Cochran's Q simulation: exact scores, the rotation and the stacked fit."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from q_heterogeneity_sim import BLOCK, draw_z, fit, jk_se, ld_scores  # noqa: E402


def panels(rng, blocks=3, n=(60, 75)):
    out = []
    for size in n:
        X = rng.normal(size=(blocks, size, BLOCK)).cumsum(axis=2)        # correlated neighbours
        X = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
        out.append(X.astype(np.float32))
    return out


def test_scores_match_direct_block_correlations():
    XE, XF = panels(np.random.default_rng(3))
    lEE, lFF, lEF = ld_scores(XE, XF, chunk=2)
    RE = np.einsum("bni,bnj->bij", XE, XE, dtype=np.float64) / XE.shape[1]
    RF = np.einsum("bni,bnj->bij", XF, XF, dtype=np.float64) / XF.shape[1]
    np.testing.assert_allclose(lEE, (RE ** 2).sum(-1).ravel(), rtol=1e-5)
    np.testing.assert_allclose(lFF, (RF ** 2).sum(-1).ravel(), rtol=1e-5)
    np.testing.assert_allclose(lEF, (RE * RF).sum(-1).ravel(), rtol=1e-5, atol=1e-5)
    corrected = ld_scores(XE, XF, correct=True)[0]
    np.testing.assert_allclose(corrected, (RE ** 2 - (1 - RE ** 2) / (XE.shape[1] - 2)).sum(-1).ravel(), rtol=1e-5)


def test_association_statistics_have_the_block_ld_as_covariance():
    rng = np.random.default_rng(5)
    X = panels(rng, blocks=1, n=(400,))[0][:, :, :6].repeat(BLOCK // 6, axis=2)[:, :, :BLOCK]
    X = np.ascontiguousarray(np.concatenate([X, X[:, :, :BLOCK - X.shape[2]]], axis=2))
    R = (X[0].T @ X[0]).astype(np.float64) / X.shape[1]
    beta = np.zeros((1, BLOCK), dtype=np.float32)
    beta[0, 0] = 0.02
    z = np.array([draw_z(X, beta, 5_000.0, rng) for _ in range(4_000)])
    np.testing.assert_allclose(z.mean(0)[:6], np.sqrt(5_000.0) * R[:6, 0] * 0.02, atol=0.06)
    np.testing.assert_allclose(np.cov(z[:, :6].T), R[:6, :6], atol=0.08)


def test_pooled_statistic_and_q_are_a_rotation_of_the_strata():
    rng = np.random.default_rng(7)
    zE, zF = rng.normal(size=(2, 500))
    seE, seF = rng.uniform(0.01, 0.2, (2, 500))
    wE, wF = 1 / seE ** 2, 1 / seF ** 2
    cE, cF = np.sqrt(wE / (wE + wF)), np.sqrt(wF / (wE + wF))
    Z, Q = cE * zE + cF * zF, (cF * zE - cE * zF) ** 2
    np.testing.assert_allclose(Z ** 2 + Q, zE ** 2 + zF ** 2)
    bE, bF = zE * seE, zF * seF                                  # the inverse-variance meta-analysis itself
    pooled = (wE * bE + wF * bF) / (wE + wF)
    np.testing.assert_allclose(Z, pooled * np.sqrt(wE + wF))
    np.testing.assert_allclose(Q, wE * (bE - pooled) ** 2 + wF * (bF - pooled) ** 2)


def test_stacked_fit_recovers_exact_moments_and_deletes_whole_groups():
    rng = np.random.default_rng(11)
    n, starts = 400, np.arange(0, 400, 40)
    cE = np.sqrt(rng.uniform(0.3, 0.9, n))
    cF = np.sqrt(1 - cE ** 2)
    xE, xF, xX = rng.uniform(1, 9, (3, n))
    one, zero = np.ones(n), np.zeros(n)
    XZ = np.column_stack((one, zero, cE ** 2 * xE, cF ** 2 * xF, 2 * cE * cF * xX))
    XQ = np.column_stack((zero, one, cF ** 2 * xE, cE ** 2 * xF, -2 * cE * cF * xX))
    truth = np.array([1.0, 1.1, 0.3, 0.2, 0.15])
    w = rng.uniform(0.5, 2, n)
    full, dele = fit([(XZ, XZ @ truth, w), (XQ, XQ @ truth, w)], starts)
    np.testing.assert_allclose(full, truth, atol=1e-9)
    np.testing.assert_allclose(dele, np.tile(truth, (starts.size, 1)), atol=1e-8)
    noisy = XZ @ truth + rng.normal(size=n)
    _, dele = fit([(XZ[:, [0, 2, 3, 4]], noisy, w)], starts)
    keep = np.arange(n) >= 40
    A = XZ[keep][:, [0, 2, 3, 4]] * np.sqrt(w[keep, None])
    direct = np.linalg.lstsq(A, noisy[keep] * np.sqrt(w[keep]), rcond=None)[0]
    np.testing.assert_allclose(dele[0], direct, atol=1e-8)
    assert jk_se(dele).shape == (4,)
