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


# ------------------------------------------------- intercept forms (q_model.py)
from q_model import weights_vary, zq_design, zq_moments  # noqa: E402

MINI = Path(__file__).resolve().parents[1] / "examples" / "mini" / "chr22_ldscores.npz"


def _expected_fit(XZ, XQ, yZ, yQ, dW):
    """Weighted least squares on exact moments, with the LDSC weights at the truth."""
    wz, wq = 1.0 / (dW * yZ ** 2), 1.0 / (dW * yQ ** 2)
    A = (XZ * wz[:, None]).T @ XZ + (XQ * wq[:, None]).T @ XQ
    b = (XZ * wz[:, None]).T @ yZ + (XQ * wq[:, None]).T @ yQ
    return np.linalg.solve(A, b)


def _chr22_asthma_like():
    """chr22 EUR/EAS example scores with the GBMI asthma sizes and strata estimates."""
    d = np.load(MINI, allow_pickle=False)
    S, V, maf, av = d["scores"], d["variances"], d["maf"], d["available"]
    e, f = 0, 2                                                    # EUR, EAS
    sel = av[e] & av[f] & (maf[e] >= 0.05) & (maf[f] >= 0.05)
    lE, lF, lX = S[e, e][sel], S[f, f][sel], S[e, f][sel]
    NE, NF = 354909.05, 61243.12                                   # report/asthma_heterogeneity.json
    M = 1664534 * lE.mean() / 23.7949                              # match the asthma per-SNP signal
    x = (NE * lE / M, NF * lF / M, np.sqrt(NE * NF) * lX / M)
    return x, (NE * V[e][sel], NF * V[f][sel]), (NE, NF), np.maximum(lE, 1.0)


def test_shared_intercepts_leak_under_inverse_variance_weights_and_three_do_not():
    (xE, xF, xX), (wE, wF), _, dW = _chr22_asthma_like()
    assert xE.size == 6415
    cE, cF = np.sqrt(wE / (wE + wF)), np.sqrt(wF / (wE + wF))
    hE, hF, g = 0.1451, 0.1976, 0.1250
    yZ, yQ = zq_moments(cE, cF, xE, xF, xX, hE, hF, g, aE=1.076, aF=0.983)
    XZ, XQ, form = zq_design(cE, cF, xE, xF, xX, "shared")
    shared = _expected_fit(XZ, XQ, yZ, yQ, dW)
    assert form == "shared"
    np.testing.assert_allclose(shared[-1] / g - 1.0, -0.0936, atol=0.001)     # -9.4% of g
    np.testing.assert_allclose(shared[:2], [1.060, 0.997], atol=0.001)       # close to the asthma fit's 1.063, 0.996
    XZ, XQ, form = zq_design(cE, cF, xE, xF, xX, "three")
    three = _expected_fit(XZ, XQ, yZ, yQ, dW)
    assert form == "three"
    np.testing.assert_allclose(three, [1.076, 0.983, 0.0, hE, hF, g], atol=1e-9)
    # Equal intercepts: no leak in either form.
    yZ, yQ = zq_moments(cE, cF, xE, xF, xX, hE, hF, g)
    XZ, XQ, _ = zq_design(cE, cF, xE, xF, xX, "shared")
    np.testing.assert_allclose(_expected_fit(XZ, XQ, yZ, yQ, dW)[-1], g, rtol=1e-9)


def test_sample_size_weights_have_no_leak_and_three_falls_back_to_shared():
    (xE, xF, xX), _, (NE, NF), dW = _chr22_asthma_like()
    n = xE.size
    cE, cF = np.full(n, np.sqrt(NE / (NE + NF))), np.full(n, np.sqrt(NF / (NE + NF)))
    assert not weights_vary(cE)
    yZ, yQ = zq_moments(cE, cF, xE, xF, xX, 0.1451, 0.1976, 0.1250, aE=1.076, aF=0.983)
    XZ, XQ, form = zq_design(cE, cF, xE, xF, xX, "three")
    assert form == "shared" and XZ.shape[1] == 5
    est = _expected_fit(XZ, XQ, yZ, yQ, dW)
    np.testing.assert_allclose(est[-3:], [0.1451, 0.1976, 0.1250], rtol=1e-9)


def test_unequal_intercept_arm_runs_and_reports_both_forms():
    from q_heterogeneity_sim import part1_cell
    rng = np.random.default_rng(13)
    XE, XF = panels(rng, blocks=200, n=(24, 26))
    P = {"XE": XE, "XF": XF, "sdE": rng.uniform(0.2, 0.7, 200 * BLOCK), "sdF": rng.uniform(0.2, 0.7, 200 * BLOCK)}
    P["scores"] = ld_scores(XE, XF)
    P["N"] = 2_000.0
    out = part1_cell(P, 0.85, 0.74, "inf", 2, np.random.default_rng(1), intE=1.08)
    assert out["ZQ3_intercept_form"] == "three" and out["stratum_intercepts"] == [1.08, 1.0]
    for key in ("strata", "ZQ", "ZQ3", "Z"):
        assert np.isfinite(out[key]["g_mean"])
