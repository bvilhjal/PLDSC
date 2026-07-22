import unittest

import numpy as np

from pldsc._lmm import (
    ld_scores,
    null_inflation,
    oracle_moment,
    residual_precision,
    score_pieces,
)


def _spd(rng, n):
    """A well-conditioned symmetric positive-definite matrix."""
    A = rng.standard_normal((n, n))
    return A @ A.T + n * np.eye(n)


class ResidualPrecisionIdentityTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(np.random.PCG64(1))
        self.N = 30
        self.V = _spd(self.rng, self.N)
        self.C = np.column_stack(
            (np.ones(self.N), self.rng.standard_normal((self.N, 2)))
        )
        self.P = residual_precision(self.V, self.C)

    def test_annihilates_the_fixed_design(self):
        # P C = 0 and C^T P = 0 (property P1)
        np.testing.assert_allclose(self.P @ self.C, 0.0, rtol=0, atol=1e-10)
        np.testing.assert_allclose(self.C.T @ self.P, 0.0, rtol=0, atol=1e-10)

    def test_is_vc_idempotent(self):
        # P V P = P (property P2) -- this is what pins the intercept to 1
        np.testing.assert_allclose(
            self.P @ self.V @ self.P, self.P, rtol=0, atol=1e-9
        )

    def test_is_symmetric(self):
        np.testing.assert_allclose(self.P, self.P.T, rtol=0, atol=1e-12)

    def test_is_positive_semidefinite(self):
        w = np.linalg.eigvalsh(self.P)
        self.assertGreaterEqual(w.min(), -1e-10)

    def test_self_weights_positive_off_design_and_zero_on_design(self):
        # q_j > 0 for a generic variant; q_j ~ 0 for a variant in span(C)
        X = self.rng.standard_normal((self.N, 4))
        _, q = score_pieces(self.P, X)
        self.assertTrue(np.all(q > 1e-8))
        _, q_design = score_pieces(self.P, self.C[:, [1]])
        np.testing.assert_allclose(q_design, 0.0, rtol=0, atol=1e-10)

    def test_handles_rank_deficient_design(self):
        # a duplicated covariate column must not break P C = 0
        C = np.column_stack((self.C, self.C[:, 1]))
        P = residual_precision(self.V, C)
        np.testing.assert_allclose(P @ C, 0.0, rtol=0, atol=1e-9)


class ReductionTests(unittest.TestCase):
    """Check the limiting geometries under their explicit scale conditions."""

    def setUp(self):
        self.rng = np.random.default_rng(np.random.PCG64(2))
        self.N = 40

    def test_reduces_to_covariate_projection_when_V_is_scalar(self):
        # V = sigma^2 I  =>  P = sigma^{-2} (I - C (C^T C)^+ C^T)
        sigma2 = 3.0
        C = np.column_stack(
            (np.ones(self.N), self.rng.standard_normal((self.N, 2)))
        )
        P = residual_precision(sigma2 * np.eye(self.N), C)
        P_C = np.eye(self.N) - C @ np.linalg.pinv(C.T @ C) @ C.T
        np.testing.assert_allclose(P, P_C / sigma2, rtol=0, atol=1e-10)

    def test_ldscore_reduces_to_sum_of_squared_correlations(self):
        # intercept-only, V = I, centered unit-norm genotypes => q_k = 1 and
        # L_j = sum_k corr(x_j, x_k)^2 exactly.
        X = self.rng.standard_normal((self.N, 6))
        X = X - X.mean(0)
        X = X / np.linalg.norm(X, axis=0)
        one = np.ones((self.N, 1))
        P = residual_precision(np.eye(self.N), one)
        L = ld_scores(P, X)
        corr = X.T @ X  # unit-norm centered => Pearson correlations
        np.testing.assert_allclose(L, (corr ** 2).sum(1), rtol=0, atol=1e-10)


class MomentTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(np.random.PCG64(3))
        self.N = 35
        self.V = _spd(self.rng, self.N)
        self.C = np.ones((self.N, 1))
        self.P = residual_precision(self.V, self.C)
        self.X = self.rng.standard_normal((self.N, 5))

    def test_oracle_null_inflation_is_one(self):
        # with the true V, the null intercept is exactly 1 for every variant
        infl = null_inflation(self.P, self.X, self.V)
        np.testing.assert_allclose(infl, 1.0, rtol=0, atol=1e-9)

    def test_misspecified_V_moves_the_intercept_off_one(self):
        V_wrong = _spd(self.rng, self.N)
        infl = null_inflation(self.P, self.X, V_wrong)
        self.assertGreater(np.max(np.abs(infl - 1.0)), 1e-3)

    def test_oracle_moment_is_one_plus_tau2_L(self):
        tau2 = 0.4
        moment = oracle_moment(self.P, self.X, tau2)
        expected = 1.0 + tau2 * ld_scores(self.P, self.X)
        np.testing.assert_allclose(moment, expected, rtol=0, atol=1e-12)

    def test_score_covariance_allows_distinct_effect_and_test_sets(self):
        # METHOD M6/M16: tested variants need not equal modeled-effect variants.
        X_effect = self.rng.standard_normal((self.N, 3))
        X_test = self.rng.standard_normal((self.N, 4))
        effect_var = np.array([0.1, 0.2, 0.4])

        _, q_test = score_pieces(self.P, X_test)
        q_inv_sqrt = np.diag(1.0 / np.sqrt(q_test))
        cross = X_test.T @ self.P @ X_effect
        B = q_inv_sqrt @ cross
        R = q_inv_sqrt @ (X_test.T @ self.P @ X_test) @ q_inv_sqrt
        expected_cov = R + B @ np.diag(effect_var) @ B.T

        phenotype_cov = (
            self.V + X_effect @ np.diag(effect_var) @ X_effect.T
        )
        direct_cov = (
            q_inv_sqrt
            @ X_test.T
            @ self.P
            @ phenotype_cov
            @ self.P
            @ X_test
            @ q_inv_sqrt
        )
        np.testing.assert_allclose(
            direct_cov, expected_cov, rtol=0, atol=1e-9
        )

        expected_moment = 1.0 + (cross * cross) @ effect_var / q_test
        np.testing.assert_allclose(
            np.diag(expected_cov), expected_moment, rtol=0, atol=1e-10
        )

    def test_windowed_ldscore_is_at_most_full_sum(self):
        M = self.X.shape[1]
        mask = np.abs(np.subtract.outer(np.arange(M), np.arange(M))) <= 1
        full = ld_scores(self.P, self.X)
        win = ld_scores(self.P, self.X, window_mask=mask.astype(float))
        self.assertTrue(np.all(win <= full + 1e-12))


if __name__ == "__main__":
    unittest.main()
