import unittest

import numpy as np

from pldsc.regression import ld_score_regression


class LdScoreRegressionTests(unittest.TestCase):
    def test_recovers_intercept_and_slope_on_clean_line(self):
        x = np.linspace(0.0, 50.0, 200)
        y = 1.3 + 0.07 * x
        r = ld_score_regression(y, x)
        self.assertAlmostEqual(r.intercept, 1.3, places=8)
        self.assertAlmostEqual(r.slope, 0.07, places=8)
        # a noiseless line has (near) zero jackknife SE
        self.assertLess(r.intercept_se, 1e-6)
        self.assertLess(r.slope_se, 1e-6)

    def test_ols_matches_numpy_polyfit(self):
        rng = np.random.default_rng(np.random.PCG64(10))
        x = rng.uniform(0, 30, 300)
        y = 1.0 + 0.5 * x + rng.standard_normal(300)
        r = ld_score_regression(y, x)
        slope, intercept = np.polyfit(x, y, 1)
        self.assertAlmostEqual(r.intercept, intercept, places=8)
        self.assertAlmostEqual(r.slope, slope, places=8)

    def test_weighted_matches_direct_wls(self):
        rng = np.random.default_rng(np.random.PCG64(11))
        x = rng.uniform(0, 10, 150)
        y = 2.0 - 0.3 * x + rng.standard_normal(150)
        w = rng.uniform(0.5, 2.0, 150)
        r = ld_score_regression(y, x, weights=w)
        X = np.column_stack((np.ones_like(x), x))
        beta = np.linalg.solve((X.T * w) @ X, (X.T * w) @ y)
        self.assertAlmostEqual(r.intercept, beta[0], places=8)
        self.assertAlmostEqual(r.slope, beta[1], places=8)

    def test_jackknife_se_tracks_analytic_ols_se(self):
        # For iid-noise OLS the delete-1 jackknife SE should match the classic
        # sqrt(sigma^2 (X^T X)^-1) SE to within a modest factor.
        rng = np.random.default_rng(np.random.PCG64(12))
        n = 400
        x = rng.uniform(0, 20, n)
        sigma = 2.0
        y = 1.0 + 0.4 * x + sigma * rng.standard_normal(n)
        r = ld_score_regression(y, x, n_blocks=n)  # delete-1
        X = np.column_stack((np.ones_like(x), x))
        resid = y - (r.intercept + r.slope * x)
        s2 = np.sum(resid ** 2) / (n - 2)
        cov = s2 * np.linalg.inv(X.T @ X)
        analytic = np.sqrt(np.diag(cov))
        self.assertTrue(0.7 < r.intercept_se / analytic[0] < 1.4)
        self.assertTrue(0.7 < r.slope_se / analytic[1] < 1.4)

    def test_two_step_weighting_runs_and_downweights_high_ld(self):
        rng = np.random.default_rng(np.random.PCG64(13))
        x = rng.uniform(0, 40, 300)
        y = 1.0 + 0.1 * x + rng.standard_normal(300) * (1 + 0.1 * x)
        r = ld_score_regression(y, x, two_step=True)
        self.assertIn("two-step", r.weighting)
        self.assertTrue(np.isfinite(r.intercept) and np.isfinite(r.slope))

    def test_rejects_bad_input(self):
        with self.assertRaisesRegex(ValueError, "equal length"):
            ld_score_regression(np.arange(5.0), np.arange(4.0))
        with self.assertRaisesRegex(ValueError, "at least 4"):
            ld_score_regression(np.arange(3.0), np.arange(3.0))


if __name__ == "__main__":
    unittest.main()
