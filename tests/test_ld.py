import unittest

import numpy as np

from pldsc import local_ld_scores, residualize_genotypes


class ResidualizationTests(unittest.TestCase):
    def test_residuals_are_orthogonal_to_intercept_and_covariate(self):
        covariate = np.arange(8, dtype=float)
        noise = np.array([1, -1, 0, 0, 1, -1, 0, 0], dtype=float)
        genotypes = np.column_stack(
            (2 + 3 * covariate + noise, -1 + covariate - noise)
        )

        residuals = residualize_genotypes(genotypes, covariate)
        design = np.column_stack((np.ones(8), covariate))

        np.testing.assert_allclose(design.T @ residuals, 0.0, atol=1e-12)


class LocalLdScoreTests(unittest.TestCase):
    def setUp(self):
        self.genotypes = np.array(
            [
                [0, 0, 1, 2],
                [1, 0, 2, 1],
                [2, 1, 1, 0],
                [0, 2, 0, 1],
                [1, 2, 1, 0],
                [2, 1, 0, 2],
            ],
            dtype=float,
        )
        self.positions = np.array([0.0, 1.0, 2.0, 3.0])

    def test_scores_match_direct_intercept_only_calculation(self):
        correlations = np.corrcoef(self.genotypes, rowvar=False)
        squared = correlations**2
        adjusted = squared - (1.0 - squared) / (self.genotypes.shape[0] - 2)
        expected = np.array(
            [
                adjusted[j, np.abs(self.positions - self.positions[j]) <= 1.0].sum()
                for j in range(self.genotypes.shape[1])
            ]
        )

        observed = local_ld_scores(
            self.genotypes, self.positions, window_cm=1.0
        )

        np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_covariate_projection_removes_structure_induced_correlation(self):
        ancestry = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)
        first_noise = np.array([1, -1, 0, 0, 1, -1, 0, 0], dtype=float)
        second_noise = np.array([0, 0, 1, -1, 0, 0, 1, -1], dtype=float)
        genotypes = np.column_stack(
            (ancestry + first_noise, ancestry + second_noise)
        )
        positions = np.array([0.0, 0.5])

        ordinary = local_ld_scores(
            genotypes,
            positions,
            finite_sample_correction=False,
            include_self=False,
        )
        adjusted = local_ld_scores(
            genotypes,
            positions,
            covariates=ancestry,
            finite_sample_correction=False,
            include_self=False,
        )

        self.assertTrue(np.all(ordinary > 0.0))
        np.testing.assert_allclose(adjusted, 0.0, atol=1e-12)

    def test_finite_sample_correction_uses_residual_rank(self):
        ancestry = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)
        first_noise = np.array([1, -1, 0, 0, 1, -1, 0, 0], dtype=float)
        second_noise = np.array([0, 0, 1, -1, 0, 0, 1, -1], dtype=float)
        genotypes = np.column_stack(
            (ancestry + first_noise, ancestry + second_noise)
        )

        observed = local_ld_scores(
            genotypes,
            [0.0, 0.5],
            covariates=ancestry,
            include_self=False,
        )

        # N=8 and rank([intercept, ancestry])=2, so nu=6 and the
        # adjustment of an r-squared value of zero is -1/(nu-1).
        np.testing.assert_allclose(observed, -1.0 / 5.0, atol=1e-12)

    def test_rejects_unsorted_positions_and_fully_explained_variants(self):
        with self.assertRaisesRegex(ValueError, "nondecreasing"):
            local_ld_scores(self.genotypes, [0.0, 2.0, 1.0, 3.0])

        ancestry = np.array([-1, -1, 1, 1], dtype=float)
        with self.assertRaisesRegex(ValueError, "zero residual variance"):
            local_ld_scores(
                np.column_stack((ancestry, [0.0, 1.0, 0.0, 1.0])),
                [0.0, 0.5],
                covariates=ancestry,
            )


if __name__ == "__main__":
    unittest.main()
