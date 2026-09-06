import unittest

import numpy as np

from pldsc import effective_pair_df, local_ld_scores, residualize_genotypes


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

        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-12)

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
        np.testing.assert_allclose(adjusted, 0.0, rtol=0, atol=1e-12)

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
        np.testing.assert_allclose(observed, -1.0 / 5.0, rtol=0, atol=1e-12)

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


class EffectivePairDfTests(unittest.TestCase):
    def test_unrelated_panel_reduces_to_n_minus_one(self):
        # A = I with an intercept only: Atilde = P, tr(P) = tr(P^2) = N - 1.
        self.assertAlmostEqual(effective_pair_df(np.eye(40)), 39.0, places=9)

    def test_extra_covariates_lower_the_df_by_their_rank(self):
        rng = np.random.default_rng(np.random.PCG64(11))
        covariates = rng.standard_normal((40, 3))
        self.assertAlmostEqual(
            effective_pair_df(np.eye(40), covariates), 36.0, places=9
        )

    def test_full_sib_families_shrink_the_df(self):
        # 10 families of 4 full sibs: A has 1 on the diagonal, 0.5 between sibs.
        block = 0.5 * np.ones((4, 4))
        np.fill_diagonal(block, 1.0)
        A = np.kron(np.eye(10), block)

        df = effective_pair_df(A)

        # Strictly below the unrelated value, and matching the closed form
        # tr(Atilde)^2 / tr(Atilde^2) computed densely.
        P = np.eye(40) - np.ones((40, 40)) / 40
        At = P @ A @ P
        expected = np.trace(At) ** 2 / np.trace(At @ At)
        self.assertLess(df, 39.0)
        self.assertAlmostEqual(df, float(expected), places=9)

    def test_rejects_non_square_and_non_finite_input(self):
        with self.assertRaisesRegex(ValueError, "square"):
            effective_pair_df(np.ones((4, 5)))
        bad = np.eye(5)
        bad[0, 0] = np.inf
        with self.assertRaisesRegex(ValueError, "finite"):
            effective_pair_df(bad)

    def test_rejects_asymmetric_relationship(self):
        bad = np.eye(6)
        bad[0, 1] = 0.5
        with self.assertRaisesRegex(ValueError, "symmetric"):
            effective_pair_df(bad)


class PairDfOverrideTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(np.random.PCG64(12))
        self.G = rng.integers(0, 3, size=(30, 6)).astype(float)
        self.pos = np.sort(rng.uniform(0.0, 5.0, 6))

    def test_override_replaces_the_residual_df_denominator(self):
        correlations = np.corrcoef(self.G, rowvar=False)
        squared = correlations**2
        adjusted = squared - (1.0 - squared) / 7.5
        expected = np.array(
            [
                adjusted[j, np.abs(self.pos - self.pos[j]) <= 2.0].sum()
                for j in range(self.G.shape[1])
            ]
        )

        observed = local_ld_scores(
            self.G, self.pos, window_cm=2.0, pair_df=7.5
        )

        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-12)

    def test_default_matches_an_explicit_residual_df_minus_one(self):
        default = local_ld_scores(self.G, self.pos, window_cm=2.0)
        explicit = local_ld_scores(
            self.G, self.pos, window_cm=2.0, pair_df=self.G.shape[0] - 2
        )
        np.testing.assert_allclose(default, explicit, rtol=0, atol=1e-12)

    def test_smaller_pair_df_subtracts_more(self):
        # Window spans the whole chromosome, so every variant has neighbours
        # and the extra subtraction is strict (the self-term is exactly 1 and
        # contributes no correction either way).
        loose = local_ld_scores(self.G, self.pos, window_cm=10.0, pair_df=100.0)
        tight = local_ld_scores(self.G, self.pos, window_cm=10.0, pair_df=10.0)
        self.assertTrue(np.all(tight < loose))

    def test_rejects_nonpositive_pair_df(self):
        with self.assertRaisesRegex(ValueError, "pair_df"):
            local_ld_scores(self.G, self.pos, pair_df=0.0)


class BaselineHardeningTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(np.random.PCG64(7))
        self.G = rng.integers(0, 3, size=(40, 12)).astype(float)
        self.pos = np.sort(rng.uniform(0.0, 10.0, 12))
        self.cov = rng.standard_normal((40, 3))

    def test_invariant_to_full_rank_covariate_reparameterization(self):
        # [1, C] and [1, C @ A] span the same subspace for invertible A, so the
        # residualized correlations -- and thus the LD scores -- must match.
        rng = np.random.default_rng(np.random.PCG64(8))
        A = rng.standard_normal((3, 3))
        A += 3 * np.eye(3)  # keep it well-conditioned / invertible
        s1 = local_ld_scores(self.G, self.pos, covariates=self.cov, window_cm=2.0)
        s2 = local_ld_scores(
            self.G, self.pos, covariates=self.cov @ A, window_cm=2.0
        )
        np.testing.assert_allclose(s1, s2, rtol=0, atol=1e-9)

    def test_accepts_multi_column_covariates(self):
        scores = local_ld_scores(
            self.G, self.pos, covariates=self.cov, window_cm=2.0
        )
        self.assertEqual(scores.shape, (12,))
        self.assertTrue(np.all(np.isfinite(scores)))

    def test_rejects_non_finite_genotypes(self):
        bad = self.G.copy()
        bad[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            local_ld_scores(bad, self.pos)

    def test_rejects_misshaped_covariates(self):
        with self.assertRaisesRegex(ValueError, "covariates"):
            local_ld_scores(self.G, self.pos, covariates=self.cov[:5])

    def test_rejects_too_few_samples(self):
        with self.assertRaisesRegex(ValueError, "at least 3 samples"):
            local_ld_scores(self.G[:2], self.pos)

    def test_rejects_position_count_mismatch(self):
        with self.assertRaisesRegex(ValueError, "one value per variant"):
            local_ld_scores(self.G, self.pos[:5])


if __name__ == "__main__":
    unittest.main()
