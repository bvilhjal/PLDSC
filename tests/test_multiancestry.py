import unittest

import numpy as np

from pldsc import local_ld_scores
from pldsc.multiancestry import (
    CrossAncestryScores,
    cross_ancestry_ld_scores,
    meta_ld_scores,
    meta_weights,
    pooled_ld_scores,
    single_reference_bias,
)


def _panel(rng, n, m, rho, freq):
    """AR(1)-latent haplotypes -> diploid genotypes with local LD."""
    z = np.empty((2 * n, m))
    z[:, 0] = rng.standard_normal(2 * n)
    for k in range(1, m):
        z[:, k] = rho * z[:, k - 1] + np.sqrt(1 - rho * rho) * rng.standard_normal(2 * n)
    from pldsc.simulate import ndtri
    h = (z < ndtri(np.tile(freq, (2 * n, 1)))).astype(float)
    return h[:n] + h[n:]


def _std(X):
    """In-sample standardization (ddof = 0): columns mean 0, mean square 1."""
    return (X - X.mean(0)) / X.std(0)


class CrossAncestryScoreTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        self.m = 60
        self.pos = np.cumsum(rng.uniform(0.05, 0.2, self.m))
        f1 = rng.uniform(0.15, 0.85, self.m)
        f2 = np.clip(f1 + rng.normal(0, 0.1, self.m), 0.05, 0.95)
        self.X1 = _panel(rng, 80, self.m, 0.8, f1)
        self.X2 = _panel(rng, 90, self.m, 0.6, f2)

    def test_single_panel_matches_local_ld_scores(self):
        got = cross_ancestry_ld_scores([self.X1], self.pos, window_cm=1.0)
        want = local_ld_scores(self.X1, self.pos, window_cm=1.0)
        np.testing.assert_allclose(got.scores[0, 0], want, rtol=1e-12, atol=1e-12)

    def test_regression_subset_and_blocking_do_not_change_scores(self):
        full = cross_ancestry_ld_scores([self.X1, self.X2], self.pos, window_cm=0.7)
        sub = np.arange(3, self.m, 4)
        part = cross_ancestry_ld_scores([self.X1, self.X2], self.pos, window_cm=0.7,
                                        regression=sub, block_size=3)
        np.testing.assert_allclose(part.scores, full.scores[:, :, sub], atol=1e-12)
        np.testing.assert_array_equal(part.index, sub)

    def test_symmetric_self_term_and_uncorrected_cross_terms(self):
        s = cross_ancestry_ld_scores([self.X1, self.X2], self.pos, window_cm=0.0)
        # A zero window keeps only the self term: 1 on and off the diagonal.
        np.testing.assert_allclose(s.scores, 1.0, atol=1e-12)
        s = cross_ancestry_ld_scores([self.X1, self.X2], self.pos, window_cm=1.0)
        np.testing.assert_allclose(s.scores[0, 1], s.scores[1, 0])
        # Duplicated panel: the cross term is the UNcorrected sum of r^2.
        d = cross_ancestry_ld_scores([self.X1, self.X1], self.pos, window_cm=1.0)
        raw = local_ld_scores(self.X1, self.pos, window_cm=1.0, finite_sample_correction=False)
        np.testing.assert_allclose(d.scores[0, 1], raw, atol=1e-12)
        self.assertTrue(np.all(d.scores[0, 0] < d.scores[0, 1]))

    def test_monomorphic_variant_contributes_nothing(self):
        X2 = self.X2.copy()
        X2[:, 7] = 2.0
        s = cross_ancestry_ld_scores([self.X1, X2], self.pos, window_cm=1.0)
        self.assertFalse(s.available[1, 7])
        self.assertTrue(s.available[0, 7])
        self.assertEqual(s.scores[1, 1, 7], 0.0)
        self.assertEqual(s.scores[0, 1, 7], 0.0)
        self.assertEqual(s.variances[1, 7], 0.0)
        # A neighbour's within-ancestry score gets no noise-floor charge for 7.
        keep = np.delete(np.arange(self.m), 7)
        ref = cross_ancestry_ld_scores([self.X2[:, keep]], self.pos[keep], window_cm=1.0)
        np.testing.assert_allclose(s.scores[1, 1, 6], ref.scores[0, 0, 6], atol=1e-12)

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            cross_ancestry_ld_scores([self.X1, self.X2[:, :-1]], self.pos)
        with self.assertRaises(ValueError):
            cross_ancestry_ld_scores([self.X1], self.pos, regression=[3, 2])
        with self.assertRaises(ValueError):
            cross_ancestry_ld_scores([self.X1], self.pos, pair_df=[0.0])


class MetaWeightTests(unittest.TestCase):
    def test_sample_size_weights_sum_to_one_and_follow_availability(self):
        avail = np.array([[True, True, False], [True, False, True]])
        om, n = meta_weights([300, 100], avail)
        np.testing.assert_allclose(n, [400, 300, 100])
        np.testing.assert_allclose(om.sum(axis=(0, 1)), 1.0)
        np.testing.assert_allclose(om[:, :, 0], [[0.5625, 0.1875], [0.1875, 0.0625]])
        np.testing.assert_allclose(om[:, :, 1], [[1, 0], [0, 0]])

    def test_ivw_equals_sample_size_at_equal_variances_and_is_smaller_otherwise(self):
        avail = np.ones((2, 2), bool)
        v = np.array([[0.3, 0.4], [0.3, 0.1]])
        ivw, _ = meta_weights([200, 200], avail, variances=v, scheme="ivw")
        ss, _ = meta_weights([200, 200], avail)
        np.testing.assert_allclose(ivw[:, :, 0], ss[:, :, 0])
        self.assertLess(ivw[:, :, 1].sum(), 1.0)       # Cauchy-Schwarz, MA4

    def test_single_ancestry_meta_reproduces_its_ld_score(self):
        S = np.random.default_rng(0).uniform(1, 5, (2, 2, 7))
        S = (S + S.transpose(1, 0, 2)) / 2
        om, _ = meta_weights([500, 0], np.ones((2, 7), bool))
        np.testing.assert_allclose(meta_ld_scores(S, om), S[0, 0])

    def test_rho_splits_within_and_cross_parts(self):
        S = np.random.default_rng(1).uniform(1, 5, (2, 2, 5))
        S = (S + S.transpose(1, 0, 2)) / 2
        om, _ = meta_weights([100, 300], np.ones((2, 5), bool))
        within = meta_ld_scores(S, om, rho=np.eye(2))
        full = meta_ld_scores(S, om)
        half = meta_ld_scores(S, om, rho=[[1, 0.5], [0.5, 1]])
        np.testing.assert_allclose(half, within + 0.5 * (full - within))
        with self.assertRaises(ValueError):
            meta_ld_scores(S, om, rho=[[1, 2], [2, 1]])


class TheoremTests(unittest.TestCase):
    """MA2 and MA12 as exact in-sample algebraic identities."""

    def setUp(self):
        rng = np.random.default_rng(5)
        self.m = 40
        self.pos = np.cumsum(rng.uniform(0.05, 0.2, self.m))
        f1 = rng.uniform(0.2, 0.8, self.m)
        f2 = np.clip(f1 + rng.normal(0, 0.15, self.m), 0.05, 0.95)
        self.G = [_panel(rng, 150, self.m, 0.8, f1), _panel(rng, 60, self.m, 0.5, f2)]
        self.N = np.array([150.0, 60.0])

    def test_meta_moment_equals_brute_force_covariance(self):
        tau = np.array([[1.0, 0.7], [0.7, 1.0]]) * 0.5 / self.m     # rho = 0.7, h2 = 0.5
        Xs = [_std(g) for g in self.G]
        # Cov(z_aj, z_bj) for genotype-fixed GWAS in each ancestry (independent samples).
        cov = np.zeros((2, 2, self.m))
        for a in range(2):
            for b in range(2):
                Ra = Xs[a].T @ Xs[a] / self.N[a]
                Rb = Xs[b].T @ Xs[b] / self.N[b]
                cov[a, b] = np.sqrt(self.N[a] * self.N[b]) * tau[a, b] * np.sum(Ra * Rb, axis=1)
        s = cross_ancestry_ld_scores(self.G, self.pos, window_cm=1e9,
                                     finite_sample_correction=False)
        for scheme in ("sample_size", "ivw"):
            om, n_tot = meta_weights(self.N, s.available, variances=s.variances, scheme=scheme)
            # z-weights from their definitions, independently of meta_weights.
            if scheme == "sample_size":
                c = np.sqrt(np.broadcast_to(self.N[:, None] / self.N.sum(), (2, self.m)))
            else:
                w = self.N[:, None] * np.stack([g.var(0) for g in self.G])
                c = np.sqrt(w / w.sum(0))
            brute = np.einsum("aj,bj,abj->j", c, c, cov)
            L = meta_ld_scores(s, om, rho=[[1, 0.7], [0.7, 1]])
            np.testing.assert_allclose(n_tot * 0.5 / self.m * L, brute, rtol=1e-10)

    def test_pooled_with_group_indicators_equals_ivw_meta(self):
        G = np.vstack(self.G)
        groups = np.repeat([0, 1], self.N.astype(int))
        indicator = (groups == 1).astype(float)
        rho = [[1, 0.6], [0.6, 1]]
        pooled = pooled_ld_scores(G, groups, self.pos, covariates=indicator,
                                  window_cm=0.8, rho=rho)
        s = cross_ancestry_ld_scores(self.G, self.pos, window_cm=0.8,
                                     finite_sample_correction=False)
        om, _ = meta_weights(self.N, s.available, variances=s.variances, scheme="ivw")
        np.testing.assert_allclose(pooled, meta_ld_scores(s, om, rho=rho), rtol=1e-10)

    def test_pooled_without_ancestry_adjustment_adds_structure(self):
        G = np.vstack(self.G)
        groups = np.repeat([0, 1], self.N.astype(int))
        adjusted = pooled_ld_scores(G, groups, self.pos, covariates=(groups == 1).astype(float),
                                    window_cm=0.8)
        naive = pooled_ld_scores(G, groups, self.pos, window_cm=0.8)
        self.assertGreater(naive.mean(), adjusted.mean())


class ReferenceBiasTests(unittest.TestCase):
    def test_exact_linear_relation(self):
        R = np.linspace(1, 9, 50)
        b = single_reference_bias(0.7 * R + 0.3, R)
        self.assertAlmostEqual(b.kappa, 0.7)
        self.assertAlmostEqual(b.offset, 0.3)
        b = single_reference_bias(R, R, weights=np.linspace(1, 2, 50))
        self.assertAlmostEqual(b.kappa, 1.0)
        self.assertAlmostEqual(b.offset, 0.0)

    def test_sample_sizes_scale_both_regressors(self):
        R = np.linspace(1, 9, 50)
        L = 0.7 * R + 0.3
        same = single_reference_bias(L, R, sample_sizes=np.full(50, 7.0))
        self.assertAlmostEqual(same.kappa, 0.7)
        n = np.where(np.arange(50) < 25, 1.0, 3.0)
        b = single_reference_bias(L, R, sample_sizes=n)
        x, y = R * n / n.mean(), L * n / n.mean()
        self.assertAlmostEqual(b.kappa, np.polyfit(x, y, 1)[0])

    def test_dataclass_passthrough(self):
        s = CrossAncestryScores(np.ones((1, 1, 3)), np.ones((1, 3), bool),
                                np.ones((1, 3)), np.arange(3))
        om, _ = meta_weights([10], s.available)
        np.testing.assert_allclose(meta_ld_scores(s, om), 1.0)


if __name__ == "__main__":
    unittest.main()
