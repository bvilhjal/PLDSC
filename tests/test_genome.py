import os
import tempfile
import unittest

import numpy as np

from pldsc.genome import (
    brute_full_scores,
    grm_eigh,
    hard_projection,
    interpolate_cm,
    pc_operator,
    read_plink_genetic_map,
    spectral_full_scores,
    standardize,
    structure_share,
    unit_columns,
    windowed_scores,
)


def _toy_genotypes(seed=0, n=60, m=40):
    rng = np.random.default_rng(np.random.PCG64(seed))
    p = rng.uniform(0.1, 0.9, size=m)
    return (rng.random((n, m)) < p).astype(float) + (rng.random((n, m)) < p)


class StandardizeTests(unittest.TestCase):
    def test_columns_are_centered_and_hwe_scaled(self):
        G = _toy_genotypes()
        X = standardize(G)
        np.testing.assert_allclose(X.mean(axis=0), 0.0, atol=1e-12)
        # Column j must equal (g_j - 2 p_j) / sqrt(2 p_j (1 - p_j)) exactly.
        p = G.mean(axis=0) / 2.0
        expected = (G - 2.0 * p) / np.sqrt(2.0 * p * (1.0 - p))
        np.testing.assert_allclose(X, expected, rtol=0, atol=1e-12)

    def test_monomorphic_column_is_left_centered_not_nan(self):
        G = _toy_genotypes(seed=99, n=30, m=5)
        G[:, 2] = 2.0  # fixed homozygote -> zero variance
        X = standardize(G)
        self.assertTrue(np.all(np.isfinite(X)))
        np.testing.assert_allclose(X[:, 2], 0.0, atol=1e-12)


class SpectralFullScoreTests(unittest.TestCase):
    def test_spectral_matches_brute_force(self):
        # The whole point of the engine: the O(M N^2) spectral score equals the
        # O(M^2) brute-force sum of squared correlations, exactly.
        Xn, _ = unit_columns(standardize(_toy_genotypes(seed=1, n=50, m=70)))
        np.testing.assert_allclose(
            spectral_full_scores(Xn), brute_full_scores(Xn), rtol=0, atol=1e-9
        )

    def test_spectral_matches_brute_force_after_whitening(self):
        # It must also hold for whitened columns (the operator-adjusted score),
        # since that is how the structure share is computed genome-wide.
        X = standardize(_toy_genotypes(seed=2, n=48, m=64))
        w, U = grm_eigh(X)
        P = pc_operator(w, U, 0.5)
        Xn, _ = unit_columns(P @ X)
        np.testing.assert_allclose(
            spectral_full_scores(Xn), brute_full_scores(Xn), rtol=0, atol=1e-9
        )

    def test_self_term_is_one_per_snp(self):
        # Each SNP's own r^2 = 1 contributes exactly 1; with independent columns
        # the score should sit just above 1.
        rng = np.random.default_rng(np.random.PCG64(3))
        Xn, _ = unit_columns(rng.standard_normal((200, 5)))
        scores = spectral_full_scores(Xn)
        self.assertTrue(np.all(scores >= 1.0 - 1e-9))


class WindowedScoreTests(unittest.TestCase):
    def test_full_window_equals_genome_wide_score(self):
        Xn, _ = unit_columns(standardize(_toy_genotypes(seed=4, n=40, m=30)))
        coord = np.arange(30, dtype=float)
        windowed = windowed_scores(Xn, coord, window=100.0)  # spans everything
        np.testing.assert_allclose(
            windowed, brute_full_scores(Xn), rtol=0, atol=1e-9
        )

    def test_window_is_a_subset_of_the_full_score(self):
        Xn, _ = unit_columns(standardize(_toy_genotypes(seed=5, n=40, m=50)))
        coord = np.sort(np.random.default_rng(6).uniform(0, 10, 50))
        narrow = windowed_scores(Xn, coord, window=0.5)
        full = brute_full_scores(Xn)
        # Squared correlations are nonnegative, so a window never exceeds the
        # full sum and always includes the self-term.
        self.assertTrue(np.all(narrow <= full + 1e-9))
        self.assertTrue(np.all(narrow >= 1.0 - 1e-9))

    def test_include_self_toggle(self):
        Xn, _ = unit_columns(standardize(_toy_genotypes(seed=7, n=40, m=20)))
        coord = np.arange(20, dtype=float)
        with_self = windowed_scores(Xn, coord, 100.0, include_self=True)
        without = windowed_scores(Xn, coord, 100.0, include_self=False)
        np.testing.assert_allclose(with_self - without, 1.0, rtol=0, atol=1e-9)

    def test_rejects_unsorted_coord(self):
        Xn, _ = unit_columns(standardize(_toy_genotypes(seed=8, n=40, m=10)))
        with self.assertRaisesRegex(ValueError, "nondecreasing"):
            windowed_scores(Xn, np.arange(10)[::-1].astype(float), 1.0)


class OperatorTests(unittest.TestCase):
    def test_pc_operator_matches_benchmark_definition(self):
        X = standardize(_toy_genotypes(seed=9, n=40, m=55))
        w, U = grm_eigh(X)
        P = pc_operator(w, U, 0.5)
        # Reproduce the benchmark's pcov closed form independently.
        d = 1.0 / (0.5 * w + 0.5)
        Vinv = (U * d) @ U.T
        v1 = Vinv @ np.ones(40)
        expected = Vinv - np.outer(v1, v1) / v1.sum()
        np.testing.assert_allclose(P, expected, rtol=0, atol=1e-10)

    def test_pc_operator_annihilates_the_intercept(self):
        X = standardize(_toy_genotypes(seed=10, n=36, m=48))
        w, U = grm_eigh(X)
        P = pc_operator(w, U, 0.7)
        np.testing.assert_allclose(P @ np.ones(36), 0.0, atol=1e-9)
        np.testing.assert_allclose(np.ones(36) @ P, 0.0, atol=1e-9)

    def test_hard_projection_removes_top_eigenvectors(self):
        X = standardize(_toy_genotypes(seed=11, n=40, m=50))
        w, U = grm_eigh(X)
        Xp = hard_projection(X, U, 3)
        top = U[:, ::-1][:, :3]
        np.testing.assert_allclose(top.T @ Xp, 0.0, atol=1e-9)

    def test_pc_operator_rejects_degenerate_h2(self):
        X = standardize(_toy_genotypes(seed=12, n=30, m=30))
        w, U = grm_eigh(X)
        for bad in (0.0, 1.0, -0.1, 1.2):
            with self.assertRaises(ValueError):
                pc_operator(w, U, bad)


class StructureShareTests(unittest.TestCase):
    def test_no_whitening_is_zero_share(self):
        s = np.array([2.0, 3.0, 4.0])
        self.assertAlmostEqual(structure_share(s, s), 0.0, places=12)

    def test_half_removed_is_half_share(self):
        o = np.array([4.0, 4.0])
        w = np.array([2.0, 2.0])
        self.assertAlmostEqual(structure_share(o, w), 0.5, places=12)


class InterpolateCmTests(unittest.TestCase):
    def test_linear_map_interpolates_exactly(self):
        map_bp = np.array([0.0, 1e6, 2e6])
        map_cm = np.array([0.0, 1.0, 2.0])  # 1 cM per Mb
        bp = np.array([5e5, 1.5e6])
        np.testing.assert_allclose(
            interpolate_cm(bp, map_bp, map_cm), [0.5, 1.5], rtol=0, atol=1e-12
        )

    def test_extrapolation_uses_flanking_rate_not_flat(self):
        map_bp = np.array([1e6, 2e6])
        map_cm = np.array([1.0, 3.0])  # 2 cM per Mb
        # Below and above the map range must keep moving, not clamp.
        cm = interpolate_cm(np.array([0.0, 3e6]), map_bp, map_cm)
        np.testing.assert_allclose(cm, [-1.0, 5.0], rtol=0, atol=1e-9)

    def test_cm_is_monotonic_across_the_boundary(self):
        map_bp = np.array([1e6, 2e6])
        map_cm = np.array([1.0, 3.0])
        bp = np.array([0.5e6, 1.0e6, 1.5e6, 2.0e6, 2.5e6])
        cm = interpolate_cm(bp, map_bp, map_cm)
        self.assertTrue(np.all(np.diff(cm) > 0))

    def test_rejects_nonincreasing_map(self):
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            interpolate_cm(np.array([1.0]), np.array([2.0, 1.0]), np.array([0.0, 1.0]))


class ReadPlinkMapTests(unittest.TestCase):
    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix=".map")
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        self.addCleanup(os.remove, path)
        return path

    def test_parses_chr_marker_cm_bp_layout(self):
        path = self._write(
            "22\t.\t0.000000\t16051347\n"
            "22\t.\t0.010291\t16052618\n"
            "22\t.\t0.018472\t16053624\n"
        )
        map_bp, map_cm = read_plink_genetic_map(path)
        np.testing.assert_allclose(map_bp, [16051347, 16052618, 16053624])
        np.testing.assert_allclose(map_cm, [0.0, 0.010291, 0.018472])

    def test_collapses_duplicate_bp_to_stay_strictly_increasing(self):
        path = self._write(
            "1 . 0.0 1000\n"
            "1 . 0.5 2000\n"
            "1 . 0.7 2000\n"   # duplicate bp
            "1 . 1.0 3000\n"
        )
        map_bp, map_cm = read_plink_genetic_map(path)
        self.assertTrue(np.all(np.diff(map_bp) > 0))
        # The reader's output must be directly usable by interpolate_cm.
        cm = interpolate_cm(np.array([1500.0, 2500.0]), map_bp, map_cm)
        self.assertEqual(cm.shape, (2,))


if __name__ == "__main__":
    unittest.main()
