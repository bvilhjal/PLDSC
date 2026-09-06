import unittest

import numpy as np

from pldsc.simulate import simulate_admixture


class AdmixtureFamilyTests(unittest.TestCase):
    def test_returns_exactly_n_indiv_when_divisible(self):
        rng = np.random.default_rng(np.random.PCG64(5))
        G, chrom, pos = simulate_admixture(rng, [10, 12], 8, families=True, sibs=4)
        self.assertEqual(G.shape, (8, 22))

    def test_rejects_indivisible_family_sizes(self):
        rng = np.random.default_rng(np.random.PCG64(6))
        with self.assertRaisesRegex(ValueError, "multiple"):
            simulate_admixture(rng, [10, 12], 10, families=True, sibs=4)


if __name__ == "__main__":
    unittest.main()
