"""Smoke tests for pldsc.plotting (skipped when matplotlib is absent)."""
import os
import tempfile
import unittest

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    HAS_MPL = True
except ImportError:                                  # pragma: no cover
    HAS_MPL = False

from pldsc.plotting import forest_plot, kappa_heatmap, kappa_profile  # noqa: E402


@unittest.skipUnless(HAS_MPL, "matplotlib not installed")
class PlottingSmokeTests(unittest.TestCase):
    def setUp(self):
        fd, self.png = tempfile.mkstemp(suffix=".png")
        os.close(fd)
        self.addCleanup(os.remove, self.png)
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_kappa_heatmap_annotates_grid(self):
        rng = np.random.default_rng(1)
        K = rng.uniform(0.6, 1.0, size=(3, 3))
        ax = kappa_heatmap(K, ["EUR 0.9+AFR", "EUR 0.75+AFR", "EUR 0.5+AFR"],
                           ["rho=1", "rho=0.8", "rho=0.6"], path=self.png)
        self.assertEqual(len(ax.texts), 9)
        self.assertGreater(os.path.getsize(self.png), 0)

    def test_kappa_profile_draws_gate_lines(self):
        K = np.full((5, 2), 0.75)
        ax = kappa_profile(K, [18, 19, 20, 21, 22], ancestries=["AFR", "EAS"],
                           path=self.png)
        # one dashed gate line + one dotted reference line among 2 data lines
        self.assertEqual(sum(line.get_linestyle() in ("--", ":")
                             for line in ax.lines), 2)

    def test_forest_plot_orders_by_relative_change(self):
        ax = forest_plot(["t1", "t2", "t3"], [0.2, 0.2, 0.2],
                         [0.1, 0.15, 0.19], [0.01, 0.01, 0.01], path=self.png)
        labels = [t.get_text() for t in ax.get_yticklabels()]
        self.assertEqual(labels, ["t3", "t2", "t1"])   # bottom-up; t1 on top

    def test_shape_mismatches_raise(self):
        with self.assertRaises(ValueError):
            kappa_heatmap(np.zeros((2, 2)), ["a"], ["b", "c"])
        with self.assertRaises(ValueError):
            kappa_profile(np.zeros((3, 2)), [1, 2])


if __name__ == "__main__":
    unittest.main()
