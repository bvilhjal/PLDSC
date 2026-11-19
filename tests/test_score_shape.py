"""Independent checks of the held-out risk and the known-truth control."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from score_shape import asthma_result, centered_products, regional_moments, residual_slope  # noqa: E402
from score_shape_sim import block_chi_square, scores  # noqa: E402


def test_held_out_predictions_and_deletion_se_match_direct_wls():
    rng = np.random.default_rng(84)
    chrom = np.repeat(np.arange(6), [19, 31, 12, 9, 35, 17])
    eur, difference, noise = rng.normal(size=(3, len(chrom)))
    matched = 0.8 * eur + 0.2 * difference
    y, w = 1.2 + matched + noise, rng.uniform(0.1, 2, len(chrom))
    result = asthma_result(chrom, eur, matched, y, w)
    def direct_gain(keep):
        losses = []
        for j in np.unique(chrom[keep]):
            train, test = keep & (chrom != j), keep & (chrom == j)
            fold = []
            for x in (eur, matched):
                a = np.column_stack((np.ones_like(x), x))
                beta = np.linalg.lstsq(a[train] * np.sqrt(w[train, None]),
                                       y[train] * np.sqrt(w[train]), rcond=None)[0]
                fold.append(np.sum(w[test] * (y[test] - a[test] @ beta)**2))
            losses.append(fold)
        return np.asarray(losses)
    losses = direct_gain(np.ones(len(chrom), dtype=bool))
    np.testing.assert_allclose(losses, [[r["eur_sse"], r["matched_sse"]] for r in result["chromosomes"]])
    deletions = np.array([direct_gain(chrom != c).sum(axis=0) for c in np.unique(chrom)])
    gains = 100 * (1 - deletions[:, 1] / deletions[:, 0])
    assert result["summary"]["gain_percent"]["se"] == pytest.approx(np.sqrt(5) * gains.std())
    beta = np.linalg.lstsq(np.column_stack((np.ones_like(y), eur, matched)) * np.sqrt(w[:, None]),
                           y * np.sqrt(w), rcond=None)[0]
    assert result["summary"]["residual_slope"]["estimate"] == pytest.approx(beta[2])


def test_residual_signal_has_correct_endpoints_and_rejects_collinearity():
    eur = np.arange(16, dtype=float)
    matched = 2 + 0.5 * eur + np.sin(eur)
    def theta(candidate, y):
        moments = regional_moments(eur, candidate, y, np.ones_like(eur), [0, 5, 10])
        return residual_slope(centered_products(moments.sum(axis=0)))
    assert theta(matched, 4 + 1.3 * eur) == pytest.approx(0, abs=1e-12)
    assert theta(matched, 4 + 1.3 * matched) == pytest.approx(1.3)
    with pytest.raises(ValueError, match="collinear"):
        theta(7 + 0.5 * eur, matched)


def test_gaussian_pair_eigenmodes_reproduce_ldsc_expectation():
    r_e = np.array([0.9, 0.25, 0.4])
    r_f = np.array([-0.2, 0.5, 0.65])
    for share in (0, 0.25, 0.5):
        # Unit means of the two independent chi-squares integrate out the GWAS.
        population_mean = block_chi_square(r_e, r_f, share, np.ones((2, 3)))
        _, matched = scores(r_e, r_f, 0, share)
        np.testing.assert_allclose(population_mean, 1 + matched)
