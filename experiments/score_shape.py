"""Compare the association patterns predicted by two LD scores.

Each regional matrix is the weighted cross-product of (1, x_EUR, x_matched,
z²). Summing these small matrices exactly reproduces SNP-level least squares.
The GBMI pilot fixes the effect correlation at 0.8 for every pair and eta=1.
Its primary endpoint is the paired reduction in weighted prediction error on
chromosomes excluded from fitting. It reuses the historical asthma sample.
"""
import numpy as np


def regional_moments(eur, matched, y, weight, starts):
    """Cross-products within consecutive groups; starts includes zero."""
    design = np.column_stack((np.ones_like(y), eur, matched, y))
    if (weight.shape != (len(design),) or not np.all(np.isfinite(design))
            or not np.all(np.isfinite(weight) & (weight > 0))):
        raise ValueError("finite predictors and responses and positive weights are required")
    if len(starts) == 0 or starts[0] != 0 or starts[-1] >= len(y) or np.any(np.diff(starts) <= 0):
        raise ValueError("groups must start at zero and contain at least one row")
    return np.stack([rows.T @ (w[:, None] * rows)
                     for rows, w in zip(np.split(design, starts[1:]), np.split(weight, starts[1:]))])


def centered_products(moment):
    """Centered sums for (EUR, matched, response), for one or more fits."""
    first = moment[..., 1:, 0]
    return moment[..., 1:, 1:] - first[..., :, None] * first[..., None, :] / moment[..., 0, 0, None, None]


def residual_slope(centered):
    """Coefficient on matched after an intercept and EUR; undefined at collinearity."""
    v_eur, v_matched = centered[..., 0, 0], centered[..., 1, 1]
    if np.any(v_eur <= 0) or np.any(v_matched <= 0):
        raise ValueError("both scores must vary within the fit")
    kappa = centered[..., 0, 1] / v_eur
    residual_var = v_matched - kappa * centered[..., 0, 1]
    if np.any(residual_var <= 1e-12 * v_matched):
        raise ValueError("the matched score is collinear with EUR")
    return (centered[..., 1, 2] - kappa * centered[..., 0, 2]) / residual_var


def chromosome_losses(moments):
    """Weighted sums of squared errors on each omitted chromosome for both models."""
    train = moments.sum(axis=0) - moments
    centered = centered_products(train)
    variances = centered[:, (0, 1), (0, 1)]
    if np.any(variances <= 0):
        raise ValueError("a training score has zero variance")
    slopes = centered[:, :2, 2] / variances
    intercepts = (train[:, 0, 3, None] - slopes * train[:, 0, 1:3]) / train[:, 0, 0, None]
    # Each contrast encodes y - intercept - slope*x on the test chromosome.
    contrasts = np.zeros((len(train), 2, 4))
    contrasts[:, :, 0] = -intercepts
    contrasts[:, (0, 1), (1, 2)] = -slopes
    contrasts[:, :, 3] = 1
    return np.einsum("nki,nij,nkj->nk", contrasts, moments, contrasts)


def gain_percent(loss):
    """Reduction in total held-out squared error as a percentage of EUR error."""
    total = loss.sum(axis=0)
    if total[0] <= 0:
        raise ValueError("relative gain requires positive EUR prediction error")
    return float(100 * (total[0] - total[1]) / total[0])


def jackknife_se(deleted):
    return float(np.std(deleted, ddof=0) * np.sqrt(len(deleted) - 1))


def asthma_result(chrom, eur, matched, y, weight):
    """Leave one chromosome out for prediction; refit all folds on each deletion."""
    labels, starts, counts = np.unique(chrom, return_index=True, return_counts=True)
    if not np.array_equal(chrom, np.repeat(labels, counts)) or len(labels) < 4:
        raise ValueError("four or more ordered chromosomes are required")
    moments = regional_moments(eur, matched, y, weight, starts)
    total = moments.sum(axis=0)
    centered = centered_products(total)
    deleted = centered_products(total - moments)
    loss = chromosome_losses(moments)
    gain_deleted = [gain_percent(chromosome_losses(np.delete(moments, j, axis=0)))
                    for j in range(len(moments))]
    def estimate(value, leave_one_out):
        return {"estimate": float(value), "se": jackknife_se(leave_one_out)}
    estimates = {}
    for k, name in enumerate(("eur_h2", "matched_h2")):
        estimates[name] = estimate(centered[k, 2] / centered[k, k], deleted[:, k, 2] / deleted[:, k, k])
    theta = residual_slope(centered)
    theta_deleted = residual_slope(deleted)
    estimates["residual_slope"] = estimate(theta, theta_deleted)
    estimates["residual_minus_matched"] = estimate(
        theta - estimates["matched_h2"]["estimate"], theta_deleted - deleted[:, 1, 2] / deleted[:, 1, 1])
    estimates["gain_percent"] = estimate(gain_percent(loss), gain_deleted)
    return {
        "rho_all_pairs": 0.8, "eta_all_groups": 1.0, "n_snps": len(chrom),
        "predictors": "N_ALL times LD score divided by M", "weight": "1/max(EUR LD score,1)",
        "residual_variance_fraction": float(1 - centered[0, 1]**2 / (centered[0, 0] * centered[1, 1])),
        "summary": estimates, "moment_order": ["intercept", "EUR", "matched", "chi_square"],
        "chromosomes": [
            {"chr": int(c), "n_snps": int(n), "weighted_crossproduct": moment.tolist(),
             "eur_sse": float(sse[0]), "matched_sse": float(sse[1])}
            for c, n, moment, sse in zip(labels, counts, moments, loss)
        ],
    }
