"""Relatedness in the LD panel: the residual-df correction is wrong.

The finite-sample adjustment in `pldsc.ld.local_ld_scores`,
`r2 - (1 - r2) / (nu - 1)` with `nu = N - rank(C)`, assumes independent
samples. This experiment tests the generalization

    E[r_hat^2] ~= tr(Atilde^2) / tr(Atilde)^2  =: 1 / N_eff,
    Atilde = P_C A P_C,

for two UNLINKED variants in a panel with relationship matrix `A`, and
quantifies what using `nu - 1` instead costs a windowed LD score.

Design notes. Genotypes come from real Mendelian transmission (founder
haplotypes, one transmitted per parent per SNP), not a Gaussian surrogate, so
the relatedness is genuine and the discrete genotype distribution is correct.
SNPs are independent by construction, so every pair is unlinked and the
observed mean squared correlation IS the noise floor. The predicted value is
computed by the shipped `pldsc.effective_pair_df`, so this validates the
released code path rather than a reimplementation.

Verdict (see `experiments/README.md`): the tr-ratio prediction is accurate to
<=1% across designs spanning a 4.7x inflation, while `nu - 1` is off by that
same 4.7x. But the resulting LD-score bias is near-uniform across variants, so
it lands on the LDSC intercept and not the slope -- which is why this is a
correctness fix and not a heritability-estimation result. See `docs/METHOD.md`
Section 10.1.

Run:  conda run -n pldsc python experiments/relatedness_effective_df.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from pldsc import effective_pair_df  # noqa: E402

SEED = 20260723
N_SNPS = 600
N_REPS = 8

# (label, families, sibs per family) -- all N = 400 so the naive df is constant.
DESIGNS = [
    ("unrelated (control)", 400, 1),
    ("100 fam x 4 sibs", 100, 4),
    ("50 fam x 8 sibs", 50, 8),
    ("25 fam x 16 sibs", 25, 16),
]


def sib_panel(rng, n_families, sibs_per_family, n_snps):
    """Genotypes for full-sib families at unlinked SNPs, by transmission."""
    freq = rng.uniform(0.05, 0.5, size=n_snps)
    # (families, parent, haplotype, snps)
    hap = (rng.random((n_families, 2, 2, n_snps)) < freq).astype(np.float64)
    out = np.empty((n_families, sibs_per_family, n_snps))
    for s in range(sibs_per_family):
        pick = rng.integers(0, 2, size=(n_families, 2, n_snps))
        transmitted = np.take_along_axis(hap, pick[:, :, None, :], axis=2)
        out[:, s, :] = transmitted[:, :, 0, :].sum(axis=1)
    return out.reshape(n_families * sibs_per_family, n_snps)


def additive_relationship(n_families, sibs_per_family):
    """1 on the diagonal, 0.5 between full sibs, 0 between families."""
    block = 0.5 * np.ones((sibs_per_family, sibs_per_family))
    np.fill_diagonal(block, 1.0)
    return np.kron(np.eye(n_families), block)


def observed_floor(G):
    """Mean squared correlation over all distinct pairs of unlinked variants."""
    Z = G - G.mean(axis=0)
    norms = np.linalg.norm(Z, axis=0)
    keep = norms > 0
    Z = Z[:, keep] / norms[keep]
    R = Z.T @ Z
    iu = np.triu_indices(R.shape[0], k=1)
    return float(np.mean(R[iu] ** 2))


def main():
    rng = np.random.default_rng(np.random.PCG64(SEED))

    print("Unlinked-pair noise floor E[r^2] (N = 400 throughout)")
    print(
        f"{'design':>22} {'N_eff':>8} {'observed':>11} {'naive 1/(nu-1)':>15}"
        f" {'tr-ratio':>10} {'obs/naive':>10} {'obs/tr':>7}"
    )
    print("-" * 92)

    rows = []
    for label, n_fam, n_sibs in DESIGNS:
        n = n_fam * n_sibs
        A = additive_relationship(n_fam, n_sibs)
        df = effective_pair_df(A)
        predicted = 1.0 / df
        naive = 1.0 / (n - 2)  # nu - 1, with nu = N - rank(intercept)

        floors = [
            observed_floor(sib_panel(rng, n_fam, n_sibs, N_SNPS))
            for _ in range(N_REPS)
        ]
        observed = float(np.mean(floors))

        print(
            f"{label:>22} {df:>8.1f} {observed:>11.6f} {naive:>15.6f}"
            f" {predicted:>10.6f} {observed / naive:>10.2f}"
            f" {observed / predicted:>7.2f}"
        )
        rows.append((label, predicted - naive))

    print()
    print("Resulting bias in a windowed LD score, by variants per window S")
    print(f"{'design':>22} {'S=100':>10} {'S=500':>10} {'S=2000':>10}")
    print("-" * 56)
    for label, gap in rows:
        print(f"{label:>22}" + "".join(f"{gap * S:>10.3f}" for S in (100, 500, 2000)))

    print()
    print("The bias is near-uniform across variants at fixed S, so it shifts the")
    print("LDSC intercept rather than the slope; h^2 is protected but any")
    print("intercept-based confounding readout on a related panel is not.")


if __name__ == "__main__":
    main()
