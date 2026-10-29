#!/usr/bin/env python3
"""Reproduce the PLDSC headline tables from the committed LD-score bundle.

No download, no genotypes, no genetic map: everything MULTI_ANCESTRY.md
Tables 3-5 need is in examples/mini/chr22_ldscores.npz (see
examples/mini/README.md for provenance). Runs in well under a second.

  python examples/reproduce_chr22_tables.py
"""
from __future__ import annotations

import os

import numpy as np

from pldsc.multiancestry import meta_ld_scores, meta_weights, single_reference_bias

HERE = os.path.dirname(os.path.abspath(__file__))
EUR_SHARES = (0.9, 0.75, 0.5)
RHOS = (1.0, 0.8, 0.6)

d = np.load(os.path.join(HERE, "mini", "chr22_ldscores.npz"))
anc, scores, maf = list(d["ancestries"]), d["scores"], d["maf"]
gwas_ok = (maf >= float(d["maf_gwas"])) & d["available"]
eur = anc.index("EUR")
eur_common = maf[eur] >= 0.05
A, lE_all = len(anc), scores[eur, eur]


def cells(X, fE, rho):
    n = np.zeros(A); n[eur], n[anc.index(X)] = fE, 1.0 - fE
    om, nj_full = meta_weights(n, gwas_ok, variances=d["variances"], scheme="sample_size")
    R = np.full((A, A), rho); np.fill_diagonal(R, 1.0)
    L = meta_ld_scores(scores, om, rho=R)[eur_common]
    nj, lE = nj_full[eur_common], lE_all[eur_common]
    x = lE * nj / nj.mean()
    w = 1.0 / (np.maximum(lE, 1.0) * (1.0 + 0.5 / x.mean() * x) ** 2)
    return (single_reference_bias(L, lE, sample_sizes=nj).kappa,
            single_reference_bias(L, lE, sample_sizes=nj, weights=w).kappa)


print("Table 3 (MULTI_ANCESTRY.md). Panel plug-in kappa for a meta regressed on the")
print(f"1000G EUR score; OLS / LDSC-style weights, EUR MAF >= 0.05 ({eur_common.sum()} SNPs)")
print(f"{'mix':>16s} " + " ".join(f"{'rho=' + str(r):>19s}" for r in RHOS))
for X in anc[1:]:
    for fE in EUR_SHARES:
        print(f"{f'EUR {fE:.2f}+{X}':>16s} "
              + " ".join(f"{c[0]:.3f} / {c[1]:.3f}".center(19) for c in (cells(X, fE, r) for r in RHOS)))

print("\nTable 4, rule cell (EUR 0.75 + X, rho = 0.8), corrected:")
rule = {X: cells(X, 0.75, 0.8) for X in anc[1:]}
print("  " + ", ".join(f"{X} {v[0]:.3f} / {v[1]:.3f}" for X, v in rule.items())
      + f"  -> KEEP: {max(abs(v[0] - 1) for v in rule.values()) >= 0.10}")
print("  recorded: AFR 0.750 / 0.750, EAS 0.841 / 0.845, "
      "SAS 0.868 / 0.870, AMR 0.907 / 0.896")

print("\nTable 5. corr(within, cross) across SNPs -- can the meta alone identify rho?")
for X in anc[1:]:
    n = np.zeros(A); n[eur], n[anc.index(X)] = 0.75, 0.25
    om, nj_full = meta_weights(n, gwas_ok, variances=d["variances"], scheme="sample_size")
    full = meta_ld_scores(scores, om)[eur_common] * nj_full[eur_common]
    within = meta_ld_scores(scores, om, rho=np.eye(A))[eur_common] * nj_full[eur_common]
    print(f"  EUR 0.75 + {X}: corr = {np.corrcoef(within, full - within)[0, 1]:.3f}")
