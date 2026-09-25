#!/usr/bin/env python3
"""Regenerate examples/mini/chr22_ldscores.npz from the local 1000G scratch data.

The bundle is everything `multiancestry_1000g_meta.py` needs to reproduce
MULTI_ANCESTRY.md Tables 3-5 *except* the raw genotypes: the (5, 5, J)
cross-ancestry LD-score tensor at the 9,654 chr22 regression SNPs, the
polymorphism mask, the per-ancestry genotype variances (inverse-variance
meta weights), and the per-ancestry MAFs (GWAS availability and the
EUR-common restriction). With it, the headline kappa tables are a
sub-second, no-download computation (examples/reproduce_chr22_tables.py).

Requires .scratch_1000g/ (chr22 VCF cache + panel + SHAPEIT4 chr22 map);
`python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22`
builds all of it. Provenance and data terms: examples/mini/README.md.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

import multiancestry_1000g_meta as M   # noqa: E402  (reuses its load path)

ANCESTRIES = M.ANCESTRIES
CHROM = M.CHROM


def main():
    panels, G, cm = M.load(10)
    reg = np.flatnonzero(M.maf(G) >= M.MAF_REGRESSION)
    s = M.cross_ancestry_ld_scores(panels, cm, window_cm=M.WINDOW_CM, regression=reg)
    maf = np.stack([M.maf(P)[reg] for P in panels])
    panel_n = np.array([P.shape[0] for P in panels], dtype=np.int64)
    assert s.scores.shape == (5, 5, reg.size)
    out = os.path.join(HERE, "mini", "chr22_ldscores.npz")
    np.savez_compressed(
        out,
        scores=s.scores, available=s.available, variances=s.variances, maf=maf,
        ancestries=np.array(ANCESTRIES), panel_n=panel_n,
        chrom=CHROM, thin=np.array(10, dtype=np.int64),
        window_cm=M.WINDOW_CM, maf_regression=M.MAF_REGRESSION, maf_gwas=M.MAF_GWAS,
    )
    print(f"wrote {out}: scores {s.scores.shape}, {os.path.getsize(out) / 1e6:.1f} MB")
    print("  mean within-ancestry scores (should match MULTI_ANCESTRY Table 2): "
          + ", ".join(f"{a} {s.scores[i, i].mean():.2f}"
                      for i, a in enumerate(ANCESTRIES)))


if __name__ == "__main__":
    main()
