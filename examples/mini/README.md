# examples/mini: the chr22 LD-score bundle

`chr22_ldscores.npz` (2.0 MB) is everything needed to reproduce the
headline tables of [docs/MULTI_ANCESTRY.md](../../docs/MULTI_ANCESTRY.md)
(Table 3, the κ grid; Table 4, the kill-rule cell; Table 5, the
identifiability correlations) **without any download**:

```bash
python examples/reproduce_chr22_tables.py    # ~0.6 s, all 36 κ cells
```

## Contents

| key | shape | meaning |
|---|---|---|
| `scores` | 5×5×9,654 | cross-ancestry LD scores ℓ^ab_j (MA1/MA7): within-ancestry on the diagonal (noise-floor corrected), cross-ancestry off it |
| `available` | 5×9,654 | variant polymorphic in panel a |
| `variances` | 5×9,654 | residual genotype variance (per-allele information for IVW weights) |
| `maf` | 5×9,654 | minor allele frequency per ancestry |
| `ancestries`, `panel_n` | 5 | EUR 503, AFR 661, EAS 504, SAS 489, AMR 347 |
| `chrom`, `thin`, `window_cm`, `maf_regression`, `maf_gwas` | 1 | the exact configuration of the recording run: chr22, 1-in-10 thinned cache, ±1 cM, regression MAF ≥ 0.05, GWAS MAF ≥ 0.01 |

Regenerate with `examples/make_mini.py` (needs `.scratch_1000g/`:
`python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22`).

## Provenance and data terms

The MIT [LICENSE](../../LICENSE) of this repository covers the **software
only**. This bundle is a derived, aggregate statistic (per-variant sums of
LD correlations) computed from:

- **1000 Genomes Project phase 3**, the 20130502 release
  (`v5a` genotypes via the S3 mirror; `samples.panel` superpopulations).
  1000G data are open with no data-use agreement; if you use or
  redistribute this bundle, please credit the consortium: *1000 Genomes
  Project Consortium, "A global reference for human genetic variation",
  Nature 526:68–74 (2015), doi:10.1038/nature15393.*
- the **SHAPEIT4 GRCh37 genetic map** (`chr22.b37.gmap.gz`) from the
  MIT-licensed [odelaneau/shapeit4](https://github.com/odelaneau/shapeit4)
  `maps/` directory.

No individual-level genotypes are included or recoverable from the bundle.
If the upstream terms matter to your use, re-derive from the sources above
with the scripts in this directory and `scripts/fetch_1000g.py`.
