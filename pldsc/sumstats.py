"""Read, QC and align GWAS summaries with allelign before LDSC."""
from __future__ import annotations

from dataclasses import dataclass
from os import PathLike

import numpy as np
import allelign

from .ldsc import LDSCResult, ldsc_chi2_keep, ldsc_h2
from .bivariate import LDSCRgResult, ldsc_rg_z
from ._common import _finite_control


@dataclass
class SumstatsLDSCResult:
    """LDSC estimate and per-trait QC and alignment records.

    ``reference_index`` contains the retained reference rows; ``n_common``
    counts aligned rows before the chi-square cap and ``n_chi2_dropped``
    counts rows removed by that cap. ``qc`` and ``alignment`` follow the
    input-trait order.
    """
    estimate: LDSCResult | LDSCRgResult
    reference_index: np.ndarray
    qc: tuple[dict, ...]
    alignment: tuple[dict, ...]
    n_common: int
    n_chi2_dropped: int


def _reference(variants, ld_scores):
    if isinstance(variants, (str, PathLike)):
        variants = allelign.read_variants(variants)
    scores = np.asarray(ld_scores, dtype=float)
    if scores.ndim != 1 or len(variants) != scores.size:
        raise ValueError("ld_scores must have one value per reference variant")
    # Genomic blocks must be contiguous. Chromosome order itself need not be
    # numerical, but a chromosome may not recur and positions must increase.
    chrom = np.asarray(variants.chrom)
    pos = np.asarray(variants.pos)
    if chrom.shape != scores.shape or pos.shape != scores.shape or not scores.size:
        raise ValueError("reference requires chromosome and position for every variant")
    if np.any(pos <= 0):
        raise ValueError("reference positions must be positive")
    normalized = np.array([allelign.normalize_chrom(str(c)) for c in np.unique(chrom)])
    _, inv = np.unique(chrom, return_inverse=True)
    cc = normalized[inv]
    starts = np.r_[True, cc[1:] != cc[:-1]]
    if np.unique(cc[starts]).size != np.count_nonzero(starts):
        raise ValueError("reference chromosomes must form contiguous blocks")
    if np.any((cc[1:] == cc[:-1]) & (pos[1:] < pos[:-1])):
        raise ValueError("reference positions must be sorted within chromosomes")
    return variants, scores


def _prepare(data, variants, read_options, qc_options, drop_ambiguous):
    if isinstance(data, (str, PathLike)):
        data = allelign.read_sumstats(data, **(read_options or {}))
    elif read_options:
        raise ValueError("read_options apply only to a summary-statistics file")
    if not isinstance(data, allelign.Sumstats):
        raise TypeError("sumstats must be an allelign.Sumstats or a path")
    keep, qc = allelign.qc_sumstats(data, **(qc_options or {}))
    h = allelign.harmonize(data.subset(keep), variants, drop_ambiguous=drop_ambiguous)
    return h, qc


def ldsc_h2_sumstats(sumstats, variants, ld_scores, *, m_snps, read_options=None,
                     qc_options=None, drop_ambiguous=True, chi2_cap=True,
                     n_blocks=100, n_iter=2, constrain_intercept=None):
    """Read/QC/align a GWAS and estimate SNP heritability.

    ``variants`` is an allelign Variants table or .bim/.pvar path; ``ld_scores``
    has one value per row in that table. ``m_snps`` is the total number of SNPs
    summed into each LD score, before restricting regression rows. Supply a
    compatible effective N in the file or ``read_options``. The standard LDSC
    chi-square cap is applied by default. Numeric ``read_options['n_eff']``
    overrides per-variant sample sizes; omit it to use the detected N source.
    See :func:`allelign.qc_sumstats` for ``qc_options`` and the returned logs.
    Reference and GWAS coordinates must use the same genome build. See
    ``docs/LDSC.md`` for inputs, examples, and interpretation.
    """
    m_snps = _finite_control("m_snps", m_snps, lower=1.0)
    variants, scores = _reference(variants, ld_scores)
    h, qc = _prepare(sumstats, variants, read_options, qc_options, drop_ambiguous)
    z = h.beta / h.se
    keep = ldsc_chi2_keep(z * z, h.n_eff)[0] if chi2_cap else np.ones(len(h), bool)
    idx = h.var_index[keep]
    estimate = ldsc_h2(z[keep] ** 2, scores[idx], h.n_eff[keep], m_snps=m_snps,
                       n_blocks=n_blocks, n_iter=n_iter,
                       constrain_intercept=constrain_intercept)
    return SumstatsLDSCResult(estimate, idx, (qc,), (h.log,), len(h), int((~keep).sum()))


def ldsc_rg_sumstats(sumstats1, sumstats2, variants, ld_scores, *, m_snps,
                     read_options1=None, read_options2=None, qc_options1=None,
                     qc_options2=None, drop_ambiguous=True, chi2_cap=True,
                     n_blocks=200, n_iter=2, constrain_intercept=None, ncores=1):
    """QC and harmonize two GWAS files to a reference, then estimate genetic correlation.

    ``m_snps`` counts the variants summed in the LD scores, not the common
    regression rows. Sign flips are performed by allelign. Both traits must
    pass the default LDSC chi-square cap. ``ncores`` accelerates independent
    jackknife fits when Numba is available. ``constrain_intercept`` fixes the
    cross-trait intercept; the marginal intercepts remain free. Reference and
    GWAS coordinates must use the same genome build. See ``docs/LDSC.md``
    for inputs, examples, and interpretation.
    """
    m_snps = _finite_control("m_snps", m_snps, lower=1.0)
    variants, scores = _reference(variants, ld_scores)
    h1, q1 = _prepare(sumstats1, variants, read_options1, qc_options1, drop_ambiguous)
    h2, q2 = _prepare(sumstats2, variants, read_options2, qc_options2, drop_ambiguous)
    idx, i, j = np.intersect1d(h1.var_index, h2.var_index, assume_unique=True, return_indices=True)
    z1, z2 = h1.beta[i] / h1.se[i], h2.beta[j] / h2.se[j]
    n1, n2 = h1.n_eff[i], h2.n_eff[j]
    keep = (ldsc_chi2_keep(z1 * z1, n1)[0] & ldsc_chi2_keep(z2 * z2, n2)[0]
            if chi2_cap else np.ones(idx.size, bool))
    estimate = ldsc_rg_z(z1[keep], z2[keep], scores[idx[keep]], n1[keep], n2[keep],
                         m_snps=m_snps, n_blocks=n_blocks, n_iter=n_iter,
                         constrain_intercept=constrain_intercept, ncores=ncores)
    return SumstatsLDSCResult(estimate, idx[keep], (q1, q2), (h1.log, h2.log),
                              idx.size, int((~keep).sum()))
