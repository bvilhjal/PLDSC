"""Read, QC and align GWAS summaries with allelign before LDSC."""
from __future__ import annotations

from dataclasses import dataclass
import gzip
import math
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


_Z_COLUMNS = ("Z", "ZSCORE", "Z_SCORE")
_MUNGED_FIELDS = {"id": ("SNP", "RSID", "ID"), "ea": ("A1",), "oa": ("A2",), "n_eff": ("N",)}


def read_ldsc_munged(path, *, n_eff=None, **columns):
    """Read an ldsc ``munge_sumstats.py`` file (``SNP A1 A2 [N] Z``).

    Z is with respect to A1. Returns an :class:`allelign.Sumstats` whose
    ``beta`` is Z and ``se`` is 1, so ``beta / se`` is the z score; there is
    no allele frequency, INFO or position, so variants are matched to the
    reference by identifier. Rows with a missing Z or N are dropped and
    counted in ``parse_log``. ``n_eff`` (a number) replaces the N column;
    a string names the N column. ``id``, ``ea`` and ``oa`` may name columns.
    """
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8") as stream:
        first = stream.readline().rstrip("\r\n")
        delim = "\t" if "\t" in first else None          # munge_sumstats writes tabs, NaN as empty
        header = first.split(delim)
        upper = [h.strip().upper() for h in header]

        def find(field, names):
            given = columns.get(field)
            if given is not None:
                if str(given).upper() not in upper:
                    raise ValueError(f"{path}: column {given!r} not found")
                return upper.index(str(given).upper())
            for name in names:
                if name in upper:
                    return upper.index(name)
            return None
        idx = {f: find(f, names) for f, names in _MUNGED_FIELDS.items()}
        if isinstance(n_eff, str):
            if n_eff.upper() not in upper:
                raise ValueError(f"{path}: N column {n_eff!r} not found")
            idx["n_eff"], n_eff = upper.index(n_eff.upper()), None
        iz = next((upper.index(n) for n in _Z_COLUMNS if n in upper), None)
        missing = [f for f in ("id", "ea", "oa") if idx[f] is None] + (["Z"] if iz is None else [])
        if n_eff is None and idx["n_eff"] is None:
            missing.append("N (or a constant n_eff)")
        if missing:
            raise ValueError(f"{path}: not an ldsc munged file; missing {', '.join(missing)}")
        ids, ea, oa, z, n = [], [], [], [], []
        dropped = 0
        for line_no, line in enumerate(stream, 2):
            fields = line.rstrip("\r\n").split(delim)
            if not any(f.strip() for f in fields):
                continue
            if len(fields) != len(header):
                raise ValueError(f"{path}:{line_no}: expected {len(header)} fields, got {len(fields)}")
            try:
                zv = float(fields[iz])
                nv = float(n_eff) if n_eff is not None else float(fields[idx["n_eff"]])
            except ValueError:
                dropped += 1
                continue
            if not (math.isfinite(zv) and math.isfinite(nv) and nv > 0):
                dropped += 1
                continue
            ids.append(fields[idx["id"]].strip())
            ea.append(fields[idx["ea"]].strip().upper())
            oa.append(fields[idx["oa"]].strip().upper())
            z.append(zv)
            n.append(nv)
    m = len(ids)
    if m == 0:
        raise ValueError(f"{path}: no rows with a finite Z and N")
    log = {"format": "ldsc-munged", "n_input": m + dropped, "n_dropped_missing_z_or_n": dropped,
           "columns": {f: header[i] for f, i in idx.items() if i is not None} | {"z": header[iz]},
           "n_eff_source": "scalar" if n_eff is not None else "column"}
    return allelign.Sumstats(np.asarray(ids, dtype=object), np.full(m, "", dtype=object),
                             np.zeros(m, dtype=np.int64), np.asarray(ea, dtype=object),
                             np.asarray(oa, dtype=object), np.asarray(z, dtype=float), np.ones(m),
                             np.asarray(n, dtype=float), np.full(m, np.nan), np.full(m, np.nan),
                             parse_log=log)


def _read(path, read_options):
    options = dict(read_options or {})
    try:
        return allelign.read_sumstats(path, **options)
    except ValueError as exc:
        if "no effect column" not in str(exc):
            raise
        munged = {k: v for k, v in options.items() if k in ("id", "ea", "oa", "n_eff")}
        if len(munged) != len(options):
            raise
        try:
            return read_ldsc_munged(path, **munged)
        except ValueError as munged_exc:
            raise ValueError(f"{exc}; nor is it an ldsc munged file ({munged_exc})") from None


def _prepare(data, variants, read_options, qc_options, drop_ambiguous):
    if isinstance(data, (str, PathLike)):
        data = _read(data, read_options)
    elif read_options:
        raise ValueError("read_options apply only to a summary-statistics file")
    if not isinstance(data, allelign.Sumstats):
        raise TypeError("sumstats must be an allelign.Sumstats or a path")
    keep, qc = allelign.qc_sumstats(data, **(qc_options or {}))
    h = allelign.harmonize(data.subset(keep), variants, drop_ambiguous=drop_ambiguous)
    try:
        qc = {**qc, "reference_coverage": allelign.reference_coverage(data, variants, h)}
    except Exception as exc:  # diagnostic only; never block a fit on it
        qc = {**qc, "reference_coverage": {"error": str(exc)}}
    return h, qc


def _w_ld(w_ld, scores):
    if w_ld is None:
        return None
    w = np.asarray(w_ld, dtype=float)
    if w.shape != scores.shape or not np.all(np.isfinite(w)):
        raise ValueError("w_ld must have one finite value per reference variant")
    return w


def ldsc_h2_sumstats(sumstats, variants, ld_scores, *, m_snps, read_options=None,
                     qc_options=None, drop_ambiguous=True, chi2_cap=None,
                     n_blocks=None, n_iter=2, constrain_intercept=None,
                     ldsc_compat=False, w_ld=None):
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

    ``ldsc_compat=True`` uses the reference ldsc estimator (see
    :mod:`pldsc.ldsc_compat`); the chi-square cap then defaults to off, as in
    ldsc for one-annotation h2, and ``w_ld`` (one value per reference row)
    supplies separate over-counting scores. ``chi2_cap=None`` means on by
    default and off with ``ldsc_compat``.
    """
    m_snps = _finite_control("m_snps", m_snps, lower=1.0)
    variants, scores = _reference(variants, ld_scores)
    w = _w_ld(w_ld, scores)
    if chi2_cap is None:
        chi2_cap = not ldsc_compat
    h, qc = _prepare(sumstats, variants, read_options, qc_options, drop_ambiguous)
    z = h.beta / h.se
    keep = ldsc_chi2_keep(z * z, h.n_eff)[0] if chi2_cap else np.ones(len(h), bool)
    idx = h.var_index[keep]
    estimate = ldsc_h2(z[keep] ** 2, scores[idx], h.n_eff[keep], m_snps=m_snps,
                       n_blocks=n_blocks, n_iter=n_iter,
                       constrain_intercept=constrain_intercept, ldsc_compat=ldsc_compat,
                       w_ld=None if w is None else w[idx])
    return SumstatsLDSCResult(estimate, idx, (qc,), (h.log,), len(h), int((~keep).sum()))


def ldsc_rg_sumstats(sumstats1, sumstats2, variants, ld_scores, *, m_snps,
                     read_options1=None, read_options2=None, qc_options1=None,
                     qc_options2=None, drop_ambiguous=True, chi2_cap=None,
                     n_blocks=200, n_iter=2, constrain_intercept=None, ncores=1,
                     ldsc_compat=False, w_ld=None):
    """QC and harmonize two GWAS files to a reference, then estimate genetic correlation.

    ``m_snps`` counts the variants summed in the LD scores, not the common
    regression rows. Sign flips are performed by allelign. Both traits must
    pass the default LDSC chi-square cap. ``ncores`` accelerates independent
    jackknife fits when Numba is available. ``constrain_intercept`` fixes the
    cross-trait intercept; the marginal intercepts remain free. Reference and
    GWAS coordinates must use the same genome build. See ``docs/LDSC.md``
    for inputs, examples, and interpretation. ``ldsc_compat``, ``w_ld`` and
    ``chi2_cap=None`` behave as in :func:`ldsc_h2_sumstats`.
    """
    m_snps = _finite_control("m_snps", m_snps, lower=1.0)
    variants, scores = _reference(variants, ld_scores)
    w = _w_ld(w_ld, scores)
    if chi2_cap is None:
        chi2_cap = not ldsc_compat
    h1, q1 = _prepare(sumstats1, variants, read_options1, qc_options1, drop_ambiguous)
    h2, q2 = _prepare(sumstats2, variants, read_options2, qc_options2, drop_ambiguous)
    idx, i, j = np.intersect1d(h1.var_index, h2.var_index, assume_unique=True, return_indices=True)
    z1, z2 = h1.beta[i] / h1.se[i], h2.beta[j] / h2.se[j]
    n1, n2 = h1.n_eff[i], h2.n_eff[j]
    keep = (ldsc_chi2_keep(z1 * z1, n1)[0] & ldsc_chi2_keep(z2 * z2, n2)[0]
            if chi2_cap else np.ones(idx.size, bool))
    estimate = ldsc_rg_z(z1[keep], z2[keep], scores[idx[keep]], n1[keep], n2[keep],
                         m_snps=m_snps, n_blocks=n_blocks, n_iter=n_iter,
                         constrain_intercept=constrain_intercept, ncores=ncores,
                         ldsc_compat=ldsc_compat, w_ld=None if w is None else w[idx[keep]])
    return SumstatsLDSCResult(estimate, idx[keep], (q1, q2), (h1.log, h2.log),
                              idx.size, int((~keep).sum()))
