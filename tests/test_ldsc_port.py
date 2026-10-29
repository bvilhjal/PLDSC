"""Known-truth, provider-parity and file-workflow checks for LDSC."""
from dataclasses import asdict

import allelign
import numpy as np
import pytest

import pldsc
from pldsc import LowRankLD, ld_scores, ldsc_h2, ldsc_rg, ldsc_rg_z


def _data(m=160):
    rng = np.random.default_rng(723)
    ell = rng.uniform(1, 5, m)
    n1 = rng.uniform(8000, 12000, m)
    n2 = rng.uniform(9000, 13000, m)
    scale1 = np.sqrt(1.1 + n1 * 0.25 * ell / m)
    scale2 = np.sqrt(1.2 + n2 * 0.31 * ell / m)
    shared, noise = rng.normal(size=(2, m))
    z1 = shared * scale1
    z2 = (0.7 * shared + np.sqrt(1 - 0.7**2) * noise) * scale2
    return ell, n1, n2, z1, z2


def _assert_results_close(a, b):
    for name, v in asdict(a).items():
        other = asdict(b)[name]
        if v is None:
            assert other is None
        else:
            np.testing.assert_allclose(v, other, rtol=2e-9, atol=2e-10, equal_nan=True)


def test_dense_and_lowrank_scores_and_reference_correction():
    rng = np.random.default_rng(91)
    u = rng.normal(size=(30, 4)) * 0.09
    R = u @ u.T
    np.fill_diagonal(R, 1.0)
    expected = (R * R).sum(axis=1)
    np.testing.assert_allclose(ld_scores(R), expected, atol=1e-12)
    np.testing.assert_allclose(ld_scores([(LowRankLD(u), np.arange(30))]), expected, atol=1e-12)
    np.testing.assert_allclose(ld_scores(R, n_ref=100), expected - (30 - expected) / 98, atol=1e-12)
    q = np.rint(R * 127).astype(np.int8)
    np.testing.assert_allclose(ld_scores(q), ((q.astype(float) / 127)**2).sum(axis=1), atol=1e-12)
    with pytest.raises(ValueError, match="duplicates"):
        ld_scores([(np.eye(2), np.array([0, 0]))])
    with pytest.raises(ValueError, match="symmetric"):
        ld_scores(np.array([[1.0, 0.5], [0.0, 1.0]]))


def test_lowrank_tiled_and_quantized_scores(monkeypatch):
    import importlib
    score_module = importlib.import_module("pldsc.ld_scores")
    rng = np.random.default_rng(61)
    q = rng.integers(-20, 21, size=(80, 9), dtype=np.int8)
    scale = 1 / 127.0
    u = q.astype(float) * scale
    R = u @ u.T
    np.fill_diagonal(R, 1)
    monkeypatch.setattr(score_module, "_LD_SCORE_WORK_BYTES", 1)
    np.testing.assert_allclose(ld_scores([(LowRankLD(q, scale), np.arange(80))]),
                               np.sum(R * R, axis=1), atol=1e-12)


def test_known_line_univariate_and_varying_n():
    ell, n1, _, _, _ = _data()
    chi = 1.2 + n1 * ell / 1000 * 0.35
    r = ldsc_h2(chi, ell, n1, m_snps=1000, n_blocks=10, n_iter=2)
    assert r.h2 == pytest.approx(0.35, abs=1e-12)
    assert r.intercept == pytest.approx(1.2, abs=1e-12)
    assert r.h2_se < 1e-12


def test_provider_parity():
    l = pytest.importorskip("ldpred3.ldsc")
    b = pytest.importorskip("bipred.ldsc")
    representation = pytest.importorskip("ldpred3.ld_repr")
    factor = np.random.default_rng(91).normal(size=(30, 4)) * .09
    reference_block = representation.LowRankLD(factor, m=30)
    blocks = [(reference_block, np.arange(30))]
    np.testing.assert_allclose(ld_scores(blocks, n_ref=100), l.ld_scores(blocks, n_ref=100),
                               rtol=1e-12, atol=1e-12)
    ell, n1, n2, z1, z2 = _data()
    beta1, beta2 = z1 / np.sqrt(n1 + z1*z1), z2 / np.sqrt(n2 + z2*z2)
    for c in (None, 0.0):
        kwargs = dict(m_snps=400, n_blocks=10, n_iter=2, constrain_intercept=c)
        _assert_results_close(ldsc_h2(z1*z1, ell, n1, **kwargs), l.ldsc_h2(z1*z1, ell, n1, **kwargs))
        ours = ldsc_rg(beta1, beta2, ell, n1, n2, **kwargs)
        _assert_results_close(ours, b.ldsc_rg(beta1, beta2, ell, n1, n2, **kwargs))
        _assert_results_close(ours, ldsc_rg_z(z1, z2, ell, n1, n2, **kwargs))
        _assert_results_close(ours, ldsc_rg_z(z1, z2, ell, n1, n2, ncores=2, **kwargs))


@pytest.mark.parametrize("intercept,expected", [
    (None, (0.71495856233018, 0.09020330907998576, 0.5328799356978597,
            -2.5706086028385933, 5.467449792050247)),
    (0.0, (0.6435873464940464, 0.09524830366994964, 0.4796848403324965, 0.0, 0.0)),
])
def test_bivariate_against_recorded_bipred(intercept, expected):
    # Reference values from bipred fa211a1; the signed-z route is independent of
    # a runtime dependency on that package.
    ell, n1, n2, z1, z2 = _data()
    r = ldsc_rg_z(z1, z2, ell, n1, n2, m_snps=400, n_blocks=10,
                  constrain_intercept=intercept)
    np.testing.assert_allclose((r.rg, r.rg_se, r.gcov, r.gcov_intercept,
                                r.gcov_intercept_se), expected, rtol=2e-9, atol=2e-10)


def test_bivariate_nearly_constant_design_falls_back_stably(monkeypatch):
    import importlib
    module = importlib.import_module("pldsc.bivariate")
    t = np.linspace(-1, 1, 80)
    ell = 1 + 1e-8*t
    z1, z2 = np.sqrt(4 + 2*t), np.sqrt(3 + t)
    args = (z1, z2, ell, 10000, 12000)
    accelerated = ldsc_rg_z(*args, m_snps=400, n_blocks=8)
    monkeypatch.setattr(module, "HAVE_NUMBA", False)
    reference = ldsc_rg_z(*args, m_snps=400, n_blocks=8)
    _assert_results_close(accelerated, reference)
    with pytest.raises(ValueError, match="singular"):
        ldsc_rg_z(z1, z2, np.ones(80), 10000, 12000, m_snps=400)


def _tables():
    ell, n1, n2, z1, z2 = _data()
    m = len(ell)
    pos = 1000 + np.arange(m)
    ids = np.array([f"rs{i}" for i in range(m)])
    ref = allelign.Variants(np.full(m, "1"), ids, pos, np.full(m, "A"), np.full(m, "C"))
    def table(z, n):
        return allelign.Sumstats(ids.copy(), np.full(m, "1"), pos.copy(),
                                 np.full(m, "A"), np.full(m, "C"), z * 0.01,
                                 np.full(m, 0.01), n, np.full(m, 0.2), np.ones(m))
    s1, s2 = table(z1, n1), table(z2, n2)
    # These are the same biological effects recorded against the other allele.
    s2.ea[::2], s2.oa[::2] = "C", "A"
    s2.beta[::2] *= -1
    return ref, ell, s1, s2


def test_allelign_file_qc_and_bivariate_allele_flips(tmp_path):
    ref, ell, s1, s2 = _tables()
    # A duplicate is removed entirely, and a second row fails MAF QC.
    s2 = s2.subset(np.r_[np.arange(len(s2)), 7])
    s2.eaf[10] = .001
    f1, f2 = tmp_path / "trait1.tsv.gz", tmp_path / "trait2.tsv.gz"
    allelign.write_sumstats(f1, s1)
    allelign.write_sumstats(f2, s2)
    out = pldsc.ldsc_rg_sumstats(f1, f2, ref, ell, m_snps=500,
                                 n_blocks=10, chi2_cap=False)
    ix = np.delete(np.arange(len(ref)), [7, 10])
    np.testing.assert_array_equal(out.reference_index, ix)
    assert out.n_common == len(ix)
    assert out.qc[1]["n_drop_duplicate"] >= 2
    z1 = s1.beta[ix] / s1.se[ix]
    z2 = s2.beta[ix] / s2.se[ix]
    z2[ix % 2 == 0] *= -1
    expected = ldsc_rg_z(z1, z2, ell[ix], s1.n_eff[ix], s2.n_eff[ix], m_snps=500, n_blocks=10)
    _assert_results_close(out.estimate, expected)
    single = pldsc.ldsc_h2_sumstats(f1, ref, ell, m_snps=500, n_blocks=10, chi2_cap=False)
    direct = ldsc_h2((s1.beta/s1.se)**2, ell, s1.n_eff, m_snps=500, n_blocks=10)
    _assert_results_close(single.estimate, direct)


def test_allelign_realistic_formats_and_sample_size_sources(tmp_path):
    ref, ell, s1, s2 = _tables()
    ref.a2[6] = s1.oa[6] = s2.oa[6] = "T"  # palindromic: dropped by default
    s2.ea[6] = "A"
    n_cases = np.linspace(3000, 3500, len(ref))
    n_controls = np.linspace(5000, 5500, len(ref))
    n2 = 4 / (1 / n_cases + 1 / n_controls)
    n1 = s1.n_eff.copy()
    n1[9] = 100  # below the default relative-N filter
    bim = tmp_path / "reference.bim"
    bim.write_text("".join(f"chr1 {ref.id[i]} 0 {ref.pos[i]} {ref.a1[i]} {ref.a2[i]}\n"
                           for i in range(len(ref))))
    f1, f2 = tmp_path / "custom.tsv", tmp_path / "case_control.csv"
    lines1 = ["SNP\tCHR\tBP\tALT\tREF\tBETA\tEFFECT_USED\tSE_USED\tACTUAL_N\tEFFECTIVE_N\tEAF\tINFO\tTEST\n"]
    lines2 = ["SNP,CHR,BP,EA,OA,OR,SE,N,COUNT_CASE,COUNT_CONTROL,EAF,INFO,TEST\n"]
    for i in range(len(ref)):
        test = "DOM" if i == 4 else ("ADD-CONDTL" if i == 8 else "ADD")
        info = 0.1 if i == 7 else 0.9
        lines1.append(f"{s1.id[i]}\t1\t{s1.pos[i]}\t{s1.ea[i]}\t{s1.oa[i]}\t5\t{s1.beta[i]}\t{s1.se[i]}\t{n1[i]}\t9999\t0.2\t{info}\t{test}\n")
        ea, oa = ("T", "G") if i == 11 else (s2.ea[i], s2.oa[i])
        pos = s2.pos[i] + (10000 if i == 10 else 0)
        lines2.append(f"{s2.id[i]},1,{pos},{ea},{oa},{np.exp(s2.beta[i])},{s2.se[i]},123456,{n_cases[i]},{n_controls[i]},0.2,0.9,ADD\n")
    lines1.append(lines1[6])  # duplicate rs5: both copies are excluded
    f1.write_text("".join(lines1))
    f2.write_text("".join(lines2))
    opts = {"beta": "EFFECT_USED", "se": "SE_USED", "n_eff": "ACTUAL_N"}
    with pytest.warns(RuntimeWarning):
        single = pldsc.ldsc_h2_sumstats(f1, bim, ell, m_snps=500, read_options=opts,
                                        n_blocks=10, chi2_cap=False)
    ix1 = np.delete(np.arange(len(ref)), [4, 5, 6, 7, 9])
    np.testing.assert_array_equal(single.reference_index, ix1)
    _assert_results_close(single.estimate, ldsc_h2((s1.beta[ix1] / s1.se[ix1]) ** 2,
                          ell[ix1], n1[ix1], m_snps=500, n_blocks=10))
    assert single.qc[0]["sumstats_parse"]["n_eff_source"].startswith("column")
    assert single.qc[0]["sumstats_parse"]["n_dropped_test"] == 1
    assert single.qc[0]["n_drop_duplicate"] == 2
    with pytest.warns(RuntimeWarning):
        out = pldsc.ldsc_rg_sumstats(f1, f2, bim, ell, m_snps=500, read_options1=opts,
                                     read_options2={"n_cases": "COUNT_CASE", "n_controls": "COUNT_CONTROL"},
                                     n_blocks=10, chi2_cap=False)
    ix = np.delete(ix1, np.flatnonzero(ix1 == 10))  # coordinate disagreement
    np.testing.assert_array_equal(out.reference_index, ix)
    assert out.alignment[1]["n_strand_flipped"] == 1
    assert out.qc[1]["sumstats_parse"]["n_eff_source"] == "case_control"
    z2 = s2.beta[ix] / s2.se[ix]
    z2[ix % 2 == 0] *= -1
    expected = ldsc_rg_z(s1.beta[ix] / s1.se[ix], z2, ell[ix], n1[ix], n2[ix],
                         m_snps=500, n_blocks=10)
    _assert_results_close(out.estimate, expected)


def test_sumstats_cap_order_and_explicit_variant_count():
    ref, ell, s1, _ = _tables()
    s1.beta[4] = 0.2
    out = pldsc.ldsc_h2_sumstats(s1, ref, ell, m_snps=500, n_blocks=10)
    assert 4 not in out.reference_index
    assert out.n_chi2_dropped >= 1
    with pytest.raises(ValueError, match="m_snps"):
        pldsc.ldsc_h2_sumstats(s1, ref, ell, m_snps=None)
    ref.pos[[3, 4]] = ref.pos[[4, 3]]
    with pytest.raises(ValueError, match="sorted"):
        pldsc.ldsc_h2_sumstats(s1, ref, ell, m_snps=500)
