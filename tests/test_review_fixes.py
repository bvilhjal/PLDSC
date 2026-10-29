"""Regressions for the numerical and inference defects in the review."""
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

from pldsc import (cross_ancestry_ld_scores, effective_pair_df, ld_score_regression,
                   meta_ld_scores, meta_weights, pooled_ld_scores, single_reference_bias)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import gbmi_asthma_casualty as G  # noqa: E402
import multiancestry_genome_kappa as GK  # noqa: E402


def test_gbmi_uses_allelign_with_original_selection_and_orientation(tmp_path, monkeypatch):
    header = ("#CHR\tPOS\tREF\tALT\trsid\tall_meta_AF\tinv_var_meta_beta\t"
              "inv_var_meta_sebeta\tN_case\tN_ctrl\tis_strand_flip\tis_diff_AF_gnomAD\n")
    def row(rs, ref="A", alt="C", beta=.2, se=.1, af=.2, flag="no", pos=100):
        return f"1\t{pos}\t{ref}\t{alt}\t{rs}\t{af}\t{beta}\t{se}\t99999\t99999\t{flag}\tno\n"
    all_rows = [row("rs1"), row("rs2", alt="T", se=.2, af=.4),
                row("rs3"), row("rs3", flag="yes"), row("rs4"),
                row("rs5", flag="yes"), row("rs6", af=0),
                row("rs7", se=1), row("rs8", ref="AA")]
    eur_rows = [row("rs1", ref="C", alt="A", beta=-.3, pos=999),
                row("rs2", ref="T", alt="A", beta=-.4),
                row("rs4", ref="T", alt="G"),  # complement only: excluded
                row("rs7", beta=.4, se=1), row("rs9"), row("other")]
    for name, rows in (("all", all_rows), ("eur", eur_rows)):
        with gzip.open(tmp_path / f"{name}.tsv.gz", "wt") as out:
            out.write(header + "".join(rows))
    monkeypatch.setattr(G, "GBMI_DIR", str(tmp_path))
    monkeypatch.setattr(G, "PATTERN", "{}.tsv.gz")
    monkeypatch.setattr(G, "STRATA", {"ALL": "all", "EUR": "eur"})
    ids = np.array([f"rs{i}" for i in range(1, 10)])
    z, n = G.load_gbmi({"rsid": ids})
    np.testing.assert_array_equal(np.flatnonzero(np.isfinite(z[0])), [0, 1, 3, 6])
    np.testing.assert_array_equal(np.flatnonzero(np.isfinite(z[1])), [0, 1, 6])
    np.testing.assert_allclose(z[0, [0, 1, 3, 6]], [2, 1, 2, .2])
    np.testing.assert_allclose(z[1, [0, 1, 6]], [3, 4, .4])
    expected_n = 4 / (np.array([.1, .2, .1, 1.0]) ** 2 * 2 *
                      np.array([.2, .4, .2, .2]) * (1 - np.array([.2, .4, .2, .2])))
    np.testing.assert_allclose(n[0, [0, 1, 3, 6]], expected_n)
    assert n[0, 6] < 0.7 * np.median(expected_n)  # historical rule has no relative-N filter


def test_rank_deficient_covariates_do_not_change_pair_df():
    n = 20
    rng = np.random.default_rng(48)
    a = rng.normal(size=(n, n))
    A = a @ a.T
    c = rng.normal(size=n)
    redundant = np.column_stack((c, c, np.ones((n, 9))))
    assert effective_pair_df(np.eye(n), redundant) == pytest.approx(n - 2)
    assert effective_pair_df(A, redundant) == pytest.approx(effective_pair_df(A, c), rel=1e-12)
    with pytest.raises(ValueError, match="residual degrees"):
        effective_pair_df(np.eye(n), np.eye(n))


def test_wls_nearly_constant_design_and_nonfinite_bias():
    x = 1 + 1e-8 * np.linspace(-1, 1, 200)
    y = 2 + .7 * (x - 1) / 1e-8
    fit = ld_score_regression(y, x, n_blocks=10)
    assert fit.slope == pytest.approx(7e7, rel=1e-10)
    np.testing.assert_allclose(fit.intercept + fit.slope * x, y, atol=2e-8)
    with pytest.raises(ValueError, match="finite"):
        single_reference_bias([1, np.nan, 3], [1, 2, 3])
    with pytest.raises(ValueError, match="finite"):
        single_reference_bias([1, 2, 3], [1, np.inf, 3])


def test_rejects_invalid_correlation_matrix_but_accepts_singular():
    scores = np.ones((3, 3, 1))
    omega = scores / 9
    rho = np.full((3, 3), -.9)
    np.fill_diagonal(rho, 1)
    with pytest.raises(ValueError, match="positive semidefinite"):
        meta_ld_scores(scores, omega, rho=rho)
    assert meta_ld_scores(scores, omega, rho=np.ones((3, 3)))[0] == pytest.approx(1)


def test_pooled_identity_with_and_without_an_extra_covariate():
    rng = np.random.default_rng(420)
    n, m = 1000, 8
    group = np.repeat([0, 1], n)
    c = rng.normal(size=2*n)
    effect = np.array([1.4, 1.2, .7, .4, 0, 0, 0, 0])
    p = 1 / (1 + np.exp(-(-.4 + c[:, None] * effect)))
    geno = rng.binomial(2, p).astype(float)
    scores = cross_ancestry_ld_scores([geno[:n], geno[n:]], np.arange(m),
                                      window_cm=10, finite_sample_correction=False)
    omega, _ = meta_weights([n, n], scores.available, variances=scores.variances, scheme="ivw")
    meta = meta_ld_scores(scores, omega)
    one = pooled_ld_scores(geno, group, np.arange(m), covariates=(group == 1), window_cm=10)
    two = pooled_ld_scores(geno, group, np.arange(m),
                           covariates=np.column_stack((group == 1, c)), window_cm=10)
    np.testing.assert_allclose(one, meta, atol=1e-12)
    assert np.max(np.abs(two - meta)) > .5


def test_joint_jackknife_sufficient_sums_match_brute_refit():
    rng = np.random.default_rng(627)
    m, B, k = 150, 10, 3
    starts = np.arange(0, m, m//B)
    x = rng.uniform(.5, 5, m)
    y = .8 + .3*x + rng.normal(size=m)
    w = rng.uniform(.2, 1, m)
    derivatives = rng.normal(size=(k, m)) * .3
    deltas = rng.normal(size=(B, k)) * .05
    got = G.varying_predictor_deletions(x, y, w, starts, derivatives, deltas)
    expected = []
    for b, lo in enumerate(starts):
        keep = np.ones(m, bool)
        keep[lo:lo+m//B] = False
        xb = x + deltas[b] @ derivatives
        X = np.column_stack([np.ones(keep.sum()), xb[keep]])
        beta = np.linalg.lstsq(X*np.sqrt(w[keep, None]), y[keep]*np.sqrt(w[keep]), rcond=None)[0]
        expected.append(beta[1])
    np.testing.assert_allclose(got, expected, atol=1e-12)


def test_bracket_respects_estimated_and_deleted_correlations():
    rho = np.eye(5)
    usable = np.eye(5, dtype=bool)
    rho[0, 2] = rho[2, 0] = .755
    usable[0, 2] = usable[2, 0] = True
    deleted = np.repeat(rho[None, :, :], 2, axis=0)
    deleted[0, 0, 2] = deleted[0, 2, 0] = .74
    upper = G.largest_psd_bracket(rho, usable, deleted)
    assert .89 < upper < .92
    for matrix in (rho, *deleted):
        assert np.linalg.eigvalsh(np.where(usable, matrix, upper))[0] >= -1e-10
        assert np.linalg.eigvalsh(np.where(usable, matrix, 1))[0] < 0


def test_cache_binding_rejects_mixed_or_changed_maps(tmp_path, monkeypatch):
    out, maps = tmp_path / "out", tmp_path / "maps"
    out.mkdir(); maps.mkdir()
    prefix = tmp_path / "panel"
    for ext in GK.CHECKSUMS:
        Path(f"{prefix}.{ext}").write_bytes(ext.encode())
    for c in range(1, 23):
        (maps / f"chr{c}.b37.gmap.gz").write_bytes(f"map{c}".encode())
    monkeypatch.setattr(GK, "OUT", str(out))
    monkeypatch.setattr(GK, "MAPS", str(maps))
    monkeypatch.setattr(GK, "BFILE", str(prefix))
    import hashlib
    hashes = {k: GK.sha256(f"{prefix}.{k}") for k in GK.CHECKSUMS}
    data = hashlib.sha256("".join(hashes[k] for k in sorted(hashes)).encode()).hexdigest()
    for c in range(1, 23):
        np.savez(out / f"chr{c}.npz", config=np.array(GK.config(data, c)),
                 scores=np.ones((5, 5, 1)), available=np.ones((5, 1), bool),
                 variances=np.ones((5, 1)), maf=np.full((5, 1), .2),
                 maf_full=np.array([.2]), pos=np.array([c]), rsid=np.array([f"rs{c}"]),
                 n_reference=np.array(1))
    assert G.load_scores()["M"] == 22
    (maps / "chr7.b37.gmap.gz").write_bytes(b"changed")
    with pytest.raises(ValueError, match="chr7"):
        G.load_scores()
    for c in range(1, 23):
        path = out / f"chr{c}.npz"
        with np.load(path, allow_pickle=False) as d:
            payload = {key: d[key] for key in d.files if key != "config"}
        np.savez(path, config=np.array(GK.legacy_config(data)), **payload)
    with pytest.raises(ValueError, match="legacy"):
        G.load_scores()
    assert G.load_scores(allow_legacy=True)["M"] == 22
