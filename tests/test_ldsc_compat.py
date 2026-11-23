"""Parity of the ldsc-compatible estimator with the reference ldsc code.

``tests/data/ldsc_parity.npz`` holds synthetic inputs and the estimates of
bulik/ldsc's own estimator code (commit 2fdeeb3) on them; it is written by
``scripts/make_ldsc_parity_fixture.py``, which documents the Python 3 port.
"""
from pathlib import Path

import numpy as np
import pytest

from pldsc import ldsc_h2, ldsc_h2_compat, ldsc_rg_compat, ldsc_rg_z

FIXTURE = Path(__file__).resolve().parent / "data" / "ldsc_parity.npz"


@pytest.fixture(scope="module")
def ref():
    with np.load(FIXTURE, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


H2_CASES = {"twostep": {}, "onestep": {"twostep": None}, "fixed": {"constrain_intercept": 1.0}}
RG_CASES = {"twostep": {}, "onestep": {"twostep": None}, "fixed": {"constrain_intercept": 0.0, "twostep": None}}


@pytest.mark.parametrize("case", sorted(H2_CASES))
def test_h2_matches_reference_ldsc(ref, case):
    r = ldsc_h2_compat(ref["z1"] ** 2, ref["ld"], ref["N1"], m_snps=float(ref["M"]), w_ld=ref["w_ld"],
                       **H2_CASES[case])
    np.testing.assert_allclose([r.h2, r.h2_se, r.intercept, r.intercept_se], ref[f"h2_{case}"],
                               rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("case", sorted(RG_CASES))
def test_rg_matches_reference_ldsc(ref, case):
    r = ldsc_rg_compat(ref["z1"], ref["z2"], ref["ld"], ref["N1"], ref["N2"], m_snps=float(ref["M"]),
                       w_ld=ref["w_ld"], **RG_CASES[case])
    got = [r.rg, r.rg_se, r.gcov, r.gcov_intercept, r.gcov_intercept_se, *r.h2, r.gcov_se, *r.h2_se]
    np.testing.assert_allclose(got, ref[f"rg_{case}"], rtol=1e-9, atol=1e-12)


def test_public_entry_points_dispatch_and_guard(ref):
    chisq, ell, N, M = ref["z1"] ** 2, ref["ld"], ref["N1"], float(ref["M"])
    a = ldsc_h2(chisq, ell, N, m_snps=M, ldsc_compat=True, w_ld=ref["w_ld"])
    assert a.h2 == pytest.approx(ref["h2_twostep"][0], rel=1e-12)
    b = ldsc_rg_z(ref["z1"], ref["z2"], ell, N, ref["N2"], m_snps=M, ldsc_compat=True, w_ld=ref["w_ld"])
    assert b.rg == pytest.approx(ref["rg_twostep"][0], rel=1e-12)
    with pytest.raises(ValueError, match="m_snps"):
        ldsc_h2(chisq, ell, N, ldsc_compat=True)
    with pytest.raises(ValueError, match="w_ld"):
        ldsc_h2(chisq, ell, N, m_snps=M, w_ld=ref["w_ld"])
    with pytest.raises(ValueError, match="n_iter"):
        ldsc_h2(chisq, ell, N, m_snps=M, ldsc_compat=True, n_iter=3)
    # Default mode is unchanged: 100 blocks, one-step, no w_ld.
    default = ldsc_h2(chisq, ell, N, m_snps=M)
    assert default.h2 == pytest.approx(ldsc_h2(chisq, ell, N, m_snps=M, n_blocks=100).h2)
    assert default.h2_se == pytest.approx(ldsc_h2(chisq, ell, N, m_snps=M, n_blocks=100).h2_se)


def test_two_step_drops_large_chi_squares_from_the_intercept_only(ref):
    chisq = ref["z1"] ** 2
    assert np.sum(chisq >= 30) > 0
    full = ldsc_h2_compat(chisq, ref["ld"], ref["N1"], m_snps=float(ref["M"]), w_ld=ref["w_ld"])
    one = ldsc_h2_compat(chisq, ref["ld"], ref["N1"], m_snps=float(ref["M"]), w_ld=ref["w_ld"], twostep=None)
    assert full.intercept != pytest.approx(one.intercept)
