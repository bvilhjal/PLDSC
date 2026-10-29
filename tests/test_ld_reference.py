"""Bundled LD references: loader invariants and CLI parity."""
import hashlib
import json
import os
import subprocess
import sys
from importlib.resources import files

import allelign
import numpy as np
import pytest

from pldsc import (available_ld_references, ldsc_h2_sumstats,
                   load_ld_reference)


def _command(*args):
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1"}
    return subprocess.run([sys.executable, "-m", "pldsc.cli", *map(str, args)],
                          text=True, capture_output=True, env=env, check=False)


def test_bundled_hm3_eur_loads_with_recorded_counts():
    assert available_ld_references() == ["hm3-eur"]
    ref = load_ld_reference("hm3-eur", include_long_range_ld=True)
    assert ref.long_range_ld is None
    assert len(ref.variants) == 1_190_268
    assert ref.scores.shape == ref.maf.shape == (len(ref.variants),)
    assert ref.m_snps == 9_997_231 and ref.m_snps_5_50 == 5_961_159
    assert (ref.ancestry, ref.genome_build) == ("EUR", "GRCh37")
    assert np.all(np.isfinite(ref.scores)) and np.all(ref.scores > 0)
    assert np.all((ref.maf > 0) & (ref.maf <= 0.5))
    chrom = np.asarray(ref.variants.chrom)
    pos = np.asarray(ref.variants.pos)
    for c in map(str, range(1, 23)):
        block = pos[chrom == c]
        assert block.size and np.all(np.diff(block) > 0), c
    bundle = files("pldsc").joinpath("data", "hm3_eur.npz").read_bytes()
    assert hashlib.sha256(bundle).hexdigest() == ref.provenance["bundle_sha256"]
    assert ref.provenance["n_dropped_nonpositive_score"] == 53


def test_bundled_rows_match_the_source_files():
    # Values read from LDscore.{1,22}.l2.ldscore.gz and 1000G.EUR.QC.{1,22}.frq.
    ref = load_ld_reference("hm3-eur")
    ids = np.asarray(ref.variants.id)
    for snp, chrom, pos, a1, a2, maf, score in (
            ("rs3094315", "1", 752566, "G", "A", 0.1636, 80.826),
            ("rs28729663", "22", 51219006, "A", "G", 0.1595, 48.435)):
        i = np.flatnonzero(ids == snp)
        assert i.size == 1, snp
        i = i[0]
        assert (ref.variants.chrom[i], ref.variants.pos[i]) == (chrom, pos)
        assert (ref.variants.a1[i], ref.variants.a2[i]) == (a1, a2)
        assert ref.maf[i] == pytest.approx(maf, abs=1e-6)
        assert ref.scores[i] == pytest.approx(score, abs=1e-4)
    assert ids[0] == "rs3094315" and ids[-1] == "rs28729663"
    # Printed with a negative score (-3.598), so not a regression variant.
    assert "rs12046770" in ref.provenance["dropped_nonpositive_score_ids"]
    assert not np.any(ids == "rs12046770")


def test_long_range_ld_regions_are_left_out_by_default():
    full = load_ld_reference("hm3-eur", include_long_range_ld=True)
    in_region = allelign.in_long_range_ld(full.variants.chrom, full.variants.pos,
                                          genome_build="GRCh37")
    ref = load_ld_reference("hm3-eur")
    assert ref.long_range_ld == {"genome_build": "GRCh37", "include_apoe": False,
                                 "n_dropped": int(in_region.sum())}
    assert in_region.sum() > 0 and len(ref.variants) == len(full.variants) - in_region.sum()
    assert np.array_equal(ref.variants.id, full.variants.id[~in_region])
    assert np.array_equal(ref.scores, full.scores[~in_region])
    assert np.array_equal(ref.maf, full.maf[~in_region])
    assert (ref.m_snps, ref.m_snps_5_50) == (full.m_snps, full.m_snps_5_50)
    mhc = (ref.variants.chrom == "6") & (ref.variants.pos >= 25_000_000) & (ref.variants.pos <= 35_000_000)
    assert not mhc.any()
    apoe = load_ld_reference("hm3-eur", include_apoe=True)
    assert apoe.long_range_ld["n_dropped"] > ref.long_range_ld["n_dropped"]
    with pytest.raises(ValueError, match="no effect"):
        load_ld_reference("hm3-eur", include_long_range_ld=True, include_apoe=True)


def test_corrupted_bundle_is_rejected(tmp_path, monkeypatch):
    import pldsc.ld_reference as module
    data = tmp_path / "data"
    data.mkdir()
    source = files("pldsc").joinpath("data")
    bundle = bytearray(source.joinpath("hm3_eur.npz").read_bytes())
    bundle[len(bundle) // 2] ^= 0xFF
    (data / "hm3_eur.npz").write_bytes(bytes(bundle))
    (data / "hm3_eur.json").write_text(source.joinpath("hm3_eur.json").read_text())
    monkeypatch.setattr(module, "files", lambda package: tmp_path)
    with pytest.raises(ValueError, match="checksum"):
        load_ld_reference("hm3-eur")


def test_unknown_reference_is_rejected():
    with pytest.raises(ValueError, match="unknown LD reference"):
        load_ld_reference("hm3-afr")
    run = _command("ld-reference", "--name", "hm3-afr")
    assert run.returncode == 2 and "unknown LD reference" in run.stderr


def test_cli_ld_reference_describes_the_bundle():
    run = _command("ld-reference")
    assert run.returncode == 0, run.stderr
    out = json.loads(run.stdout)
    assert out["available"] == ["hm3-eur"]
    assert out["reference"]["name"] == "hm3-eur"
    assert out["reference"]["m_snps"] == 9_997_231
    assert out["long_range_ld"]["n_dropped"] > 0


@pytest.fixture
def bundled_gwas(tmp_path):
    reference = load_ld_reference("hm3-eur")
    rng = np.random.default_rng(417)
    variants_all = reference.variants
    a1, a2 = np.asarray(variants_all.a1), np.asarray(variants_all.a2)
    palindromic = (((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A"))
                   | ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C")))
    eligible = ((np.asarray(variants_all.chrom) == "1") & ~palindromic
                & (reference.maf >= 0.05))
    keep = np.flatnonzero(eligible)[:400]
    assert keep.size == 400
    variants = reference.variants
    ell = reference.scores[keep]
    z = rng.standard_normal(keep.size) * np.sqrt(1 + 50_000 * 0.25 * ell / reference.m_snps)
    ss = allelign.Sumstats(
        variants.id[keep], variants.chrom[keep], variants.pos[keep],
        variants.a1[keep], variants.a2[keep], z * 0.01, np.full(keep.size, 0.01),
        np.full(keep.size, 50_000.0), reference.maf[keep], np.ones(keep.size))
    path = tmp_path / "trait.tsv.gz"
    allelign.write_sumstats(path, ss)
    return reference, path


def test_cli_ld_reference_fit_matches_python(bundled_gwas):
    reference, path = bundled_gwas
    run = _command("h2", "--sumstats", path, "--ld-reference", "hm3-eur",
                   "--n-blocks", "8", "--no-chi2-cap")
    assert run.returncode == 0, run.stderr
    out = json.loads(run.stdout)
    expected = ldsc_h2_sumstats(path, reference.variants, reference.scores,
                                m_snps=reference.m_snps, n_blocks=8, chi2_cap=False)
    assert out["estimate"]["h2"] == pytest.approx(expected.estimate.h2)
    assert out["n_retained"] == len(expected.reference_index) == 400
    assert out["inputs"]["ld_reference"]["name"] == "hm3-eur"
    assert out["inputs"]["reference"] is None


def test_cli_ld_reference_conflicts_with_explicit_inputs(bundled_gwas, tmp_path):
    _, path = bundled_gwas
    run = _command("h2", "--sumstats", path, "--ld-reference", "hm3-eur",
                   "--reference", tmp_path / "unused.bim")
    assert run.returncode == 2 and "not both" in run.stderr
    run = _command("h2", "--sumstats", path)
    assert run.returncode == 2 and "give --reference" in run.stderr
    run = _command("h2", "--sumstats", path, "--ld-reference", "")
    assert run.returncode == 2 and "unknown LD reference" in run.stderr


def test_cli_m_snps_overrides_the_bundled_count(bundled_gwas):
    reference, path = bundled_gwas
    runs = [_command("h2", "--sumstats", path, "--ld-reference", "hm3-eur",
                     "--n-blocks", "8", "--no-chi2-cap", *extra)
            for extra in ((), ("--m-snps", reference.m_snps_5_50))]
    assert all(run.returncode == 0 for run in runs), [run.stderr for run in runs]
    bundled, common = (json.loads(run.stdout) for run in runs)
    assert bundled["inputs"]["ld_reference"]["m_snps"] == "bundled"
    assert common["inputs"]["ld_reference"]["m_snps"] == "--m-snps"
    assert common["options"]["m_snps"] == reference.m_snps_5_50
    assert common["estimate"]["h2"] == pytest.approx(
        bundled["estimate"]["h2"] * reference.m_snps_5_50 / reference.m_snps)
    assert common["estimate"]["intercept"] == pytest.approx(bundled["estimate"]["intercept"])
