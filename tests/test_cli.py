"""End-to-end checks of command-line LDSC against the Python API."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import allelign
import numpy as np
import pytest

from pldsc import ldsc_h2_sumstats, ldsc_rg_sumstats


@pytest.fixture
def gwas_files(tmp_path):
    rng = np.random.default_rng(901)
    m = 120
    ids = np.array([f"rs{i}" for i in range(m)])
    pos = np.r_[1_000_000 + 1000 * np.arange(60), 30_000_000 + 1000 * np.arange(20),
                36_000_000 + 1000 * np.arange(40)]
    ell = rng.uniform(1, 4, m)
    n1, n2 = 8000., 7000.
    z1 = rng.normal(size=m) * np.sqrt(1 + n1 * .12 * ell / m)
    z2 = .6 * z1 + rng.normal(size=m) * np.sqrt(1 + n2 * .08 * ell / m)
    reference = tmp_path / "reference.bim"
    reference.write_text("".join(f"6\t{ids[i]}\t0\t{pos[i]}\tA\tC\n" for i in range(m)))
    scores = tmp_path / "scores.npy"
    np.save(scores, ell)
    def write(path, z, n, swap):
        ss = allelign.Sumstats(ids, np.full(m, "6"), pos, np.where(swap, "C", "A"),
                                np.where(swap, "A", "C"), np.where(swap, -z, z) * .01,
                                np.full(m, .01), np.full(m, n), np.full(m, .2), np.ones(m))
        ss.eaf[3] = .001
        allelign.write_sumstats(path, ss)
    one, two = tmp_path / "one.tsv.gz", tmp_path / "two.tsv.gz"
    write(one, z1, n1, np.zeros(m, bool))
    write(two, z2, n2, np.arange(m) % 2 == 0)
    return reference, scores, one, two, ell


def _command(*args):
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    return subprocess.run([sys.executable, "-m", "pldsc.cli", *map(str, args)],
                          text=True, capture_output=True, env=env, check=False)


def test_cli_fits_match_python_and_record_processing(gwas_files, tmp_path):
    reference, scores, one, two, ell = gwas_files
    common = ["--reference", reference, "--ld-scores", scores, "--m-snps", "120",
              "--n-blocks", "8", "--no-chi2-cap"]
    run = _command("h2", "--sumstats", one, "--column", "n_eff=N_EFF", *common)
    assert run.returncode == 0, run.stderr
    out = json.loads(run.stdout)
    expected = ldsc_h2_sumstats(one, reference, ell, m_snps=120, n_blocks=8,
                                 read_options={"n_eff": "N_EFF"}, chi2_cap=False)
    assert out["estimate"]["h2"] == pytest.approx(expected.estimate.h2)
    assert out["estimate"]["h2_se"] == pytest.approx(expected.estimate.h2_se)
    assert out["n_retained"] == len(expected.reference_index)
    assert out["qc"][0]["n_drop_lowMAF"] == 1
    assert out["options"]["read_options"][0]["n_eff"] == "N_EFF"

    target = tmp_path / "rg.json"
    run = _command("rg", "--sumstats1", one, "--sumstats2", two, "--n-eff2", "7000",
                   "--ncores", "2", "--out", target, *common)
    assert run.returncode == 0, run.stderr
    assert run.stdout == ""
    out = json.loads(target.read_text())
    expected = ldsc_rg_sumstats(one, two, reference, ell, m_snps=120, n_blocks=8,
                                 read_options2={"n_eff": 7000}, ncores=2, chi2_cap=False)
    assert out["estimate"]["rg"] == pytest.approx(expected.estimate.rg)
    assert out["estimate"]["rg_se"] == pytest.approx(expected.estimate.rg_se)
    assert out["estimate"]["gcov"] == pytest.approx(expected.estimate.gcov)
    assert out["alignment"][1]["n_flipped"] == expected.alignment[1]["n_flipped"]
    assert out["n_common"] == expected.n_common
    assert out["qc"][1]["sumstats_parse"]["n_eff_source"] == "scalar"
    previous = target.read_bytes()
    again = _command("rg", "--sumstats1", one, "--sumstats2", two, "--out", target, *common)
    assert again.returncode != 0
    assert target.read_bytes() == previous


def test_cli_optional_allelign_region_filter(gwas_files):
    reference, scores, one, two, ell = gwas_files
    args = ["--reference", reference, "--ld-scores", scores, "--m-snps", "120", "--n-blocks", "8",
            "--exclude-long-range-ld-grch37", "--include-apoe", "--no-chi2-cap"]
    run = _command("rg", "--sumstats1", one, "--sumstats2", two, *args)
    assert run.returncode == 0, run.stderr
    out = json.loads(run.stdout)
    assert [entry["n_dropped"] for entry in out["long_range_ld"]] == [20, 20]
    ss1, ss2 = allelign.read_sumstats(one), allelign.read_sumstats(two)
    mask = ~allelign.in_long_range_ld(ss1.chrom, ss1.pos, genome_build="GRCh37", include_apoe=True)
    expected = ldsc_rg_sumstats(ss1.subset(mask), ss2.subset(mask), reference, ell,
                                 m_snps=120, n_blocks=8, chi2_cap=False)
    assert out["estimate"]["rg"] == pytest.approx(expected.estimate.rg)
    assert out["n_retained"] == len(expected.reference_index)


def test_cli_rejects_conflicting_inputs_and_writes_strict_json(gwas_files, tmp_path):
    reference, scores, one, _, _ = gwas_files
    common = ["h2", "--sumstats", one, "--reference", reference, "--ld-scores", scores, "--m-snps", "120"]
    run = _command(*common, "--n-eff", "9000", "--column", "n_eff=N_EFF")
    assert run.returncode == 2 and "one sample-size source" in run.stderr
    run = _command(*common, "--column", "beta=BETA", "--column", "or=BETA")
    assert run.returncode == 2 and "single effect source" in run.stderr
    run = _command(*common, "--include-apoe")
    assert run.returncode == 2 and "requires" in run.stderr
    run = _command(*common, "--n-blocks", "1")
    assert run.returncode == 0, run.stderr
    out = json.loads(run.stdout, parse_constant=lambda value: pytest.fail(f"nonstandard JSON: {value}"))
    assert out["estimate"]["h2_se"] is None
    bad = tmp_path / "bad.npz"
    np.savez(bad, values=np.ones(120))
    run = _command("h2", "--sumstats", one, "--reference", reference, "--ld-scores", bad, "--m-snps", "120")
    assert run.returncode == 2 and ".npy" in run.stderr


def test_cli_region_filter_rejects_missing_coordinates(gwas_files, tmp_path):
    reference, scores, one, _, _ = gwas_files
    ss = allelign.read_sumstats(one)
    ss.pos[0] = 0
    missing = tmp_path / "missing.tsv"
    allelign.write_sumstats(missing, ss)
    run = _command("h2", "--sumstats", missing, "--reference", reference, "--ld-scores", scores,
                   "--m-snps", "120", "--exclude-long-range-ld-grch37")
    assert run.returncode == 2 and "requires valid chromosome and position" in run.stderr


def test_prepare_standard_ldsc_scores_aligns_and_rejects_bad_joins(gwas_files, tmp_path):
    reference, _, _, _, ell = gwas_files
    variants = allelign.read_variants(reference)
    order = np.random.default_rng(37).permutation(len(variants))
    def lines(indices):
        return "CHR SNP BP L2\n" + "".join(
            f"chr6 {variants.id[i]} {variants.pos[i]} {ell[i]}\n" for i in indices)
    first, second = tmp_path / "chr-a.ldscore.gz", tmp_path / "chr-b.ldscore"
    with gzip.open(first, "wt") as stream:
        stream.write(lines(order[:60]))
    second.write_text(lines(order[60:]) + "6 rs_other 42 1.5\n")
    out = tmp_path / "aligned.npy"
    args = ["prepare-scores", "--reference", reference, "--score-files", first, second, "--out", out]
    run = _command(*args)
    assert run.returncode == 0, run.stderr
    np.testing.assert_allclose(np.load(out), ell)
    assert json.loads(run.stdout)["n_ignored_score_rows"] == 1
    assert _command(*args).returncode == 2  # do not overwrite a score vector

    duplicate = tmp_path / "duplicate.ldscore"
    duplicate.write_text(lines([order[0]]))
    bad_out = tmp_path / "bad.npy"
    run = _command("prepare-scores", "--reference", reference, "--score-files", first, second, duplicate,
                   "--out", bad_out)
    assert run.returncode == 2 and "duplicate LD score" in run.stderr
    assert not bad_out.exists()
    run = _command("prepare-scores", "--reference", reference, "--score-files", first, "--out", bad_out)
    assert run.returncode == 2 and "missing for 60" in run.stderr
    mismatch = tmp_path / "mismatch.ldscore"
    i = order[0]
    mismatch.write_text(f"CHR SNP BP L2\n6 {variants.id[i]} {variants.pos[i] + 1} 1.0\n")
    run = _command("prepare-scores", "--reference", reference, "--score-files", mismatch, "--out", bad_out)
    assert run.returncode == 2 and "coordinates disagree" in run.stderr


def test_notebook_code_runs_in_a_fresh_python_process():
    project = Path(__file__).resolve().parents[1]
    notebook = project / "examples" / "ldsc_notebook.ipynb"
    program = """
import json, math, sys
notebook = json.load(open(sys.argv[1], encoding='utf-8'))
scope = {'__name__': '__main__'}
for cell in notebook['cells']:
    if cell['cell_type'] == 'code':
        exec(compile(''.join(cell['source']), '<notebook>', 'exec'), scope)
assert math.isfinite(scope['h2'].estimate.h2)
assert math.isfinite(scope['rg'].estimate.rg)
assert scope['rg'].alignment[1]['n_flipped'] > 0
assert len(scope['rg'].reference_index) > 3000
"""
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    run = subprocess.run([sys.executable, "-c", program, str(notebook)], cwd=project / "examples",
                         text=True, capture_output=True, env=env, check=False)
    assert run.returncode == 0, run.stderr
