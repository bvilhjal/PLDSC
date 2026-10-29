"""Warm paired timings and numerical parity against LDpred3 and bipred.

Run from a stable checkout on AC power with BLAS threads pinned:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python experiments/ldsc_provider_parity_bench.py
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pldsc  # noqa: E402


def check_power():
    if sys.platform != "darwin":
        print("Power state is not checked on this platform.")
        return
    battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
    settings = subprocess.check_output(["pmset", "-g"], text=True)
    if "Battery Power" in battery or any(
            line.strip().split()[:2] == ["lowpowermode", "1"] for line in settings.splitlines()):
        raise RuntimeError("benchmark requires AC power and Low Power Mode off")


def compare(a, b):
    for name, value in asdict(a).items():
        other = asdict(b)[name]
        if value is None:
            assert other is None
        else:
            np.testing.assert_allclose(value, other, rtol=2e-9, atol=2e-10, equal_nan=True)


def main():
    check_power()
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        if os.environ.get(name) != "1":
            raise RuntimeError(f"set {name}=1 before running the benchmark")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--variants", type=int, default=100_000)
    ap.add_argument("--blocks", type=int, default=100)
    ap.add_argument("--repeats", type=int, default=3)
    args = ap.parse_args()
    if min(args.variants, args.blocks, args.repeats) < 2:
        ap.error("variants, blocks and repeats must be at least 2")
    import ldpred3
    from ldpred3.ldsc import ldsc_h2 as reference_h2
    import bipred
    from bipred.ldsc import ldsc_rg as reference_rg
    from pldsc._numba import HAVE_NUMBA
    if not HAVE_NUMBA:
        raise RuntimeError("Numba is required for this accelerated comparison")
    import numba
    print(f"host={platform.platform()} Python={platform.python_version()} NumPy={np.__version__} "
          f"Numba={numba.__version__} PLDSC={pldsc.__version__} "
          f"LDpred3={ldpred3.__version__} bipred={bipred.__version__}")
    for label, path in (("pldsc univariate", Path(pldsc.__file__).parent / "ldsc.py"),
                        ("pldsc bivariate kernels", Path(pldsc.__file__).parent / "_ldsc_kernels.py"),
                        ("ldpred3 univariate", Path(ldpred3.__file__).parent / "ldsc.py"),
                        ("bipred bivariate", Path(bipred.__file__).parent / "ldsc.py"),
                        ("bipred bivariate kernels", Path(bipred.__file__).parent / "_ldsc_kernels.py")):
        print(f"{label}: sha256 {hashlib.sha256(path.read_bytes()).hexdigest()}")
    print(f"m={args.variants}, blocks={args.blocks}, repeats={args.repeats}, seed=723")
    rng = np.random.default_rng(723)
    m = args.variants
    ell = rng.uniform(1, 8, m)
    n1, n2 = rng.uniform(60_000, 100_000, (2, m))
    a, b = rng.normal(size=(2, m))
    z1 = a * np.sqrt(1.1 + n1 * .25 * ell / m)
    z2 = (.7*a + np.sqrt(1-.7**2)*b) * np.sqrt(1.2 + n2 * .3 * ell / m)
    beta1, beta2 = z1 / np.sqrt(n1 + z1*z1), z2 / np.sqrt(n2 + z2*z2)
    cases = [
        ("univariate", lambda: pldsc.ldsc_h2(z1*z1, ell, n1, n_blocks=args.blocks),
         lambda: reference_h2(z1*z1, ell, n1, n_blocks=args.blocks)),
        ("bivariate serial", lambda: pldsc.ldsc_rg(beta1, beta2, ell, n1, n2, n_blocks=args.blocks),
         lambda: reference_rg(beta1, beta2, ell, n1, n2, n_blocks=args.blocks)),
        ("bivariate 4 threads", lambda: pldsc.ldsc_rg(beta1, beta2, ell, n1, n2, n_blocks=args.blocks, ncores=4),
         lambda: reference_rg(beta1, beta2, ell, n1, n2, n_blocks=args.blocks, ncores=4)),
    ]
    for name, ours, reference in cases:
        compare(ours(), reference())  # warm both and verify
        elapsed = [[], []]
        for i in range(args.repeats):
            for j in ((0, 1) if i % 2 == 0 else (1, 0)):
                fn = (ours, reference)[j]
                t = time.perf_counter(); result = fn(); elapsed[j].append(time.perf_counter() - t)
                if j == 0:
                    last_ours = result
                else:
                    last_reference = result
            compare(last_ours, last_reference)
        print(f"{name}: PLDSC {np.median(elapsed[0]):.3f}s; "
              f"provider {np.median(elapsed[1]):.3f}s; parity passed")


if __name__ == "__main__":
    main()
