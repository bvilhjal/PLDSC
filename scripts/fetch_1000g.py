#!/usr/bin/env python3
"""Fetch the 1000 Genomes phase-3 data the real-data experiments need.

A fresh clone of this repository is 100% data-free: `.scratch_1000g/` is
gitignored and no downloader used to exist, so the three data-carrying
experiments (`multiancestry_1000g_meta.py`, `realdata_genome_loco.py`,
`realdata_1000g_benchmarks.py`) all died instantly on missing files. This
script is the documented path back (see experiments/README.md and
docs/MULTI_ANCESTRY.md "Data").

What it fetches, and from where:

  VCFs         1000G phase 3 (20130502 release), one per autosome. The S3
               mirror is primary because the EBI FTP was unreachable from
               the 2026-09-25 environment. THE TWO MIRRORS SERVE DIFFERENT
               RELEASES OF THE SAME CALLSET: S3 has the original `v5a`
               files, EBI has `v5b` (a later re-release). Every recorded
               result in this repo was computed from the S3 `v5a` bytes, so
               v5a is the default and its sha256 is verified; the EBI `v5b`
               fallback is accepted with a printed warning that numbers may
               differ in the last digits. Both are written to the expected
               local name `ALL.chrN.vcf.gz` (the docstring recipe in
               multiancestry_1000g_meta.py used a manual symlink for this).
  panel        integrated_call_samples_v3.20130502.ALL.panel -> samples.panel.
  maps         the SHAPEIT4 GRCh37 genetic map tarball (odelaneau/shapeit4),
               extracted per chromosome to `chrN.b37.gmap.gz` (the format
               multiancestry_1000g_meta.py reads), AND converted to plink
               format `plink.chrN.GRCh37.map` (the format
               realdata_genome_loco.py reads). The conversion is needed
               because the historical Beagle host for the plink-format maps
               (bochet.gcc.biostat.washington.edu) returned 404 on
               2026-09-25; the generated files carry a provenance header
               and are a marker-for-marker reformatting of the same map, so
               cM windows agree to interpolation resolution.
  baselineLD   v2.2 annotation tarball (Alkes group; the host moved from
               data.broadinstitute.org to alkesgroup.broadinstitute.org --
               both are tried). Only needed for the D1 annotation gate work.
  cache        --build-cache re-runs the thinned VCF parser in
               experiments/realdata_1000g_benchmarks.py (pure Python; budget
               minutes per chromosome).

Checksums below are sha256 recorded from the local copies that produced the
results in docs/MULTI_ANCESTRY.md and the 2026-08-11 gate runs; files that
match are skipped, so the script is resumable.

Usage:
  python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel
  python scripts/fetch_1000g.py --chrom 1 ... 22 --maps all --build-cache 1 ... 22
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRATCH = os.path.join(ROOT, ".scratch_1000g")
MAPS = os.path.join(SCRATCH, "maps")
CACHE = os.path.join(SCRATCH, "cache")

S3 = "https://1000genomes.s3.amazonaws.com/release/20130502"
EBI = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_A = ("ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a."
         "20130502.genotypes.vcf.gz")
VCF_B = ("ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b."
         "20130502.genotypes.vcf.gz")
SHAPEIT4_MAPS = ("https://raw.githubusercontent.com/odelaneau/shapeit4/"
                 "master/maps/genetic_maps.b37.tar.gz")
BASELINELD = ("https://alkesgroup.broadinstitute.org/LDSCORE/"
              "baselineLD_v2.2.tgz")
BASELINELD_OLD = ("https://data.broadinstitute.org/alkesgroup/LDSCORE/"
                  "baselineLD_v2.2.tgz")

# sha256 of the copies that produced the recorded results (v5a mirror).
CHECKSUMS = {
    "ALL.chr10.vcf.gz": "98ac7174805ad0f7e07a632a673ad9e6b902a330ab3f17bdffe9ecd4e47408e1",
    "ALL.chr21.vcf.gz": "1942e0709625d77928b191e5a1bfbf523aad2fee1764994a599ae4623305cace",
    "ALL.chr22.vcf.gz": "a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5",
    "samples.panel": "b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858",
    "baselineLD_v2.2.tgz": "f7f4596eedf7e7598dc153203abbbb59539e664d0ae5caa1c195bc3a6031cac4",
    "genetic_maps.b37.tar.gz": "25a5eab7b4a9b47203cd2889fc77b7941242bdbae934aa7a7f3e1ad3184346e0",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, expect_sha=None, label=""):
    """Stream url to dest (resumable at file granularity, not byte ranges)."""
    if os.path.exists(dest) and (expect_sha is None or sha256(dest) == expect_sha):
        print(f"  [ok] {os.path.basename(dest)} already present")
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    print(f"  [get] {label or os.path.basename(dest)} <- {url}")
    with urllib.request.urlopen(url, timeout=60) as r, open(tmp, "wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    if expect_sha is not None:
        got = sha256(tmp)
        if got != expect_sha:
            os.remove(tmp)
            raise SystemExit(f"checksum mismatch for {dest}: {got} != {expect_sha}")
    os.replace(tmp, dest)
    print(f"  [done] {dest} ({os.path.getsize(dest) / 1e6:.0f} MB)")
    return dest


def fetch_vcf(chrom):
    """One autosome's phase-3 VCF, S3 v5a primary, EBI v5b fallback."""
    dest = os.path.join(SCRATCH, f"ALL.chr{chrom}.vcf.gz")
    try:
        fetch(f"{S3}/{VCF_A.format(c=chrom)}", dest,
              expect_sha=CHECKSUMS.get(f"ALL.chr{chrom}.vcf.gz"),
              label=f"chr{chrom} v5a (S3)")
        return
    except Exception as e:                      # noqa: BLE001 - fall through
        print(f"  [warn] S3 v5a failed ({e}); trying the EBI v5b fallback")
    fetch(f"{EBI}/{VCF_B.format(c=chrom)}", dest, label=f"chr{chrom} v5b (EBI)")
    print(f"  [warn] chr{chrom} is the v5b re-release; recorded numbers used "
          "v5a -- expect agreement only to the last digits")


def fetch_panel():
    fetch(f"{S3}/integrated_call_samples_v3.20130502.ALL.panel",
          os.path.join(SCRATCH, "samples.panel"),
          expect_sha=CHECKSUMS["samples.panel"], label="samples.panel")


def fetch_maps(chroms):
    """SHAPEIT4 gmap files, plus plink-format maps derived from them."""
    tar_path = fetch(SHAPEIT4_MAPS, os.path.join(MAPS, "genetic_maps.b37.tar.gz"),
                     expect_sha=CHECKSUMS["genetic_maps.b37.tar.gz"],
                     label="SHAPEIT4 b37 maps")
    os.makedirs(MAPS, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        for c in chroms:
            member = f"chr{c}.b37.gmap.gz"
            tf.extract(member, MAPS, filter="data")
    # plink-format maps for realdata_genome_loco.py: `chrom . cM bp` per line.
    # The Beagle host for these died (404, 2026-09-25), so reformat the same
    # SHAPEIT4 markers; read_plink_genetic_map skips short/header lines.
    import gzip
    for c in chroms:
        src = os.path.join(MAPS, f"chr{c}.b37.gmap.gz")
        dst = os.path.join(MAPS, f"plink.chr{c}.GRCh37.map")
        with gzip.open(src, "rt") as inp, open(dst, "w") as out:
            # keep the header to <= 3 whitespace fields; the reader skips
            # short lines and '#'-leading lines
            out.write("# SHAPEIT4-GRCh37-derived (scripts/fetch_1000g.py)\n")
            for i, line in enumerate(inp):
                if i == 0:                       # `pos chr cM` header
                    continue
                pos, _, cm = line.split()
                out.write(f"{c}\t.\t{cm}\t{pos}\n")
        print(f"  [ok] plink.chr{c}.GRCh37.map (converted from SHAPEIT4)")


def fetch_baselineld():
    try:
        fetch(BASELINELD, os.path.join(SCRATCH, "baselineLD_v2.2.tgz"),
              expect_sha=CHECKSUMS["baselineLD_v2.2.tgz"], label="baselineLD v2.2")
    except Exception as e:                      # noqa: BLE001 - fall through
        print(f"  [warn] alkesgroup host failed ({e}); trying the old data. host")
        fetch(BASELINELD_OLD, os.path.join(SCRATCH, "baselineLD_v2.2.tgz"),
              expect_sha=CHECKSUMS["baselineLD_v2.2.tgz"], label="baselineLD v2.2")


def build_cache(chroms):
    """Run the experiments' own thinned VCF parser (pure Python; minutes/chr)."""
    sys.path.insert(0, os.path.join(ROOT, "experiments"))
    import realdata_1000g_benchmarks as B      # noqa: E402
    for c in chroms:
        B.THIN.setdefault(c, 25)                # 1-in-10 only for chr22
        B.parse_and_cache(c)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chrom", type=int, nargs="+", default=[],
                    help="autosomes to fetch (e.g. 21 22, or 1 ... 22)")
    ap.add_argument("--panel", action="store_true", help="fetch samples.panel")
    ap.add_argument("--maps", nargs="*", default=None, metavar="CHR",
                    help="fetch genetic maps; 'all' = 1..22, no value = 21 22 10")
    ap.add_argument("--baselineld", action="store_true",
                    help="fetch the baselineLD v2.2 annotation tarball")
    ap.add_argument("--build-cache", type=int, nargs="+", default=[], metavar="CHR",
                    help="parse fetched VCFs into the thinned .scratch_1000g/cache")
    args = ap.parse_args()

    for c in args.chrom:
        fetch_vcf(c)
    if args.panel:
        fetch_panel()
    if args.maps is not None:
        chroms = list(range(1, 23)) if "all" in args.maps else \
            [int(x) for x in args.maps] if args.maps else [21, 22, 10]
        fetch_maps(chroms)
    if args.baselineld:
        fetch_baselineld()
    if args.build_cache:
        build_cache(args.build_cache)


if __name__ == "__main__":
    main()
