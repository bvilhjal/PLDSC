#!/usr/bin/env python3
"""Build the bundled ``hm3-eur`` LD reference in ``pldsc/data/``.

Inputs are the two archives of the LDSC reference deposit (Bulik-Sullivan
et al., Zenodo record 7768714): the 1000 Genomes Phase 3 EUR LD scores
printed on HapMap3 SNPs (1 cM windows, summed over the 9,997,231 variants
of the EUR QC panel) and the matching panel allele frequencies, which
supply the alleles the score files lack. Downloads are MD5-checked
against the record and cached under ``.scratch_ld_reference/``.

The output is one compressed ``.npz`` (variant table, minor allele
frequencies and scores in reference order) and a provenance ``.json``.
Printed SNPs whose score is not positive (53 in this deposit, an artefact
of the unbiased score estimator on a 489-person panel) are dropped from
the regression set and listed in the provenance; ``m_snps`` is unchanged,
because it counts the variants summed into the scores, not the rows kept.

Run from the repository root, on the maintainer's machine only::

    python scripts/build_ld_reference.py
"""
import argparse
import gzip
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

RECORD = "https://zenodo.org/api/records/7768714/files/{name}/content"
ARCHIVES = {
    "1000G_Phase3_ldscores.tgz": "e4352ccf778f296835d73985350a863b",
    "1000G_Phase3_frq.tgz": "ac29686ffd5b6378789857a522ebca77",
}
USER_AGENT = "pldsc-build-ld-reference (+https://github.com/bvilhjal/PLDSC)"


def fetch(name, md5_expected, cache_dir):
    """Return the cached archive path, downloading and MD5-checking if needed."""
    path = cache_dir / name
    if path.exists() and hashlib.md5(path.read_bytes()).hexdigest() == md5_expected:
        return path
    cache_dir.mkdir(parents=True, exist_ok=True)
    url = RECORD.format(name=name)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    print(f"downloading {url}")
    with urllib.request.urlopen(request, timeout=900) as response, path.open("wb") as out:
        while chunk := response.read(1 << 20):
            out.write(chunk)
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    if md5 != md5_expected:
        raise ValueError(f"{name}: md5 {md5} does not match the Zenodo record")
    return path


def member_text(tar, member):
    handle = tar.extractfile(member)
    if handle is None:
        raise ValueError(f"archive member {member} is not a regular file")
    return handle.read().decode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", default=".scratch_ld_reference",
                        help="download cache (default: .scratch_ld_reference)")
    parser.add_argument("--out-dir", default="pldsc/data",
                        help="output directory (default: pldsc/data)")
    args = parser.parse_args()
    cache_dir, out_dir = Path(args.cache_dir), Path(args.out_dir)
    paths = {name: fetch(name, md5, cache_dir) for name, md5 in ARCHIVES.items()}

    chrom, ids, pos, a1, a2, maf, scores = [], [], [], [], [], [], []
    dropped_nonpositive, missing_alleles = [], 0
    m_snps = m_snps_5_50 = 0
    with tarfile.open(paths["1000G_Phase3_frq.tgz"]) as frq_tar, \
            tarfile.open(paths["1000G_Phase3_ldscores.tgz"]) as ld_tar:
        for c in range(1, 23):
            alleles = {}
            frq_lines = member_text(
                frq_tar, f"1000G_Phase3_frq/1000G.EUR.QC.{c}.frq").splitlines()
            header = frq_lines[0].split()
            i = {name: header.index(name) for name in ("SNP", "A1", "A2", "MAF")}
            for line in frq_lines[1:]:
                f = line.split()
                if f[i["SNP"]] in alleles:
                    raise ValueError(f"chr{c}: duplicate frequency row for {f[i['SNP']]}")
                alleles[f[i["SNP"]]] = (f[i["A1"]], f[i["A2"]], float(f[i["MAF"]]))
            m_snps += int(member_text(ld_tar, f"LDscore/LDscore.{c}.l2.M"))
            m_snps_5_50 += int(member_text(ld_tar, f"LDscore/LDscore.{c}.l2.M_5_50"))
            raw = ld_tar.extractfile(f"LDscore/LDscore.{c}.l2.ldscore.gz").read()
            lines = gzip.decompress(raw).decode("utf-8").splitlines()
            header = lines[0].split()
            j = {name: header.index(name) for name in ("CHR", "SNP", "BP", "L2")}
            previous_pos = 0
            for line in lines[1:]:
                f = line.split()
                if f[j["CHR"]] != str(c):
                    raise ValueError(f"chr{c}: score file row on chromosome {f[j['CHR']]}")
                snp, bp, score = f[j["SNP"]], int(f[j["BP"]]), float(f[j["L2"]])
                if bp < previous_pos:
                    raise ValueError(f"chr{c}: score file positions are not sorted")
                previous_pos = bp
                if snp not in alleles:
                    missing_alleles += 1
                    continue
                if not np.isfinite(score) or score <= 0:
                    dropped_nonpositive.append(snp)
                    continue
                ref_a1, ref_a2, ref_maf = alleles[snp]
                if len(ref_a1) != 1 or len(ref_a2) != 1 or ref_a1 == ref_a2:
                    raise ValueError(f"{snp}: unexpected alleles {ref_a1}/{ref_a2}")
                chrom.append(str(c)); ids.append(snp); pos.append(bp)
                a1.append(ref_a1); a2.append(ref_a2); maf.append(ref_maf)
                scores.append(score)
    if missing_alleles:
        raise ValueError(f"{missing_alleles} printed SNPs have no frequency-row alleles")

    if max(map(len, ids)) > 16:
        raise ValueError("variant IDs longer than 16 characters would be truncated")
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir / "hm3_eur.npz"
    np.savez_compressed(
        npz_path,
        chrom=np.asarray(chrom, dtype="S2"), id=np.asarray(ids, dtype="S16"),
        pos=np.asarray(pos, dtype=np.int32), a1=np.asarray(a1, dtype="S1"),
        a2=np.asarray(a2, dtype="S1"), maf=np.asarray(maf, dtype=np.float32),
        scores=np.asarray(scores, dtype=np.float32))
    provenance = {
        "name": "hm3-eur",
        "description": ("HapMap3 regression SNPs with 1000 Genomes Phase 3 EUR LD "
                        "scores (1 cM windows), as published for LD Score regression."),
        "ancestry": "EUR",
        "genome_build": "GRCh37",
        "panel": "1000G EUR QC (489 individuals)",
        "n_variants": len(ids),
        "m_snps": m_snps,
        "m_snps_5_50": m_snps_5_50,
        "n_dropped_nonpositive_score": len(dropped_nonpositive),
        "dropped_nonpositive_score_ids": dropped_nonpositive,
        "allele_convention": "a1 is the panel minor allele (PLINK .frq convention); maf is its frequency",
        "source": {
            "record": "https://doi.org/10.5281/zenodo.7768714",
            "files": {name: {"md5": md5, "url": RECORD.format(name=name)}
                      for name, md5 in ARCHIVES.items()},
        },
        "built_by": "scripts/build_ld_reference.py",
        "bundle_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
    }
    (out_dir / "hm3_eur.json").write_text(json.dumps(provenance, indent=1) + "\n")
    scores_arr = np.asarray(scores)
    print(f"wrote {npz_path} ({npz_path.stat().st_size / 1e6:.1f} MB): "
          f"{len(ids)} variants, m_snps={m_snps}, m_snps_5_50={m_snps_5_50}, "
          f"{len(dropped_nonpositive)} dropped (non-positive score)")
    print(f"LD scores: mean {scores_arr.mean():.2f}, "
          f"range {scores_arr.min():.3f} to {scores_arr.max():.1f}")


if __name__ == "__main__":
    main()
