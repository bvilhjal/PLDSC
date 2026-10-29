"""Bundled LD-score references, ready for the sumstats fits.

A bundled reference is a variant table, its LD scores and the variant
count the scores were summed over, packaged so that

    ref = load_ld_reference("hm3-eur")
    ldsc_rg_sumstats(file1, file2, ref.variants, ref.scores, m_snps=ref.m_snps)

needs no downloads or file hunting. Bundles are built once from public
deposits by ``scripts/build_ld_reference.py``; each carries a provenance
record naming its source files and checksums. Scores are stored as
float32 and returned as float64; the regression SNP set, panel, genome
build and both published variant counts (``m_snps`` for every variant
summed into the scores, ``m_snps_5_50`` for the MAF 5-50% subset the
original LDSC software uses by default) are part of the record, because
the choice rescales h2 while leaving rg unchanged when applied
consistently.
"""
from __future__ import annotations

import hashlib
import io
import json
from dataclasses import dataclass
from importlib.resources import files

import allelign
import numpy as np

#: Bundled references by name, as (scores bundle, provenance) file names.
_REFERENCES = {"hm3-eur": ("hm3_eur.npz", "hm3_eur.json")}


@dataclass
class LDReference:
    """One bundled LD-score reference (see :func:`load_ld_reference`)."""

    name: str
    variants: allelign.Variants   # regression variants, in score order
    scores: np.ndarray            # one LD score per variant row
    maf: np.ndarray               # reference-panel minor allele frequency
    m_snps: int                    # variants summed into each score
    m_snps_5_50: int               # MAF 5-50% subset count (original LDSC default)
    ancestry: str
    genome_build: str
    provenance: dict
    long_range_ld: dict | None = None  # regions left out, or None if kept

    def __repr__(self):
        return (f"LDReference({self.name!r}: {len(self.variants)} variants, "
                f"m_snps={self.m_snps}, {self.ancestry}, {self.genome_build})")


def available_ld_references():
    """Names of the bundled LD references, sorted."""
    return sorted(_REFERENCES)


def _strings(values):
    return np.asarray(values).astype(str)


def load_ld_reference(name="hm3-eur", *, include_long_range_ld=False, include_apoe=False):
    """Load a bundled LD reference by name (see :func:`available_ld_references`).

    The returned scores are all finite and positive and align row-by-row
    with ``variants``. ``m_snps`` counts every variant summed into them;
    ``m_snps_5_50`` counts the MAF 5-50% subset that the original LDSC
    software uses by default. Either can be passed as ``m_snps`` to the
    fits; it sets which variants h2 refers to. Loading reads only package
    data, checks the recorded checksum and needs no network.

    By default, variants in allelign's long-range LD regions for the
    bundle's genome build (including the MHC) are left out of the
    regression set, and ``long_range_ld`` records how many; the variant
    counts are unchanged. ``include_long_range_ld=True`` keeps them.
    ``include_apoe=True`` adds the APOE region to the regions left out.
    """
    if include_long_range_ld and include_apoe:
        raise ValueError("include_apoe has no effect when long-range LD regions are kept")
    if name not in _REFERENCES:
        raise ValueError(
            f"unknown LD reference {name!r}; available: "
            + ", ".join(available_ld_references()))
    bundle_name, provenance_name = _REFERENCES[name]
    data = files("pldsc").joinpath("data")
    provenance = json.loads(data.joinpath(provenance_name).read_text())
    raw = data.joinpath(bundle_name).read_bytes()
    if hashlib.sha256(raw).hexdigest() != provenance["bundle_sha256"]:
        raise ValueError(f"{bundle_name}: checksum differs from its provenance record")
    with np.load(io.BytesIO(raw), allow_pickle=False) as bundle:
        arrays = {key: bundle[key] for key in bundle.files}
    n = len(arrays["id"])
    if any(len(arrays[key]) != n for key in ("chrom", "pos", "a1", "a2", "maf", "scores")):
        raise ValueError(f"{bundle_name}: bundled arrays disagree in length")
    if n != provenance["n_variants"]:
        raise ValueError(f"{bundle_name}: variant count disagrees with its provenance")
    scores = np.asarray(arrays["scores"], dtype=np.float64)
    if not np.all(np.isfinite(scores)) or np.any(scores <= 0):
        raise ValueError(f"{bundle_name}: scores must be finite and positive")
    variants = allelign.Variants(
        chrom=_strings(arrays["chrom"]), id=_strings(arrays["id"]),
        pos=np.asarray(arrays["pos"], dtype=np.int64),
        a1=_strings(arrays["a1"]), a2=_strings(arrays["a2"]))
    maf = np.asarray(arrays["maf"], dtype=np.float64)
    long_range_ld = None
    if not include_long_range_ld:
        build = provenance["genome_build"]
        in_region = allelign.in_long_range_ld(variants.chrom, variants.pos, genome_build=build,
                                              include_apoe=include_apoe)
        variants, scores, maf = variants.subset(~in_region), scores[~in_region], maf[~in_region]
        long_range_ld = {"genome_build": build, "include_apoe": include_apoe,
                         "n_dropped": int(in_region.sum())}
    return LDReference(
        name=name, variants=variants, scores=scores, maf=maf,
        m_snps=int(provenance["m_snps"]), m_snps_5_50=int(provenance["m_snps_5_50"]),
        ancestry=provenance["ancestry"], genome_build=provenance["genome_build"],
        provenance=provenance, long_range_ld=long_range_ld)
