"""Command-line LDSC with allelign-backed summary-statistics processing."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
import json
import math
import os
from pathlib import Path
import sys

import allelign
from allelign.sumstats import COLUMN_ALIASES
import numpy as np

from . import __version__
from .ld_reference import available_ld_references, load_ld_reference
from .sumstats import ldsc_h2_sumstats, ldsc_rg_sumstats


def _column(value):
    field, sep, name = value.partition("=")
    field, name = field.strip(), name.strip()
    if not sep or field not in COLUMN_ALIASES or not name:
        raise argparse.ArgumentTypeError(
            "column must be FIELD=COLUMN; FIELD is one of " + ", ".join(COLUMN_ALIASES))
    return field, name


def _read_options(columns, n_eff):
    options = {}
    for field, name in columns:
        if field in options:
            raise ValueError(f"column {field!r} was specified more than once")
        options[field] = name
    if "beta" in options and "or" in options:
        raise ValueError("map a single effect source: beta or odds ratio")
    counts = {"n_cases", "n_controls"} & options.keys()
    if counts and len(counts) != 2:
        raise ValueError("map both n_cases and n_controls columns")
    if (n_eff is not None and ({"n_eff", "n_cases", "n_controls"} & options.keys())
            or "n_eff" in options and counts):
        raise ValueError("choose one sample-size source per trait")
    if n_eff is not None:
        if not math.isfinite(n_eff) or n_eff <= 0:
            raise ValueError("constant effective N must be positive and finite")
        options["n_eff"] = n_eff
    return options


def _parser():
    parser = argparse.ArgumentParser(prog="pldsc", description="LD Score regression (h2, rg) using allelign for GWAS loading and QC.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="analysis", required=True)
    score_parser = sub.add_parser("prepare-scores", help="align LDSC text scores to a reference and write .npy")
    score_parser.add_argument("--reference", required=True, help="sorted .bim or .pvar variant table")
    score_parser.add_argument("--score-files", nargs="+", required=True, help="LDSC text files, optionally gzip compressed")
    score_parser.add_argument("--out", required=True, help="new .npy output file")
    score_parser.add_argument("--id-column", default="SNP", help="variant ID column (default: SNP)")
    score_parser.add_argument("--chrom-column", default="CHR", help="chromosome column (default: CHR)")
    score_parser.add_argument("--pos-column", default="BP", help="position column (default: BP)")
    score_parser.add_argument("--score-column", default="L2", help="LD score column (default: L2)")
    ref_parser = sub.add_parser("ld-reference", help="describe a bundled LD reference")
    ref_parser.add_argument("--name", default="hm3-eur",
                            help="bundled reference name (default: hm3-eur; one of "
                                 + ", ".join(available_ld_references()) + ")")
    for analysis, blocks in (("h2", 100), ("rg", 200)):
        p = sub.add_parser(analysis, help="SNP heritability" if analysis == "h2" else "genetic correlation")
        if analysis == "h2":
            p.add_argument("--sumstats", required=True, help="GWAS file, optionally gzip compressed")
            p.add_argument("--column", action="append", type=_column, default=[], metavar="FIELD=COLUMN")
            p.add_argument("--n-eff", type=float, help="constant effective N; overrides all file rows")
        else:
            for i in (1, 2):
                p.add_argument(f"--sumstats{i}", required=True, help=f"GWAS file for trait {i}")
                p.add_argument(f"--column{i}", action="append", type=_column, default=[], metavar="FIELD=COLUMN")
                p.add_argument(f"--n-eff{i}", type=float, help=f"constant effective N for trait {i}")
            p.add_argument("--ncores", type=int, default=1, help="parallel jackknife workers (default: 1)")
        p.add_argument("--reference", help="sorted .bim or .pvar variant table")
        p.add_argument("--ld-scores", help=".npy vector in reference row order")
        p.add_argument("--m-snps", type=int,
                       help="number of SNPs summed into each LD score; with --ld-reference, "
                            "overrides its count (e.g. its MAF 5-50%% count)")
        p.add_argument("--ld-reference", metavar="NAME",
                       help="bundled LD reference (one of " + ", ".join(available_ld_references())
                            + ") instead of --reference/--ld-scores")
        p.add_argument("--n-blocks", type=int, help=f"jackknife blocks (default: {blocks}; 200 with --ldsc-compat)")
        p.add_argument("--n-iter", type=int, default=2, help="weight updates (default: 2)")
        p.add_argument("--constrain-intercept", type=float, help="fixed intercept (cross-trait for rg)")
        p.add_argument("--min-n-ratio", type=float, default=0.7, help="minimum N relative to median (default: 0.7)")
        p.add_argument("--min-maf", type=float, default=0.01, help="minimum available MAF (default: 0.01)")
        p.add_argument("--min-info", type=float, default=0.7, help="minimum available INFO (default: 0.7)")
        p.add_argument("--allow-ambiguous", action="store_true", help="retain palindromic alleles when orientation is known")
        p.add_argument("--no-chi2-cap", action="store_true", help="disable the LDSC high-chi-square filter")
        p.add_argument("--ldsc-compat", action="store_true",
                       help="use the reference ldsc estimator (two-step at chi2 30, fixed-weight jackknife, "
                            "200 blocks, no chi-square cap); requires an explicit --m-snps, normally the "
                            "MAF 5-50%% count (5961159 for hm3-eur)")
        p.add_argument("--w-ld-scores",
                       help="with --ldsc-compat and --reference: .npy over-counting scores (ldsc --w-ld) in "
                            "reference row order (default: the regression scores)")
        p.add_argument("--include-long-range-ld-grch37", action="store_true",
                       help="keep reference variants in allelign's GRCh37 long-range LD regions, "
                            "including the MHC (left out by default)")
        p.add_argument("--include-apoe", action="store_true",
                       help="also leave out the APOE region (GRCh37)")
        p.add_argument("--min-snps", type=int, default=MIN_SNPS,
                       help=f"refuse to report a fit on fewer regression SNPs (default: {MIN_SNPS:,}, the "
                            "count below which ldsc warns); see --allow-low-coverage")
        p.add_argument("--allow-low-coverage", action="store_true",
                       help="report the fit even when fewer than --min-snps SNPs remain")
        p.add_argument("--sumstats-build", metavar="BUILD",
                       help="genome build of the GWAS coordinates (e.g. GRCh37, GRCh38); checked against the "
                            "reference build and refused on a mismatch")
        p.add_argument("--reference-build", metavar="BUILD",
                       help="genome build of an explicit --reference (the bundled reference is GRCh37); "
                            "the long-range LD filter needs GRCh37")
        p.add_argument("--out", help="write JSON to a new file (default: stdout; existing files are not overwritten)")
    return parser


def _prepare_scores(args):
    if not args.out.endswith(".npy"):
        raise ValueError("--out must end in .npy")
    variants = allelign.read_variants(args.reference)
    ids = {}
    for i, rsid in enumerate(variants.id):
        rsid = str(rsid).strip()
        if not rsid or rsid == "." or rsid in ids:
            raise ValueError("reference variant IDs must be present and unique")
        ids[rsid] = i
    if not ids:
        raise ValueError("reference contains no variants")
    reference_chrom = {str(c): allelign.normalize_chrom(str(c)) for c in np.unique(variants.chrom)}
    file_chrom = {}
    scores = np.empty(len(variants), dtype=float)
    found = np.zeros(len(variants), dtype=bool)
    ignored = 0
    for path in args.score_files:
        opener = gzip.open if path.endswith((".gz", ".bgz", ".bgzip")) else open
        with opener(path, "rt", encoding="utf-8") as stream:
            header = stream.readline().split()
            names = [name.lstrip("#").upper() for name in header]
            if len(set(names)) != len(names):
                raise ValueError(f"{path}: duplicate score-file column names")
            wanted = [args.id_column, args.chrom_column, args.pos_column, args.score_column]
            try:
                id_col, chr_col, pos_col, score_col = [names.index(name.lstrip("#").upper()) for name in wanted]
            except ValueError:
                raise ValueError(f"{path}: score file must have the requested ID, chromosome, position and score columns") from None
            for line_no, line in enumerate(stream, 2):
                fields = line.split()
                if not fields:
                    continue
                if len(fields) != len(names):
                    raise ValueError(f"{path}:{line_no}: expected {len(names)} fields, got {len(fields)}")
                rsid = fields[id_col]
                i = ids.get(rsid)
                if i is None:
                    ignored += 1
                    continue
                if found[i]:
                    raise ValueError(f"duplicate LD score for reference variant {rsid!r}")
                raw_chrom = fields[chr_col]
                if raw_chrom not in file_chrom:
                    file_chrom[raw_chrom] = allelign.normalize_chrom(raw_chrom)
                try:
                    pos, score = int(fields[pos_col]), float(fields[score_col])
                except ValueError:
                    raise ValueError(f"{path}:{line_no}: invalid position or LD score") from None
                if (file_chrom[raw_chrom] != reference_chrom[str(variants.chrom[i])]
                        or pos != variants.pos[i]):
                    raise ValueError(f"reference and LD-score coordinates disagree for {rsid!r}")
                if not math.isfinite(score) or score <= 0:
                    raise ValueError(f"LD score for {rsid!r} must be positive and finite")
                scores[i], found[i] = score, True
    if not found.all():
        missing = np.flatnonzero(~found)
        raise ValueError(f"LD scores missing for {len(missing)} reference variants; first: {variants.id[missing[0]]}")
    with open(args.out, "xb") as stream:
        np.save(stream, scores, allow_pickle=False)
    return {"reference": args.reference, "score_files": args.score_files, "score_column": args.score_column,
            "output": args.out, "n_matched": len(scores), "n_ignored_score_rows": ignored}


def _json_safe(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


MIN_SNPS = 200_000


def _check_out(path):
    """Fail before fitting if --out cannot be created."""
    if path is None:
        return
    target = Path(path)
    if target.exists():
        raise ValueError(f"--out {path} already exists; existing files are not overwritten")
    parent = target.parent if str(target.parent) else Path(".")
    if not parent.is_dir():
        raise ValueError(f"--out directory {parent} does not exist")
    if not os.access(parent, os.W_OK):
        raise ValueError(f"--out directory {parent} is not writable")


def _m_convention(m_snps, bundled):
    if bundled is None:
        return "user-supplied (--m-snps)"
    if m_snps == bundled.m_snps_5_50:
        return "M_5_50: MAF 5-50% variants (the published ldsc convention)"
    if m_snps == bundled.m_snps:
        return (f"all {bundled.m_snps:,} variants summed into the scores; h2 is "
                f"{bundled.m_snps / bundled.m_snps_5_50:.3f} times the M_5_50 value at the same slope")
    return "user-supplied (--m-snps)"


def _run(args):
    _check_out(args.out)
    if args.min_snps < 1:
        raise ValueError("--min-snps must be positive")
    if args.include_apoe and args.include_long_range_ld_grch37:
        raise ValueError("--include-apoe has no effect with --include-long-range-ld-grch37")
    if args.analysis == "h2":
        paths = [args.sumstats]
        options = [_read_options(args.column, args.n_eff)]
    else:
        paths = [args.sumstats1, args.sumstats2]
        options = [_read_options(args.column1, args.n_eff1), _read_options(args.column2, args.n_eff2)]
    if args.ldsc_compat and args.m_snps is None:
        raise ValueError("--ldsc-compat requires an explicit --m-snps; the reference ldsc results use "
                         "the MAF 5-50% count (5961159 for hm3-eur)")
    if args.w_ld_scores is not None and not (args.ldsc_compat and args.ld_reference is None):
        raise ValueError("--w-ld-scores needs --ldsc-compat and the explicit --reference/--ld-scores inputs")
    if args.ldsc_compat and args.no_chi2_cap:
        raise ValueError("--ldsc-compat already leaves out the chi-square cap")
    if args.n_blocks is None:
        args.n_blocks = 200 if (args.ldsc_compat or args.analysis == "rg") else 100
    w_ld = None
    bundled = None
    explicit = (args.reference, args.ld_scores, args.m_snps)
    if args.ld_reference is not None:
        if args.reference is not None or args.ld_scores is not None:
            raise ValueError("give either --ld-reference or the explicit "
                             "--reference/--ld-scores inputs, not both")
        bundled = load_ld_reference(args.ld_reference,
                                    include_long_range_ld=args.include_long_range_ld_grch37,
                                    include_apoe=args.include_apoe)
        if args.reference_build is not None:
            raise ValueError("--reference-build applies to an explicit --reference; the bundle records its own")
        variants_input, scores, long_range_ld = bundled.variants, bundled.scores, bundled.long_range_ld
        m_snps = bundled.m_snps if args.m_snps is None else args.m_snps
        reference_build = bundled.genome_build
        reference_record = {"name": bundled.name, "ancestry": bundled.ancestry,
                            "genome_build": bundled.genome_build,
                            "bundle_sha256": bundled.provenance["bundle_sha256"],
                            "m_snps": "bundled" if args.m_snps is None else "--m-snps"}
    else:
        if any(value is None for value in explicit):
            raise ValueError("give --reference, --ld-scores and --m-snps, "
                             "or name a bundled --ld-reference")
        variants_input, m_snps, reference_record = args.reference, args.m_snps, None
        reference_build = (None if args.reference_build is None
                           else allelign.normalize_build(args.reference_build, name="--reference-build"))
        scores = np.load(args.ld_scores, allow_pickle=False, mmap_mode="r")
        if not isinstance(scores, np.ndarray):
            scores.close()
            raise ValueError("--ld-scores must be a single .npy vector, not an archive")
        if args.w_ld_scores is not None:
            w_ld = np.load(args.w_ld_scores, allow_pickle=False)
            if not isinstance(w_ld, np.ndarray) or w_ld.shape != scores.shape:
                raise ValueError("--w-ld-scores must be a .npy vector with one value per reference variant")
        long_range_ld = None
        if not args.include_long_range_ld_grch37 and reference_build not in (None, "GRCh37"):
            raise ValueError(f"the long-range LD regions are GRCh37 coordinates but --reference-build is "
                             f"{reference_build}; pass --include-long-range-ld-grch37 to skip the filter")
        if not args.include_long_range_ld_grch37:
            # Leave out reference rows, so GWAS files matched by ID need no positions.
            variants_input = allelign.read_variants(args.reference)
            if len(variants_input) != scores.size:
                raise ValueError("ld_scores must have one value per reference variant")
            in_region = allelign.in_long_range_ld(variants_input.chrom, variants_input.pos,
                                                  genome_build="GRCh37", include_apoe=args.include_apoe)
            variants_input = variants_input.subset(~in_region)
            scores = np.asarray(scores)[~in_region]
            if w_ld is not None:
                w_ld = np.asarray(w_ld)[~in_region]
            long_range_ld = {"genome_build": "GRCh37", "include_apoe": args.include_apoe,
                             "n_dropped": int(in_region.sum()),
                             "reference_build": reference_build or "not declared; GRCh37 assumed"}
            if reference_build is None:
                print("pldsc: warning: --reference-build not given; the long-range LD filter assumes "
                      "GRCh37 reference coordinates", file=sys.stderr)
    declared = (None if args.sumstats_build is None
                else allelign.normalize_build(args.sumstats_build, name="--sumstats-build"))
    build_check = allelign.check_build(declared, reference_build)
    qc_options = {"min_n_ratio": args.min_n_ratio, "min_maf": args.min_maf, "min_info": args.min_info}
    common = dict(m_snps=m_snps, n_blocks=args.n_blocks, n_iter=args.n_iter,
                  constrain_intercept=args.constrain_intercept, drop_ambiguous=not args.allow_ambiguous,
                  chi2_cap=not (args.no_chi2_cap or args.ldsc_compat), ldsc_compat=args.ldsc_compat)
    if args.analysis == "h2":
        result = ldsc_h2_sumstats(paths[0], variants_input, scores, read_options=options[0],
                                 qc_options=qc_options, w_ld=w_ld, **common)
    else:
        result = ldsc_rg_sumstats(paths[0], paths[1], variants_input, scores, read_options1=options[0],
                                  read_options2=options[1], qc_options1=qc_options, qc_options2=qc_options,
                                  ncores=args.ncores, w_ld=w_ld, **common)
    n_retained, n_reference = len(result.reference_index), int(np.shape(scores)[0])
    if n_retained < args.min_snps and not args.allow_low_coverage:
        raise ValueError(
            f"only {n_retained:,} of {n_reference:,} reference variants ({n_retained / n_reference:.2%}) "
            f"remain for the regression, fewer than --min-snps {args.min_snps:,}. An LDSC fit on so few "
            "variants is unreliable. Check identifiers, alleles and genome build (the QC log's "
            "reference_coverage names the likely cause), or pass --allow-low-coverage")
    return _json_safe({
        "pldsc_version": __version__, "allelign_version": allelign.__version__, "analysis": args.analysis,
        "inputs": {"sumstats": paths, "reference": args.reference, "ld_scores": args.ld_scores,
                   "ld_reference": reference_record},
        "options": {**common, "read_options": options, "qc_options": qc_options,
                    "estimator": "ldsc-compatible (two-step, fixed-weight jackknife)" if args.ldsc_compat
                                 else "pldsc one-step (weights refit in every jackknife deletion)",
                    "w_ld": (args.w_ld_scores or "regression LD scores") if args.ldsc_compat else None,
                    "ncores": args.ncores if args.analysis == "rg" else None},
        "long_range_ld": long_range_ld,
        "estimate": asdict(result.estimate), "qc": result.qc, "alignment": result.alignment,
        "n_common": result.n_common, "n_chi2_dropped": result.n_chi2_dropped,
        "n_retained": n_retained,
        "coverage": {"n_regression_snps": n_retained, "n_reference_variants": n_reference,
                     "fraction": n_retained / n_reference, "min_snps": args.min_snps,
                     "low_coverage_allowed": bool(args.allow_low_coverage and n_retained < args.min_snps)},
        "m_snps": {"value": m_snps, "convention": _m_convention(m_snps, bundled)},
        "build_check": build_check,
    })


def _ld_reference(args):
    reference = load_ld_reference(args.name)
    record = {key: value for key, value in reference.provenance.items()
              if key != "dropped_nonpositive_score_ids"}
    return {"available": available_ld_references(), "reference": record,
            "long_range_ld": reference.long_range_ld}


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    if args.analysis in ("h2", "rg") and args.out:
        try:
            _check_out(args.out)
        except ValueError as exc:
            print(f"pldsc: error: {exc}", file=sys.stderr)
            return 2
    try:
        if args.analysis == "prepare-scores":
            result = _prepare_scores(args)
        elif args.analysis == "ld-reference":
            result = _ld_reference(args)
        else:
            result = _run(args)
        output = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.analysis in ("h2", "rg") and args.out:
            with open(args.out, "x", encoding="utf-8") as stream:
                stream.write(output)
        else:
            sys.stdout.write(output)
    except (OSError, ValueError) as exc:
        print(f"pldsc: error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
