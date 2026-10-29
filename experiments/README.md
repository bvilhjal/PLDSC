# Experiments

The scripts here support the multi-ancestry research and the earlier,
concluded estimator investigation. Detailed protocols, results, corrections
and negative findings are retained in [RESULTS.md](RESULTS.md). See
[MULTI_ANCESTRY.md](../docs/MULTI_ANCESTRY.md) for the current model and
its limitations. Run commands from the repository root with PLDSC installed.

## Multi-ancestry studies

- [Monte Carlo validation](multiancestry_meta_sim.py):
  `python experiments/multiancestry_meta_sim.py --n-ref 6000`. The default
  400-person reference-panel run failed its original calibration gate; the
  larger-panel rerun and its separately stated gates are documented in
  MULTI_ANCESTRY section 6.
- [Chromosome 22 panel analysis](multiancestry_1000g_meta.py): requires the
  1000 Genomes scratch data described below. The aggregate bundle reproduces
  selected tables without a download: `python examples/reproduce_chr22_tables.py`.
- [Genome-wide analysis](multiancestry_genome_kappa.py): requires the
  bigsnpr 1000 Genomes genotype files and maps described in the script.
- [GBMI asthma](gbmi_asthma_casualty.py): requires six GBMI summary-statistics
  files and the genome-wide score cache. The script records their sources,
  checksums and analysis rules.
- [Reference-score noise](reference_score_noise.py): a known-truth simulation
  illustrating the limitations of panel plug-in slope comparisons.

For the chromosome 22 genotype-based analysis, fetch the sources using
`python scripts/fetch_1000g.py --maps 22 --chrom 22 --panel --build-cache 22`.
The data and caches live in gitignored scratch directories. See
[CONTRIBUTING.md](../CONTRIBUTING.md) and each script for prerequisites;
the committed aggregate bundle's provenance is in [examples/mini](../examples/mini).

## Other investigations

[RESULTS.md](RESULTS.md) describes the structure decomposition, annotation
and relatedness studies and the concluded mixed-model estimator line. The
scripts remain here so that the recorded configurations can be inspected
and rerun. The [reading guide](../docs/READING_ORDER.md) distinguishes this
historical record from current research.

For regression performance and numerical parity against installed LDpred3 and
bipred checkouts, see [ldsc_provider_parity_bench.py](ldsc_provider_parity_bench.py)
and the [LDSC guide](../docs/LDSC.md#interpretation-and-performance).
