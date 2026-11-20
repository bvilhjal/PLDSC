# Which LD score matches a GWAS?

Read the [PDF](methods.pdf) or the [Markdown edition](methods.md).
Sections 1–6 cover ordinary LDSC, multi-ancestry meta-analysis, known-truth
simulations, 1000 Genomes panels and the GBMI asthma application. Sections 7–9
cover mixed models, relatives, enrichment and family mechanisms. Section 6.1
tests prediction on held-out asthma chromosomes. Section 10 measures power
against known truth and specifies an independent-cohort test. Section 11 asks
what Cochran's Q adds, in simulation and in asthma. Section 12 links the experiments and inputs.
The [reproduction record](reproduction.txt) contains the reviewed runs'
numerical output and the SHA-256 digests of the empirical inputs.

## Rebuild the report

Edit [methods.md](methods.md). From the `report` directory, generate the
[standalone LaTeX source](methods.tex) with Pandoc 3.8 and build the PDF
with Tectonic. The second command also works on its own with the distributed
source and the four figure images.

```bash
pandoc --from=markdown-implicit_figures+gfm_auto_identifiers --to=latex --template=layout.tex \
  --wrap=preserve methods.md -o methods.tex
tectonic methods.tex
```

## Reproduce the figures

From the repository root, install PLDSC with its `plot` extra and run
`python report/reference_figure.py`. The four figures use distributed inputs.

**Table 1. Inputs distributed with each figure.**

| Figure | Inputs |
|:--|:--|
| [1. Reference composition](reference_composition.png) | The [chr22 score bundle](../examples/mini) contains the two scores for every plotted SNP. |
| [2. Asthma prediction](asthma_prediction.png) | [22 chromosome cross-products and held-out errors](asthma_score_shape.json) reproduce both fits and their uncertainty. |
| [3. False positives and power](score_shape_power.png) | [Simulation counts](score_shape_simulation.json) cover all 18 cells, with 1,000 replicates per cell. |
| [4. Cochran's Q](heterogeneity.png) | [Simulation summaries](q_heterogeneity_simulation.json) for every cell and the [asthma estimates](asthma_heterogeneity.json) with their standard errors. |

To regenerate the simulation counts, run `python experiments/score_shape_sim.py --out report/score_shape_simulation.json`.
To regenerate the chromosome summaries, first build the genome-wide scores as described in the
[experiment index](../experiments/README.md), then run
`python experiments/gbmi_asthma_casualty.py --shape-out report/asthma_score_shape.json`.
Use `--score-dir` for a separate directory of source- and map-bound chromosome scores.
For Figure 4, run `python experiments/q_heterogeneity_sim.py --out report/q_heterogeneity_simulation.json`
and `python experiments/gbmi_asthma_heterogeneity.py --out report/asthma_heterogeneity.json`.
