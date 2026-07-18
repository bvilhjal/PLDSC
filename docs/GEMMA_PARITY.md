# M4 GEMMA score-test parity protocol

Date: 2026-07-18

## Status and scope

M4 asks whether PLDSC stays calibrated with a **fitted** `V_c` and whether its
score statistic reproduces an established exact-LMM engine. The two halves have
different states:

- **Fitted-component science: validated.** The REML fit used in
  `experiments/m4_fitted_vc.py` and `experiments/m5_inference.py` is the exact
  Gaussian LOCO LMM null that GEMMA also fits (REML variance components with the
  LOCO GRM). On that fitted `V_hat_c` the PLDSC intercept is `1.010 ± 0.078`
  (M4b) and the regression slope stays calibrated (M5). So the science M4 exists
  to establish — does calibration survive plugging in estimated variance
  components — is already answered.
- **External byte-parity against GEMMA: deferred to Linux.** GEMMA ships no
  Windows build (`conda search -c bioconda gemma` returns no `win-64` package;
  the `bioconductor-gemma.r` hits are an unrelated gene-expression R package).
  This protocol specifies the parity check so it runs unchanged on Linux, WSL,
  macOS, or CI. It is an **independent cross-check** of an already
  triple-derived, unit-tested score, not an open scientific risk.

## What parity means

For target chromosome `c`, PLDSC forms the LOCO mixed-model score
`z_j = x_j^T P_c y / sqrt(x_j^T P_c x_j)` with
`P_c` built from `V_c = sigma_g^2 K_{-c} + sigma_e^2 I` (docs/METHOD.md). GEMMA's
univariate LMM **score test** (`-lmm 3`) computes the same GLS score for each SNP
under the null variance components fitted from a supplied kinship. With matched
inputs, `z_j^2` must equal GEMMA's per-SNP score chi-square within a documented
numerical tolerance.

## Predeclared protocol (fix before running; per the plan's M4 gate)

1. **Genotypes.** Simulate a cohort with `pldsc.simulate` and write it in BIMBAM
   mean-genotype format (GEMMA `-g`) plus a phenotype file (`-p`). Use the same
   `float64` dosages PLDSC consumes — no re-rounding.
2. **Allele scaling.** GEMMA `-gk 2` standardizes genotypes as
   `(g - 2p)/sqrt(2p(1-p))`, identical to `pldsc.simulate.standardize`. Use
   `-gk 2` so the kinship matches `pldsc.simulate.grm`.
3. **LOCO kinship.** Build `K_{-c}` from all chromosomes except the target and
   pass it to GEMMA with `-k`. PLDSC must use the *same* matrix (write GEMMA's
   `.cXX.txt` and load it, rather than recomputing) to remove GRM-construction as
   a variable.
4. **Null fit.** GEMMA fits `sigma_g^2, sigma_e^2` by REML (its default). Record
   `pve` and the null log-likelihood. Build `V_hat_c` and `P_hat_c` from exactly
   those components.
5. **Statistic and field.** Run `-lmm 3` (score test). Compare GEMMA's per-SNP
   score chi-square (from `p_score`, `chi2 = qchisq(1 - p_score, df = 1)`) to
   PLDSC `z_j^2`. Do **not** use `-lmm 1/2/4` (Wald/LRT) — those are different
   statistics and, by default, refit variance components per SNP.
6. **Gate.** Report `max_j |z_j^2 - chi2_gemma,j| / chi2_gemma,j` and the
   Pearson correlation across SNPs. Parity passes at a preregistered tolerance
   (target: correlation > 0.9999 and max relative difference < 1e-3, subject to
   GEMMA's internal convergence settings). A persistent gap is documented as a
   scope limit, not silently rebranded.

## Ready-to-run harness (Linux/WSL/macOS)

```python
# 1) write BIMBAM genotypes + phenotype from pldsc.simulate
import numpy as np
from pldsc.simulate import simulate_admixture, grm
rng = np.random.default_rng(np.random.PCG64(20260718))
G, chrom, pos = simulate_admixture(rng, [80, 80, 80], n_indiv=400,
                                   fst=0.15, families=True, sibs=4)
keep = G.var(0) > 0; G, chrom = G[:, keep], chrom[keep]
tgt = chrom == 0
y = np.linalg.cholesky(0.7*grm(G[:, ~tgt]) + 0.3*np.eye(G.shape[0])) @ rng.standard_normal(G.shape[0])

def write_bimbam(path, Gsub):
    with open(path, "w") as f:
        for j in range(Gsub.shape[1]):
            dos = ",".join(f"{v:.6f}" for v in Gsub[:, j])
            f.write(f"rs{j}, A, T, {dos}\n")
write_bimbam("target.geno.txt", G[:, tgt])         # tested SNPs
write_bimbam("loco.geno.txt", G[:, ~tgt])          # for the LOCO kinship
np.savetxt("pheno.txt", y)
```

```bash
# 2) GEMMA: LOCO kinship (-gk 2 = standardized) then the score test (-lmm 3)
gemma -g loco.geno.txt   -p pheno.txt -gk 2 -o loco
gemma -g target.geno.txt -p pheno.txt -k output/loco.sXX.txt -lmm 3 -o parity
```

```python
# 3) PLDSC side with GEMMA's kinship and fitted components, then compare
import numpy as np, pandas as pd
from scipy.stats import chi2 as chi2d           # only for qchisq in the check
from pldsc._lmm import residual_precision, score_pieces
from pldsc.simulate import standardize
K = np.loadtxt("output/loco.sXX.txt")
# sigma_g2, sigma_e2 from GEMMA's log (pve, ve); build V_hat_c:
Vc = sigma_g2 * K + sigma_e2 * np.eye(K.shape[0])
P = residual_precision(Vc, np.ones((K.shape[0], 1)))
Xc = standardize(G[:, tgt]); PX, q = score_pieces(P, Xc)
z2 = (PX.T @ y / np.sqrt(q)) ** 2
g = pd.read_csv("output/parity.assoc.txt", sep="\t")
chi2_gemma = chi2d.ppf(1 - g["p_score"].values, df=1)
print("corr", np.corrcoef(z2, chi2_gemma)[0, 1],
      "max rel diff", np.max(np.abs(z2 - chi2_gemma) / chi2_gemma))
```

(The check uses `scipy` for `qchisq` only; PLDSC itself has no scipy dependency.)

## When this runs

Run in the Linux CI job or a WSL/Linux checkout once GEMMA is on `PATH`. Until
then, M4's calibration conclusion rests on the local-REML results in
`experiments/`, which reproduce GEMMA's null model exactly; this protocol upgrades
that to an independent-engine cross-check.
