# PLDSC methods report

[Methods report (PDF)](methods.pdf) describes univariate and bivariate LDSC,
QC, computation, interpretation and the connection to the multi-ancestry
research. The [LaTeX source](methods.tex) is maintained with the PDF. For
runnable examples and file formats, see the [LDSC guide](../docs/LDSC.md).

From the PLDSC root, build with Tectonic:

```bash
tectonic --outdir report report/methods.tex
```

The report was reviewed against version 0.1.0.dev3. Version 0.1.0.dev4 leaves
the regression implementation unchanged.
Tests and parity checks do not establish real-data calibration or a measured
speedup.
