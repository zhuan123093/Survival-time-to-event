# ACTG175 Data Sources

The `ACTG175.txt` file comes from the R package `speff2trial` and contains 2,139 participants and 27 variables. The accompanying `ACTG175.Rd` file documents the variables and endpoint. Both files come from commit `9fa954e262daa0dd8d078d8a843537467f7603df` in the CRAN mirror.

- [ACTG175 data file](https://github.com/cran/speff2trial/blob/9fa954e262daa0dd8d078d8a843537467f7603df/data/ACTG175.txt)
- [ACTG175 variable documentation](https://github.com/cran/speff2trial/blob/9fa954e262daa0dd8d078d8a843537467f7603df/man/ACTG175.Rd)

The `sources.json` file records the download URLs, file sizes, and SHA-256 checksums. The analysis script selects the two monotherapy arms in Z3 from the full dataset and writes the selected source rows and participant IDs to [`sample_index.csv`](../results/sample_index.csv).

Hammer et al. (1996) reported the trial in the *New England Journal of Medicine*, 335(15), 1081–1090. [Trial report](https://doi.org/10.1056/NEJM199610103351501)
