# ACTG175 Z3 Application

The analysis compares didanosine (`A = 1`) with zidovudine (`A = 0`) among participants who received more than 52 weeks of prior antiretroviral treatment. The selection rule is `strat == 3 & arms %in% c(0, 3)`. There are 434 participants and 148 events, with 56 events among 221 participants in the didanosine group and 92 among 213 in the zidovudine group.

The endpoint is time to the first occurrence of a decline of at least 50% in CD4 count, an AIDS-defining event, or death. The observed time is stored in `days` and the event indicator in `cens`. The data source and documentation for the 27 variables are in the [data README](data/README.md).

## Adjustment sets

The three nested sets use the same 434 participants. The dimension `p` excludes the intercept.

| `p` | Covariates |
| --- | --- |
| 5 | `cd40`, `preanti`, `cd80`, `age`, `wtkg` |
| 8 | Adds `karnof`, `symptom`, `gender` |
| 14 | Adds `hemo`, `homo`, `drugs`, `race`, `oprior`, `z30` |

Continuous covariates enter as linear terms. Binary covariates use their original 0/1 coding before scaling. All covariates are centered and scaled across the analysis sample.

## Reproduce the application

Run the following command from the repository root with the R packages `survival` and `MASS` installed.

~~~sh
Rscript --vanilla application/code/run.R
~~~

The script uses `code/data.R`, `code/methods.R`, and `code/chiang_loo.R` to compute `lr`, `ols`, `dof`, `loo`, `cloo`, and `ploo`. The `ols` test is the Ye et al. covariate-adjusted log-rank test. The `loo` test uses a cross-fitted numerator based on within-arm regressions with intercepts and the Chiang mHC3 variance estimator. Events tied at the same time are pooled. The `lr` test uses the tie-adjusted variance from `survival::survdiff`, and the covariate adjustments use the observed treatment proportion, 221/434.

Results are written to `application/results/`. Use `--output=path` to specify a different output directory.

| File | Contents |
| --- | --- |
| `test_statistics.csv` | Numerators, variances, standard deviations, and test statistics for the 18 method-by-dimension settings, plus components I and II for PLOO and CLOO |
| `table.csv`, `table.tex` | The application table values in CSV and LaTeX. The LaTeX table uses `booktabs`. |
| `sample_index.csv`, `sample_summary.csv` | Selected source rows and participant IDs, along with participant and event counts by arm |
| `design_diagnostics.csv` | Covariate sets, within-arm design matrix ranks, and maximum leverage |
| `formula_checks.csv` | Numerical checks of the derived-outcome score identity and the dense PLOO and CLOO formulas |

The table reports the numerator `sqrt(n) U`, its estimated standard deviation `sigma`, and the test statistic `sqrt(n) U / sigma`.

The trial was reported by [Hammer et al. (1996)](https://doi.org/10.1056/NEJM199610103351501). The BibTeX entry is in [`reference.bib`](reference.bib).
