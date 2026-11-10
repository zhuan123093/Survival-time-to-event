# Covariate-Adjusted Log-Rank Tests

This repository contains the code and results for a simulation study of covariate-adjusted log-rank tests and a reanalysis of the ACTG175 Z3 subgroup. The simulations compare seven tests in completely randomized experiments (CRE), Bernoulli randomized experiments (BRE), and stratified covariate-adaptive randomization (CAR). The ACTG175 application compares six tests and reproduces the application table in the paper.

## Repository layout

| Path | Contents |
| --- | --- |
| `INDEX.csv` | Links each of the nine figures to its plotting data and simulation results |
| `MANIFEST.csv` | File sizes and SHA-256 checksums |
| `code/settings.csv` | The 24 simulation settings, including sample sizes, designs, hypotheses, allocation ratios, adjustment dimensions, and replicate offsets |
| `code/data.R` | Event censoring, risk-set derived outcomes, and the unadjusted log-rank test |
| `code/methods.R` | OLS, degrees-of-freedom correction, CLOO, PLOO, oracle, and shared matrix calculations |
| `code/chiang_loo.R` | The cross-fitted numerator with intercepts and the Chiang mHC3 variance estimator |
| `code/experiment.R`, `code/run.R` | Data generation, random seeds, and execution or replay of simulation replicates |
| `code/summarize.R`, `code/plot.py` | Simulation summaries, figures, and plotting data |
| `results/raw/` | 24 RDS files with replicate-level results, with 10,000 replicates per setting |
| `results/summary.csv` | 700 rows of simulation summaries |
| `results/variance_components.csv`, `results/pooled_variance_components.csv` | Variance components for LOO and for CLOO, PLOO, and oracle |
| `results/invalid_trials.csv` | Invalid replicates and their test numerators and variances |
| `results/figures/`, `results/plot_data/` | Nine figures in PDF and PNG and the CSV data used to plot them |
| `results/code_check.csv`, `results/summary_check.csv`, `results/plot_check.csv` | Comparisons of replayed replicates, regenerated summaries, and plotted values with the saved results |
| `application/` | ACTG175 data and variable documentation, analysis code, results, and the original trial citation |

## Simulation design

The seven tests are `lr`, `ols`, `dof`, `loo`, `cloo`, `ploo`, and `oracle`. The standard LOO test uses a cross-fitted numerator based on within-arm regressions with intercepts and the Chiang mHC3 variance estimator. The adjustment dimension `p` excludes the intercept.

In the stratified design, the first ten covariate coefficients are −0.4, 0.5, and 0.2 in the three strata, respectively. The tests are computed within strata and retain the between-stratum variance term. Sample Gram matrices and regression cross-moments use `n` for normalization. PLOO uses the projection weights `H`. The oracle uses the known population covariance with weights `H_ij^oracle = X_i' Sigma^(-1) X_j / n`. The CLOO, PLOO, and oracle variance calculations separate components I and II. The oracle uses the off-diagonal factor `(1 + H_ij^oracle) H_ij^oracle` in component I.

Under the alternative, the treatment hazard ratio is `exp(-0.5)`. Follow-up ends at the 70% event cutoff. CRE and BRE use replicate IDs 1–10,000 and CAR uses 20,001–30,000. The complete settings are in [`code/settings.csv`](code/settings.csv).

## Simulation results

In `results/raw/`, `n100_CRE_1to1_H0.rds` contains the Gaussian CRE results for `n = 100` with 1:1 allocation under the null. The `n200`, `n300`, `n500`, and `n1000` files cover the sample-size experiment. The `gamma_n100_CRE_*` files cover the Gamma covariate experiment, and the `sign_low_*` files cover the stratified design.

Use `readRDS()` to load a replicate-level results file. The `raw[replicate, p, field]` array stores the test numerators, variances, and selected standardized statistics. The `parts[replicate, p, component]` array stores the LOO variance components. The seed rules in `code/experiment.R` reproduce the covariates and survival data.

The `pooled_components[replicate, p, method, component]` array stores I, II, and the between-stratum term for CLOO, PLOO, and oracle. Their sum is the variance estimate. The between-stratum term is zero for CRE and BRE. The component averages are in [`results/pooled_variance_components.csv`](results/pooled_variance_components.csv).

Each row of `results/summary.csv` represents one simulation setting, adjustment dimension, and method. The columns `rate`, `valid`, and `invalid` give the rejection proportion and the counts of valid and invalid replicates. In the higher-dimensional BRE settings with 2:1 allocation, 1,536 mHC3 variances are undefined. Rejection rates use each method's valid replicates. The `value` column in each plotting CSV gives the rejection rate shown in the figure.

| Figure stem in `results/figures/` | Content |
| --- | --- |
| `introduction_CRE` | Seven-method CRE comparison in the introduction |
| `CRE_main` | CRE simulations in the main text |
| `CAR_main` | Stratified randomization simulations in the main text |
| `BRE_appendix` | BRE simulations in the appendix |
| `CRE_with_ols` | CRE simulations including OLS |
| `BRE_with_ols` | BRE simulations including OLS |
| `CAR_with_ols` | Stratified randomization simulations including OLS |
| `sample_size_appendix` | Sample-size experiment at `p = 20` |
| `gamma_CRE_appendix` | CRE simulations with Gamma covariates |

## Reproduce the results

Run the following commands from the repository root. The R scripts use `MASS` and `survival`. The plotting script uses `numpy`, `pandas`, and `matplotlib`.

~~~sh
# Regenerate the nine simulation figures.
python code/plot.py

# Check the summary against the saved replicate-level results.
Rscript --vanilla code/summarize.R

# Replay selected saved replicates.
Rscript --vanilla code/run.R --mode=replay

# Run setting 1 from code/settings.csv. The runner skips existing output.
Rscript --vanilla code/run.R --mode=run --job=1 --B=10000

# Rebuild the summary after all 24 simulation settings are available.
Rscript --vanilla code/summarize.R --collect

# Reproduce the ACTG175 Z3 application table.
Rscript --vanilla application/code/run.R
~~~

The ACTG175 data are included in the repository. The subgroup definition, adjustment sets, and output files are described in the [application README](application/README.md).
