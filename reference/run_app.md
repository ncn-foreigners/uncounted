# Launch the uncounted Shiny Application

An interactive Shiny application for estimating unauthorised population
sizes using the `uncounted` package. The app provides a multi-tab
interface for data upload, model fitting, population size estimation,
model comparison, and diagnostics.

## Usage

``` r
run_app(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Details

The application has five tabs:

- Data:

  Upload CSV data or use the built-in `irregular_migration` dataset.
  Select which columns map to observed counts (m), auxiliary counts (n),
  reference population (N), and country identifier.

- Model Fitting:

  Specify and fit multiple models with different methods (Poisson, NB,
  OLS, NLS) and covariate structures. View coefficient summaries.

- Population Size:

  Compute population size estimates stratified by user-specified
  grouping variables. Download results as CSV.

- Model Comparison:

  Compare fitted models using AIC, BIC, and pseudo-R\\^2\\ measures. Run
  likelihood ratio tests.

- Diagnostics:

  Residual plots, rootograms, gamma profiles, and leave-one-out
  sensitivity analysis.

## Examples

``` r
if (interactive()) {
  run_app()
}
```
