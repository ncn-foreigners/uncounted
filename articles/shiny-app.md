# Interactive Shiny Application

## Overview

The `uncounted` package includes an interactive Shiny application that
provides a point-and-click interface for estimating unauthorised
population sizes. The app wraps the package’s core functions in a
multi-tab interface suitable for exploring data, fitting models, and
assessing model quality.

## Installation

The app requires three additional packages:

``` r
install.packages(c("shiny", "DT", "bslib"))
```

## Launching the App

``` r
library(uncounted)
run_app()
```

This opens the application in your default browser.

## Tabs

### 1. Data

Upload a CSV file or use the built-in `irregular_migration` dataset
(1,382 observations of border apprehensions across 145 countries,
2019–2024).

Map your columns to the model variables:

- **Observed (m)**: count of detected unauthorised migrants
- **Auxiliary (n)**: auxiliary detection count (e.g., police
  identifications)
- **Reference pop (N)**: size of the lawful reference population
- **Countries**: country identifier for clustering

### 2. Model Fitting

Specify and fit multiple models. For each model, choose:

- **Method**: Poisson (recommended), Negative Binomial, OLS, or NLS
- **cov_alpha**: formula for community anchor covariates (e.g.,
  `~ year * ukr + sex`)
- **cov_beta**: formula for detection elasticity covariates (e.g.,
  `~ year`)
- **Gamma**: estimate jointly, fix at a value, or exclude
- **Constrained**: enforce alpha in (0,1) and beta \> 0

Click “Fit Model” to estimate. The model is stored and can be compared
with others. The summary panel shows coefficients with robust standard
errors.

### 3. Population Size

Select a fitted model and specify a grouping formula (e.g., `~ year` or
`~ year + sex`). The app computes population size estimates with
delta-method confidence intervals and optional bias correction.

Results can be downloaded as CSV.

### 4. Model Comparison

When two or more models are fitted, this tab shows a comparison table
with:

- AIC, BIC, log-likelihood
- Deviance, Pearson chi-squared
- Pseudo R-squared measures (correlation, explained deviance,
  Cameron-Windmeijer)

A likelihood ratio test can be run for any pair of nested models.

### 5. Diagnostics

Select a fitted model and choose a diagnostic:

- **Residual plots**: 4-panel Zhang (2008) diagnostics
- **Rootogram**: count distribution fit (hanging, suspended, or
  standing)
- **Gamma profile**: sensitivity of population size to gamma
- **LOO (by country)**: leave-one-out influence analysis
- **Exploratory**: log-log diagnostic scatter plots

## Example Workflow

``` r
library(uncounted)
run_app()
```

1.  Select “Built-in (irregular_migration)” and verify variable mapping
2.  In Model Fitting, set method = “poisson”, cov_alpha = “~ year +
    sex”, cov_beta = “~ year”, name = “Poisson S3”. Click Fit Model.
3.  Change method to “nb”, name = “NB S3”. Click Fit Model.
4.  Change cov_alpha to “~ year \* ukr + sex”, method = “poisson”, name
    = “Poisson S8”. Click Fit Model.
5.  Go to Model Comparison to see AIC differences.
6.  Go to Population Size, select “Poisson S8”, type “~ year” in the
    stratify field, click Compute.
7.  Go to Diagnostics, select “Poisson S8”, and explore residuals,
    rootograms, and gamma profiles.
