# Changelog

## uncounted 0.5.1

- Added Aniela Czerniawska as package author (original GitHub version).

- New
  [`dfbeta.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/dfbeta.uncounted.md)
  method for leave-one-out influence on coefficients. Wrapper around
  [`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md)
  returning the familiar n x p matrix of coefficient changes.

- New
  [`dfpopsize()`](https://ncn-foreigners.github.io/uncounted/reference/dfpopsize.md)
  generic and
  [`dfpopsize.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/dfpopsize.md)
  method for leave-one-out influence on population size estimates.
  Returns a named vector of per-observation (or per-country) changes in
  total xi.

## uncounted 0.5.0

### Breaking changes

- The `vcov_type` parameter in
  [`estimate_hidden_pop()`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md)
  has been renamed to `vcov`. It now accepts either a character string
  (`"HC0"` through `"HC5"`) or a **function** that takes the fitted
  object and returns a variance-covariance matrix (e.g.,
  [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
  or a custom wrapper around
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)).

- The `countries` parameter no longer triggers cluster-robust standard
  errors. It is now used **only** for grouping in
  [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  and
  [`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md).
  Use the new `cluster` parameter for cluster-robust variance
  estimation.

### New features

- New `cluster` parameter in
  [`estimate_hidden_pop()`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md)
  accepts a one-sided formula (e.g., `~ country_code`) to request
  cluster-robust variance estimation via
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html).

- Full compatibility with the **sandwich** package through new S3
  methods: `bread()`, `estfun()`,
  [`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html), and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). This
  allows direct use of
  [`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html),
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html),
  and
  [`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)
  on `uncounted` objects.

- New `type = "working"` option in
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) returns the
  score residuals used internally for sandwich variance computation.

- New
  [`update.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/update.uncounted.md)
  method enables re-fitting with modified arguments (e.g., different
  weights or vcov). Supports `evaluate = FALSE` to return the
  unevaluated call.

- New `weights.uncounted()` method returns observation weights from the
  fit.

- Compatibility with
  [`fwb::vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.html)
  for fractional weighted bootstrap variance estimation, thanks to
  [`update()`](https://rdrr.io/r/stats/update.html) and
  [`weights()`](https://rdrr.io/r/stats/weights.html) methods.

- The `vcov` parameter accepts a function, enabling custom variance
  estimators:

  ``` r
  # Cluster-robust via function
  fit <- estimate_hidden_pop(...,
    vcov = function(x) sandwich::vcovCL(x, cluster = data$country))

  # Observation-level via string (default)
  fit <- estimate_hidden_pop(..., vcov = "HC3")

  # Cluster-robust via string + cluster formula
  fit <- estimate_hidden_pop(..., vcov = "HC1", cluster = ~ country_code)
  ```

### Bug fixes

- Fixed a bug where `HC0` and `HC3` (and all other HC types) produced
  **identical** standard errors when clustering was active. The old
  implementation ignored `vcov_type` entirely in the cluster branch of
  the internal sandwich computation. Variance estimation is now
  delegated to the **sandwich** package, which correctly handles all HC
  types for both observation-level and cluster-robust cases.

- Fixed a bug in the cluster-robust meat matrix computation where an
  extra $\sqrt{\mu_{i}}$ factor was incorrectly included in the Poisson
  and NB score contributions. This is now resolved by using
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html).

### Internal

- Variance-covariance computation is now performed post-estimation using
  the **sandwich** package rather than the internal
  `.compute_sandwich_vcov()` function. The internal function is retained
  for the estimator-level fallback `vcov_model` (used in bias
  correction).

- Each estimator now returns `model_matrix_full`, `bread_weights`, and
  `score_residuals` to support the `bread()` and `estfun()` S3 methods.

- Added **sandwich** to `Imports` in `DESCRIPTION`.

- New file `R/sandwich_methods.R` containing
  [`bread.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/bread.uncounted.md),
  [`estfun.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/estfun.uncounted.md),
  [`hatvalues.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/hatvalues.uncounted.md),
  and
  [`model.matrix.uncounted()`](https://ncn-foreigners.github.io/uncounted/reference/model.matrix.uncounted.md).

- Expanded test suite to 245 tests covering sandwich methods,
  cluster-robust variance, function-based vcov, and the HC0 != HC3
  regression test.

## uncounted 0.4.0

### New features

- [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  and
  [`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
  gain a `by` parameter for stratified population size estimates (e.g.,
  `by = ~ year`, `by = ~ year + sex`). Previously, stratification was
  only available through covariate-varying alpha.

- New pseudo $R^{2}$ measures in
  [`compare_models()`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md):
  explained deviance ($R_{D}^{2}$) and Cameron-Windmeijer ($R_{CW}^{2}$)
  for Poisson and Negative Binomial models.

- Interactive **Shiny application** via
  [`run_app()`](https://ncn-foreigners.github.io/uncounted/reference/run_app.md)
  with five tabs: Data, Model Fitting, Population Size, Model
  Comparison, and Diagnostics.

### Improvements

- Rootogram: white background bars for hanging style, light-blue bars
  for standing and suspended styles for clearer visual contrast.

- [`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md)
  call arguments are now force-evaluated and stored in the call object,
  fixing “object not found” errors when refitting from within Shiny or
  other non-standard evaluation contexts.

- Shiny app: automatic conversion of year-like columns to factors to
  prevent continuous covariate issues.

### Bug fixes

- Fixed Cameron-Windmeijer reference year: 1997 to 1996 (JBES 14(2),
  209-220).

- Fixed non-ASCII characters in country names (`Cote d'Ivoire`,
  `Sao Tome e Principe`) causing encoding issues on some platforms.

## uncounted 0.3.1

- Renamed internal `sigmoid` to `logit` / `inv_logit` throughout the
  package for consistency with standard statistical terminology.

## uncounted 0.3.0

### New features

- Added `irregular_migration` built-in dataset (Poland 2019-2024) with
  documentation.

- Added package vignettes: `getting-started` and `model-workflow`.

- Full roxygen2 documentation for all exported functions.

## uncounted 0.2.0

### New features

- [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  replaces
  [`xi()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  as the primary interface for population size estimation.
  [`xi()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  is retained as a deprecated alias.

- Bias-corrected population size estimates using second-order Taylor
  expansion (Jensen’s inequality correction).

- Delta-method confidence intervals for total population size.

- [`profile_gamma()`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md)
  for profiling the gamma parameter over a grid.

- [`compare_models()`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md)
  for side-by-side model comparison (AIC, BIC, log-likelihood).

- [`lrtest()`](https://ncn-foreigners.github.io/uncounted/reference/lrtest.md)
  for likelihood-ratio tests between nested models.

- [`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md)
  for leave-one-out diagnostics by observation or by country.

- [`compare_loo()`](https://ncn-foreigners.github.io/uncounted/reference/compare_loo.md)
  for comparing LOO influence across models.

- Constrained estimation (`constrained = TRUE`): logit link for
  $\alpha \in (0,1)$ and log link for $\beta > 0$.

- Fractional weighted bootstrap via
  [`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
  using the **fwb** package (Xu et al. 2020), with cluster support.

- Observation weights via `weights` parameter.

- HC4, HC4m, and HC5 robust standard error types.

- [`rootogram()`](https://ncn-foreigners.github.io/uncounted/reference/rootogram.md)
  for assessing distributional fit (hanging, standing, suspended).

### Improvements

- Default `gamma = "estimate"` (was `NULL`).

- Model-based (non-sandwich) variance stored as `vcov_model` for bias
  correction, separate from HC-robust `vcov`.

- Improved numerical stability for extreme parameter values.

- `~0` formula treated as intercept-only (not empty model).

## uncounted 0.1.0

- Initial release.

- [`estimate_hidden_pop()`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md)
  for fitting the power-law model
  $E\left( m_{i} \right) = N_{i}^{\alpha_{i}}\left( \gamma + n_{i}/N_{i} \right)^{\beta_{i}}$
  using OLS, NLS, Poisson PML, or Negative Binomial MLE.

- Covariate-varying $\alpha$ and $\beta$ via formula interface
  (`cov_alpha`, `cov_beta`).

- Gamma parameter: estimated, fixed, or excluded.

- HC0-HC3 robust standard errors.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) S3 methods.
