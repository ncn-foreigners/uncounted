# uncounted 1.4.0

## New features

* **`frailty_sensitivity()`**: New omitted-frailty identification
  sensitivity analysis for Poisson and NB MLE fits. The method linearizes the
  fitted count model with its IRLS working response, indexes omitted shared
  frailty with weighted partial-`R^2` parameters, and reports how the grouped
  hidden-population estimand `Xi` changes under first-order bias formulas.

* **Grouped `Xi_t` support in omitted-frailty sensitivity**:
  `frailty_sensitivity()` mirrors `popsize(by = ...)` group construction so
  users can study yearly or subgroup hidden-population targets with the same
  aggregation logic as the main population-size output.

## Documentation

* **Theory-rich manual page for omitted-frailty sensitivity**: The new
  `frailty_sensitivity()` help page documents the full model setup
  `mu = xi * rho`, the omitted-frailty relaxation, IRLS linearization,
  weighted partial-`R^2` sensitivity parameters, bias formulas, delta-method
  conversion to `Xi`, robustness values, interpretation, and limitations.

* **Diagnostics and theory vignette updates**: The diagnostics vignette now
  includes a compact omitted-frailty workflow, and the theory vignette
  distinguishes fixed-center dependence bounds, moving-`Xi` dependence
  profiles, bootstrap threshold summaries, and the new IRLS-based
  omitted-frailty diagnostic.

# uncounted 1.3.0

## New features

* **`dependence_bounds()`**: New fixed-center identification envelope for the
  separability assumption underlying `E(m_i | N_i, n_i) = xi_i * rho_i`.
  Reports lower and upper total-population bounds over user-supplied `Gamma`
  values and optionally computes a tipping-point `gamma_star`.

* **`sensitivity_dependence()` deprecation**: Kept as a backward-compatible
  wrapper to `dependence_bounds()` and now emits a deprecation warning.

* **`profile_dependence()`**: New moving-`xi` dependence profile for Poisson
  and NB MLE fits. Refits the model over a grid of fixed dependence offsets
  `delta`, returning `delta`, `kappa = exp(delta)`, the plug-in total
  population-size estimate, and the corresponding log-likelihood.

* **`robustness_dependence()`**: New one-dimensional robustness-value summary
  built on top of the dependence profile. Reports the smallest dependence
  strength needed for the profiled total population-size estimate to cross a
  user-specified threshold or relative target.

* **`exceedance_popsize()`**: New bootstrap threshold helper that computes the
  empirical bootstrap tail area `P*(xi > c)` or `P*(xi < c)` from an
  `uncounted_boot` object, using per-replicate totals when available.

## Documentation

* **Dependence diagnostics workflow**: The README, diagnostics vignette, and
  theory vignette now distinguish fixed-center bounds, moving-`xi` profiles,
  robustness-value summaries, and bootstrap exceedance probabilities. The new
  examples use simulated data where auxiliary counts scale with the reference
  population, matching the real-data structure more closely.

# uncounted 1.1.0

## New features

* **Bounded detection links** (`link_rho`): Count-model fits now support
  `link_rho = "power"` (paper baseline), `"cloglog"`, `"logit"`, and
  `"probit"`.
  The bounded links constrain the detection component to `(0, 1)` while
  preserving the population-size identity
  `xi_i = mu_i / rho_i = N_i^alpha_i`. The legacy spelling `"logistic"`
  remains accepted as a backward-compatible alias and is normalized to
  `"logit"`.

* **Moment-based count estimation** (`estimator = "gmm"` / `"el"`):
  Poisson and Negative Binomial models can now be estimated with the same
  score equations via the `momentfit` package in addition to the default
  MLE path. The fitted object stores the estimator type and propagates it
  through `update()`, `loo()`, `bootstrap_popsize()`, `predict()`, and the
  Shiny app.

* **Robust covariance support for non-MLE count fits**: `sandwich::vcovHC()`,
  `sandwich::vcovCL()`, and `fwb::vcovFWB()` now work for `estimator = "gmm"`
  and `estimator = "el"`. HC0/HC1 are supported directly; HC2+ requests are
  downgraded to HC1 with a message.

* **`profile_alpha()` and `profile_beta()`**: Profile likelihood functions
  analogous to `profile_gamma()`. Evaluate log-likelihood and population size
  over a grid of values for a chosen alpha or beta coefficient.
  `reoptimize = FALSE` (default) holds other parameters at their MLE
  (concentrated profile); `reoptimize = TRUE` re-optimizes all nuisance
  parameters at each grid point (true profile). Both produce two-panel
  diagnostic plots.

* **`profile.uncounted()` S3 method**: Dispatches `stats::profile()` to
  `profile_gamma()`, `profile_alpha()`, or `profile_beta()` via the `param`
  argument. Usage: `profile(fit, param = "alpha")`.

* **`bootstrap_popsize()` per-replicate parameters**: New `$boot_params` field
  in the output — an R x n_params matrix storing coefficient estimates (plus
  gamma and theta when applicable) from each bootstrap replicate. Enables
  post-hoc analysis of parameter correlations and bootstrap distributions
  without writing manual FWB loops.

## Bug fixes

* **Non-MLE bias correction**: `popsize(..., bias_correction = TRUE)` now uses
  the same Fisher-style model-based variance for Poisson/NB mean parameters at
  the fitted `gmm`/`el` estimates as in the MLE case, matching the paper’s
  lognormal bias-correction derivation.

* **`predict()` gamma stabilization**: Predictions now apply the same cap for
  covariate-varying gamma linear predictors as the fitting code, avoiding
  overflow and keeping `predict(newdata = ...)` aligned with the fitted model.

* **Cluster-robust variance guard**: Clustered HC/CR covariance paths now stop
  cleanly when fewer than two clusters are supplied instead of failing with a
  low-level matrix error.

* **Shiny data/model-comparison fixes**: Uploaded count columns selected as
  `m`, `n`, or `N` are no longer silently coerced to factors by the
  auto-factor heuristic, and the comparison plot now falls back gracefully
  when likelihood-based metrics such as AIC are unavailable for non-MLE fits.

* **`bootstrap_popsize(..., total = TRUE)`**: The Total row now appears in
  `$popsize` and `$popsize_full` data frames (previously only in the separate
  `$total` list). The bias-corrected total now uses the correct full-gradient
  delta-method computation from `popsize(total = TRUE)` instead of naively
  summing independent per-group bias corrections.

## Documentation

* **Theory and API documentation refresh**: The README, theory vignette, and
  `estimate_hidden_pop()` help page now distinguish the paper’s theoretical
  decomposition `mu = xi * rho` from the baseline empirical power-law
  specification and the package’s bounded-link extensions. The constrained
  parameterization is documented explicitly as
  `alpha_values = logit^{-1}(X_alpha a)` and `beta_values = exp(X_beta b)`.

* **README and examples**: The README now reflects the current package surface
  (`iols`, bounded `link_rho`, count-model `gmm`/`el`), uses the `ukr`
  interaction example in `cov_alpha`, suppresses bootstrap progress via
  `verbose = FALSE`, and updates the funding line to NCN OPUS 20.

# uncounted 1.0.0

## New features

* **Covariate-varying gamma** (`cov_gamma`): The gamma offset parameter can now
  vary across observations via a formula interface, analogous to `cov_alpha` and
  `cov_beta`. The parameterization is
  $\gamma_i = \exp(\mathbf{X}_{\gamma,i}' \boldsymbol{\delta})$, ensuring
  positivity. Available for `method = "poisson"` and `method = "nb"`.

    ```r
    fit <- estimate_hidden_pop(
      data, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
      method = "poisson", gamma = "estimate",
      cov_gamma = ~ sex
    )
    ```

  Gamma coefficients are estimated on the log scale jointly with alpha, beta,
  and theta. `summary()` displays response-scale gamma values with delta-method
  standard errors per covariate group. `predict()`, `tidy()`, `loo()`, and
  `bootstrap_popsize()` all support the new parameter.

## Bug fixes

* **`popsize()` NA crash in `loo()`**: When LOO drops an observation that makes
  the refitted model unstable (e.g., rank-deficient design matrix), `popsize()`
  could receive NA coefficients. Bare `> 0` comparisons now use `isTRUE()` to
  handle NA values safely. The `loo()` loop also wraps `popsize()` in
  `tryCatch()` so individual failed iterations are skipped rather than killing
  the entire analysis.

* **`tidy()` gamma standard error**: The scalar gamma row now reports the
  correct response-scale SE via the delta method
  (`SE(gamma) = gamma * SE(log gamma)`). The Wald statistic and p-value are
  set to `NA` because there is no natural null hypothesis for a positive
  nuisance parameter. Confidence intervals use the log-scale transformation
  (`gamma * exp(+/- z * SE_log)`) ensuring positivity.

* **`lrtest()` gamma nesting**: Now validates gamma submodel nesting via
  column-space inclusion when both models use `cov_gamma`. Also correctly
  detects that a fixed-gamma model is not nested in a `cov_gamma` design
  without an intercept (e.g., `~ 0 + z`).

* **`profile_gamma()` on `cov_gamma` fits**: Now stops with a clear error
  instead of silently profiling a scalar gamma grid on a varying-gamma model.
  The Shiny diagnostics tab also handles this case gracefully.

## Internal

* Version bump to 1.0.0 reflecting the new `cov_gamma` feature and stabilized
  API.

* 52 new tests for `cov_gamma` including simulation-based coefficient recovery,
  one-column formula edge cases, `lrtest()` nesting validation, and
  `profile_gamma()` error handling. 650 tests total.

* `.count_params()` updated to count multi-parameter gamma for correct AIC/BIC.

# uncounted 0.7.3

## Bug fixes

* **`predict()` factor encoding**: Now uses the training data's `terms` object
  and `xlev` for consistent dummy encoding in `newdata`. Reordered factor
  levels or single-level subsets produce correct predictions.

* **`tidy()` confidence intervals**: OLS/NLS models now use `qt()` (matching
  the t-based p-values) instead of `qnorm()` for confidence intervals.

# uncounted 0.7.2

## New features

* **`plot(popsize(fit))`**: New S3 plot method for population size estimates.
  `type = "estimate"` (default) shows bias-corrected estimates with CIs;
  `type = "compare"` shows plug-in and bias-corrected side by side.

* **`compare_popsize(fit1, fit2, ...)`**: Compare population size estimates
  across models, with print and plot methods. Supports `by` parameter for
  grouped comparison and custom labels.

* **`predict(fit, newdata)`**: Prediction method supporting new data, with
  `type = "response"` (counts) or `type = "link"` (log scale). Enables
  compatibility with the `marginaleffects` package.

* **`tidy(fit)`** and **`glance(fit)`**: broom-compatible methods enabling
  `modelsummary::modelsummary()` for side-by-side model comparison tables.
  Requires `generics` package (in Suggests).

# uncounted 0.7.1

## Bug fixes

* **Multiplicative bias correction**: Replaced the subtractive Taylor
  approximation `xi_BC = xi - bias` with the exact multiplicative lognormal
  form `xi_BC = sum N^alpha * exp(-0.5 * (log N)^2 * x'Vx)`. This is exact
  under normality of alpha-hat, always positive, and reduces bias by ~94%
  in Monte Carlo simulations (vs ~88% for the old subtractive form).

* **Constrained BC direction**: Under the logit constraint, the second
  derivative `h''(eta)` can be negative when alpha is near 1, so the bias
  correction can go either direction. The code no longer assumes upward bias
  for constrained models.

* **iOLS model-based variance**: Uses `(Z'Z)^{-1}` (GPML Fisher information)
  without sigma^2 scaling, preventing inflated bias correction on zero-heavy
  data.

* **Removed nonexistent reference**: Replaced Beresewicz, Gudaszewski & Walsh
  (2025) with the correct Zhang (2008) and Beresewicz & Pawlukiewicz (2020)
  references.

## Enhancements

* **Extended iOLS tests**: 39 new tests covering bootstrap, LOO, dfbeta,
  cluster-robust SEs, HC types, the year*ukr migration specification, and
  iOLS vs Poisson numerical comparisons. 493 tests total.

# uncounted 0.7.0

## New features

* **iOLS estimator** (`method = "iols"`): Iterated OLS targeting the Gamma PML
  (GPML) score equations, following Benatia, Bellego & Pape (2024). Two-phase
  algorithm: Phase 1 warms up with increasing delta and empirical centering;
  Phase 2 uses the exact limiting transform `y_tilde = log(mu) + (m/mu - 1)/(1+rho)`
  to solve `Z'(m/mu - 1) = 0`. Converges on complex specifications including
  `year * UKR` interactions with 50% zeros. Supports HC0-HC5 sandwich SEs
  via GPML score residuals. Currently requires fixed gamma (`gamma = "estimate"`
  not yet supported).

* **LOO rank-deficiency detection**: `loo(by = "country")` now checks for
  rank-deficient design matrices before refitting. When dropping a country
  removes a covariate level (e.g., Ukraine with `year * UKR`), the refit is
  skipped and marked `converged = FALSE` with a warning.

## Enhancements

* **`lrtest()` nesting checks**: Uses column-space inclusion via QR projection
  to detect non-nested model pairs. Also warns for cross-method comparisons
  (except Poisson vs NB), different constrained settings, and incompatible
  gamma specifications.

* **`compare_models()` warns for mixed likelihood types**: Now includes iOLS
  in the pseudo-loglik warning alongside OLS/NLS, preventing silent comparison
  of GPML deviance with Poisson/NB count likelihoods.

* **Test coverage improvements**: Added tests for HC4/HC4m/HC5, LOO
  `summary()` and `plot(type = "coef")`, bootstrap `point_estimate`/`total`/`by`
  parameters, constrained summary output, and `profile_gamma(plot = TRUE)`.
  440 tests total.

# uncounted 0.6.0

## Bug fixes

* **Total CI in printer**: `.print_popsize_table()` now uses the delta-method
  total CI from `attr(ps, "total")` instead of naively summing subgroup
  confidence bounds, which ignored between-group correlation.

* **Constrained bias correction**: The second-order Taylor bias correction for
  constrained models now includes the logit curvature term
  `alpha''(eta) = alpha(1-alpha)(1-2*alpha)`. Previously only `[alpha'(eta)]^2`
  was used, producing incomplete bias correction when alpha was far from 0.5.

* **`lrtest()` nesting check**: Now uses column-space inclusion (QR projection)
  to detect non-nested model pairs, replacing the previous heuristic. Also
  checks gamma specification compatibility.

## Enhancements

* **NB sandwich with theta**: The Negative Binomial estimator now includes
  theta (dispersion) in the sandwich variance-covariance computation via a
  dedicated NB covariance path using the full per-observation score vector and
  the observed Hessian from `optim(hessian = TRUE)`. This properly accounts for
  theta estimation uncertainty in coefficient standard errors. The fitted object
  now includes `theta_se`, `score_full`, and `hessian_nll`. Supports HC0, HC1,
  and cluster-robust variance; HC2+ falls back to HC1 with a message.

* **CI workflows**: Added GitHub Actions for R-CMD-check (5 platforms),
  CRAN-like check, test coverage (Codecov), and pkgdown site deployment.

* **Unit test redesign**: 340+ tests using a shared 200-row simulated fixture,
  split into CRAN (fast) and extended (GitHub-only) suites. Includes base R
  oracle tests against `glm()`, `lm()`, `sandwich::vcovHC()`, and
  `MASS::glm.nb()`.

## Documentation

* Fixed references: replaced incorrect Zhang (2008) citation, added URLs,
  added Beresewicz & Pawlukiewicz (2020).

* OLS `log(m+1)` fitted-value scale limitation documented.

* Added `URL` and `BugReports` fields to DESCRIPTION.

* `year` column in `irregular_migration` dataset converted to factor.

# uncounted 0.5.1

* Added Aniela Czerniawska as package author (original GitHub version).

* New `dfbeta.uncounted()` method for leave-one-out influence on coefficients.
  Wrapper around `loo()` returning the familiar n x p matrix of coefficient
  changes.

* New `dfpopsize()` generic and `dfpopsize.uncounted()` method for leave-one-out
  influence on population size estimates. Returns a named vector of
  per-observation (or per-country) changes in total xi.

# uncounted 0.5.0

## Breaking changes

* The `vcov_type` parameter in `estimate_hidden_pop()` has been renamed to
  `vcov`. It now accepts either a character string (`"HC0"` through `"HC5"`) or
  a **function** that takes the fitted object and returns a variance-covariance
  matrix (e.g., `sandwich::vcovHC` or a custom wrapper around
  `sandwich::vcovCL()`).

* The `countries` parameter no longer triggers cluster-robust standard errors.
  It is now used **only** for grouping in `popsize()` and `loo()`. Use the new
  `cluster` parameter for cluster-robust variance estimation.

## New features

* New `cluster` parameter in `estimate_hidden_pop()` accepts a one-sided formula
 (e.g., `~ country_code`) to request cluster-robust variance estimation via
 `sandwich::vcovCL()`.

* Full compatibility with the **sandwich** package through new S3 methods:
  `bread()`, `estfun()`, `hatvalues()`, and `model.matrix()`. This allows
  direct use of `sandwich::vcovHC()`, `sandwich::vcovCL()`, and
  `sandwich::sandwich()` on `uncounted` objects.

* New `type = "working"` option in `residuals()` returns the score residuals
  used internally for sandwich variance computation.

* New `update.uncounted()` method enables re-fitting with modified arguments
  (e.g., different weights or vcov). Supports `evaluate = FALSE` to return
  the unevaluated call.

* New `weights.uncounted()` method returns observation weights from the fit.

* Compatibility with `fwb::vcovFWB()` for fractional weighted bootstrap
  variance estimation, thanks to `update()` and `weights()` methods.

* The `vcov` parameter accepts a function, enabling custom variance estimators:

    ```r
    # Cluster-robust via function
    fit <- estimate_hidden_pop(...,
      vcov = function(x) sandwich::vcovCL(x, cluster = data$country))

    # Observation-level via string (default)
    fit <- estimate_hidden_pop(..., vcov = "HC3")

    # Cluster-robust via string + cluster formula
    fit <- estimate_hidden_pop(..., vcov = "HC1", cluster = ~ country_code)
    ```

## Bug fixes

* Fixed a bug where `HC0` and `HC3` (and all other HC types) produced
  **identical** standard errors when clustering was active. The old
  implementation ignored `vcov_type` entirely in the cluster branch of the
  internal sandwich computation. Variance estimation is now delegated to the
  **sandwich** package, which correctly handles all HC types for both
  observation-level and cluster-robust cases.

* Fixed a bug in the cluster-robust meat matrix computation where an extra
  $\sqrt{\mu_i}$ factor was incorrectly included in the Poisson and NB score
  contributions. This is now resolved by using `sandwich::vcovCL()`.

## Internal

* Variance-covariance computation is now performed post-estimation using
  the **sandwich** package rather than the internal `.compute_sandwich_vcov()`
  function. The internal function is retained for the estimator-level fallback
  `vcov_model` (used in bias correction).

* Each estimator now returns `model_matrix_full`, `bread_weights`, and
  `score_residuals` to support the `bread()` and `estfun()` S3 methods.

* Added **sandwich** to `Imports` in `DESCRIPTION`.

* New file `R/sandwich_methods.R` containing `bread.uncounted()`,
  `estfun.uncounted()`, `hatvalues.uncounted()`, and
  `model.matrix.uncounted()`.

* Expanded test suite to 245 tests covering sandwich methods, cluster-robust
  variance, function-based vcov, and the HC0 != HC3 regression test.

# uncounted 0.4.0

## New features

* `popsize()` and `bootstrap_popsize()` gain a `by` parameter for stratified
  population size estimates (e.g., `by = ~ year`, `by = ~ year + sex`).
  Previously, stratification was only available through covariate-varying alpha.

* New pseudo $R^2$ measures in `compare_models()`: explained deviance ($R^2_D$)
  and Cameron-Windmeijer ($R^2_{CW}$) for Poisson and Negative Binomial models.

* Interactive **Shiny application** via `run_app()` with five tabs: Data,
  Model Fitting, Population Size, Model Comparison, and Diagnostics.

## Improvements

* Rootogram: white background bars for hanging style, light-blue bars for
  standing and suspended styles for clearer visual contrast.

* `loo()` call arguments are now force-evaluated and stored in the call object,
  fixing "object not found" errors when refitting from within Shiny or other
  non-standard evaluation contexts.

* Shiny app: automatic conversion of year-like columns to factors to prevent
  continuous covariate issues.

## Bug fixes

* Fixed Cameron-Windmeijer reference year: 1997 to 1996 (JBES 14(2), 209-220).

* Fixed non-ASCII characters in country names (`Cote d'Ivoire`,
  `Sao Tome e Principe`) causing encoding issues on some platforms.

# uncounted 0.3.1

* Renamed internal `sigmoid` to `logit` / `inv_logit` throughout the package
  for consistency with standard statistical terminology.

# uncounted 0.3.0

## New features

* Added `irregular_migration` built-in dataset (Poland 2019-2024) with
  documentation.

* Added package vignettes: `getting-started` and `model-workflow`.

* Full roxygen2 documentation for all exported functions.

# uncounted 0.2.0

## New features

* `popsize()` replaces `xi()` as the primary interface for population size
  estimation. `xi()` is retained as a deprecated alias.

* Bias-corrected population size estimates using second-order Taylor expansion
  (Jensen's inequality correction).

* Delta-method confidence intervals for total population size.

* `profile_gamma()` for profiling the gamma parameter over a grid.

* `compare_models()` for side-by-side model comparison (AIC, BIC, log-likelihood).

* `lrtest()` for likelihood-ratio tests between nested models.

* `loo()` for leave-one-out diagnostics by observation or by country.

* `compare_loo()` for comparing LOO influence across models.

* Constrained estimation (`constrained = TRUE`): logit link for
  $\alpha \in (0, 1)$ and log link for $\beta > 0$.

* Fractional weighted bootstrap via `bootstrap_popsize()` using the **fwb**
  package (Xu et al. 2020), with cluster support.

* Observation weights via `weights` parameter.

* HC4, HC4m, and HC5 robust standard error types.

* `rootogram()` for assessing distributional fit (hanging, standing, suspended).

## Improvements

* Default `gamma = "estimate"` (was `NULL`).

* Model-based (non-sandwich) variance stored as `vcov_model` for bias
  correction, separate from HC-robust `vcov`.

* Improved numerical stability for extreme parameter values.

* `~0` formula treated as intercept-only (not empty model).

# uncounted 0.1.0

* Initial release.

* `estimate_hidden_pop()` for fitting the power-law model
  $E(m_i) = N_i^{\alpha_i} (\gamma + n_i/N_i)^{\beta_i}$ using OLS, NLS,
  Poisson PML, or Negative Binomial MLE.

* Covariate-varying $\alpha$ and $\beta$ via formula interface (`cov_alpha`,
  `cov_beta`).

* Gamma parameter: estimated, fixed, or excluded.

* HC0-HC3 robust standard errors.

* `print()` and `summary()` S3 methods.
