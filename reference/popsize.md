# Estimated Population Size

Computes the estimated total unauthorized population \\\hat{\xi} =
\sum\_{i} N_i^{\hat{\alpha}}\\ for each group defined by the covariates
in alpha. Includes bias correction via a multiplicative lognormal
adjustment for unconstrained models (and a second-order Taylor
approximation for constrained models), plus log-normal delta-method
confidence intervals on the population-size scale.

## Usage

``` r
popsize(object, ...)

# S3 method for class 'uncounted'
popsize(
  object,
  by = NULL,
  level = 0.95,
  bias_correction = TRUE,
  total = FALSE,
  ...
)
```

## Arguments

- object:

  An `"uncounted"` object.

- ...:

  Additional arguments (ignored).

- by:

  Optional formula specifying grouping variables for stratified
  population size estimation (e.g., `~ year`, `~ country`,
  `~ year + sex`). The variables must exist in the data used to fit the
  model. When `NULL` (default), groups are defined by the alpha
  covariate pattern from `cov_alpha`. When provided, the same fitted
  alpha values are used but summed over the `by`-groups instead. This is
  analogous to `singleRcapture::stratifyPopsize()`.

- level:

  Confidence level for intervals (default 0.95).

- bias_correction:

  Logical; apply analytical bias correction? Default TRUE. Uses
  model-based variance (not HC-robust) to avoid overcorrection from
  inflated leverage-driven standard errors.

- total:

  Logical; if `TRUE` and multiple groups exist, compute a delta-method
  total with SE and CI, stored in `attr(result, "total")`. Default
  `FALSE`. **Warning**: for panel data where groups are defined by time
  periods (e.g., `by = ~ year`), the total sums population estimates
  across years. This is generally not meaningful because the same
  individuals may appear in multiple years. The total is only
  appropriate when groups represent non-overlapping subpopulations
  (e.g., `by = ~ sex` within a single year).

## Value

A data frame with columns:

- group:

  Group label derived from alpha covariates, or `"(all)"` when no alpha
  covariates are specified.

- estimate:

  Plug-in estimate \\\hat{\xi}\_g = \sum N_i^{\hat{\alpha}\_g}\\.

- estimate_bc:

  Bias-corrected estimate \\\hat{\xi}^{BC}\_g\\; `NA` when
  `bias_correction = FALSE`.

- lower:

  Lower bound of the (bias-corrected) confidence interval.

- upper:

  Upper bound of the (bias-corrected) confidence interval.

- share_pct:

  Group share as percentage of total \\\hat{\xi}\\.

When multiple groups exist, `attr(result, "total")` contains the total
estimate, bias-corrected estimate, standard error, and CI.

## Details

**Point estimate.** For a group \\g\\ with reference populations \\N_1,
\ldots, N_n\\ and estimated exponent \\\hat{\alpha}\_g\\, the plug-in
population size is \$\$\hat{\xi}\_g = \sum\_{i=1}^{n}
N_i^{\hat{\alpha}\_g}.\$\$

**Bias correction.** Because \\\xi\\ is a nonlinear function of
\\\hat{\alpha}\\, the plug-in estimate is biased upward (by Jensen's
inequality, \\E\[N^{\hat{\alpha}}\] \geq N^{E\[\hat{\alpha}\]}\\ when
\\N \> 1\\). A second-order Taylor expansion of \\h(\alpha) = \sum_i
N_i^{\alpha}\\ around \\\alpha_0\\ gives the approximate bias \$\$
\mathrm{Bias}(\hat{\xi}\_g) \approx \frac{1}{2} \sum\_{i=1}^{n}
N_i^{\alpha_g} (\log N_i)^2 \\ \mathbf{x}\_g' \mathbf{V} \mathbf{x}\_g,
\$\$ where \\\mathbf{x}\_g\\ is the design vector for group \\g\\ and
\\\mathbf{V}\\ is the variance-covariance matrix of \\\hat{\alpha}\\.
For unconstrained models, the exact multiplicative correction is used:
\$\$ \hat{\xi}^{BC}\_g = \sum\_{i=1}^{n} N_i^{\hat{\alpha}\_g}
\exp\\\left(-\frac{1}{2} (\log N_i)^2 \mathbf{x}\_g' \mathbf{V}
\mathbf{x}\_g\right), \$\$ which is exact under normality of
\\\hat{\alpha}\\ and always positive. For constrained models the
subtractive Taylor correction is used instead (the logistic-normal
integral has no closed form), and the bias can be positive or negative
depending on \\\alpha\\. Model-based variance is used for bias
correction rather than HC-robust variance. For Poisson and NB count
models this means the Fisher-style inverse information for the
mean-model parameters, evaluated at the fitted coefficients. This
remains true when the coefficients were obtained by `estimator = "gmm"`
or `estimator = "el"`: the robust HC or FWB covariance is still used for
confidence intervals, but the bias correction uses the same Fisher-style
model variance as in the MLE case, evaluated at the non-MLE estimate.
For iOLS/GPML, the model-based variance is
\\(\mathbf{Z}'\mathbf{Z})^{-1}\\ (Gamma Fisher information, no
dispersion scaling).

When `constrained = TRUE`, the delta method accounts for the logit link:
\\\mathrm{Var}(\alpha) = \mathrm{Var}(\eta) \cdot \[\sigma'(\eta)\]^2\\
where \\\sigma'(\eta) = \alpha(1-\alpha)\\.

**Confidence intervals via the delta method.** Let \\\mathbf{g}\_g =
\partial \hat{\xi}\_g / \partial \boldsymbol{\alpha}\\ denote the
gradient of the plug-in estimator with respect to the alpha
coefficients. For unconstrained models, \\\mathbf{g}\_g = \sum_i
N_i^{\hat{\alpha}\_i} \log(N_i)\mathbf{x}\_i;\\ for constrained models
the same expression is multiplied by the derivative of the inverse-logit
map. Using the HC-robust covariance \\\mathbf{V}\\, the package computes
\$\$ \widehat{\mathrm{Var}}(\hat{\xi}\_g) = \mathbf{g}\_g^\top
\mathbf{V}\mathbf{g}\_g. \$\$ To preserve positivity, the subgroup
interval is then reported on a log-normal scale: \$\$ \hat{\xi}\_g
\exp\\\left( \pm z\_{\alpha/2}
\frac{\widehat{\mathrm{se}}(\hat{\xi}\_g)}{\hat{\xi}\_g} \right). \$\$
When `bias_correction = TRUE`, the lower and upper bounds are rescaled
by \\\hat{\xi}^{BC}\_g / \hat{\xi}\_g\\, matching the returned
bias-corrected point estimate.

**Total across groups.** When multiple groups exist, the total
\\\hat{\xi} = \sum_g \hat{\xi}\_g\\ has its own delta-method CI computed
via the gradient \\\nabla\_\alpha \xi\\ and a log-normal approximation
for positivity. This is stored in `attr(result, "total")`.

## Examples

``` r
# Simulate synthetic data for 5 countries, 3 years each
set.seed(42)
n_obs <- 15
sim_data <- data.frame(
  country = rep(paste0("C", 1:5), each = 3),
  year    = rep(2018:2020, 5),
  N       = rpois(n_obs, lambda = 500000),
  n       = rpois(n_obs, lambda = 1000),
  m       = rpois(n_obs, lambda = 50)
)

# Fit a Poisson model
fit <- estimate_hidden_pop(
  data = sim_data, observed = ~m, auxiliary = ~n,
  reference_pop = ~N, method = "poisson"
)
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.

# Population size with bias correction and 95% CI
popsize(fit)
#>   group observed   estimate estimate_bc lower upper share_pct
#> 1 (all)      752 0.03262174           0     0     0       100

# Without bias correction
popsize(fit, bias_correction = FALSE)
#>   group observed   estimate estimate_bc         lower         upper share_pct
#> 1 (all)      752 0.03262174          NA 3.457442e-219 3.077935e+215       100

# 90% confidence interval
popsize(fit, level = 0.90)
#>   group observed   estimate estimate_bc lower upper share_pct
#> 1 (all)      752 0.03262174           0     0     0       100
```
