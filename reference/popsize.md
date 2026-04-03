# Estimated Population Size

Computes the estimated total unauthorized population \\\hat{\xi} =
\sum\_{i} N_i^{\hat{\alpha}}\\ for each group defined by the covariates
in alpha. Includes bias correction via second-order Taylor expansion and
confidence intervals via monotone transformation of the Wald interval on
the link scale.

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

  Logical; apply Taylor-expansion bias correction? Default TRUE. Uses
  model-based variance (not HC-robust) to avoid overcorrection from
  inflated leverage-driven standard errors.

## Value

A data frame with columns:

- group:

  Group label derived from alpha covariates, or `"(all)"` when no alpha
  covariates are specified.

- estimate:

  Plug-in estimate \\\hat{\xi}\_g = \sum N_i^{\hat{\alpha}\_g}\\.

- estimate_bc:

  Bias-corrected estimate \\\hat{\xi}^{BC}\_g\\.

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

**Bias correction via Taylor expansion.** Because \\\xi\\ is a nonlinear
function of \\\hat{\alpha}\\, the plug-in estimate is biased upward (by
Jensen's inequality, \\E\[N^{\hat{\alpha}}\] \geq
N^{E\[\hat{\alpha}\]}\\ when \\N \> 1\\). A second-order Taylor
expansion of \\h(\alpha) = \sum_i N_i^{\alpha}\\ around \\\alpha_0\\
gives the approximate bias \$\$ \mathrm{Bias}(\hat{\xi}\_g) \approx
\frac{1}{2} \sum\_{i=1}^{n} N_i^{\alpha_g} (\log N_i)^2 \\
\mathbf{x}\_g' \mathbf{V} \mathbf{x}\_g, \$\$ where \\\mathbf{x}\_g\\ is
the design vector for group \\g\\ and \\\mathbf{V}\\ is the
variance-covariance matrix of \\\hat{\alpha}\\. The bias-corrected
estimate is \\\hat{\xi}^{BC}\_g = \hat{\xi}\_g -
\widehat{\mathrm{Bias}}\\. Model-based (homoscedastic) variance is used
for bias correction rather than HC-robust variance, because HC3 can be
inflated by high-leverage observations in skewed data, leading to
overcorrection.

When `constrained = TRUE`, the delta method accounts for the logit link:
\\\mathrm{Var}(\alpha) = \mathrm{Var}(\eta) \cdot \[\sigma'(\eta)\]^2\\
where \\\sigma'(\eta) = \alpha(1-\alpha)\\.

**Confidence intervals via monotone transformation.** A Wald interval is
first constructed on the link scale for the linear predictor:
\$\$\hat{\eta}\_g \pm z\_{\alpha/2} \cdot
\mathrm{se}(\hat{\eta}\_g),\$\$ where the standard error uses the
HC-robust variance. The interval endpoints are then mapped through the
monotone transformation \\g(\alpha) = \sum_i N_i^{\alpha}\\ (increasing
for \\N_i \geq 1\\) to obtain \\\[\hat{\xi}\_L, \hat{\xi}\_U\]\\. When
`constrained = TRUE`, the logit link is applied before exponentiation.
Bias correction is also applied to the CI bounds at their respective
alpha values.

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
#>   group observed   estimate estimate_bc          lower          upper share_pct
#> 1 (all)      752 0.03262174    -1977.66 -2.096793e-214 -1.865297e+220       100

# Without bias correction
popsize(fit, bias_correction = FALSE)
#>   group observed   estimate estimate_bc         lower         upper share_pct
#> 1 (all)      752 0.03262174  0.03262174 3.458684e-219 3.076829e+215       100

# 90% confidence interval
popsize(fit, level = 0.90)
#>   group observed   estimate estimate_bc          lower          upper share_pct
#> 1 (all)      752 0.03262174    -1977.66 -1.604466e-179 -2.437659e+185       100
```
