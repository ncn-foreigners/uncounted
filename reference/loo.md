# Leave-One-Out Sensitivity Analysis

Refits the model dropping one observation (or one country/group) at a
time. Returns the impact on coefficient estimates, fitted values, and
\\\xi\\ (total population size estimate).

## Usage

``` r
loo(object, ...)

# S3 method for class 'uncounted'
loo(object, by = c("obs", "country"), verbose = FALSE, ...)
```

## Arguments

- object:

  An `"uncounted"` object (fitted model).

- ...:

  Additional arguments (ignored).

- by:

  What to leave out: `"obs"` drops one row at a time, `"country"` drops
  all rows for one country at a time (requires `countries` to have been
  specified in the original fit).

- verbose:

  Logical; print progress?

## Value

An object of class `"uncounted_loo"` with components:

- coefficients:

  Matrix (`n_drops` x `p`): coefficients from each refit.

- xi:

  Data frame with \\\xi\\ estimates from each refit.

- dropped:

  Character vector: label of dropped observation or country.

- dfbeta:

  Matrix (`n_drops` x `p`): influence on the coefficients,
  \\\hat{\beta} - \hat{\beta}\_{(-i)}\\, the same sign as
  [`stats::dfbeta()`](https://rdrr.io/r/stats/influence.measures.html).
  Positive means the unit pulls the coefficient up.

- dxi:

  Numeric vector: influence on the population size for a fixed target,
  using the refit's coefficients on all rows (see Details).

- dxi_subset:

  Numeric vector: change in population size when the unit also leaves
  the target (see Details).

- full_coefs:

  Full-model coefficients.

- full_ps:

  Full-model population size data frame.

- full_ps_total:

  Full-model total \\\hat{\xi}\\.

- by:

  The `by` argument used.

- converged:

  Logical vector: whether each refit converged.

- n_drops:

  Number of leave-one-out iterations.

## Details

**Purpose.** LOO sensitivity analysis assesses how much each observation
(or country) influences the estimated population size. An observation
with a large \\\|\Delta\xi\|\\ relative to the full-model \\\hat{\xi}\\
is "influential" and may warrant closer inspection.

**`by = "obs"` vs `by = "country"`.** With `by = "obs"`, each row is
dropped individually (n refits). This identifies individual data points
that drive the estimate, useful for detecting outliers or data errors.
With `by = "country"`, all rows belonging to one country are dropped
simultaneously (one refit per country). This measures each country's
overall contribution to \\\hat{\xi}\\ and is more relevant for assessing
structural sensitivity: would the conclusion change if a country were
excluded?

**Interpretation of `dxi`.** `dxi` measures influence on the estimate
for a fixed target population: \\\Delta\xi_i = \sum_j
N_j^{\hat\alpha\_{j(-i)}} - \sum_j N_j^{\hat\alpha_j}\\, where both sums
run over *all* rows and \\\hat\alpha\_{(-i)}\\ comes from the refit
without unit \\i\\. It therefore reflects only how the coefficients
move.

- \\\Delta\xi_i \< 0\\: without \\i\\ the estimate is lower (the unit
  was pulling the estimate up).

- \\\Delta\xi_i \> 0\\: without \\i\\ the estimate is higher (the unit
  was pulling the estimate down).

`dxi_subset` is the change when the unit also leaves the target: the
total of
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
on the refit, over the remaining rows only, minus \\\hat\xi\\. It is
mostly negative because it includes the unit's own share of \\\hat\xi\\,
so it is not an influence measure, but it answers “what would we
estimate if this unit were not in the data at all?”

## Examples

``` r
# Simulate synthetic data
set.seed(42)
n_obs <- 15
sim_data <- data.frame(
  country = rep(paste0("C", 1:5), each = 3),
  year    = rep(2018:2020, 5),
  N       = rpois(n_obs, lambda = 500000),
  n       = rpois(n_obs, lambda = 1000),
  m       = rpois(n_obs, lambda = 50)
)

fit <- estimate_hidden_pop(
  data = sim_data, observed = ~m, auxiliary = ~n,
  reference_pop = ~N, method = "poisson",
  countries = ~country
)
#> Warning: Some alpha values < 0 (min = -0.046). Consider using constrained = TRUE.
#> Estimated gamma is at its lower bound (1e-10); the optimum may lie outside gamma_bounds = c(1e-10, 0.5).

# LOO by observation (drops one row at a time)
loo_obs <- loo(fit, by = "obs")
#> Warning: Some alpha values < 0 (min = -0.064). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.211). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.057). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.025). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.451). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.051). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.087). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.064). Consider using constrained = TRUE.
print(loo_obs)       # top 10 most influential observations
#> Leave-one-out sensitivity analysis
#> Dropped by: obs 
#> N iterations: 15 
#> Converged: 15 / 15 
#> 
#> Full model xi: 8.2 
#> LOO xi range: 0 to 49.6 
#> 
#> Most influential (by |delta xi|):
#>  dropped  dxi pct_change
#>       13 41.4     503.14
#>       12 29.6     360.64
#>        6 12.9     157.44
#>        9  9.6     116.51
#>        4  8.4     102.78
#>        1  8.2      99.69
#>        8 -8.2     -99.51
#>        3 -7.3     -88.56
#>       14  7.0      84.99
#>       11 -3.4     -41.85
summary(loo_obs)     # coefficient and xi stability
#> Leave-one-out sensitivity analysis
#> Dropped by: obs 
#> N iterations: 15 (converged: 15 )
#> 
#> === Coefficient stability ===
#>            full  loo_mean   loo_sd   loo_min   loo_max
#> alpha -0.045833 -0.052940 0.131550 -0.450711  0.091103
#> beta  -0.726555 -0.741571 0.278221 -1.583679 -0.435062
#> 
#> === Xi stability ===
#> Full xi: 8.2 
#> LOO mean: 14.6 
#> LOO range: 0 to 49.6 
#> Max |%change|: 503.14 %
plot(loo_obs)        # bar plot of dxi

plot(loo_obs, type = "coef")  # DFBETA plots per coefficient


# LOO by country (drops all rows for one country)
loo_ctry <- loo(fit, by = "country")
#> Warning: Some alpha values < 0 (min = -0.171). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.406). Consider using constrained = TRUE.
print(loo_ctry)
#> Leave-one-out sensitivity analysis
#> Dropped by: country 
#> N iterations: 5 
#> Converged: 5 / 5 
#> 
#> Full model xi: 8.2 
#> LOO xi range: 0.1 to 152.3 
#> 
#> Most influential (by |delta xi|):
#>  dropped   dxi pct_change
#>       C5 144.0    1752.38
#>       C2  73.9     898.60
#>       C4  14.5     176.66
#>       C3  -8.1     -99.11
#>       C1  -6.6     -80.74
```
