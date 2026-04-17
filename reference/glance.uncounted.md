# Glance at an uncounted model

Returns a one-row data frame summarizing model fit, compatible with
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html).

## Usage

``` r
# S3 method for class 'uncounted'
glance(x, ...)
```

## Arguments

- x:

  An `uncounted` object.

- ...:

  Additional arguments (ignored).

## Value

A one-row data frame with columns: `method`, `nobs`, `logLik`, `AIC`,
`BIC`, `deviance`, `df.residual`, `gamma`, `theta`.

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == "2019", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.005)
glance(fit)
#>    method estimator link_rho nobs    logLik      AIC      BIC deviance
#> 1 POISSON       MLE    power  203 -1011.373 2026.747 2033.373 1560.452
#>   df.residual
#> 1         201
```
