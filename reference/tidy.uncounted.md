# Tidy an uncounted model

Extracts a tidy coefficient table from a fitted `uncounted` object.
Compatible with
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html).

## Usage

``` r
# S3 method for class 'uncounted'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  An `uncounted` object.

- conf.int:

  Logical; include confidence intervals? Default FALSE.

- conf.level:

  Confidence level for intervals (default 0.95).

- ...:

  Additional arguments (ignored).

## Value

A data frame with columns: `term`, `estimate`, `std.error`, `statistic`,
`p.value`, and optionally `conf.low`, `conf.high`.

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == "2019", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.005)
tidy(fit)
#>    term  estimate  std.error statistic      p.value
#> 1 alpha 0.8426007 0.05247666 16.056676 5.133340e-58
#> 2  beta 0.7116506 0.12053883  5.903912 3.549823e-09
tidy(fit, conf.int = TRUE)
#>    term  estimate  std.error statistic      p.value  conf.low conf.high
#> 1 alpha 0.8426007 0.05247666 16.056676 5.133340e-58 0.7397483 0.9454531
#> 2  beta 0.7116506 0.12053883  5.903912 3.549823e-09 0.4753988 0.9479024
```
