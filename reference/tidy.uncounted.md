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
#> Error in tidy(fit): could not find function "tidy"
tidy(fit, conf.int = TRUE)
#> Error in tidy(fit, conf.int = TRUE): could not find function "tidy"
```
