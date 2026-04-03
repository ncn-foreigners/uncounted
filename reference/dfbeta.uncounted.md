# Leave-One-Out Influence on Coefficients

Computes the change in estimated coefficients when each observation (or
country) is dropped. This is a convenience wrapper around
[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md).

## Usage

``` r
# S3 method for class 'uncounted'
dfbeta(model, by = c("obs", "country"), ...)
```

## Arguments

- model:

  A fitted `uncounted` object.

- by:

  Character: `"obs"` (default) drops one observation at a time,
  `"country"` drops all observations for one country at a time (requires
  `countries` in the original fit).

- ...:

  Additional arguments passed to
  [`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md).

## Value

A matrix with one row per dropped unit and one column per coefficient.
Each entry is the change in the coefficient estimate relative to the
full model (\\\hat\beta\_{(-i)} - \hat\beta\\).

## See also

[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md),
[`dfpopsize`](https://ncn-foreigners.github.io/uncounted/reference/dfpopsize.md)

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == 2019 & irregular_migration$sex == "m", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.001, countries = ~ country_code)
#> Error in lm.fit(Z_start, y_start): 0 (non-NA) cases
db <- dfbeta(fit)
#> Error: object 'fit' not found
head(db)
#> Error: object 'db' not found
```
