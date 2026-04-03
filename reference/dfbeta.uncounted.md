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
d <- irregular_migration[irregular_migration$year == "2019" & irregular_migration$sex == "Male", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.001, countries = ~ country_code)
db <- dfbeta(fit)
head(db)
#>           alpha          beta
#> 1  4.618242e-03  1.448704e-02
#> 2 -1.379703e-04 -4.409443e-04
#> 3  3.970809e-03  1.264425e-02
#> 4  3.043501e-04  9.516993e-04
#> 5  1.259357e-03  9.422888e-04
#> 6 -2.749349e-05 -8.519565e-05
```
