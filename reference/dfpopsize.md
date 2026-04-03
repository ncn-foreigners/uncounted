# Leave-One-Out Influence on Population Size

Computes the change in the total estimated population size \\\hat\xi =
\sum_i N_i^{\hat\alpha_i}\\ when each observation (or country) is
dropped. This is a convenience wrapper around
[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md).

## Usage

``` r
dfpopsize(model, ...)

# S3 method for class 'uncounted'
dfpopsize(model, by = c("obs", "country"), ...)
```

## Arguments

- model:

  A fitted `uncounted` object.

- ...:

  Additional arguments passed to
  [`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md).

- by:

  Character: `"obs"` (default) or `"country"`.

## Value

A named numeric vector. Each element is the change in \\\hat\xi\\ when
that unit is dropped (\\\hat\xi\_{(-i)} - \hat\xi\\).

## See also

[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md),
[`dfbeta.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/dfbeta.uncounted.md)

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == 2019 & irregular_migration$sex == "m", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.001, countries = ~ country_code)
#> Error in lm.fit(Z_start, y_start): 0 (non-NA) cases
dp <- dfpopsize(fit)
#> Error: object 'fit' not found
head(dp)
#> Error: object 'dp' not found
```
