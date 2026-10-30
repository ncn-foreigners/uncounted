# Leave-One-Out Influence on Population Size

Computes the influence of each observation (or country) on the total
estimated population size \\\hat\xi = \sum_i N_i^{\hat\alpha_i}\\,
holding the target population fixed. This is a convenience wrapper
around
[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md);
see there for the definition.

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

A named numeric vector. Each element is the change in \\\hat\xi\\,
summed over all rows, when the coefficients are re-estimated without
that unit (`dxi` in
[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md)).

## See also

[`loo`](https://ncn-foreigners.github.io/uncounted/reference/loo.md),
[`dfbeta.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/dfbeta.uncounted.md)

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == "2019" & irregular_migration$sex == "Male", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.001, countries = ~ country_code)
dp <- dfpopsize(fit)
head(dp)
#>          1          2          3          4          5          6 
#> 4156.06841 -120.79823 3560.03090  267.14848 1111.51823  -24.08684 
```
