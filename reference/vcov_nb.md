# Theta-Aware Full Variance-Covariance for NB Models

Returns the full sandwich variance-covariance matrix for NB models,
including the theta (dispersion) parameter. For non-NB models, returns
`vcov(object)` unchanged.

## Usage

``` r
vcov_nb(object, vcov_type = "HC1", cluster = NULL)
```

## Arguments

- object:

  A fitted `uncounted` object.

- vcov_type:

  HC type for the sandwich (`"HC0"` or `"HC1"`). Default `"HC1"`. HC2+
  are not available for the theta-aware path.

- cluster:

  Optional cluster vector for cluster-robust variance.

## Value

A square variance-covariance matrix. For NB models, dimensions are
`p + 1` (or `p + 2` if gamma is estimated), where the extra row/column
corresponds to `log(theta)`. For non-NB models, returns `vcov(object)`.

## Details

This is the explicit interface for the theta-aware covariance that
[`vcov()`](https://rdrr.io/r/stats/vcov.html) uses internally for NB
fits. Use this when you need the full matrix including theta, or when
you want to be explicit about which covariance you are requesting.

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == "2019", ]
fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.005)
vcov_nb(fit)  # includes theta row/column
#>              alpha        beta            
#> alpha 0.0008780469 0.001285950 0.001801256
#> beta  0.0012859499 0.002185531 0.002006601
#>       0.0018012560 0.002006601 0.044811875
vcov(fit)     # alpha/beta only (same as vcov_nb submatrix)
#>              alpha        beta
#> alpha 0.0008780469 0.001285950
#> beta  0.0012859499 0.002185531
```
