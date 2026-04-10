# Profile Likelihood for Beta Coefficients

Evaluates the log-likelihood and population size over a grid of values
for a chosen beta coefficient. When `reoptimize = FALSE`, the population
size \\\hat\xi\\ is constant because beta does not enter \\\xi = \sum
N^{\alpha}\\ directly.

## Usage

``` r
profile_beta(
  object,
  coef_index = 1,
  grid = NULL,
  reoptimize = FALSE,
  plot = TRUE,
  ...
)
```

## Arguments

- object:

  An `"uncounted"` object.

- coef_index:

  Integer: which alpha coefficient to profile (1 = intercept).

- grid:

  Numeric vector of values to evaluate. If `NULL` (default),
  auto-generates 30 points spanning \\\pm 3\\ SE around the MLE.

- reoptimize:

  Logical. If `TRUE`, optimize all other parameters at each grid point
  (true profile likelihood). If `FALSE` (default), hold everything else
  at the MLE (concentrated profile, faster).

- plot:

  Logical; produce the plot? Default `TRUE`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a data frame with columns: `value`, `xi`, `loglik`.

## See also

[`profile_alpha`](https://ncn-foreigners.github.io/uncounted/reference/profile_alpha.md),
[`profile_gamma`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md),
[`profile.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/profile.uncounted.md)
