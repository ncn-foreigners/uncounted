# Extract Posterior Population-Size Draws

Extract posterior draws of grouped population-size estimates from
Bayesian uncounted fits or Bayesian
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
results.

## Usage

``` r
posterior_popsize_draws(object, ...)

# S3 method for class 'uncounted_bayes'
posterior_popsize_draws(
  object,
  by = NULL,
  include_total = FALSE,
  format = c("matrix", "data.frame"),
  ...
)

# S3 method for class 'uncounted_popsize_bayes'
posterior_popsize_draws(
  object,
  include_total = FALSE,
  format = c("matrix", "data.frame"),
  ...
)
```

## Arguments

- object:

  A Bayesian fitted model from
  [`estimate_hidden_pop_bayes`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop_bayes.md)
  or a Bayesian population-size object from
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

- ...:

  Additional arguments passed to
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  when `object` is an `"uncounted_bayes"` fit.

- by:

  Optional grouping formula when `object` is an `"uncounted_bayes"` fit.

- include_total:

  Logical; if `TRUE`, append a `"Total"` draw column.

- format:

  Output format: `"matrix"` or `"data.frame"`.

## Value

A draw-by-group matrix or wide data frame.
