# Bayesian Residuals

Bayesian Residuals

## Usage

``` r
# S3 method for class 'uncounted_bayes'
residuals(
  object,
  type = c("response", "pearson"),
  summary = c("mean", "median"),
  ...
)
```

## Arguments

- object:

  An `"uncounted_bayes"` object.

- type:

  Residual type: response or Pearson.

- summary:

  Posterior summary used for fitted values.

- ...:

  Additional arguments passed to
  [`fitted.uncounted_bayes`](https://ncn-foreigners.github.io/uncounted/reference/fitted.uncounted_bayes.md).

## Value

A numeric vector of residuals.
