# Predict from a Bayesian uncounted Model

Predict from a Bayesian uncounted Model

## Usage

``` r
# S3 method for class 'uncounted_bayes'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "counts", "draws"),
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  An `"uncounted_bayes"` object.

- newdata:

  Optional data frame for prediction.

- type:

  Prediction scale. `"response"` returns posterior expected observed
  counts, `"link"` returns the log mean, `"counts"` returns posterior
  predictive replicated counts, and `"draws"` returns response-scale
  posterior expected-count draws.

- summary:

  Logical or character. Use `FALSE` to return draws; `TRUE` or `"mean"`
  to return posterior means; `"median"` to return posterior medians.

- ...:

  Additional arguments passed to the underlying brms posterior
  prediction function.

## Value

A numeric vector when summarized, otherwise a draw-by-observation
matrix.
