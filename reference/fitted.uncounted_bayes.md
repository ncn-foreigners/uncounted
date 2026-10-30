# Fitted Values from a Bayesian uncounted Model

Posterior expected observed counts for the data used to fit the model.

## Usage

``` r
# S3 method for class 'uncounted_bayes'
fitted(object, summary = c("mean", "median", "draws"), ...)
```

## Arguments

- object:

  An `"uncounted_bayes"` object.

- summary:

  Posterior summary to return: `"mean"` (the default), `"median"`, or
  `"draws"` for the full draw-by-observation matrix.

- ...:

  Additional arguments passed to the underlying brms expected-prediction
  function.

## Value

A numeric vector of posterior summaries, or a draw-by-observation matrix
when `summary = "draws"`.
