# Extract the brmsfit from a Bayesian uncounted Model

Extract the `brmsfit` stored inside an `"uncounted_bayes"` object.

## Usage

``` r
as_brmsfit(object, ...)

# S3 method for class 'uncounted_bayes'
as_brmsfit(object, ...)
```

## Arguments

- object:

  An object.

- ...:

  Additional arguments, currently ignored.

## Value

The underlying `brmsfit`.
