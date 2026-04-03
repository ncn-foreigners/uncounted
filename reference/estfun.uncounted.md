# Empirical estimating functions (score contributions)

Returns an n x p matrix of per-observation score contributions for use
with the sandwich package. Row i is the gradient of the
(quasi-)log-likelihood contribution of observation i with respect to the
parameter vector.

## Usage

``` r
# S3 method for class 'uncounted'
estfun(x, ...)
```

## Arguments

- x:

  An uncounted object

- ...:

  Ignored
