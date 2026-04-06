# Bread matrix for sandwich estimator

Returns the bread component B such that the model-based variance is B/n.
Specifically: \\B = n (Z' W Z)^{-1}\\ where Z is the model matrix and
\\W = \mathrm{diag}(\mathrm{bread\\weights})\\.

## Usage

``` r
# S3 method for class 'uncounted'
bread(x, ...)
```

## Arguments

- x:

  An uncounted object

- ...:

  Ignored

## Note

For NB models, these methods operate on the mean-model parameters
(alpha, beta, and optionally gamma) only, excluding theta. The stored
[`vcov()`](https://rdrr.io/r/stats/vcov.html) on NB objects uses a
dedicated theta-aware path instead. Calling
[`sandwich::vcovHC()`](https://rdrr.io/pkg/sandwich/man/vcovHC.html)
directly on an NB object gives theta-conditional standard errors, which
differ from [`vcov()`](https://rdrr.io/r/stats/vcov.html).
