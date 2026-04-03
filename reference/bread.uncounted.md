# Bread matrix for sandwich estimator

Returns the bread component B such that the model-based variance is B/n.
Specifically: B = n \* (Z' W Z)^-1 where Z is the model matrix and W =
diag(bread_weights).

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
