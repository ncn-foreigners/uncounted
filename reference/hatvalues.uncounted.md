# Hat values (leverage)

Computes diagonal of the hat matrix \\H = W^{1/2} Z (Z'WZ)^{-1} Z'
W^{1/2}\\ where Z is the model matrix and W is the diagonal weight
matrix.

## Usage

``` r
# S3 method for class 'uncounted'
hatvalues(model, ...)
```

## Arguments

- model:

  An uncounted object

- ...:

  Ignored
