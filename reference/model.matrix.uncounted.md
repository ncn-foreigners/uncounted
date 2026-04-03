# Extract model matrix

Returns the Jacobian/design matrix used in sandwich computation. When
gamma was estimated, the returned matrix includes the gamma column.

## Usage

``` r
# S3 method for class 'uncounted'
model.matrix(object, ...)
```

## Arguments

- object:

  An uncounted object

- ...:

  Ignored
