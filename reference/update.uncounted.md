# Update and re-fit an uncounted model

Update and re-fit an uncounted model

## Usage

``` r
# S3 method for class 'uncounted'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  An `uncounted` object.

- ...:

  Arguments to override in the original call (e.g., `weights = w`,
  `vcov = "HC0"`).

- evaluate:

  Logical. If `TRUE` (default), evaluate the updated call and return the
  new fit. If `FALSE`, return the unevaluated call.

## Value

An `uncounted` object (if `evaluate = TRUE`) or an unevaluated call (if
`evaluate = FALSE`).
