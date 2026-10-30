# Cluster-robust covariance for uncounted models

Cluster-robust covariance for uncounted models

## Usage

``` r
# S3 method for class 'uncounted'
vcovCL(x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, ...)
```

## Arguments

- x:

  An uncounted object

- cluster:

  Clustering variable

- type:

  HC type

- sandwich:

  Return the full sandwich when `TRUE`

- fix:

  Passed to
  [`sandwich::meatCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)

- ...:

  Ignored
