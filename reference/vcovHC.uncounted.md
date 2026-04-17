# Heteroskedasticity-consistent covariance for uncounted models

Heteroskedasticity-consistent covariance for uncounted models

## Usage

``` r
# S3 method for class 'uncounted'
vcovHC(
  x,
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL,
  ...
)
```

## Arguments

- x:

  An uncounted object

- type:

  HC type

- omega:

  Optional omega function passed to
  [`sandwich::meatHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)

- ...:

  Ignored
