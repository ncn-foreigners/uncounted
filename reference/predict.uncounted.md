# Predict from an uncounted model

Returns predicted values from a fitted `uncounted` model, optionally for
new data. This method enables compatibility with the `marginaleffects`
package.

## Usage

``` r
# S3 method for class 'uncounted'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  An `uncounted` object.

- newdata:

  Optional data frame for predictions. If `NULL`, returns fitted values
  from the original data.

- type:

  Character: `"response"` (default) returns predicted counts
  \\\hat{\mu}\_i = \exp(\mathbf{z}\_i'\hat{\boldsymbol{\theta}})\\;
  `"link"` returns the linear predictor
  \\\mathbf{z}\_i'\hat{\boldsymbol{\theta}}\\.

- ...:

  Additional arguments (ignored).

## Value

A numeric vector of predicted values.

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year == "2019", ]
fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
                           gamma = 0.005)
head(predict(fit))
#> [1]  2.5691032  0.2256555  2.1005157  1.9776050 22.8677881  0.5587853
head(predict(fit, newdata = d[1:5, ]))
#> [1]  2.5691032  0.2256555  2.1005157  1.9776050 22.8677881
head(predict(fit, type = "link"))
#> [1]  0.9435569 -1.4887458  0.7421829  0.6818865  3.1297293 -0.5819900
```
