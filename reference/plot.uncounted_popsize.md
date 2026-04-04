# Plot Population Size Estimates

Visualizes population size estimates by group with confidence intervals.

## Usage

``` r
# S3 method for class 'uncounted_popsize'
plot(x, type = c("estimate", "compare"), ...)
```

## Arguments

- x:

  An `"uncounted_popsize"` object from
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

- type:

  `"estimate"` (default) shows bias-corrected estimates with CIs;
  `"compare"` shows plug-in and bias-corrected side by side.

- ...:

  Additional graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).
