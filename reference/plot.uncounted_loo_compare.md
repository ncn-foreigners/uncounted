# Plot LOO Comparison: Scatter and Barplot

Produces diagnostic plots comparing observation influence across two
models.

## Usage

``` r
# S3 method for class 'uncounted_loo_compare'
plot(x, type = c("scatter", "bar"), n = 20, label_top = 10, ...)
```

## Arguments

- x:

  An `"uncounted_loo_compare"` object.

- type:

  `"scatter"` (default) or `"bar"`.

- n:

  Number of top observations to show in barplot (default 20).

- label_top:

  Number of points to label in scatter (default 10).

- ...:

  Additional arguments passed to `plot`.

## Details

**Scatter plot (`type = "scatter"`).** Each point is one
observation/country. The x-axis is \\\\\Delta\xi\\ for Model 1, the
y-axis for Model 2. The diagonal line marks equal influence. Points near
the origin are non-influential under both models. Points far from the
diagonal have divergent influence – they affect the two models
differently. The top `label_top` most influential points are labeled.
See
[`compare_loo`](https://ncn-foreigners.github.io/uncounted/reference/compare_loo.md)
for quadrant interpretation.

**Bar plot (`type = "bar"`).** Horizontal side-by-side bar chart of the
top `n` most influential observations (by maximum absolute percentage
change across both models). Red bars are Model 1, blue bars are Model 2.
Useful for quickly seeing which observations matter most and whether the
models agree on direction.
