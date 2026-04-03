# Plot LOO Sensitivity Results

Produces diagnostic plots for the leave-one-out analysis.

- `type = "xi"`:

  Horizontal bar plot of \\\Delta\xi\\ for each dropped
  observation/country, sorted by value. Blue bars indicate positive
  change (dropping the obs increases the estimate); red bars indicate
  negative change (dropping the obs decreases the estimate).

- `type = "coef"`:

  DFBETA plots: one panel per coefficient showing \\\Delta\beta_j\\ for
  each observation. Points beyond \\\pm 2 \\ \mathrm{SD}\\ are
  highlighted in red.

## Usage

``` r
# S3 method for class 'uncounted_loo'
plot(x, type = c("xi", "coef"), ...)
```

## Arguments

- x:

  An `"uncounted_loo"` object.

- type:

  `"xi"` (default) or `"coef"`.

- ...:

  Additional arguments passed to `plot`/`barplot`.
