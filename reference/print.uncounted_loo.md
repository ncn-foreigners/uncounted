# Print LOO Sensitivity Results

Displays a summary header (drop type, convergence count, full-model
\\\hat{\xi}\\, and the LOO range) followed by a table of the top 10 most
influential observations/countries ranked by \\\|\Delta\xi\|\\. The
table shows:

- dropped:

  The observation index or country label that was removed.

- dxi:

  \\\Delta\xi\\: change in total population size estimate.

- pct_change:

  Percentage change relative to full-model \\\hat{\xi}\\.

## Usage

``` r
# S3 method for class 'uncounted_loo'
print(x, ...)
```

## Arguments

- x:

  An `"uncounted_loo"` object.

- ...:

  Additional arguments (ignored).
