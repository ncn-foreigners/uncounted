# Print LOO Comparison Results

Displays the full-model \\\hat{\xi}\\ for each model and a table of the
top `n` most influential observations/countries ranked by the maximum
absolute percentage change across both models.

## Usage

``` r
# S3 method for class 'uncounted_loo_compare'
print(x, n = 15, ...)
```

## Arguments

- x:

  An `"uncounted_loo_compare"` object.

- n:

  Number of top observations to display (default 15).

- ...:

  Additional arguments (ignored).
