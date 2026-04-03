# Summary of LOO Sensitivity Results

Prints coefficient stability (full-model value, LOO mean, SD, min, max
for each coefficient) and xi stability (full-model value, LOO mean, LOO
range, maximum absolute percentage change). Useful for assessing whether
the model is driven by a few observations or is broadly stable.

## Usage

``` r
# S3 method for class 'uncounted_loo'
summary(object, ...)
```

## Arguments

- object:

  An `"uncounted_loo"` object.

- ...:

  Additional arguments (ignored).
