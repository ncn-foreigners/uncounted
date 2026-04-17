# Print Bootstrap Population Size Results

Prints a formatted table showing all point estimate types and the
bootstrap confidence interval. The output columns are:

- Plugin:

  Plug-in estimate \\\hat{\xi} = \sum N^{\hat{\alpha}}\\.

- Plugin (BC):

  Analytical bias-corrected plug-in estimate from
  [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

- Boot median:

  Median of the bootstrap distribution (recommended).

- Boot mean:

  Mean of the bootstrap distribution.

- CI lower / CI upper:

  Bootstrap confidence interval bounds.

A Total row is added when there are multiple groups. The header
indicates which point estimate type and CI method were selected.

## Usage

``` r
# S3 method for class 'uncounted_boot'
print(x, ...)
```

## Arguments

- x:

  An `"uncounted_boot"` object.

- ...:

  Additional arguments (ignored).
