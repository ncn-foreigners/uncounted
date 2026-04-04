# Compare Population Size Across Models

Produces a side-by-side comparison of population size estimates from
multiple fitted models.

## Usage

``` r
compare_popsize(..., labels = NULL, by = NULL, bias_correction = TRUE)
```

## Arguments

- ...:

  Fitted `uncounted` objects to compare.

- labels:

  Character vector of model labels. Defaults to model method.

- by:

  Optional formula for grouping (passed to
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)).

- bias_correction:

  Logical; apply bias correction? Default TRUE.

## Value

An object of class `"uncounted_popsize_compare"` with components `table`
(data frame) and `labels` (model labels).
