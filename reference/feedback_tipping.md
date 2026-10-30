# Tipping Points under Detection Feedback

Finds how much detection feedback it takes to overturn a statement about
the population size.

## Usage

``` r
feedback_tipping(
  x,
  threshold = NULL,
  group = NULL,
  hypothesis = NULL,
  level = 0.95
)

# S3 method for class 'uncounted_feedback_tipping'
print(x, ...)
```

## Arguments

- x:

  An `"uncounted_feedback"` object.

- threshold:

  A single value for a group total. Use with `group`. The tested
  direction is the side of the threshold on which the estimate lies at
  \\\delta = 0\\.

- group:

  Group label for `threshold`. Defaults to `"Total"` when available, or
  the only group.

- hypothesis:

  A comparison in the `xi[...]` syntax of
  [`hypotheses_popsize`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md),
  such as `"xi[year == 2024] - xi[year == 2019] > 0"`. The filters refer
  to the `by` variables of `x`.

- level:

  Level of the one-sided bootstrap test: the statement is supported at
  \\\delta\\ when the `1 - level` quantile of the replicate contrasts is
  above 0 (for `>`), or the `level` quantile is below 0 (for `<`). For
  `=`, both tails use `(1 - level) / 2`.

- ...:

  Ignored.

## Value

An object of class `"uncounted_feedback_tipping"`, a list with
`statement`, `alternative`, `baseline` (the contrast at \\\delta = 0\\),
`point` and `test` (each `c(negative, positive)`, `NA` when not
reached), `supported_at_zero` (`NA` without a bootstrap),
`admissible_range`, and `curve`, a data frame with `delta`, `contrast`,
`lower`, `upper` and `supported`.

## Details

Two kinds of tipping point are reported on each side of 0, within the
`delta_grid` of `x`:

- `point`: the \\\delta\\ nearest 0 at which the point estimate of the
  contrast changes sign, found exactly between grid points.

- `test`: with a bootstrap, the grid value nearest 0 at which the test
  stops supporting the statement.

If the test does not support the statement at \\\delta = 0\\, there is
nothing to overturn: `test` is `NA` and `supported_at_zero` is `FALSE`.
Tipping points are reported whether or not they are admissible, with the
admissible range shown alongside.

## See also

[`feedback_sensitivity`](https://ncn-foreigners.github.io/uncounted/reference/feedback_sensitivity.md),
[`hypotheses_popsize`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md)
