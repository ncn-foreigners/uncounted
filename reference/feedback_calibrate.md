# Calibrate the Feedback Parameter

Summarizes \\s_i = \xi\_{0i} / (N_i + \xi\_{0i})\\, the share of the
unauthorized population in the foreign population of a cell. If
inspections are spread over the whole foreign population, the implied
mechanical congestion elasticity of detection is \\-s\\, which helps to
choose a `delta_grid` for
[`feedback_sensitivity`](https://ncn-foreigners.github.io/uncounted/reference/feedback_sensitivity.md).

## Usage

``` r
feedback_calibrate(fit)

# S3 method for class 'uncounted_feedback_calibration'
print(x, ...)
```

## Arguments

- fit:

  An `"uncounted"` object.

- x:

  Object to print.

- ...:

  Ignored.

## Value

An object of class `"uncounted_feedback_calibration"`, a list with `s`
(one value per cell), `mean_s`, `weighted_mean_s` (weighted by
\\\xi\_{0i}\\) and `congestion_elasticity`, the negative of
`weighted_mean_s`.

## See also

[`feedback_sensitivity`](https://ncn-foreigners.github.io/uncounted/reference/feedback_sensitivity.md)
