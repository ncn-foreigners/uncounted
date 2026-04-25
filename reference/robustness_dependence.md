# Dependence Robustness Value for Population Size

Summarises the smallest dependence strength needed for the profiled
population-size estimate to cross a target value.

## Usage

``` r
robustness_dependence(
  object,
  q = 0.25,
  direction = c("decrease", "increase"),
  threshold = NULL,
  delta_grid = seq(-1, 1, length.out = 81)
)

# S3 method for class 'uncounted_dependence_robustness'
print(x, ...)
```

## Arguments

- object:

  Either an `"uncounted"` object or a data frame returned by
  [`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md)
  containing at least `delta`, `kappa`, and `xi`.

- q:

  Relative change target used when `threshold = NULL`. For
  `direction = "decrease"`, the target is `baseline_xi * (1 - q)`. For
  `direction = "increase"`, the target is `baseline_xi * (1 + q)`.

- direction:

  Which target crossing to report: `"decrease"` or `"increase"`.

- threshold:

  Optional positive numeric target for \\\xi\\. When supplied, `q` is
  ignored for target construction.

- delta_grid:

  Numeric grid passed to
  [`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md)
  when `object` is an `"uncounted"` fit.

- x:

  Object to print.

- ...:

  Additional arguments passed to print methods.

## Value

An object of class `"uncounted_dependence_robustness"` with components
`baseline_xi`, `target_xi`, `rv_delta`, `rv_kappa`, `rv_Gamma`,
`direction`, `q`, `threshold`, and `reached`.

## Details

This helper adapts the robustness-value idea to the one-dimensional
dependence profile. It searches the profiled \\\hat\xi(\delta)\\ values
and returns the smallest absolute `delta` whose profile crosses the
requested target. The corresponding multiplicative scales are
`rv_kappa = exp(rv_delta)` and `rv_Gamma = exp(abs(rv_delta))`.

## See also

[`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md),
[`dependence_bounds`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md)

## Examples

``` r
set.seed(123)
d <- data.frame(
  N = round(exp(rnorm(20, mean = 7, sd = 0.35)))
)
d$n <- rpois(20, lambda = pmax(1, 0.03 * d$N))
d$m <- rpois(20, lambda = d$N^0.5 * (0.01 + d$n / d$N)^0.8)

fit <- estimate_hidden_pop(
  data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
  method = "poisson", gamma = 0.01
)

robustness_dependence(fit, q = 0.10, delta_grid = seq(-0.5, 0.5, length.out = 9))
#> Dependence robustness analysis
#> Baseline xi: 87 
#> Target xi: 78 
#> rv_delta: 0.375 
#> rv_kappa: 1.455 
#> rv_Gamma: 1.455 
```
