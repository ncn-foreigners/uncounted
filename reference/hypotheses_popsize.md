# Hypothesis Tests for Population Size Estimates

Tests scalar hypotheses and contrasts involving population size
estimates returned by
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
or
[`bootstrap_popsize`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md).

## Usage

``` r
hypotheses_popsize(
  object,
  hypothesis,
  estimate,
  level = 0.95,
  df = Inf,
  hypothesis_side = c("alternative", "null"),
  include_total = FALSE,
  rope = NULL,
  ...
)
```

## Arguments

- object:

  An `"uncounted_popsize"` object from
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md),
  an `"uncounted_popsize_bayes"` object from
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md),
  or an `"uncounted_boot"` object from
  [`bootstrap_popsize`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md).

- hypothesis:

  Numeric null value(s), a character hypothesis, or a function.
  Character hypotheses can refer to population size estimates with
  `xi[...]` filters, such as `"xi[year == 2024] > 200000"` or
  `"xi[year == 2024] - xi[year == 2019] > 0"`. Positional aliases `b1`,
  `b2`, ... are also available.

- estimate:

  Character string naming the estimate to test. For frequentist
  [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  objects, use `"estimate"` or `"estimate_bc"`. For Bayesian
  [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  objects, use `"estimate"`, `"median"`, or `"mean"`; when omitted, the
  estimate type stored by `popsize.uncounted_bayes()` is used. For
  [`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
  objects, use `"plugin"`, `"boot_median"`, or `"boot_mean"`.

- level:

  Confidence level for Wald intervals on the tested contrast.

- df:

  Degrees of freedom for p values and confidence intervals. The default
  `Inf` uses the standard normal distribution.

- hypothesis_side:

  Character string. `"alternative"` (default) treats directional
  character hypotheses such as `"xi[...] < 15000"` as the alternative
  hypothesis to support. `"null"` treats the directional character
  hypothesis as the null hypothesis to test against. Equality hypotheses
  and numeric/function hypotheses are always interpreted as
  null-equality tests.

- include_total:

  Logical; if `TRUE`, append a `"Total"` xi target when multiple groups
  are available.

- rope:

  Optional numeric region of practical equivalence for Bayesian
  contrasts. Supply a positive scalar for `[-rope, rope]` or a
  two-element vector `c(lower, upper)`.

- ...:

  Additional arguments, currently ignored.

## Value

A data frame with one row per hypothesis and columns `hypothesis`,
`null_hypothesis`, `alternative_hypothesis`, `estimate`, `null.value`,
`contrast`, `std.error`, `statistic`, `p.value`, `s.value`, `conf.low`,
`conf.high`, `alternative`, `method`, `estimate_type`, and `n_draws`.
Bayesian results also include `p_greater`, `p_less`, `p_h1`, `p_h0`,
`posterior_odds`, `evidence_ratio`, and optional ROPE columns.

## Details

For analytical
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
results, tests use Wald inference with delta-method standard errors for
nonlinear functions of the fitted alpha coefficients. For bootstrap
results, the same contrast is evaluated in each fractional weighted
bootstrap draw, and the standard error is the bootstrap standard
deviation of those contrasts. The p value is a Wald p value, not an
empirical exceedance probability; use
[`exceedance_popsize`](https://ncn-foreigners.github.io/uncounted/reference/exceedance_popsize.md)
for bootstrap tail areas.

For Bayesian
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
results, character hypotheses are interpreted as posterior propositions.
A directional hypothesis such as `"xi[...] < 15000"` reports
`p_h1 = Pr(H1 | data)` and `p_h0 = Pr(H0 | data)` from the posterior
contrast draws. Equality hypotheses are summarized by posterior contrast
intervals; exact point-null probabilities are not assigned.

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[irregular_migration$year %in% c("2019", "2024"), ]
fit <- estimate_hidden_pop(
  data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
  method = "poisson", gamma = 0.005
)
ps_year <- popsize(fit, by = ~year)
hypotheses_popsize(ps_year, "xi[year == 2024] > 200000",
                   estimate = "estimate_bc")
#> Population-size hypothesis tests
#>                 hypothesis            null_hypothesis    alternative_hypothesis
#>  xi[year == 2024] > 200000 xi[year == 2024] <= 200000 xi[year == 2024] > 200000
#>  estimate null.value  contrast std.error statistic   p.value    s.value
#>  64014.08      2e+05 -135985.9  61634.58 -2.206325 0.9863194 0.01987322
#>   conf.low conf.high alternative                 method estimate_type n_draws
#>  -256787.5 -15184.36     greater Delta method Wald test   estimate_bc      NA
hypotheses_popsize(ps_year, "xi[year == 2024] - xi[year == 2019] > 0",
                   estimate = "estimate_bc")
#> Population-size hypothesis tests
#>                               hypothesis
#>  xi[year == 2024] - xi[year == 2019] > 0
#>                           null_hypothesis
#>  xi[year == 2024] - xi[year == 2019] <= 0
#>                   alternative_hypothesis estimate null.value contrast std.error
#>  xi[year == 2024] - xi[year == 2019] > 0 27188.37          0 27188.37  27307.62
#>  statistic   p.value  s.value  conf.low conf.high alternative
#>  0.9956331 0.1597142 2.646435 -26333.59  80710.33     greater
#>                  method estimate_type n_draws
#>  Delta method Wald test   estimate_bc      NA
```
