# Hypothesis Tests for Population Size

## Goal

The
[`hypotheses_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md)
function tests statements about estimated population sizes, denoted here
by $\xi$. It is designed for questions such as:

- Is the estimated population size in 2024 larger than 200,000?
- Did the estimated population size increase between 2019 and 2024?
- Is a custom contrast of several grouped population-size estimates
  equal to zero?

The function uses an explicit quantity-indexed syntax:

``` r
xi[year == 2024]
xi[year == 2024] - xi[year == 2019]
```

This avoids ambiguous bare labels and leaves room for future extensions
to other quantities, such as `alpha[...]`.

``` r
library(uncounted)
data(irregular_migration)
```

## Fit A Model

For a compact example, fit a Poisson model with a fixed gamma offset and
then compute population size by year.

``` r
fit <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  gamma = 0.005,
  countries = ~ country
)

ps_year <- popsize(fit, by = ~ year)
ps_year
#>   group observed estimate estimate_bc    lower     upper share_pct
#> 1  2019     6604 40906.90    40897.65 19160.12  87296.86  11.89222
#> 2  2020     3398 44720.59    44710.31 20803.07  96092.18  13.00091
#> 3  2021     3105 52686.31    52673.92 24277.41 114284.94  15.31666
#> 4  2022     2949 65068.18    65051.86 29171.43 145064.68  18.91625
#> 5  2023     4362 69246.09    69229.04 31266.13 153285.99  20.13083
#> 6  2024     6687 71352.31    71335.30 32617.45 156012.32  20.74313
```

The examples below use the bias-corrected estimate explicitly. The
`estimate` argument is required so that the tested estimand is never
chosen silently.

By default, directional character hypotheses are interpreted as the
alternative hypothesis. For example, `"xi[year == 2024] > 200000"` tests

$$H_{0}:\xi_{2024} \leq 200000\quad\text{against}\quad H_{1}:\xi_{2024} > 200000.$$

The returned table includes explicit `null_hypothesis` and
`alternative_hypothesis` columns. If you want the character expression
to be the null hypothesis instead, use `hypothesis_side = "null"`.

## Threshold Tests

To test whether the 2024 population size is larger than 200,000:

``` r
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] > 200000",
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>                 hypothesis            null_hypothesis    alternative_hypothesis
#>  xi[year == 2024] > 200000 xi[year == 2024] <= 200000 xi[year == 2024] > 200000
#>  estimate null.value  contrast std.error statistic   p.value      s.value
#>   71335.3      2e+05 -128664.7  28481.64 -4.517461 0.9999969 4.514601e-06
#>   conf.low conf.high alternative                 method estimate_type n_draws
#>  -184487.7 -72841.71     greater Delta method Wald test   estimate_bc      NA
```

This is a one-sided Wald test. The contrast is

$${\widehat{\xi}}_{2024} - 200000.$$

A positive statistic supports the alternative that $\xi_{2024}$ is
larger than the threshold.

To treat the same directional statement as the null hypothesis, write:

``` r
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] > 200000",
  estimate = "estimate_bc",
  hypothesis_side = "null"
)
#> Population-size hypothesis tests
#>                 hypothesis            null_hypothesis    alternative_hypothesis
#>  xi[year == 2024] > 200000 xi[year == 2024] >= 200000 xi[year == 2024] < 200000
#>  estimate null.value  contrast std.error statistic      p.value  s.value
#>   71335.3      2e+05 -128664.7  28481.64 -4.517461 3.129278e-06 18.28574
#>   conf.low conf.high alternative                 method estimate_type n_draws
#>  -184487.7 -72841.71        less Delta method Wald test   estimate_bc      NA
```

This asks whether there is evidence against
$H_{0}:\xi_{2024} \geq 200000$ in favor of $H_{1}:\xi_{2024} < 200000$.

The compact form also works when there is exactly one grouping variable:

``` r
hypotheses_popsize(
  ps_year,
  "xi[2024] > 200000",
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>         hypothesis    null_hypothesis alternative_hypothesis estimate
#>  xi[2024] > 200000 xi[2024] <= 200000      xi[2024] > 200000  71335.3
#>  null.value  contrast std.error statistic   p.value      s.value  conf.low
#>       2e+05 -128664.7  28481.64 -4.517461 0.9999969 4.514601e-06 -184487.7
#>  conf.high alternative                 method estimate_type n_draws
#>  -72841.71     greater Delta method Wald test   estimate_bc      NA
```

## Change Over Time

To test whether the estimated population size increased between 2019 and
2024:

``` r
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] - xi[year == 2019] > 0",
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>                               hypothesis
#>  xi[year == 2024] - xi[year == 2019] > 0
#>                           null_hypothesis
#>  xi[year == 2024] - xi[year == 2019] <= 0
#>                   alternative_hypothesis estimate null.value contrast std.error
#>  xi[year == 2024] - xi[year == 2019] > 0 30437.64          0 30437.64  12659.76
#>  statistic     p.value  s.value conf.low conf.high alternative
#>   2.404282 0.008102137 6.947482 5624.959  55250.32     greater
#>                  method estimate_type n_draws
#>  Delta method Wald test   estimate_bc      NA
```

The standard error uses the full covariance matrix of the two estimates.
This matters because year-specific estimates are functions of the same
fitted model parameters and are generally dependent.

Two-sided equality tests use `=` or `==`:

``` r
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] = xi[year == 2019]",
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>                           hypothesis                     null_hypothesis
#>  xi[year == 2024] = xi[year == 2019] xi[year == 2024] = xi[year == 2019]
#>                alternative_hypothesis estimate null.value contrast std.error
#>  xi[year == 2024] != xi[year == 2019]  71335.3   40897.65 30437.64  12659.76
#>  statistic    p.value  s.value conf.low conf.high alternative
#>   2.404282 0.01620427 5.947482 5624.959  55250.32   two.sided
#>                  method estimate_type n_draws
#>  Delta method Wald test   estimate_bc      NA
```

## Multiple Groups

Filters can match more than one row. In that case, `xi[...]` is the sum
of the matching population-size estimates, with covariance handled
before testing.

``` r
ps_year_sex <- popsize(fit, by = ~ year + sex)

hypotheses_popsize(
  ps_year_sex,
  "xi[year == 2024] - xi[year == 2019] > 0",
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>                               hypothesis
#>  xi[year == 2024] - xi[year == 2019] > 0
#>                           null_hypothesis
#>  xi[year == 2024] - xi[year == 2019] <= 0
#>                   alternative_hypothesis estimate null.value contrast std.error
#>  xi[year == 2024] - xi[year == 2019] > 0 30437.64          0 30437.64  12659.76
#>  statistic     p.value  s.value conf.low conf.high alternative
#>   2.404282 0.008102128 6.947483 5624.964  55250.32     greater
#>                  method estimate_type n_draws
#>  Delta method Wald test   estimate_bc      NA
```

Custom functions are useful when the contrast is easier to write in R
code. The function receives a named vector of xi estimates and,
optionally, the group metadata.

``` r
hypotheses_popsize(
  ps_year,
  function(x, groups) {
    c(change_2024_2019 = unname(x["2024"] - x["2019"]))
  },
  estimate = "estimate_bc"
)
#> Population-size hypothesis tests
#>        hypothesis      null_hypothesis alternative_hypothesis estimate
#>  change_2024_2019 change_2024_2019 = 0  change_2024_2019 != 0 30437.64
#>  null.value contrast std.error statistic    p.value  s.value conf.low conf.high
#>           0 30437.64  12659.76  2.404282 0.01620427 5.947482 5624.959  55250.32
#>  alternative                 method estimate_type n_draws
#>    two.sided Delta method Wald test   estimate_bc      NA
```

## Statistical Theory

For analytical
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
results,
[`hypotheses_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md)
uses Wald inference for nonlinear functions of fitted model parameters.

Let $\theta$ denote the fitted alpha coefficients. Under standard
large-sample regularity conditions,

$$\widehat{\theta} \approx N\left( \theta,{\widehat{V}}_{\theta} \right).$$

Population size is a nonlinear function of those parameters:

$$\xi = h(\theta) = \sum\limits_{i}N_{i}^{\alpha_{i}}.$$

The delta method gives the approximate sampling variance:

$$\widehat{Var}\left( \widehat{\xi} \right) = \nabla h\left( \widehat{\theta} \right)\prime{\widehat{V}}_{\theta}\nabla h\left( \widehat{\theta} \right).$$

For a threshold hypothesis such as

$$H_{0}:\xi_{2024} \leq c\quad\text{against}\quad H_{1}:\xi_{2024} > c,$$

the reported statistic is

$$z = \frac{{\widehat{\xi}}_{2024} - c}{SE\left( {\widehat{\xi}}_{2024} \right)}.$$

For an increase hypothesis,

$$H_{0}:\xi_{2024} - \xi_{2019} \leq 0\quad\text{against}\quad H_{1}:\xi_{2024} - \xi_{2019} > 0,$$

the statistic is

$$z = \frac{{\widehat{\xi}}_{2024} - {\widehat{\xi}}_{2019}}{SE\left( {\widehat{\xi}}_{2024} - {\widehat{\xi}}_{2019} \right)}.$$

The default reference distribution is standard normal. If you supply a
finite `df`, the function uses a t distribution instead.

## Bootstrap Tests

For objects returned by
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md),
the same contrast is evaluated inside every fractional weighted
bootstrap draw. For example, for an increase test:

$$d^{*} = \xi_{2024}^{*} - \xi_{2019}^{*}.$$

The bootstrap standard error is `sd(d*)`, and the function reports a
Wald statistic using that standard error.

``` r
boot_year <- bootstrap_popsize(
  fit,
  R = 19,
  cluster = ~ country_code,
  seed = 2025,
  by = ~ year,
  verbose = FALSE
)

hypotheses_popsize(
  boot_year,
  "xi[year == 2024] - xi[year == 2019] > 0",
  estimate = "boot_median"
)
#> Population-size hypothesis tests
#>                               hypothesis
#>  xi[year == 2024] - xi[year == 2019] > 0
#>                           null_hypothesis
#>  xi[year == 2024] - xi[year == 2019] <= 0
#>                   alternative_hypothesis estimate null.value contrast std.error
#>  xi[year == 2024] - xi[year == 2019] > 0 35586.71          0 35586.71  17107.26
#>  statistic   p.value  s.value conf.low conf.high alternative        method
#>   2.080211 0.0187531 5.736727 2057.093  69116.33     greater FWB Wald test
#>  estimate_type n_draws
#>    boot_median      19
```

This Wald p value is different from an empirical bootstrap exceedance
probability such as `mean(xi_star > 200000)`. Use
[`exceedance_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/exceedance_popsize.md)
for those direct bootstrap tail summaries.

## Power And Limitations

The power of these tests depends on the size of the contrast, the
standard error, the number and quality of observations, clustering or
panel dependence, model specification, and bootstrap stability. A small
p value is easier to obtain when the estimated contrast is large
relative to its standard error.

The helper reports the estimate, contrast, standard error, Wald
statistic, p value, S value, and confidence interval. It does not
guarantee that a scientifically important change will be detected with
high probability.

Formal power analysis should be simulation-based: choose a
data-generating process, simulate many datasets under the alternative of
interest, refit the model, apply the hypothesis test, and estimate the
rejection rate.
