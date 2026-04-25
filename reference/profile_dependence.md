# Profile Dependence Sensitivity for Population Size

Refits the model over a grid of fixed dependence offsets and records how
the plug-in total population-size estimate changes.

## Usage

``` r
profile_dependence(
  object,
  delta_grid = seq(-1, 1, length.out = 41),
  plot = TRUE,
  ...
)
```

## Arguments

- object:

  An `"uncounted"` object fitted with `estimator = "mle"` and
  `method = "poisson"` or `method = "nb"`.

- delta_grid:

  Numeric vector of dependence offsets to evaluate. The baseline
  no-dependence case is `delta = 0`.

- plot:

  Logical; produce the two-panel profile plot? Default `TRUE`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a data frame with columns:

- `delta`:

  Fixed dependence offset used in the refit.

- `kappa`:

  Multiplicative distortion `exp(delta)`.

- `xi`:

  Plug-in total population-size estimate `sum(popsize(fit_i)$estimate)`.

- `loglik`:

  Refitted log-likelihood, or `NA` when the refit failed.

## Details

The bounded helper
[`dependence_bounds`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md)
keeps the baseline point estimate fixed and only widens the
identification envelope. In contrast, `profile_dependence()` imposes a
parametric dependence model and refits the mean function over a fixed
grid:

\$\$\mu_i(\delta) = \exp(\delta)\\\xi_i\rho_i,\$\$

or equivalently

\$\$\log \mu_i(\delta) = \delta + \alpha_i \log N_i + \log(\rho_i).\$\$

Because the model is refitted for each `delta`, the reported
population-size estimate \\\hat\xi(\delta)\\ can move with the
sensitivity parameter. This is a model-based sensitivity analysis, not a
partial-identification bound.

## See also

[`dependence_bounds`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md),
[`robustness_dependence`](https://ncn-foreigners.github.io/uncounted/reference/robustness_dependence.md)

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

profile_dependence(fit, delta_grid = c(-0.25, 0, 0.25), plot = FALSE)
```
