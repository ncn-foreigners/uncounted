# Dependence Bounds for Total Population Size

Quantifies how the estimated total hidden population changes when the
identifying separability assumption behind \\E(m_i \mid N_i, n_i) =
\xi_i \rho_i\\ is relaxed.

## Usage

``` r
dependence_bounds(
  object,
  Gamma = c(1, 1.05, 1.1, 1.25, 1.5, 2),
  bias_correction = TRUE,
  threshold = NULL,
  threshold_side = c("lower", "upper")
)

sensitivity_dependence(
  object,
  Gamma = c(1, 1.05, 1.1, 1.25, 1.5, 2),
  bias_correction = TRUE,
  threshold = NULL,
  threshold_side = c("lower", "upper")
)

# S3 method for class 'uncounted_dependence_bounds'
print(x, ...)

# S3 method for class 'uncounted_sensitivity'
print(x, ...)
```

## Arguments

- object:

  An `"uncounted"` object.

- Gamma:

  Numeric vector of sensitivity parameters. Values must be finite and at
  least 1. Duplicates are removed and the remaining values are sorted
  increasingly before the bounds are computed.

- bias_correction:

  Logical; use the bias-corrected total from
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)?
  Default `TRUE`. When `FALSE`, the plug-in total is used instead.

- threshold:

  Optional positive numeric threshold used to compute a tipping-point
  `gamma_star`. Default `NULL`.

- threshold_side:

  Which bound should be compared against `threshold`? `"lower"`
  (default) finds the first `Gamma` value with `lower <= threshold`;
  `"upper"` finds the first value with `upper >= threshold`.

- x:

  Object to print.

- ...:

  Additional arguments passed to print methods.

## Value

An object of class `"uncounted_dependence_bounds"`, a list with:

- `table`:

  Data frame with columns `Gamma`, `estimate`, `lower`, `upper`,
  `pct_change_lower`, and `pct_change_upper`.

- `baseline`:

  One-row data frame containing the baseline total columns `estimate`,
  `estimate_bc`, `se`, `lower`, and `upper`.

- `bias_correction`:

  Logical; whether the baseline used the analytical bias correction.

- `threshold`:

  The supplied threshold, or `NULL`.

- `threshold_side`:

  Which bound was compared against the threshold.

- `gamma_star`:

  `NULL` when no threshold was supplied; otherwise the first `Gamma` on
  the supplied grid that crosses the threshold, or `NA_real_` if the
  grid never crosses it.

## Details

The fitted baseline model writes the expected observed count as

\$\$\mu_i = \xi_i \rho_i,\$\$

where \\\xi_i\\ is the latent population component and \\\rho_i\\ is the
detection component. This function relaxes the separability assumption
by introducing a multiplicative distortion term:

\$\$\mu_i = \xi_i \rho_i \kappa_i.\$\$

The sensitivity model imposes the bound

\$\$1 / \Gamma \le \kappa_i \le \Gamma, \qquad \Gamma \ge 1.\$\$

When \\\Gamma = 1\\, the distortion disappears and the returned bounds
collapse to the baseline total. Larger values of \\\Gamma\\ widen the
identification envelope around the baseline estimate.

The reported bounds are a sensitivity analysis for an untestable
identifying assumption. They are not confidence intervals and should not
be merged with the sampling uncertainty returned by
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

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

dependence_bounds(fit, Gamma = c(1, 1.1, 1.25))
#> Dependence bounds analysis
#> Baseline total: 4 
#> Interpretation: bounds reflect identification sensitivity, not sampling uncertainty.
#> 
#>  Gamma estimate    lower    upper pct_change_lower pct_change_upper
#>   1.00 4.269743 4.269743 4.269743         0.000000                0
#>   1.10 4.269743 3.881585 4.696718        -9.090909               10
#>   1.25 4.269743 3.415795 5.337179       -20.000000               25
dependence_bounds(
  fit,
  Gamma = c(1, 1.1, 1.25),
  threshold = 900,
  threshold_side = "lower"
)
#> Dependence bounds analysis
#> Baseline total: 4 
#> Interpretation: bounds reflect identification sensitivity, not sampling uncertainty.
#> 
#>  Gamma estimate    lower    upper pct_change_lower pct_change_upper
#>   1.00 4.269743 4.269743 4.269743         0.000000                0
#>   1.10 4.269743 3.881585 4.696718        -9.090909               10
#>   1.25 4.269743 3.415795 5.337179       -20.000000               25
#> 
#> Threshold (lower bound): 900
#> gamma_star: 1
```
