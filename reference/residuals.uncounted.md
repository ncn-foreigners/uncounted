# Residuals for uncounted models

Extracts residuals of the specified type from a fitted `uncounted`
model.

## Usage

``` r
# S3 method for class 'uncounted'
residuals(
  object,
  type = c("response", "pearson", "deviance", "anscombe", "working"),
  ...
)
```

## Arguments

- object:

  An object of class `"uncounted"`.

- type:

  Type of residual: `"response"` (default), `"pearson"`, `"deviance"`,
  or `"anscombe"`.

- ...:

  Additional arguments (ignored).

## Value

A numeric vector of residuals with the same length as the number of
observations.

## Details

Four residual types are supported:

- **Response residuals**:

  \$\$r_i = m_i - \hat\mu_i\$\$

- **Pearson residuals**:

  \$\$r_i = \frac{m_i - \hat\mu_i}{\sqrt{V(\hat\mu_i)}}\$\$ where the
  variance function \\V(\mu)\\ depends on the model:

  - Poisson: \\V(\mu) = \mu\\

  - Negative Binomial: \\V(\mu) = \mu + \mu^2/\theta\\

  - OLS/NLS: \\V(\mu) = \mu^2 \hat\sigma^2\\ (delta-method approximation
    on the log scale)

- **Deviance residuals**:

  Signed square root of the individual deviance contributions: \$\$r_i =
  \mathrm{sign}(m_i - \hat\mu_i)\\\sqrt{d_i}\$\$ where \\d_i\\ is the
  \\i\\-th term in the deviance sum (see
  [`deviance.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/deviance.uncounted.md)).
  For OLS/NLS models, standardised log-scale residuals are returned
  instead.

- **Anscombe residuals**:

  Variance-stabilising residuals.

      For Poisson (McCullagh & Nelder, 1989):
      \deqn{r_i = \frac{3\bigl(m_i^{2/3} - \hat\mu_i^{2/3}\bigr)}
            {2\,\hat\mu_i^{1/6}}}

      For Negative Binomial (Hilbe, 2011):
      \deqn{r_i = \frac{3\bigl(m_i^{2/3} - \hat\mu_i^{2/3}\bigr)}
            {2\,\hat\mu_i^{1/6}\,(1 + \hat\mu_i/\theta)^{1/6}}}

      For OLS/NLS: standardised log-scale residuals.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman & Hall.
[doi:10.1007/978-1-4899-3242-6](https://doi.org/10.1007/978-1-4899-3242-6)

Hilbe, J. M. (2011). *Negative Binomial Regression*, 2nd ed. Cambridge
University Press.

Pierce, D. A. and Schafer, D. W. (1986). Residuals in generalized linear
models. *Journal of the American Statistical Association*, 81(396),
977–986.

## Examples

``` r
set.seed(123)
df <- data.frame(
  N = rep(1000, 30),
  n = rpois(30, lambda = 50)
)
df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)

fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                           reference_pop = ~N, method = "poisson")
#> Warning: Some alpha values > 1 (max = 3.258). Population size estimates may be unreliable. Consider using constrained = TRUE or simplifying cov_alpha.
#> Estimated gamma is at its upper bound (0.5); the optimum may lie outside gamma_bounds = c(1e-10, 0.5).

# Response residuals (default)
head(residuals(fit))
#> [1] -1.9390196 -1.2498869 -0.1381474 -0.5235144 -1.4998266 -1.0709919

# Pearson residuals
head(residuals(fit, type = "pearson"))
#> [1] -1.3924869 -0.6062923 -0.1294922 -0.3295533 -0.6395383 -0.6111486

# Deviance and Anscombe residuals
head(residuals(fit, type = "deviance"))
#> [1] -1.9692738 -0.6403853 -0.1322539 -0.3420634 -0.6725144 -0.6531233
head(residuals(fit, type = "anscombe"))
#> [1] -2.0887303 -0.6407280 -0.1322640 -0.3421464 -0.6728165 -0.6536473
```
