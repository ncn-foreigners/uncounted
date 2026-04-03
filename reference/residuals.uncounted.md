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
#> Warning: Some alpha values > 1 (max = 3.271). Population size estimates may be unreliable. Consider using constrained = TRUE or simplifying cov_alpha.

# Response residuals (default)
head(residuals(fit))
#> [1] -1.9426446 -1.2439567 -0.1250438 -0.5329768 -1.4553729 -1.0807200

# Pearson residuals
head(residuals(fit, type = "pearson"))
#> [1] -1.3937879 -0.6038371 -0.1178902 -0.3348827 -0.6231062 -0.6157254

# Deviance and Anscombe residuals
head(residuals(fit, type = "deviance"))
#> [1] -1.9711137 -0.6376606 -0.1201818 -0.3477951 -0.6544409 -0.6583125
head(residuals(fit, type = "anscombe"))
#> [1] -2.0906818 -0.6379992 -0.1201894 -0.3478822 -0.6547205 -0.6588482
```
