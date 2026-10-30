# Log-Likelihood for uncounted models

Extracts the log-likelihood from a fitted `uncounted` model.

## Usage

``` r
# S3 method for class 'uncounted'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"uncounted"`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `"logLik"` with attributes `df` (number of estimated
parameters) and `nobs` (number of observations).

## Details

For Poisson and Negative Binomial models, the true log-likelihood is
returned as computed during maximum-likelihood estimation. For OLS and
NLS models, a Gaussian pseudo-log-likelihood is returned, on the scale
each method optimises:

\$\$\ell = -\frac{n}{2}\bigl(\log(2\pi\hat\sigma^2) + 1\bigr)\$\$

where \\\hat\sigma^2 = n^{-1}\sum(\log m_i - \log\hat\mu_i)^2\\ for OLS
and \\\hat\sigma^2 = n^{-1}\sum(m_i - \hat\mu_i)^2\\ for NLS. iOLS
reports a Gamma pseudo-log-likelihood without normalising constants.

These values, and AIC/BIC based on them, are comparable only among
models of the same family: OLS with OLS, NLS with NLS, iOLS with iOLS,
and Poisson with NB. They live on different scales across families. See
[`compare_models`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md).
With `weights`, every family reports the weighted sum \\\sum_i w_i \log
f(m_i)\\, so weights of 2 double the value.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman & Hall.
[doi:10.1007/978-1-4899-3242-6](https://doi.org/10.1007/978-1-4899-3242-6)

## Examples

``` r
# Simulate synthetic data
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
logLik(fit)
#> 'log Lik.' -48.35317 (df=3)
AIC(fit)
#> [1] 102.7063
BIC(fit)
#> [1] 106.9099
```
