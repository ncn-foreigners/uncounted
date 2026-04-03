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
NLS models, a Gaussian pseudo-log-likelihood is returned, computed on
the log scale as

\$\$\ell = -\frac{n}{2}\bigl(\log(2\pi\hat\sigma^2) + 1\bigr)\$\$

where \\\hat\sigma^2 = n^{-1}\sum(\log m_i - \log\hat\mu_i)^2\\. This
pseudo-log-likelihood enables AIC/BIC comparison among OLS/NLS variants
but is not comparable to the true likelihood of count models.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman & Hall.

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
#> Warning: Some alpha values > 1 (max = 3.271). Population size estimates may be unreliable. Consider using constrained = TRUE or simplifying cov_alpha.
logLik(fit)
#> 'log Lik.' -48.35565 (df=3)
AIC(fit)
#> [1] 102.7113
BIC(fit)
#> [1] 106.9149
```
