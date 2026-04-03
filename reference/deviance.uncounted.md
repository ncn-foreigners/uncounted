# Deviance for uncounted models

Computes the deviance (twice the difference between the saturated and
fitted log-likelihoods) for a fitted `uncounted` model.

## Usage

``` r
# S3 method for class 'uncounted'
deviance(object, ...)
```

## Arguments

- object:

  An object of class `"uncounted"`.

- ...:

  Additional arguments (ignored).

## Value

A numeric scalar giving the deviance.

## Details

The formula depends on the estimation method:

**Poisson:** \$\$D = 2\sum\_{i}\bigl\[m_i\log(m_i/\hat\mu_i) - (m_i -
\hat\mu_i)\bigr\]\$\$ with the convention \\0 \log 0 = 0\\.

**Negative Binomial:** \$\$D = 2\sum\_{i}\bigl\[m_i\log(m_i/\hat\mu_i) -
(m_i+\theta)\log\bigl((m_i+\theta)/(\hat\mu_i+\theta)\bigr)\bigr\]\$\$

**OLS / NLS:** \$\$D = \sum\_{i}(\log m_i - \log\hat\mu_i)^2\$\$ i.e.,
the residual sum of squares on the log scale.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman & Hall.

## Examples

``` r
set.seed(123)
df <- data.frame(
  N = rep(1000, 30),
  n = rpois(30, lambda = 50)
)
df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)

fit_po <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                              reference_pop = ~N, method = "poisson")
#> Warning: Some alpha values > 1 (max = 3.271). Population size estimates may be unreliable. Consider using constrained = TRUE or simplifying cov_alpha.
deviance(fit_po)
#> [1] 17.61431
```
