# Likelihood ratio test for nested uncounted models

Performs a likelihood ratio test between two nested `uncounted` models
fitted by maximum likelihood (Poisson or Negative Binomial). The test
statistic is

## Usage

``` r
lrtest(object1, object2)
```

## Arguments

- object1:

  An `"uncounted"` model (Poisson or NB).

- object2:

  An `"uncounted"` model (Poisson or NB). The function automatically
  identifies the simpler model regardless of argument order.

## Value

An object of class `"uncounted_lrtest"` with components `statistic`,
`df`, `p_value`, `boundary`, model descriptions, and log-likelihoods.

## Details

\$\$LR = 2\bigl(\ell_2 - \ell_1\bigr)\$\$

where \\\ell_1\\ and \\\ell_2\\ are the maximised log-likelihoods of the
simpler and more complex models, respectively. Under \\H_0\\ the
statistic follows a \\\chi^2\\ distribution with degrees of freedom
equal to the difference in the number of estimated parameters.

**Boundary correction for Poisson vs NB.** When comparing a Poisson
model (H0) against a Negative Binomial model (H1), the dispersion
parameter \\\theta\\ lies on the boundary of its parameter space under
\\H_0\\ (\\\theta \to \infty\\). The standard \\\chi^2\\ approximation
is invalid in this setting. Following Self & Liang (1987), the p-value
is corrected as

\$\$p = 0.5\\\Pr(\chi^2_1 \> LR)\$\$

which corresponds to a 50:50 mixture of a point mass at zero and a
\\\chi^2_1\\.

This function is intended for comparing count models only (Poisson and
NB). OLS/NLS models are not supported because they use a
pseudo-log-likelihood. The models must be fitted to the same data
(identical number of observations).

## References

Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605–610.

## Examples

``` r
set.seed(42)
df <- data.frame(
  N = rep(1000, 30),
  n = rpois(30, lambda = 50)
)
df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)

fit_po <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                              reference_pop = ~N, method = "poisson")
#> Warning: Some alpha values < 0 (min = -1.28). Consider using constrained = TRUE.
fit_nb <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                              reference_pop = ~N, method = "nb")
#> Warning: Some alpha values < 0 (min = -1.244). Consider using constrained = TRUE.

# Test for overdispersion (Poisson vs NB)
lrtest(fit_po, fit_nb)
#> Likelihood ratio test
#> ---------------------------------------- 
#> Model 1: POISSON   (logLik = -52.51 )
#> Model 2: NB   (logLik = -52.52 )
#> LR statistic: 0 on 1 df
#> (Boundary-corrected: 0.5 * P(chi2 > LR), Self & Liang 1987)
#> p-value: 0.5 
```
