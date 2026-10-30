# Likelihood ratio test for nested uncounted models

Performs a likelihood ratio test between two nested `uncounted` models
fitted by maximum likelihood (Poisson or Negative Binomial). The test
statistic is

## Usage

``` r
lrtest(object1, object2, test = c("LR", "Wald"))
```

## Arguments

- object1:

  An `"uncounted"` model (Poisson or NB).

- object2:

  An `"uncounted"` model (Poisson or NB). The function automatically
  identifies the simpler model regardless of argument order.

- test:

  `"LR"` (default) for the likelihood ratio test, or `"Wald"` for a
  robust Wald test (see Details).

## Value

An object of class `"uncounted_lrtest"` with components `test`,
`statistic`, `df`, `p_value`, `boundary`, `dispersion`, model
descriptions, and log-likelihoods. For `test = "Wald"` also `vcov_type`.

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
is invalid in this setting. Following Chernoff (1954) and Self & Liang
(1987), with \\q\\ further parameters tested at interior points
(covariates in the NB model only), the p-value is

\$\$p = 0.5\\\Pr(\chi^2_q \> LR) + 0.5\\\Pr(\chi^2\_{q+1} \> LR)\$\$

where \\\chi^2_0\\ is a point mass at zero. With the same mean model
(\\q = 0\\) this is \\0.5\Pr(\chi^2_1 \> LR)\\ for \\LR \> 0\\ and 1 for
\\LR = 0\\, which happens with probability one half under \\H_0\\.

**Nesting.** A Poisson model is nested in an NB model only if its alpha,
beta and gamma specifications are contained in the NB model's. An NB
model is never nested in a Poisson model, however many covariates the
Poisson model has, so such pairs are an error: use
[`compare_models`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md)
for them. For two models of the same family, specifications that are not
nested give a warning.

**Overdispersion.** The \\\chi^2\\ reference also needs the larger model
to be correctly specified. Otherwise the statistic converges to \\\sum_j
\lambda_j Z_j^2\\, where \\\lambda_j\\ are eigenvalues of \\A^{-1}B\\ in
the tested directions, \\A\\ is the Hessian and \\B\\ the outer product
of the scores (Kent, 1982; White, 1982). With overdispersion the
\\\lambda_j\\ are roughly the dispersion, so comparing two Poisson
models rejects far too often: at a nominal 5\\ NB size of 2 and an
irrelevant covariate. The result therefore includes the Pearson
dispersion of the larger model (`NA` for Poisson against NB, which tests
overdispersion itself), and printing adds a note when it exceeds 1.5.
The same applies to `weights` that are not frequency weights, since the
weighted log-likelihood is then a pseudo-likelihood.

**Robust alternative.** `test = "Wald"` tests the larger model's extra
alpha and beta coefficients with its stored covariance (the robust
sandwich by default), \\(Rb)'(RVR')^{-1}(Rb) \sim \chi^2_r\\, and stays
near its nominal level under overdispersion (about 8\\ with n = 200). It
needs two models of the same method, with the same gamma specification
and constraint setting, so it does not cover Poisson against NB or
comparisons of gamma specifications.

This function is intended for comparing count models only (Poisson and
NB). OLS/NLS models are not supported because they use a
pseudo-log-likelihood. The models must be fitted to the same
observations and weights: the same `observed`, `auxiliary`,
`reference_pop` and `weights`.

## References

Chernoff, H. (1954). On the distribution of the likelihood ratio.
*Annals of Mathematical Statistics*, 25(3), 573–578.
[doi:10.1214/aoms/1177728725](https://doi.org/10.1214/aoms/1177728725)

Kent, J. T. (1982). Robust properties of likelihood ratio tests.
*Biometrika*, 69(1), 19–27.
[doi:10.1093/biomet/69.1.19](https://doi.org/10.1093/biomet/69.1.19)

Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605–610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

White, H. (1982). Maximum likelihood estimation of misspecified models.
*Econometrica*, 50(1), 1–25.
[doi:10.2307/1912526](https://doi.org/10.2307/1912526)

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
#> Estimated gamma is at its upper bound (0.5); the optimum may lie outside gamma_bounds = c(1e-10, 0.5).
fit_nb <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                              reference_pop = ~N, method = "nb")
#> Warning: Some alpha values < 0 (min = -1.28). Consider using constrained = TRUE.
#> Estimated gamma is at its upper bound (0.5); the optimum may lie outside gamma_bounds = c(1e-10, 0.5).

# Test for overdispersion (Poisson vs NB)
lrtest(fit_po, fit_nb)
#> Likelihood ratio test
#> ---------------------------------------- 
#> Model 1: POISSON   (logLik = -52.51 )
#> Model 2: NB   (logLik = -52.51 )
#> LR statistic: 0 on 1 df
#> (Boundary-corrected: 0.5 * chi2_0 + 0.5 * chi2_1, Self & Liang 1987)
#> p-value: 1 
```
