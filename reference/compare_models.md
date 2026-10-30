# Compare multiple uncounted models

Produces a side-by-side comparison table for two or more fitted
`uncounted` models. The table includes the following metrics for each
model:

## Usage

``` r
compare_models(..., sort_by = c("AIC", "BIC", "loglik"))
```

## Arguments

- ...:

  Two or more objects of class `"uncounted"`. Named arguments are used
  as model labels.

- sort_by:

  Column to sort by: `"AIC"` (default), `"BIC"`, or `"loglik"`.

## Value

An object of class `"uncounted_comparison"` containing a `data.frame` in
`$table` and the model objects in `$models`.

## Details

- `Method`:

  Estimation method (OLS, NLS, POISSON, NB).

- `Constrained`:

  Whether the model was fitted with parameter constraints.

- `n_par`:

  Number of estimated parameters.

- `logLik`:

  Log-likelihood (true for count models, Gaussian pseudo-log-likelihood
  for OLS/NLS on the scale each one optimises). See
  [`logLik.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/logLik.uncounted.md).

- `AIC, BIC`:

  Akaike and Bayesian information criteria.

- `Deviance`:

  Model deviance (see
  [`deviance.uncounted`](https://ncn-foreigners.github.io/uncounted/reference/deviance.uncounted.md)).

- `Pearson_X2`:

  Pearson chi-squared statistic \\\sum(m_i - \hat\mu_i)^2 /
  V(\hat\mu_i)\\.

- `RMSE`:

  Root mean squared error of response residuals.

- `R2_cor`:

  Squared correlation between observed and fitted values:
  \\\mathrm{cor}(m, \hat\mu)^2\\.

- `R2_D`:

  Explained deviance pseudo \\R^2\\: \\1 - D(\text{model}) /
  D(\text{null})\\, where the null model is the same specification
  without covariates (single \\\alpha\\, single \\\beta\\). Only for
  Poisson/NB. See Cameron and Trivedi (2013).

- `R2_CW`:

  Cameron–Windmeijer (1996) pseudo \\R^2\\: \\1 - \sum (m_i -
  \hat\mu_i)^2 / V(\hat\mu_i) \big/ \sum (m_i - \bar m)^2 / V(\bar m)\\.
  Uses the model-implied variance function \\V(\mu)\\. Only for
  Poisson/NB.

**Comparability.** `logLik`, `AIC`, `BIC` and `Deviance` are comparable
only within one likelihood family: OLS (log-scale Gaussian), NLS
(raw-scale Gaussian), iOLS (Gamma pseudo-likelihood without constants),
or Poisson/NB (count likelihoods). A table can still mix families, for
example to line up `RMSE` and `R2_cor`, which are on the response scale
for every method. In that case a warning is issued and the columns are
kept, but their values must not be compared across families.

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

comp <- compare_models(Poisson = fit_po, NB = fit_nb, sort_by = "AIC")
#> Warning: Some alpha values < 0 (min = -1.28). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -1.28). Consider using constrained = TRUE.
comp
#> Model comparison
#> ------------------------------------------------------------ 
#>    Model  Method Estimator  Link Constrained n_par logLik    AIC    BIC
#>  Poisson POISSON       MLE power       FALSE     3 -52.51 111.02 115.22
#>       NB      NB       MLE power       FALSE     4 -52.51 113.02 118.62
#>  Deviance Pearson_X2 RMSE R2_cor R2_D  R2_CW
#>     33.82      26.32 1.43 0.1019    0 0.1208
#>     33.82      26.32 1.43 0.1019    0 0.1208
```
