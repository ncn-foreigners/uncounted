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
  for OLS/NLS).

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

A warning is issued when comparing OLS/NLS pseudo-log-likelihoods with
true count-model likelihoods, since AIC/BIC values are not comparable
across these families.

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
#> Warning: Some alpha values < 0 (min = -1.274). Consider using constrained = TRUE.
fit_nb <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
                              reference_pop = ~N, method = "nb")
#> Warning: Some alpha values < 0 (min = -1.244). Consider using constrained = TRUE.
#> NB with theta in sandwich: HC3 not available, using HC1 correction.

comp <- compare_models(Poisson = fit_po, NB = fit_nb, sort_by = "AIC")
#> Warning: Some alpha values < 0 (min = -1.274). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -1.244). Consider using constrained = TRUE.
#> NB with theta in sandwich: HC3 not available, using HC1 correction.
comp
#> Model comparison
#> ------------------------------------------------------------ 
#>    Model  Method Constrained n_par logLik    AIC    BIC Deviance Pearson_X2
#>  Poisson POISSON       FALSE     3 -52.51 111.02 115.22    33.82      26.32
#>       NB      NB       FALSE     4 -52.52 113.05 118.65    33.85      26.35
#>  RMSE R2_cor R2_D  R2_CW
#>  1.43 0.1019    0 0.1209
#>  1.43 0.1018    0 0.1198
```
