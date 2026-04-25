# Methods for Bayesian Hidden-Population Models

S3 methods for Bayesian hidden-population models fitted with
[`estimate_hidden_pop_bayes`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop_bayes.md).

## Usage

``` r
# S3 method for class 'uncounted_bayes'
summary(
  object,
  total = FALSE,
  level = 0.95,
  diagnostics = TRUE,
  ...
)

# S3 method for class 'uncounted_bayes'
coef(object, summary = c("median", "mean", "draws"), ...)

# S3 method for class 'uncounted_bayes'
fitted(object, summary = c("mean", "median", "draws"), ...)

# S3 method for class 'uncounted_bayes'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "counts", "draws"),
  summary = TRUE,
  ...
)

# S3 method for class 'uncounted_bayes'
residuals(
  object,
  type = c("response", "pearson"),
  summary = c("mean", "median"),
  ...
)

# S3 method for class 'uncounted_bayes'
loo(object, ...)

# S3 method for class 'uncounted_bayes'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'uncounted_bayes'
glance(x, loo = FALSE, ...)
```

## Arguments

- object:

  An `"uncounted_bayes"` object.

- total:

  Logical; include a total population-size summary in
  [`summary()`](https://rdrr.io/r/base/summary.html).

- level:

  Credible interval level for summaries.

- diagnostics:

  Logical; include MCMC diagnostics in
  [`summary()`](https://rdrr.io/r/base/summary.html).

- summary:

  Posterior summary to return. For
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), use
  `"mean"`, `"median"`, or `"draws"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html), use
  `TRUE`/`"mean"`, `"median"`, or `FALSE` to return draws. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), use `"mean"`
  or `"median"` fitted values.

- newdata:

  Optional data frame for prediction. If `NULL`, predictions are for the
  original data.

- type:

  Prediction or residual type. For
  [`predict()`](https://rdrr.io/r/stats/predict.html), `"response"`
  returns posterior expected observed counts, `"link"` returns the log
  mean, `"counts"` returns replicated posterior predictive counts, and
  `"draws"` returns response-scale expected-count draws. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), `"response"`
  and `"pearson"` are supported.

- x:

  An `"uncounted_bayes"` object for broom-style methods.

- conf.int:

  Logical; include credible intervals in
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html).

- conf.level:

  Credible interval level for
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html).

- loo:

  Logical; if `TRUE`,
  [`glance()`](https://generics.r-lib.org/reference/glance.html) also
  computes [`brms::loo()`](https://mc-stan.org/loo/reference/loo.html)
  summaries.

- ...:

  Additional arguments passed to the underlying brms posterior
  prediction or LOO functions where relevant.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) invisibly returns a
list containing coefficient, population-size, and diagnostic summaries.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns posterior
coefficient summaries or draws.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) return numeric
summaries or draw matrices.
[`residuals()`](https://rdrr.io/r/stats/residuals.html) returns numeric
residuals.
[`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md)
returns the object from
[`brms::loo()`](https://mc-stan.org/loo/reference/loo.html).
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) return
data frames.

## Details

These methods report posterior quantities rather than Wald standard
errors and p values. Coefficient tables include posterior means,
medians, posterior standard deviations, credible intervals, R-hat, bulk
ESS, and tail ESS. Prediction methods delegate to the underlying brms
posterior prediction functions.

## Examples

``` r
if (FALSE) { # \dontrun{
summary(fit_bayes, total = TRUE)
coef(fit_bayes)
head(fitted(fit_bayes))
head(predict(fit_bayes, type = "counts", summary = FALSE))
head(residuals(fit_bayes))
tidy(fit_bayes, conf.int = TRUE)
glance(fit_bayes)
loo(fit_bayes)
} # }
```
