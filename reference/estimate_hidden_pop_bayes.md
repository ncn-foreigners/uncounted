# Bayesian Hidden-Population Models via brms

Fit Bayesian Poisson or negative-binomial hidden-population models using
brms. The model mirrors the count-model part of
[`estimate_hidden_pop`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md):
\$\$\log(\mu_i) = \alpha_i \log(N_i) + \log\rho_i,\$\$ where \\\rho_i =
h(\beta_i \log(\gamma_i + n_i / N_i))\\. Posterior draws of \\\alpha_i\\
are propagated to population-size draws \\\xi = \sum_i N_i^{\alpha_i}\\
by
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

## Usage

``` r
estimate_hidden_pop_bayes(
  data,
  observed,
  auxiliary,
  reference_pop,
  method = c("poisson", "nb"),
  cov_alpha = NULL,
  cov_beta = NULL,
  gamma = "estimate",
  cov_gamma = NULL,
  gamma_bounds = c(1e-10, 0.5),
  link_rho = c("power", "cloglog", "logit", "probit"),
  constrained = FALSE,
  weights = NULL,
  countries = NULL,
  seed = NULL,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  cores = getOption("mc.cores", 1L),
  control = NULL,
  prior = NULL,
  backend = c("rstan", "cmdstanr"),
  ...
)
```

## Arguments

- data:

  Data frame.

- observed:

  One-sided formula naming the observed count.

- auxiliary:

  One-sided formula naming the auxiliary count.

- reference_pop:

  One-sided formula naming the reference population.

- method:

  Count likelihood. Currently `"poisson"` and `"nb"` are supported.

- cov_alpha:

  Optional one-sided formula for alpha covariates.

- cov_beta:

  Optional one-sided formula for beta covariates.

- gamma:

  Numeric fixed gamma, `"estimate"` for one bounded scalar gamma, or
  `NULL` to use \\n_i / N_i\\ directly.

- cov_gamma:

  Reserved for a future release. It must be `NULL`.

- gamma_bounds:

  Length-two numeric bounds used when `gamma = "estimate"`.

- link_rho:

  Detection link: `"power"`, `"cloglog"`, `"logit"`, or `"probit"`.

- constrained:

  Logical. If `TRUE`, alpha is mapped through an inverse-logit transform
  and beta through an exponential transform.

- weights:

  Optional numeric observation weights.

- countries:

  Optional one-sided formula for a country/group label stored on the
  fitted object.

- seed, chains, iter, warmup, cores, control:

  Passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- prior:

  Optional brms prior specification. If `NULL`, weakly informative
  defaults from
  [`default_priors_uncounted_bayes()`](https://ncn-foreigners.github.io/uncounted/reference/default_priors_uncounted_bayes.md)
  are used.

- backend:

  Stan backend passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Defaults to `"rstan"`; `"cmdstanr"` is available when installed.

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

An object of class `"uncounted_bayes"` containing the `brmsfit` and the
design metadata needed by
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data(irregular_migration)
d <- irregular_migration
d$ukr <- as.integer(d$country_code == "UKR")

fit_bayes <- estimate_hidden_pop_bayes(
  data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
  method = "poisson",
  cov_alpha = ~ year * ukr + sex,
  cov_beta = ~ year,
  gamma = 0.005,
  chains = 2, iter = 1000, seed = 1
)

ps <- popsize(fit_bayes, by = ~ year + country_code, total = TRUE)
hypotheses_popsize(
  ps,
  "xi[year == 2024 & country_code == 'UKR'] < 15000"
)
} # }
```
