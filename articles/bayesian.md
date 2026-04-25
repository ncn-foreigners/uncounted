# Bayesian hidden-population models

The frequentist workflow in `uncounted` estimates a count model and then
uses delta-method or bootstrap uncertainty for population-size
functionals. The Bayesian workflow uses the same mean structure but
samples from the posterior distribution of the model parameters.

For observation `i`,

$$\log\left( \mu_{i} \right) = \alpha_{i}\log\left( N_{i} \right) + \log h\{\beta_{i}\log\left( \gamma_{i} + n_{i}/N_{i} \right)\},$$

where `m_i` is the observed count, `N_i` is the reference population,
`n_i` is the auxiliary count, and `h()` is chosen by `link_rho`. The
posterior draws of `alpha_i` imply posterior draws of the
hidden-population estimand

$$\xi_{g} = \sum\limits_{i \in g}N_{i}^{\alpha_{i}}.$$

No Wald approximation is needed for `xi`:
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
evaluates this expression in every posterior draw.

## Fitting a Bayesian model

The Bayesian fitting function mirrors the count-model arguments of
[`estimate_hidden_pop()`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md).
The first release supports Poisson and negative binomial likelihoods
through `brms`, with `rstan` as the default backend.

``` r
library(uncounted)
data(irregular_migration)

d <- irregular_migration
d$ukr <- as.integer(d$country_code == "UKR")

fit_bayes <- estimate_hidden_pop_bayes(
  data = d,
  observed = ~m,
  auxiliary = ~n,
  reference_pop = ~N,
  method = "poisson",
  cov_alpha = ~ year * ukr + sex,
  cov_beta = ~ year,
  gamma = 0.005,
  chains = 2,
  iter = 1000,
  seed = 20260425,
  refresh = 0
)
```

The example is not evaluated during package checks because compiling and
sampling Stan models can take several minutes. To run the vignette code,
set `UNCOUNTED_RUN_STAN_VIGNETTE=true`.

## Posterior population size

``` r
ps_year_country <- popsize(
  fit_bayes,
  by = ~ year + country_code,
  total = TRUE
)

head(ps_year_country)
posterior_popsize_draws(ps_year_country, include_total = TRUE)[1:3, 1:4]
```

The Bayesian
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
method reports posterior medians, means, posterior standard deviations,
and credible intervals. The `estimate` column is the chosen point
summary; by default it is the posterior median. The full posterior draw
matrix is stored on the object and can be extracted with
[`posterior_popsize_draws()`](https://ncn-foreigners.github.io/uncounted/reference/posterior_popsize_draws.md).

## Posterior hypotheses

Bayesian hypotheses are posterior propositions. For example:

``` r
hypotheses_popsize(
  ps_year_country,
  "xi[year == 2024 & country_code == 'UKR'] < 15000"
)

hypotheses_popsize(
  ps_year_country,
  "xi[year == 2024 & country_code == 'UKR'] -
     xi[year == 2019 & country_code == 'UKR'] < 0"
)
```

With the default `hypothesis_side = "alternative"`, the expression is
treated as `H1`. Thus `xi[...] < 15000` means:

- `H0: xi[...] >= 15000`
- `H1: xi[...] < 15000`

The output reports `p_h1 = Pr(H1 | data)` and `p_h0 = Pr(H0 | data)`
from the posterior contrast draws. This is not a frequentist p value. It
is a direct posterior probability under the model and priors.

If you want the written expression to be the null hypothesis, use
`hypothesis_side = "null"`:

``` r
hypotheses_popsize(
  ps_year_country,
  "xi[year == 2024 & country_code == 'UKR'] < 15000",
  hypothesis_side = "null"
)
```

Equality hypotheses are different. A continuous posterior assigns zero
probability to an exact point null, so
[`hypotheses_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md)
summarizes equality-style contrasts with posterior means, posterior
standard deviations, and credible intervals rather than a posterior
point-null probability. For practical equivalence, use a ROPE:

``` r
hypotheses_popsize(
  ps_year_country,
  "xi[year == 2024 & country_code == 'UKR'] -
     xi[year == 2019 & country_code == 'UKR'] = 0",
  rope = c(-5000, 5000)
)
```

## Priors and diagnostics

[`estimate_hidden_pop_bayes()`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop_bayes.md)
supplies weakly informative defaults through
[`default_priors_uncounted_bayes()`](https://ncn-foreigners.github.io/uncounted/reference/default_priors_uncounted_bayes.md),
but these are only starting points. The posterior for `xi` can be
sensitive to priors on alpha because `N_i^alpha_i` is nonlinear and
`N_i` can be large.

Recommended workflow:

- inspect default priors with
  [`default_priors_uncounted_bayes()`](https://ncn-foreigners.github.io/uncounted/reference/default_priors_uncounted_bayes.md);
- run prior predictive checks in `brms` when changing priors;
- examine R-hat, bulk/tail effective sample sizes, divergent
  transitions, and treedepth warnings;
- use posterior predictive checks to verify that the count model can
  reproduce observed `m`;
- compare Poisson and negative-binomial fits with `loo` when both are
  scientifically plausible.

The underlying `brmsfit` is available with
[`as_brmsfit()`](https://ncn-foreigners.github.io/uncounted/reference/as_brmsfit.md):

``` r
brms_fit <- as_brmsfit(fit_bayes)
summary(brms_fit)
```

## Power and design analysis

In the frequentist Wald setting, power is tied to the contrast size,
standard error, test level, sample size, dependence, and model
specification. In the Bayesian setting, there is no single Wald-test
power number attached to
[`hypotheses_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md).

For formal design analysis, use simulation:

1.  Specify data-generating assumptions and priors.
2.  Simulate many datasets under alternatives of interest.
3.  Fit the Bayesian model to each dataset.
4.  Compute posterior propositions, such as
    `Pr(xi_2024 > 200000 | data)`.
5.  Estimate the fraction of simulations meeting your decision rule, for
    example `Pr(H1 | data) > 0.95`.

That simulation-based quantity is the Bayesian analogue of power for a
chosen decision rule.
