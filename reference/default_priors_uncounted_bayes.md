# Default Priors for Bayesian Hidden-Population Models

Default Priors for Bayesian Hidden-Population Models

## Usage

``` r
default_priors_uncounted_bayes(
  method = c("poisson", "nb"),
  constrained = FALSE,
  gamma = "estimate",
  gamma_bounds = c(1e-10, 0.5)
)
```

## Arguments

- method:

  Count likelihood. Currently `"poisson"` and `"nb"` are supported.

- constrained:

  Logical. If `TRUE`, alpha is mapped through an inverse-logit transform
  and beta through an exponential transform.

- gamma:

  Numeric fixed gamma, `"estimate"` for one bounded scalar gamma, or
  `NULL` to use \\n_i / N_i\\ directly.

- gamma_bounds:

  Length-two numeric bounds used when `gamma = "estimate"`.

## Value

A brms prior specification.
