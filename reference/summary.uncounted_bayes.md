# Summarize a Bayesian Hidden-Population Model

Summarize a Bayesian Hidden-Population Model

## Usage

``` r
# S3 method for class 'uncounted_bayes'
summary(object, total = FALSE, level = 0.95, diagnostics = TRUE, ...)
```

## Arguments

- object:

  An `"uncounted_bayes"` object.

- total:

  Logical; include a total population-size summary?

- level:

  Credible interval level.

- diagnostics:

  Logical; include HMC and posterior diagnostic summaries?

- ...:

  Additional arguments, currently ignored.

## Value

Invisibly, a list with coefficient, population-size, and diagnostic
summaries.
