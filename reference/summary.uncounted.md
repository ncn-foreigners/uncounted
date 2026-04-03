# Summary of an Uncounted Population Model

Prints a detailed summary of the fitted model, including coefficient
estimates with robust standard errors, z- or t-values, p-values, and
goodness-of-fit statistics (log-likelihood, AIC, BIC, deviance when
available). The output concludes with estimated population sizes
\\\hat{\xi} = \sum_i N_i^{\hat{\alpha}\_i}\\ for each covariate group
defined by `cov_alpha`, along with bias-corrected estimates and
confidence intervals computed by
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

## Usage

``` r
# S3 method for class 'uncounted'
summary(object, total = FALSE, ...)
```

## Arguments

- object:

  An `"uncounted"` object returned by
  [`estimate_hidden_pop`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md).

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the coefficient summary matrix (estimates, standard
errors, test statistics, and p-values).

## Details

For OLS and NLS fits, inference uses t-statistics with \\n - p\\ degrees
of freedom. For Poisson and NB fits, z-statistics (normal approximation)
are used. Standard errors are always HC-robust (type controlled by
`vcov_type` in the original call).

## See also

[`estimate_hidden_pop`](https://ncn-foreigners.github.io/uncounted/reference/estimate_hidden_pop.md),
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
