# Profile Likelihood for Uncounted Models

S3 method for [`profile`](https://rdrr.io/r/stats/profile.html) that
dispatches to
[`profile_gamma`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md),
[`profile_alpha`](https://ncn-foreigners.github.io/uncounted/reference/profile_alpha.md),
[`profile_beta`](https://ncn-foreigners.github.io/uncounted/reference/profile_beta.md),
or
[`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md)
depending on `param`.

## Usage

``` r
# S3 method for class 'uncounted'
profile(fitted, param = c("gamma", "alpha", "beta", "dependence"), ...)
```

## Arguments

- fitted:

  An `"uncounted"` object.

- param:

  Character: which parameter to profile. One of `"gamma"`, `"alpha"`,
  `"beta"`, or `"dependence"`.

- ...:

  Additional arguments passed to the specific profile function.

## Value

Invisibly, a data frame from the dispatched function.
