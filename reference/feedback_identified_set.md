# Identified Set for Population Size under Detection Feedback

Bounds the population size over a range of \\\delta\\, with a confidence
interval for the true value or for the whole set.

## Usage

``` r
feedback_identified_set(
  x,
  delta_range = NULL,
  level = 0.95,
  method = c("imbens_manski", "union"),
  se = c("sd", "iqr")
)

# S3 method for class 'uncounted_feedback_set'
print(x, ...)
```

## Arguments

- x:

  An `"uncounted_feedback"` object.

- delta_range:

  Length-2 range of \\\delta\\. Defaults to `x$admissible_range`.

- level:

  Confidence level.

- method:

  `"imbens_manski"` covers the true \\\xi_g\\ (Imbens and Manski, 2004).
  `"union"` covers the whole set, so it is wider.

- se:

  How to measure the spread of the replicate bounds on the log scale:
  `"sd"`, or `"iqr"` for \\\mathrm{IQR} / 1.349\\.

- ...:

  Ignored.

## Value

A data frame of class `"uncounted_feedback_set"` with one row per group
and columns `group`, `L`, `U`, `delta_at_min`, `delta_at_max`,
`ci_lower`, `ci_upper`, `c`, `Delta`, `se_L` and `se_U`. Without a
bootstrap in `x` the interval columns are `NA`.

## Details

With \\D\\ = `delta_range`, the bounds are \\L_g = \min\_{\delta \in D}
\xi_g(\delta)\\ and \\U_g = \max\_{\delta \in D} \xi_g(\delta)\\. They
are found on a grid and refined with
[`optimize()`](https://rdrr.io/r/stats/optimize.html), not read from the
ends of \\D\\, because group totals need not be monotone in \\\delta\\.
Replicate bounds are replicate-wise minima and maxima over \\D\\.

`"imbens_manski"`, on the log scale with \\\Delta = \log U_g - \log
L_g\\, solves for \\c\\ in \$\$\Phi\left(c + \frac{\Delta}{\max(s_L,
s_U)}\right) - \Phi(-c) = \mathrm{level},\$\$ which lies between
`qnorm(level)` and `qnorm((1 + level) / 2)`, and returns \\\exp(\log
L_g - c s_L)\\ and \\\exp(\log U_g + c s_U)\\. `"union"` returns the
\\(1 - \mathrm{level}) / 2\\ quantile of the replicate \\L_g\\ and the
\\(1 + \mathrm{level}) / 2\\ quantile of the replicate \\U_g\\.

**Validity.** The Imbens–Manski interval requires superefficient
estimation of \\\Delta\\ (Stoye, 2009). Stoye's Lemma 3 gives this when
the bounds are jointly asymptotically normal and \\U_g \ge L_g\\ by
construction, as here. Normality can fail when the minimum or maximum
over \\D\\ is attained at more than one \\\delta\\, so treat the
intervals with caution if `delta_at_min` or `delta_at_max` is interior
to \\D\\. \\D\\ itself is treated as fixed.

## References

Imbens, G. W. and Manski, C. F. (2004). Confidence intervals for
partially identified parameters. *Econometrica*, 72(6), 1845–1857.
[doi:10.1111/j.1468-0262.2004.00555.x](https://doi.org/10.1111/j.1468-0262.2004.00555.x)

Stoye, J. (2009). More on confidence intervals for partially identified
parameters. *Econometrica*, 77(4), 1299–1315.
[doi:10.3982/ECTA7347](https://doi.org/10.3982/ECTA7347)

## See also

[`feedback_sensitivity`](https://ncn-foreigners.github.io/uncounted/reference/feedback_sensitivity.md)
