# Detection-Feedback Sensitivity Analysis for Population Size

Relaxes the "no feedback" assumption of the main model, under which
detection does not depend on the size of the unauthorized population,
and reports the population size consistent with each degree of feedback.
The analysis is a deterministic transform of the fitted cell-level
estimates: nothing is refitted.

## Usage

``` r
feedback_sensitivity(
  fit,
  anchor,
  delta_grid = seq(-0.3, 0.3, by = 0.01),
  K = NULL,
  anchor_by = NULL,
  by = NULL,
  total = FALSE,
  boot = NULL,
  level = 0.95,
  sublinear = TRUE,
  p_max = 1,
  keep_cells = FALSE,
  delta_prior = NULL
)

# S3 method for class 'uncounted_feedback'
print(x, ...)

# S3 method for class 'uncounted_feedback'
plot(x, group = NULL, plausible = NULL, ...)
```

## Arguments

- fit:

  An `"uncounted"` object.

- anchor:

  How to choose \\K_i\\, the size at which the main model's detection
  function is assumed correct. There is no default, because the anchor
  can reverse the direction of the result (see Details). `"geomean"`
  uses the geometric mean of \\\xi\_{0i}\\ over all cells (or within
  `anchor_by`), so "large" means large in absolute terms. `"share"` uses
  \\K_i = \bar{s} N_i\\ with \\\bar{s}\\ the geometric mean of
  \\\xi\_{0i} / N_i\\, so "large" means large relative to the reference
  population. `"value"` uses `K`.

- delta_grid:

  Numeric vector of \\\delta\\ values. Must contain 0, and every value
  must exceed \\-1\\.

- K:

  Anchor for `anchor = "value"`: a positive scalar, or a positive vector
  with one value per observation.

- anchor_by:

  Optional one-sided formula: compute the `"geomean"` or `"share"`
  anchor within groups. With `~ country_code` the levels of countries
  barely move, so the analysis concerns changes within countries.

- by:

  Optional one-sided formula for the reported groups, with the same
  meaning as in
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

- total:

  Logical; add a `"Total"` row when there are several groups? With
  `by = ~ year` this adds up the years and is not a population size,
  because the same people can be present in several years (see
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)).
  Use it only for non-overlapping groups.

- boot:

  Optional
  [`bootstrap_popsize`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
  result for the same fit. Its replicate coefficients are transformed;
  the bootstrap is not re-run. The anchor stays at its point-estimate
  value in every replicate, because it is part of the assumption rather
  than an estimate.

- level:

  Level of the percentile bootstrap intervals.

- sublinear:

  Logical; apply constraint A1 (see Details)?

- p_max:

  Upper limit for the implied detection in constraint A2.

- keep_cells:

  Logical; return cell-level matrices in `$cells`?

- delta_prior:

  Optional Monte Carlo sensitivity analysis (Greenland, 2005): a
  distribution for \\\delta\\ instead of a range. Either
  `list(mean = , sd = )` for a normal distribution truncated to the
  admissible range, a function of `R` returning `R` draws, or `TRUE` for
  a uniform distribution on the admissible range. One draw is used per
  bootstrap replicate, so `boot` is required. The data do not update
  this distribution, because \\\delta\\ is not identified: the result
  depends completely on the distribution chosen, so report it and try
  others.

- x:

  Object to print or plot.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- group:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): labels
  of the groups to draw. Defaults to `"Total"` when available, otherwise
  the first six groups.

- plausible:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): optional
  length-2 range of \\\delta\\ to mark as plausible.

## Value

An object of class `"uncounted_feedback"`, a list with:

- `profile`:

  Data frame with columns `delta`, `group`, `estimate`, `boot_median`,
  `lower`, `upper` and `admissible`.

- `admissible_range`:

  `c(delta_L, delta_U)`, or `NA`s when \\\delta = 0\\ is inadmissible.

- `anchor`:

  List with the `type`, the anchor `K` per cell, and `by`.

- `diagnostics`:

  Data frame per \\\delta\\: `A1`, `A2`, `A3`, `binding`, `max_p`,
  `max_p_row`, and `d_total`, the derivative of the total with respect
  to \\\delta\\.

- `a2_excluded`:

  Data frame of the cells left out of A2.

- `cells`:

  With `keep_cells = TRUE`: matrices `xi`, `detection` and `alpha`
  (cells by \\\delta\\). Implied `alpha` is \\\hat\alpha_i / (1 +
  \delta)\\ for constant anchors, \\1 - (1 - \hat\alpha_i) / (1 +
  \delta)\\ for `"share"`, and `NA` for `"value"`.

- `monte_carlo`:

  With `delta_prior`: data frame with `group`, `median`, `lower` and
  `upper`. The \\\delta\\ draws are in `attr(, "delta")` and the totals
  in `attr(, "draws")` (groups by replicates), from which changes
  between years can be computed.

## Details

**Model.** The main model assumes \\E(p_i \mid M_i, N_i, n_i) =
\rho_i\\. Here \\\delta\\ is the elasticity of detection with respect to
the unauthorized population \\M_i\\:

\$\$E(p_i \mid M_i, N_i, n_i) = \rho_i (M_i / K_i)^\delta, \qquad \delta
\> -1,\$\$

which implies \\\mu_i \approx \rho_i \xi_i (\xi_i / K_i)^\delta\\. The
approximation treats \\M_i\\ as equal to its mean. Holding the fitted
\\\hat\mu_i\\ and \\\hat\rho_i\\ fixed, the population size consistent
with \\\delta\\ is

\$\$\xi_i(\delta) = K_i (\xi\_{0i} / K_i)^{1 / (1 + \delta)}, \qquad
\xi\_{0i} = N_i^{\hat\alpha_i}.\$\$

**Fit preservation.** Fitted values, likelihood and residuals are the
same for every \\\delta\\: only the split of \\\hat\mu_i\\ into size and
detection changes. So \\\delta\\ is not identified and the data cannot
choose among its values (compare Link, 2003; Barker et al., 2018).
Refitting with an offset would not help. The offset model is \\\log\mu_i
= (1+\delta)\alpha_i \log N_i + \beta_i \log r_i - \delta \log K\\, and
because the model has no intercept, the likelihood varies with
\\\delta\\ only through the normalization \\-\delta \log K\\: with \\K =
1\\ it does not vary at all.

**The anchor is an assumption.** \\K_i\\ is a pivot: for \\\delta \>
0\\, cells above it shrink and cells below it grow. Totals are dominated
by large cells, so the direction depends on which side of the pivot they
lie. Under `"geomean"` large countries are far above \\K\\ and the total
falls with \\\delta\\. Under `"share"` they usually have a below-typical
\\\xi\_{0i} / N_i\\, so they lie below \\K_i\\ and the total can rise.
With `"geomean"` computed over exactly the summed cells, the total is
non-increasing in \\\delta\\. Report results under more than one anchor.

**Admissibility**, from point estimates only:

- A1 (`sublinear = TRUE`, constant anchors only): \\\delta \> \max_i
  \hat\alpha_i - 1\\, so that size stays sublinear in \\N\\. It does not
  apply to `"share"` or to a `K` that varies across cells. With
  `anchor_by`, the anchor is constant within groups and A1 is applied.

- A2: implied detection \\\hat\mu_i / \xi_i(\delta) \le\\ `p_max`. Cells
  that already exceed `p_max` at \\\delta = 0\\ (possible because the
  power link does not bound \\\rho_i\\) would make every \\\delta\\
  inadmissible, so they are left out of A2 with a warning and listed in
  `$a2_excluded`.

- A3: for every group, \\\sum\_{i \in g} \xi_i(\delta) \ge \sum\_{i \in
  g} m_i\\.

`$admissible_range` is the largest run of consecutive admissible grid
points containing 0, so its ends are limited by `delta_grid`. The range
can hinge on very small cells; `$diagnostics` shows which constraint
binds and where.

**Point estimates.** `estimate` uses the plug-in \\\xi\_{0i}\\, because
the analytical bias correction in
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
was derived for \\\delta = 0\\. With a bootstrap, `boot_median` is the
recommended point estimate, as in
[`bootstrap_popsize`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md).

## References

Barker, R. J., Schofield, M. R., Link, W. A., and Sauer, J. R. (2018).
On the reliability of N-mixture models for count data. *Biometrics*,
74(1), 369–377.
[doi:10.1111/biom.12734](https://doi.org/10.1111/biom.12734)

Greenland, S. (2005). Multiple-bias modelling for analysis of
observational data. *Journal of the Royal Statistical Society: Series
A*, 168(2), 267–306.
[doi:10.1111/j.1467-985X.2004.00349.x](https://doi.org/10.1111/j.1467-985X.2004.00349.x)

Link, W. A. (2003). Nonidentifiability of population size from
capture-recapture data with heterogeneous detection probabilities.
*Biometrics*, 59(4), 1123–1130.
[doi:10.1111/j.0006-341X.2003.00129.x](https://doi.org/10.1111/j.0006-341X.2003.00129.x)

## See also

[`feedback_identified_set`](https://ncn-foreigners.github.io/uncounted/reference/feedback_identified_set.md),
[`feedback_tipping`](https://ncn-foreigners.github.io/uncounted/reference/feedback_tipping.md),
[`feedback_calibrate`](https://ncn-foreigners.github.io/uncounted/reference/feedback_calibrate.md),
[`dependence_bounds`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md),
[`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md)
(where `delta` is a different quantity, a constant shift of the mean)

## Examples

``` r
set.seed(123)
d <- data.frame(N = round(exp(rnorm(40, mean = 7, sd = 1))))
d$n <- rpois(40, lambda = pmax(1, 0.03 * d$N))
d$m <- rpois(40, lambda = d$N^0.6 * (0.01 + d$n / d$N)^0.6)
fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.01)

feedback_sensitivity(fit, anchor = "geomean")
#> Detection-feedback sensitivity analysis
#> Anchor: geomean (K = 40.17)
#> The anchor is an assumption and can change the direction of the result.
#> Admissible range: [-0.3, 0.3]
#> 
#> Plug-in estimates:
#>       d=-0.2 d=-0.1 d=0.0 d=0.1 d=0.2
#> (all) 1,894  1,831  1,787 1,754 1,730
feedback_sensitivity(fit, anchor = "share")
#> Detection-feedback sensitivity analysis
#> Anchor: share (K from 5.358 to 229.3)
#> The anchor is an assumption and can change the direction of the result.
#> Admissible range: [-0.3, 0.3]
#> 
#> Plug-in estimates:
#>       d=-0.2 d=-0.1 d=0.0 d=0.1 d=0.2
#> (all) 1,713  1,751  1,787 1,819 1,848
```
