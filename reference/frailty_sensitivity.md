# Omitted-Frailty Sensitivity Analysis for Hidden Population Size

Computes a `sensemakr`-style omitted-frailty sensitivity analysis for
the hidden-population estimand \\\Xi\\. The method linearizes the fitted
Poisson or Negative Binomial mean model with one IRLS step, treats the
omitted shared frailty as an unobserved weighted confounder, and reports
how large the resulting first-order bias would have to be to change the
plug-in hidden-population estimate.

## Usage

``` r
frailty_sensitivity(
  object,
  by = NULL,
  r2_d = seq(0, 0.25, by = 0.01),
  r2_y = seq(0, 0.25, by = 0.01),
  q = c(0.1, 0.25, 0.5),
  threshold = NULL,
  direction = c("decrease", "increase"),
  df_method = c("auto", "model", "cluster"),
  plot = TRUE,
  ...
)

# S3 method for class 'uncounted_frailty_sensitivity'
print(x, ...)

# S3 method for class 'uncounted_frailty_sensitivity'
plot(x, group = NULL, type = c("ratio", "xi"), side = c("lower", "upper"), ...)
```

## Arguments

- object:

  An `"uncounted"` object.

- by:

  Optional one-sided formula defining the \\\Xi\\ targets. This uses the
  same grouping logic as
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md).

- r2_d:

  Numeric vector of weighted partial-\\R^2\\ values describing how
  strongly the omitted frailty aligns with the targeted alpha contrast.
  Values must satisfy `0 <= r2_d < 1`. Duplicates are removed and the
  grid is sorted increasingly.

- r2_y:

  Numeric vector of weighted partial-\\R^2\\ values describing how
  strongly the omitted frailty predicts the IRLS working outcome after
  conditioning on the fitted model. Values must satisfy
  `0 <= r2_y <= 1`. Duplicates are removed and the grid is sorted
  increasingly.

- q:

  Numeric vector of relative target changes used when
  `threshold = NULL`. For `direction = "decrease"`, each target is
  `Xi_hat * (1 - q)` and values must lie in `[0, 1]`. For
  `direction = "increase"`, the targets are `Xi_hat * (1 + q)` and
  values must be non-negative.

- threshold:

  Optional positive numeric threshold for \\\Xi\\. When supplied, `q` is
  ignored for target construction and the same threshold is applied to
  every reported group.

- direction:

  Which target crossing to summarize: `"decrease"` or `"increase"`.

- df_method:

  Degrees-of-freedom rule for the omitted-frailty formulas. `"model"`
  uses `n_obs - rank(X_work)`. `"cluster"` uses the number of observed
  clusters minus the same rank when cluster metadata is present on the
  fitted object. `"auto"` prefers cluster degrees of freedom when
  available and valid, otherwise falls back to the model-based rule.

- plot:

  Logical; if `TRUE` (default), draw a contour plot for the first
  reported group.

- ...:

  Additional graphical arguments passed to
  `plot.uncounted_frailty_sensitivity` when `plot = TRUE`.

- x:

  Object to print or plot.

- group:

  Optional group label for the plot method. When `NULL`, the first
  available group is used.

- type:

  Plot scale: `"ratio"` for values relative to baseline \\\Xi\\, or
  `"xi"` for raw \\\Xi\\ values.

- side:

  Which sensitivity surface to draw: `"lower"` or `"upper"`.

## Value

An object of class `"uncounted_frailty_sensitivity"` with:

- `baseline`:

  A
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  table computed with `bias_correction = FALSE` and the requested `by =`
  grouping.

- `working`:

  Data frame with one row per reported target and columns including
  `group`, `xi_hat`, `theta_hat`, `se_theta`, `d_xi_dtheta`, and `df`.

- `surface`:

  Tidy grid over `group`, `r2_d`, and `r2_y` containing `bias_theta`,
  `xi_lower`, `xi_upper`, `xi_ratio_lower`, and `xi_ratio_upper`.

- `robustness`:

  Data frame with robustness values for each requested target.

- `settings`:

  List of resolved inputs and defaults.

## Details

The baseline uncounted mean model factorizes the observed-count mean as

\$\$\mu_i = E(m_i \mid N_i, n_i) = \xi_i \rho_i,\$\$

where \\\xi_i\\ is the latent population component and \\\rho_i\\ is the
detection component. Under the default power-link specification,

\$\$ \log \mu_i = \alpha_i \log N_i + \beta_i \log(\gamma + n_i / N_i).
\$\$

This helper studies violations of the identifying assumption by
introducing an omitted shared frailty \\U_i\\ into the log-mean:

\$\$ \log \mu_i^\dagger = \alpha_i \log N_i + \beta_i \log(\gamma + n_i
/ N_i) + \delta U_i. \$\$

The frailty is not observed or estimated directly. Instead, its strength
is indexed by two weighted partial-\\R^2\\ quantities:

- `r2_d`:

  The weighted partial association between the omitted frailty and the
  targeted alpha contrast after conditioning on the fitted nuisance
  regressors.

- `r2_y`:

  The weighted partial association between the omitted frailty and the
  working IRLS outcome after conditioning on the targeted contrast and
  the same nuisance regressors.

For a fitted Poisson or NB model, the method linearizes the mean model
with the working response

\$\$ Y_i^\* = \hat\eta_i + \frac{m_i - \hat\mu_i}{\hat\mu_i}, \$\$

using working weights

\$\$ w_i = \hat\mu_i \qquad \text{(Poisson)} \$\$

and

\$\$ w_i = \frac{\hat\mu_i}{1 + \hat\mu_i / \hat\theta} \qquad
\text{(NB2)}. \$\$

The alpha block of the weighted linearized design is \\X\_{\alpha,i}
\log N_i\\; the beta block is \\X\_{\beta,i} \log(\gamma + n_i/N_i)\\.
When scalar gamma is estimated, the nuisance design also includes the
derivative column

\$\$ \frac{\partial \eta_i}{\partial \gamma} =
\frac{\hat\beta_i}{\hat\gamma + n_i / N_i}. \$\$

The primary estimand is not a raw alpha coefficient but the grouped
hidden-population estimate

\$\$ \hat\Xi_g = \sum\_{i \in g} N_i^{\hat\alpha_i}. \$\$

For each requested group, the helper computes the gradient of
\\\hat\Xi_g\\ with respect to the alpha-coefficient vector and uses that
gradient to define a one-dimensional alpha contrast. This makes the
omitted frailty analysis applicable to vector-valued `cov_alpha` models
such as `~ year * ukr + sex`. The reported contrast is an internal
device: it aligns the sensitivity calculation with the direction in
coefficient space that matters most for the chosen \\\Xi_g\\.

Let \\\theta_g\\ denote the targeted alpha contrast, with estimated
standard error \\\widehat{\mathrm{se}}(\hat\theta_g)\\ and degrees of
freedom \\\nu\\. The first-order omitted-frailty bias is

\$\$ B\_{\theta,g}(r2_d, r2_y) = \widehat{\mathrm{se}}(\hat\theta_g)
\sqrt{ \nu \frac{r2_y \\ r2_d}{1 - r2_d} }. \$\$

The hidden-population effect is then approximated by a delta method:

\$\$ \hat\Xi\_{g,\mathrm{lower}} \approx \hat\Xi_g - \frac{\partial
\Xi_g}{\partial \theta_g} B\_{\theta,g}(r2_d, r2_y), \$\$

\$\$ \hat\Xi\_{g,\mathrm{upper}} \approx \hat\Xi_g + \frac{\partial
\Xi_g}{\partial \theta_g} B\_{\theta,g}(r2_d, r2_y). \$\$

The returned `surface` table stores these lower/upper approximations as
both raw \\\Xi\\ values and ratios relative to the baseline estimate.

The `robustness` table summarizes two robustness values for threshold
questions:

- `rv_equal`:

  The minimum equal-strength partial-\\R^2\\ required for the omitted
  frailty to move the targeted \\\Xi\\ to the requested threshold under
  \\r2_d = r2_y\\.

- `rv_extreme`:

  The corresponding treatment-side strength when the omitted frailty is
  given the extreme outcome-side scenario \\r2_y = 1\\.

These robustness values are derived from the same first-order bias
mapping. They quantify how much residual alignment between the omitted
frailty and the targeted alpha contrast would be required to overturn a
substantive claim.

**Interpretation.** This is an identification-sensitivity analysis,
distinct from the package's existing distributional sensitivity
(`method = "nb"`) and bootstrap uncertainty summaries. It is inspired by
omitted-variable-bias diagnostics such as sensemakr, but it is not a
literal OLS port: the count model is first linearized with its IRLS
working representation, and the hidden-population estimand is handled
through a grouped delta-method contrast.

**Limitations.** The method is first-order and deliberately restricted
to `link_rho = "power"`, unconstrained Poisson/NB MLE fits without
`cov_gamma`. It does not yet implement benchmark covariates or a
simulation-based calibration layer.

## References

Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity:
Extending omitted variable bias. *Journal of the Royal Statistical
Society: Series B*, 82(1), 39–67.
[doi:10.1111/rssb.12348](https://doi.org/10.1111/rssb.12348)

Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
irregular foreigners in Poland using non-linear count regression models.
*arXiv preprint* arXiv:2008.09407.

Zhang, L.-C. (2008). Developing methods for determining the number of
unauthorized foreigners in Norway. *Documents* 2008/11, Statistics
Norway.

## See also

[`dependence_bounds`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md),
[`profile_dependence`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md),
[`robustness_dependence`](https://ncn-foreigners.github.io/uncounted/reference/robustness_dependence.md),
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)

## Examples

``` r
data(irregular_migration)
d <- irregular_migration[
  irregular_migration$m > 0 & irregular_migration$n > 0,
]
keep <- unique(d$country)[1:8]
d <- droplevels(d[d$country %in% keep, ])

fit <- estimate_hidden_pop(
  data = d,
  observed = ~m,
  auxiliary = ~n,
  reference_pop = ~N,
  method = "poisson",
  cov_alpha = ~ sex,
  gamma = 0.005
)

frailty_sensitivity(
  fit,
  by = ~ year,
  r2_d = c(0, 0.05, 0.10),
  r2_y = c(0, 0.05, 0.10),
  plot = FALSE
)
#> Omitted-frailty sensitivity analysis
#> Targets: 7 | DF method: model 
#> 
#>  group   xi_hat
#>   2019 2073.808
#>   2020 2228.864
#>   2021 2575.444
#>   2022 3340.510
#>   2023 3721.804
#>   2024 3866.374
#> 
#> Total xi: 17,807 
#> 
#> Working contrasts:
#>  group theta_hat   se_theta d_xi_dtheta df
#>   2019 0.5176578 0.05241501    23139.03 68
#>   2020 0.5170375 0.05229633    25443.33 68
#>   2021 0.5142046 0.05175732    31247.05 68
#>   2022 0.5139312 0.05170557    42729.20 68
#>   2023 0.5115794 0.05126221    48485.58 68
#>   2024 0.5146515 0.05184202    49637.53 68
#>  Total 0.5143723 0.05178909   220678.21 68
#> 
#> Robustness values:
#>  group direction    q threshold    xi_hat target_xi   rv_equal   rv_extreme
#>   2019  decrease 0.10        NA  2073.808  1866.427 0.02052158 0.0004297737
#>   2019  decrease 0.25        NA  2073.808  1555.356 0.05051239 0.0026800386
#>   2019  decrease 0.50        NA  2073.808  1036.904 0.09844193 0.0106346506
#>   2020  decrease 0.10        NA  2228.864  2005.978 0.02010820 0.0004124667
#>   2020  decrease 0.25        NA  2228.864  1671.648 0.04951054 0.0025723466
#>   2020  decrease 0.50        NA  2228.864  1114.432 0.09654025 0.0102105909
#>   2021  decrease 0.10        NA  2575.444  2317.900 0.01912597 0.0003727965
#>   2021  decrease 0.25        NA  2575.444  1931.583 0.04712749 0.0023254270
#>   2021  decrease 0.50        NA  2575.444  1287.722 0.09200847 0.0092372662
#>   2022  decrease 0.10        NA  3340.510  3006.459 0.01816832 0.0003360830
#>   2022  decrease 0.25        NA  3340.510  2505.382 0.04480055 0.0020968191
#>   2022  decrease 0.50        NA  3340.510  1670.255 0.08757206 0.0083348462
#>   2023  decrease 0.10        NA  3721.804  3349.624 0.01799476 0.0003296363
#>   2023  decrease 0.25        NA  3721.804  2791.353 0.04437845 0.0020566676
#>   2023  decrease 0.50        NA  3721.804  1860.902 0.08676611 0.0081762232
#>   2024  decrease 0.10        NA  3866.374  3479.736 0.01805514 0.0003318718
#>   2024  decrease 0.25        NA  3866.374  2899.780 0.04452530 0.0020705910
#>   2024  decrease 0.50        NA  3866.374  1933.187 0.08704655 0.0082312336
#>  Total  decrease 0.10        NA 17806.803 16026.123 0.01871677 0.0003568721
#>  Total  decrease 0.25        NA 17806.803 13355.103 0.04613363 0.0022262795
#>  Total  decrease 0.50        NA 17806.803  8903.402 0.09011499 0.0088460367
#>  status
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
```
