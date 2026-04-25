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

- ...:

  Additional graphical arguments passed to
  `plot.uncounted_frailty_sensitivity` when `plot = TRUE`.

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
#>   2019 2073.815
#>   2020 2228.872
#>   2021 2575.455
#>   2022 3340.524
#>   2023 3721.820
#>   2024 3866.390
#> 
#> Total xi: 17,807 
#> 
#> Working contrasts:
#>  group theta_hat   se_theta d_xi_dtheta df
#>   2019 0.5176582 0.05241505    23139.11 68
#>   2020 0.5170379 0.05229637    25443.42 68
#>   2021 0.5142050 0.05175736    31247.17 68
#>   2022 0.5139316 0.05170561    42729.37 68
#>   2023 0.5115798 0.05126226    48485.78 68
#>   2024 0.5146519 0.05184207    49637.74 68
#>  Total 0.5143727 0.05178913   220679.09 68
#> 
#> Robustness values:
#>  group direction    q threshold    xi_hat target_xi   rv_equal   rv_extreme
#>   2019  decrease 0.10        NA  2073.815  1866.434 0.02052156 0.0004297730
#>   2019  decrease 0.25        NA  2073.815  1555.362 0.05051235 0.0026800345
#>   2019  decrease 0.50        NA  2073.815  1036.908 0.09844186 0.0106346342
#>   2020  decrease 0.10        NA  2228.872  2005.985 0.02010818 0.0004124661
#>   2020  decrease 0.25        NA  2228.872  1671.654 0.04951051 0.0025723426
#>   2020  decrease 0.50        NA  2228.872  1114.436 0.09654018 0.0102105752
#>   2021  decrease 0.10        NA  2575.455  2317.909 0.01912596 0.0003727959
#>   2021  decrease 0.25        NA  2575.455  1931.591 0.04712745 0.0023254234
#>   2021  decrease 0.50        NA  2575.455  1287.727 0.09200840 0.0092372519
#>   2022  decrease 0.10        NA  3340.524  3006.472 0.01816831 0.0003360825
#>   2022  decrease 0.25        NA  3340.524  2505.393 0.04480051 0.0020968158
#>   2022  decrease 0.50        NA  3340.524  1670.262 0.08757199 0.0083348333
#>   2023  decrease 0.10        NA  3721.820  3349.638 0.01799474 0.0003296358
#>   2023  decrease 0.25        NA  3721.820  2791.365 0.04437841 0.0020566643
#>   2023  decrease 0.50        NA  3721.820  1860.910 0.08676604 0.0081762101
#>   2024  decrease 0.10        NA  3866.390  3479.751 0.01805512 0.0003318713
#>   2024  decrease 0.25        NA  3866.390  2899.792 0.04452527 0.0020705878
#>   2024  decrease 0.50        NA  3866.390  1933.195 0.08704649 0.0082312207
#>  Total  decrease 0.10        NA 17806.876 16026.188 0.01871676 0.0003568715
#>  Total  decrease 0.25        NA 17806.876 13355.157 0.04613359 0.0022262759
#>  Total  decrease 0.50        NA 17806.876  8903.438 0.09011493 0.0088460227
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
